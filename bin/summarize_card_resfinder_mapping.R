#!/usr/bin/env Rscript
# ============================================================================
# summarize_mapping.R
# ----------------------------------------------------------------------------
# Summary visualizations for the CARD <-> ResFinder merged gene-name dictionary
# produced by build_mapping.py.
#
# WHAT THIS PRODUCES
# -------------------
#   01_match_status_overview.<fmt>    How many of CARD's 2,713 ARO terms were
#                                      matched, left unmatched, or are
#                                      categorically out of ResFinder's
#                                      acquired-gene scope (chromosomal
#                                      mutation / regulatory models).
#   02_match_tier_breakdown.<fmt>     Which matching rule found each match
#                                      (exact, bla-prefix stripping, van-cluster
#                                      structural rule, historical-name alias
#                                      bridge, ...) and its confidence.
#   03_resfinder_class_distribution.<fmt>
#                                      Which ResFinder resistance classes
#                                      (Beta-Lactamase, Aminoglycoside, ...)
#                                      account for the matched genes.
#   04_drugclass_coverage.<fmt>       Per CARD "ARO Categories" drug/mechanism
#                                      term, how much of CARD's catalogue has a
#                                      ResFinder counterpart vs. not.
#   05_prevalence_by_match_status.<fmt>
#                                      Does CARD's NCBI whole-genome-sequence
#                                      prevalence differ between genes that do
#                                      and don't have a ResFinder counterpart?
#   06_summary_dashboard.<fmt>        Panels 1, 2, 3, and 5 combined into one
#                                      at-a-glance figure (via patchwork).
#
# INPUT
# -----
#   --input      card_resfinder_merged_dictionary.csv (required; produced by
#                build_mapping.py)
#   --ambiguous  ambiguous_cases_review.csv (optional; only used for a short
#                console summary of the 21 multi-candidate cases)
#
# REQUIREMENTS
# ------------
#   R packages: optparse, readr, dplyr, tidyr, stringr, forcats, ggplot2,
#   scales, patchwork. The script will try to install any that are missing
#   via install.packages(). If that fails because your environment has no
#   CRAN access, install the Debian/Ubuntu binary packages instead, e.g.:
#     sudo apt-get install r-cran-ggplot2 r-cran-dplyr r-cran-tidyr \
#          r-cran-readr r-cran-stringr r-cran-forcats r-cran-scales \
#          r-cran-optparse r-cran-patchwork
#
# USAGE EXAMPLES
# --------------
#   # Defaults: reads from /mnt/user-data/outputs, writes plots alongside it
#   Rscript summarize_mapping.R
#
#   # Explicit paths, PDF output, larger figures
#   Rscript summarize_mapping.R --input results/merged.csv \
#       --outdir results/plots --format pdf --width 10 --height 6
#
#   Run `Rscript summarize_mapping.R --help` for the full option list.
# ============================================================================

# ---- 0. Package management -------------------------------------------------

required_pkgs <- c(
  "optparse", "readr", "dplyr", "tidyr", "stringr",
  "forcats", "ggplot2", "scales", "patchwork"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing R packages: ", paste(missing_pkgs, collapse = ", "))
  tryCatch(
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org"),
    error = function(e) {
      stop(
        "Could not install package(s) automatically (", conditionMessage(e), ").\n",
        "If this machine has no CRAN access, install the Debian/Ubuntu binaries instead:\n",
        "  sudo apt-get install ",
        paste0("r-cran-", missing_pkgs, collapse = " "),
        call. = FALSE
      )
    }
  )
}

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ---- 1. Command-line interface ----------------------------------------------

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "/mnt/user-data/outputs/card_resfinder_merged_dictionary.csv",
              dest = "input_csv",
              help = "Path to card_resfinder_merged_dictionary.csv, produced by build_mapping.py. [default: %default]"),
  make_option(c("-a", "--ambiguous"), type = "character",
              default = "/mnt/user-data/outputs/ambiguous_cases_review.csv",
              dest = "ambiguous_csv",
              help = "Path to ambiguous_cases_review.csv (optional -- used only for a short console summary; pass an invalid path or omit the file to skip it). [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "/mnt/user-data/outputs/plots", dest = "outdir",
              help = "Directory to write plot files to (created if missing). [default: %default]"),
  make_option(c("--format"), type = "character", default = "png", dest = "format",
              help = "Output image format: 'png' or 'pdf'. [default: %default]"),
  make_option(c("--width"), type = "double", default = 8.5, dest = "width",
              help = "Plot width in inches (applies to all individual panels). [default: %default]"),
  make_option(c("--height"), type = "double", default = 5.5, dest = "height",
              help = "Plot height in inches (applies to all individual panels). [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300, dest = "dpi",
              help = "Resolution in dots per inch (PNG output only). [default: %default]"),
  make_option(c("--top-n"), type = "integer", default = 12, dest = "top_n",
              help = "Number of top drug/mechanism categories to show in the coverage plot. [default: %default]")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Generate summary plots from the CARD<->ResFinder merged gene-name dictionary (build_mapping.py output).",
  epilogue = "\nExample: Rscript summarize_mapping.R --input results/merged.csv --outdir results/plots --format pdf\n"
))

if (!file.exists(opt$input_csv)) {
  stop(
    "Input CSV not found: '", opt$input_csv, "'\n",
    "  Run build_mapping.py first, or pass --input <path to card_resfinder_merged_dictionary.csv>.",
    call. = FALSE
  )
}
if (!opt$format %in% c("png", "pdf")) {
  stop("--format must be 'png' or 'pdf', got: '", opt$format, "'", call. = FALSE)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load & clean ---------------------------------------------------------

#' Read the merged dictionary CSV and tidy it for plotting.
#'
#' Converts build_mapping.py's empty-string "no value" convention to proper
#' NA, and turns Match_Status / Match_Tier / Match_Confidence into ordered
#' factors so plots come out in a sensible, consistent order rather than
#' alphabetical.
#'
#' @param path Path to card_resfinder_merged_dictionary.csv.
#' @return A tidied tibble.
load_merged_dictionary <- function(path) {
  # Explicit column types, rather than readr's row-sampling auto-detection:
  # several text columns (e.g. Excluded_Candidates_ClassMismatch) are non-empty
  # for only a handful of rows out of 2,713, so relying on readr's default
  # guess_max risks it sampling an all-blank prefix and mistyping the column
  # (observed in testing: it was inferred as logical, not character).
  col_spec <- cols(
    .default = col_character(),
    N_Pathogens_With_Hit = col_double(),
    Max_NCBI_WGS_Prevalence_Pct = col_double(),
    Has_Perfect_Strict_Hit = col_logical()
  )
  df <- readr::read_csv(path, col_types = col_spec, progress = FALSE)

  status_levels <- c(
    "matched", "no_reliable_match",
    "out_of_scope_chromosomal_point_mutation",
    "out_of_scope_rRNA_mutation",
    "out_of_scope_regulatory_overexpression"
  )
  status_labels <- c(
    "Matched", "No reliable match",
    "Out of scope: chromosomal mutation",
    "Out of scope: rRNA mutation",
    "Out of scope: regulatory/overexpression"
  )
  tier_levels <- c("exact", "case_insensitive", "van_structural", "loose", "debla", "alias_bridge")
  tier_labels <- c(
    "Exact string match", "Case-insensitive match", "van-cluster structural rule",
    "Punctuation-normalized match", "bla-prefix stripped", "Historical-name alias bridge"
  )

  df %>%
    mutate(
      across(c(ResFinder_Gene_Name, ResFinder_Class, Match_Tier, Match_Confidence,
               Alternate_ResFinder_Synonyms, Excluded_Candidates_ClassMismatch),
             ~ na_if(., "")),
      Match_Status = factor(Match_Status, levels = status_levels, labels = status_labels),
      Match_Tier_Label = factor(Match_Tier, levels = tier_levels, labels = tier_labels),
      Match_Confidence = factor(Match_Confidence, levels = c("High", "Medium"))
    )
}

df <- load_merged_dictionary(opt$input_csv)
cat(sprintf("Loaded %s (%d rows) from %s\n", basename(opt$input_csv), nrow(df), opt$input_csv))

# ---- 3. Shared plot theme ----------------------------------------------------

#' Consistent minimal theme used across every plot in this script.
theme_report <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle = element_text(color = "grey35", margin = margin(b = 10)),
      plot.caption = element_text(color = "grey50", size = rel(0.75), hjust = 0),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "grey25"),
      legend.position = "bottom",
      legend.title = element_text(size = rel(0.85))
    )
}

#: Consistent palette: greys for "not matched" states, a warm accent for
#: "matched", used wherever match-vs-not is the key contrast.
PAL_MATCHED <- c("Matched" = "#1b7837", "No reliable match" = "#b2182b",
                  "Out of scope: chromosomal mutation" = "#878787",
                  "Out of scope: rRNA mutation" = "#b0b0b0",
                  "Out of scope: regulatory/overexpression" = "#d0d0d0")
PAL_CONFIDENCE <- c("High" = "#2166ac", "Medium" = "#92c5de")

n_total <- nrow(df)
caption_source <- stringr::str_wrap(
  sprintf(
    "Source: CARD RGI wildcard prevalence x ResFinder DB, reconciled by build_mapping.py  |  n = %s CARD ARO terms",
    scales::comma(n_total)
  ),
  width = 100
)

#' Save a ggplot (or patchwork composite) to opt$outdir in the requested format.
#'
#' @param plot A ggplot/patchwork object.
#' @param name File basename (without extension).
#' @param width,height Override the global --width/--height for this plot.
save_plot <- function(plot, name, width = opt$width, height = opt$height) {
  path <- file.path(opt$outdir, paste0(name, ".", opt$format))
  ggsave(
    filename = path, plot = plot, width = width, height = height,
    dpi = opt$dpi, bg = "white"
  )
  cat("  wrote", path, "\n")
  invisible(path)
}

# ---- 4. Plot 1: match status overview ----------------------------------------

plot_match_status <- function(df) {
  counts <- df %>%
    count(Match_Status, .drop = FALSE) %>%
    mutate(pct = n / sum(n))

  ggplot(counts, aes(x = fct_rev(Match_Status), y = n, fill = Match_Status)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%s  (%s)", scales::comma(n), scales::percent(pct, accuracy = 0.1))),
              hjust = -0.05, size = 3.4, color = "grey20") +
    scale_fill_manual(values = PAL_MATCHED) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
    coord_flip(clip = "off") +
    labs(
      title = "Where CARD's 2,713 ARO terms landed",
      subtitle = "Match status against ResFinder DB's acquired-gene catalogue",
      x = NULL, y = "Number of CARD ARO terms",
      caption = caption_source
    ) +
    theme_report()
}

# ---- 5. Plot 2: match tier breakdown (matched entries only) -----------------

plot_match_tier <- function(df) {
  matched <- df %>% filter(Match_Status == "Matched")
  counts <- matched %>% count(Match_Tier_Label, Match_Confidence)

  ggplot(counts, aes(x = fct_reorder(Match_Tier_Label, n, .fun = sum), y = n, fill = Match_Confidence)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = scales::comma(n)), hjust = -0.15, size = 3.3, color = "grey20") +
    scale_fill_manual(values = PAL_CONFIDENCE, name = "Match confidence") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_flip(clip = "off") +
    labs(
      title = "Which rule found each match",
      subtitle = sprintf("%s matched CARD ARO terms, by matching tier", scales::comma(nrow(matched))),
      x = NULL, y = "Number of matched ARO terms",
      caption = caption_source
    ) +
    theme_report()
}

# ---- 6. Plot 3: ResFinder class distribution among matched genes -----------

plot_resfinder_class <- function(df, top_n = 15) {
  counts <- df %>%
    filter(Match_Status == "Matched", !is.na(ResFinder_Class)) %>%
    count(ResFinder_Class, sort = TRUE) %>%
    slice_max(n, n = top_n)

  ggplot(counts, aes(x = fct_reorder(ResFinder_Class, n), y = n)) +
    geom_col(fill = "#2166ac", width = 0.7) +
    geom_text(aes(label = scales::comma(n)), hjust = -0.15, size = 3.3, color = "grey20") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    coord_flip(clip = "off") +
    labs(
      title = "Matched genes by ResFinder resistance class",
      subtitle = "Which classes account for the CARD<->ResFinder overlap",
      x = NULL, y = "Number of matched CARD ARO terms",
      caption = caption_source
    ) +
    theme_report()
}

# ---- 7. Plot 4: drug/mechanism-category coverage ----------------------------

plot_drugclass_coverage <- function(df, top_n = 12) {
  cat_df <- df %>%
    filter(!is.na(ARO_Categories), ARO_Categories != "") %>%
    tidyr::separate_rows(ARO_Categories, sep = ";\\s*") %>%
    mutate(matched_flag = ifelse(Match_Status == "Matched", "Matched to ResFinder", "No ResFinder match"))

  top_cats <- cat_df %>%
    count(ARO_Categories, sort = TRUE) %>%
    slice_max(n, n = top_n) %>%
    pull(ARO_Categories)

  plot_df <- cat_df %>%
    filter(ARO_Categories %in% top_cats) %>%
    count(ARO_Categories, matched_flag) %>%
    group_by(ARO_Categories) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    mutate(matched_flag = factor(matched_flag, levels = c("No ResFinder match", "Matched to ResFinder")))

  ggplot(plot_df, aes(x = fct_reorder(ARO_Categories, total), y = n, fill = matched_flag)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Matched to ResFinder" = "#1b7837", "No ResFinder match" = "#d6d6d6"),
                       name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(
      title = sprintf("ResFinder coverage of CARD's top %d drug/mechanism categories", top_n),
      subtitle = "A CARD ARO term can carry more than one category, so totals overlap across bars",
      x = NULL, y = "Number of CARD ARO terms",
      caption = caption_source
    ) +
    theme_report()
}

# ---- 8. Plot 5: NCBI prevalence by match status -----------------------------

plot_prevalence_by_status <- function(df) {
  plot_df <- df %>% filter(!is.na(Max_NCBI_WGS_Prevalence_Pct), Max_NCBI_WGS_Prevalence_Pct > 0)
  n_dropped <- sum(is.na(df$Max_NCBI_WGS_Prevalence_Pct) | df$Max_NCBI_WGS_Prevalence_Pct <= 0)

  ggplot(plot_df, aes(x = fct_rev(Match_Status), y = Max_NCBI_WGS_Prevalence_Pct, fill = Match_Status)) +
    geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.35, width = 0.55, show.legend = FALSE) +
    scale_fill_manual(values = PAL_MATCHED) +
    scale_y_log10(labels = scales::label_percent(scale = 1, accuracy = 0.01)) +
    coord_flip() +
    labs(
      title = "Does ResFinder coverage track with real-world prevalence?",
      subtitle = stringr::str_wrap(
        "CARD's max NCBI whole-genome-sequence detection %, by match status (log scale)",
        width = 75
      ),
      x = NULL, y = "Max NCBI WGS prevalence (%, log scale)",
      caption = stringr::str_wrap(
        sprintf("%s  |  %d ARO terms with zero/NA prevalence omitted from this panel",
                caption_source, n_dropped),
        width = 100
      )
    ) +
    theme_report()
}

# ---- 9. Build, save, and combine ---------------------------------------------

cat("\nBuilding plots ...\n")
p1 <- plot_match_status(df)
p2 <- plot_match_tier(df)
p3 <- plot_resfinder_class(df)
p4 <- plot_drugclass_coverage(df, top_n = opt$top_n)
p5 <- plot_prevalence_by_status(df)

save_plot(p1, "01_match_status_overview")
save_plot(p2, "02_match_tier_breakdown")
save_plot(p3, "03_resfinder_class_distribution")
save_plot(p4, "04_drugclass_coverage", height = opt$height * 1.15)
save_plot(p5, "05_prevalence_by_match_status")

dashboard <- (p1 + labs(subtitle = NULL, caption = NULL)) +
  (p2 + labs(subtitle = NULL, caption = NULL)) +
  (p3 + labs(subtitle = NULL, caption = NULL)) +
  (p5 + labs(title = "Prevalence by match status", subtitle = NULL, caption = NULL)) +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(
    title = "CARD \u2194 ResFinder gene-name reconciliation: summary",
    subtitle = caption_source,
    theme = theme(
      plot.title = element_text(face = "bold", size = rel(1.3)),
      plot.subtitle = element_text(color = "grey35")
    )
  )
save_plot(dashboard, "06_summary_dashboard", width = opt$width * 1.7, height = opt$height * 1.7)

# ---- 10. Console summary (including the ambiguous-cases file, if present) ---

cat("\n==================== Summary ====================\n")
status_tbl <- df %>% count(Match_Status, name = "n") %>% mutate(pct = scales::percent(n / sum(n), accuracy = 0.1))
print(as.data.frame(status_tbl), row.names = FALSE)

if (!is.null(opt$ambiguous_csv) && file.exists(opt$ambiguous_csv)) {
  amb <- readr::read_csv(opt$ambiguous_csv, col_types = cols(.default = col_character()), progress = FALSE)
  cat(sprintf(
    "\n%d CARD entries had more than one raw ResFinder candidate before re-ranking\n(see %s for full detail).\n",
    nrow(amb), opt$ambiguous_csv
  ))
} else {
  cat("\n(No ambiguous-cases file found at '", opt$ambiguous_csv, "' -- skipped.)\n", sep = "")
}

cat("\nAll plots written to:", opt$outdir, "\n")
