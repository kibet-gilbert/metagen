// modules/cardrgi_db/main.nf
process CARDRGI_DB {
    tag "CARDRGI_DB+${cversion}+${wversion}"

    input:
    val cversion
    val wversion

    output:
    path "./card.json", emit: cardJson
    path "./wildcard", emit: wildcardDir
    path "./localDB", emit: localDBDir

    script:
    """
    set -euo pipefail

    # echo "[\$(date '+%F %T')] Downloading CARD data (v${cversion})..."
    wget -qO data.tgz https://card.mcmaster.ca/latest/data
    tar -xzf data.tgz card.json

    # echo "[\$(date '+%F %T')] Downloading WildCARD data (v${wversion})..."
    wget -qO wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
    tar -xjf wildcard_data.tar.bz2 -C wildcard
    gunzip wildcard/*.gz || true

    # echo "[\$(date '+%F %T')] Creating annotations..."
    rgi card_annotation -i card.json > card_annotation_v${cversion}.log 2>&1
    rgi wildcard_annotation \
        -i wildcard \
        --card_json card.json \
        -v ${wversion} > wildcard_annotation_v${wversion}.log 2>&1

    # echo "[\$(date '+%F %T')] Loading CARD and WildCARD into RGI..."
    rgi load \
        --card_json card.json \
        --local --debug \
        --card_annotation card_database_v${cversion}.fasta \
        --card_annotation_all_models card_database_v${cversion}_all.fasta \
        --wildcard_annotation wildcard_database_v${wversion}.fasta \
        --wildcard_annotation_all_models wildcard_database_v${wversion}_all.fasta \
        --wildcard_index wildcard/index-for-model-sequences.txt \
        --wildcard_version ${wversion} \
        --amr_kmers wildcard/all_amr_61mers.txt \
        --kmer_database wildcard/61_kmer_db.json \
        --kmer_size 61

    # echo "[\$(date '+%F %T')] CARDRGI_DB is ready."
    """
}
