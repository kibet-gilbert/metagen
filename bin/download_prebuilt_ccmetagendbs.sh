#!/usr/bin/env bash

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

wget --post-data="browser=false" \
     --continue --timeout=30 --tries=20 --waitretry=5 \
     -O RefSeq_bf.zip \
     "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=Lqaic1pBmpDdqX8ofv1C1128247855&browser=true&filename=RefSeq_bf.zip"

wget --post-data="browser=false" \
     --continue --timeout=30 --tries=20 --waitretry=5 \
     -O ncbi_nt_kma.zip \
     "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=i8yedNiYfdjrBfGJ8Y5z1128247857&browser=true&filename=ncbi_nt_kma.zip"

# wget --post-data="browser=false" \
#      --continue --timeout=30 --tries=20 --waitretry=5 \
#      -O ncbi_nt_no_env_11jun2019.zip \
#      "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=ko6MbZXl7FWjAS3jsItV1128247851&browser=true&filename=ncbi_nt_no_env_11jun2019.zip"

