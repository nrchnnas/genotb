#!/bin/bash
cd /sessions/gifted-happy-davinci/mnt/outputs
mkdir -p cryptic/vcf
dl1(){ run="$1"; url="$2"; out="cryptic/vcf/${run}.vcf.gz";
  [ -s "$out" ] && exit 0
  curl -s --max-time 20 "$url" -o "$out.tmp" && mv "$out.tmp" "$out" 2>/dev/null || rm -f "$out.tmp" 2>/dev/null; }
export -f dl1
awk -F'\t' '{print $1" "$2}' cryptic/manifest_full.tsv | \
  xargs -P 40 -n 2 bash -c 'dl1 "$0" "$1"'
