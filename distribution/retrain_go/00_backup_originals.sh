#!/usr/bin/env bash
# Back up the 3,930 own-scaler GO CAV directories before 01_retrain_go_cavs.sbatch
# overwrites them in place.  ~400 MB compressed; this is the undo button.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/backup_own_scaler_cavs_$(date +%Y%m%d).tar.zst"

[ -f "$HERE/retrain_spans.txt" ] || { echo "retrain_spans.txt missing — run classify_scalers.py first"; exit 1; }
[ -e "$OUT" ] && { echo "Backup already exists: $OUT"; exit 0; }

sed 's/\.span$/_cav/' "$HERE/retrain_spans.txt" > "$HERE/.backup_list.txt"
n=$(wc -l < "$HERE/.backup_list.txt")
echo "Backing up $n CAV directories -> $OUT"
tar --zstd -cf "$OUT" -T "$HERE/.backup_list.txt"
rm -f "$HERE/.backup_list.txt"
ls -lh "$OUT"
echo "Restore with:  tar --zstd -xf $OUT -C /"
