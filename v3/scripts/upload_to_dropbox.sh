#!/bin/bash
# Upload the v3 staging area to Dropbox via rclone.
# Idempotent — re-runs skip already-uploaded files.
set -u

SRC=/media/anton/scratch/Pv4_dropbox_staging
DST=dropbox:Pv4_v3
LOG=/tmp/rclone_pv4_upload.log

if [ ! -x ~/.local/bin/rclone ]; then echo "rclone not installed"; exit 1; fi
if ! ~/.local/bin/rclone listremotes 2>/dev/null | grep -q "^dropbox:"; then
  echo "rclone remote 'dropbox:' not configured. Run: rclone config"
  exit 1
fi

echo "Uploading $(du -sh $SRC | cut -f1) → $DST"
echo "Log: $LOG"
echo

~/.local/bin/rclone copy "$SRC" "$DST" \
  --progress \
  --transfers 4 \
  --checkers 8 \
  --tpslimit 12 \
  --retries 5 \
  --low-level-retries 10 \
  --log-file "$LOG" \
  --log-level INFO \
  --stats 30s

echo
echo "=== Upload complete. Verifying ==="
~/.local/bin/rclone check "$SRC" "$DST" --one-way --combined /tmp/rclone_check.txt 2>&1 | tail -10
echo
echo "Combined check report: /tmp/rclone_check.txt"
grep -c "^= " /tmp/rclone_check.txt 2>/dev/null && echo "files OK"
grep -c "^[!*+-]" /tmp/rclone_check.txt 2>/dev/null && echo "files with problems"
