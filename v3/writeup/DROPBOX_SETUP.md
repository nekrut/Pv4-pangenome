# Dropbox setup for the 342 GB v3 archive upload

Both rclone and the official Dropbox client are installed to user-space (no sudo).

## Option A: rclone (recommended for this one-shot upload)

### One-time OAuth config
```bash
rclone config
```
Interactive prompts:
1. `n` (new remote)
2. name: `dropbox`
3. storage: `dropbox` (autocompletes)
4. Leave `client_id` / `client_secret` blank (rclone uses its own app credentials)
5. Edit advanced config: `n`
6. Use auto config: `y` → opens browser → sign in → grant → returns
7. Confirm: `y`, then `q`

Test:
```bash
rclone lsd dropbox:
rclone mkdir dropbox:Pv4_v3
```

### Upload command
```bash
rclone copy /media/anton/scratch/Pv4_dropbox_staging dropbox:Pv4_v3 \
  --progress --transfers 4 --checkers 8 \
  --log-file /tmp/rclone_pv4_upload.log
```

For ~342 GB on a 50 Mbps upload link: ~16 hours wall. On a 1 Gbps link: ~1 hour.

### Verify after upload
```bash
rclone check /media/anton/scratch/Pv4_dropbox_staging dropbox:Pv4_v3 \
  --one-way --combined /tmp/rclone_check.txt
```

## Option B: Official Dropbox client (daemon-based, for ongoing sync)

### First run
```bash
~/.dropbox-dist/dropboxd
```
Prints a URL to a sign-in page. Open it in a browser, sign in, grant access. Daemon then runs in the foreground (Ctrl-C to detach is NOT supported — leave the terminal open, or use systemd / nohup).

The daemon creates `~/Dropbox/` and syncs everything bidirectionally.

### Headless control via dropbox.py
```bash
~/.local/bin/dropbox.py start              # start daemon in background
~/.local/bin/dropbox.py status             # syncing status
~/.local/bin/dropbox.py filestatus PATH    # per-file sync state
~/.local/bin/dropbox.py exclude add PATH   # don't sync this
~/.local/bin/dropbox.py stop               # stop daemon
```

### Selective sync caveat
The default `~/Dropbox/` folder syncs to ALL devices on your account. For a one-time 342 GB archive that doesn't need to live on your laptop:
1. After OAuth, use `dropbox.py exclude add` to stop syncing folders you don't want locally
2. Or use rclone instead (no auto-sync)

## Why rclone for THIS upload
- 342 GB into ~/Dropbox would replicate to every device on your account
- The official client uses ~250 MB RAM continuously
- rclone does one-shot uploads then exits

## File inventory
- Local staging: `/media/anton/scratch/Pv4_dropbox_staging/` (~342 GB)
- MD5 manifest: `writeup/LARGE_FILES_DROPBOX.tsv` (132 entries)
- Suggested Dropbox path: `Pv4_v3/`
