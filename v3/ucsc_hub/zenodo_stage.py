#!/usr/bin/env python3
"""
Stage a UCSC track hub on Zenodo (or Zenodo Sandbox) as a byte-range file host.

Zenodo has a flat file namespace, so every hub file is uploaded under a flattened
key (path separators -> '__') and every relative reference inside hub.txt,
genomes.txt and each trackDb.txt is rewritten to the absolute Zenodo file URL
(https://<host>/api/records/<id>/files/<key>/content), which serves HTTP range
requests.

Auth: set ZENODO_TOKEN (personal access token, scopes deposit:write +
deposit:actions). Base host via ZENODO_BASE (default https://zenodo.org;
use https://sandbox.zenodo.org to test).

Usage:
    ZENODO_TOKEN=xxx ZENODO_BASE=https://sandbox.zenodo.org \
        python3 zenodo_stage.py --hub /path/to/ucsc_hub [--publish]

Without --publish it creates the deposition and uploads everything but leaves it
as a draft for review (draft files are NOT anonymously readable — publish to make
the hub live). Prints the hub.txt URL to give to UCSC / hubCheck.
"""

import argparse
import os
import sys
import requests

# trackDb / genomes.txt / hub.txt settings whose value is a path to another file.
URL_SETTINGS = {
    "genomesFile", "trackDb", "groups", "twoBitPath", "htmlPath",
    "bigDataUrl", "linkDataUrl", "summary", "frames", "bigDataIndex", "searchTrix",
}
# extensions that count as hub data files to upload
DATA_EXT = (".bb", ".2bit", ".txt", ".html")


def collect_files(hub):
    files = []
    for root, _dirs, names in os.walk(hub):
        for n in names:
            if n.endswith(DATA_EXT):
                full = os.path.join(root, n)
                rel = os.path.relpath(full, hub)
                files.append(rel)
    return sorted(files)


def key_of(relpath):
    "flattened, collision-free Zenodo file key"
    return relpath.replace("/", "__")


def url_of(base, recid, relpath):
    return f"{base}/api/records/{recid}/files/{key_of(relpath)}/content"


def rewrite_text(relpath, text, base, recid, valid_rel):
    """Rewrite path-valued settings in a hub text file to absolute Zenodo URLs.
    Paths in a file are resolved relative to that file's own directory."""
    basedir = os.path.dirname(relpath)
    out = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            parts = stripped.split(None, 1)
            if len(parts) == 2 and parts[0] in URL_SETTINGS:
                key, val = parts
                target = os.path.normpath(os.path.join(basedir, val)) if basedir else val
                target = target.replace(os.sep, "/")
                if target in valid_rel:
                    indent = line[: len(line) - len(line.lstrip())]
                    out.append(f"{indent}{key} {url_of(base, recid, target)}")
                    continue
        out.append(line)
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hub", required=True)
    ap.add_argument("--publish", action="store_true")
    ap.add_argument("--title", default="P. vivax 8-strain pangenome — UCSC track hub")
    args = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    base = os.environ.get("ZENODO_BASE", "https://zenodo.org").rstrip("/")
    if not token:
        sys.exit("ZENODO_TOKEN not set")
    hub = os.path.abspath(args.hub)
    s = requests.Session()
    p = {"access_token": token}

    files = collect_files(hub)
    valid_rel = set(files)
    print(f"{len(files)} files to stage from {hub}", file=sys.stderr)

    # 1. create deposition
    r = s.post(f"{base}/api/deposit/depositions", params=p, json={})
    r.raise_for_status()
    dep = r.json()
    recid = dep["id"]
    bucket = dep["links"]["bucket"]
    print(f"deposition {recid}  bucket {bucket}", file=sys.stderr)

    text_exts = (".txt", ".html")
    # 2. upload binaries as-is; rewrite text files first
    for rel in files:
        key = key_of(rel)
        full = os.path.join(hub, rel)
        if rel.endswith(text_exts):
            with open(full, encoding="utf-8", errors="replace") as fh:
                body = rewrite_text(rel, fh.read(), base, recid, valid_rel)
            data = body.encode()
            r = s.put(f"{bucket}/{key}", data=data, params=p)
        else:
            with open(full, "rb") as fh:
                r = s.put(f"{bucket}/{key}", data=fh, params=p)
        r.raise_for_status()
        print(f"  uploaded {key}", file=sys.stderr)

    # 3. minimal metadata
    meta = {"metadata": {
        "title": args.title,
        "upload_type": "dataset",
        "description": "UCSC track hub for the BRC P. vivax 8-strain pangenome "
                       "(alignments, projected annotations, selection scans). "
                       "Load hub.txt in the UCSC Genome Browser.",
        "creators": [{"name": "Nekrutenko, Anton"}],
    }}
    r = s.put(f"{base}/api/deposit/depositions/{recid}", params=p, json=meta)
    r.raise_for_status()

    hub_url = url_of(base, recid, "hub.txt")
    if args.publish:
        r = s.post(f"{base}/api/deposit/depositions/{recid}/actions/publish", params=p)
        r.raise_for_status()
        print("PUBLISHED", file=sys.stderr)
    else:
        print("DRAFT (not published — files not yet anonymously readable)", file=sys.stderr)
    print(hub_url)


if __name__ == "__main__":
    main()
