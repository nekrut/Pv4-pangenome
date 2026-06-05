# Staging the hub on UCSC hubSpace (working runbook)

**Status:** working as of 2026-06-04. The full *P. vivax* 8-strain hub
(`BRC_Pangenome_Pv_v1`, 4.97 GB, 253 files) is staged on UCSC hubSpace via
`hubtools up`. This replaces the Dropbox hosting, which cannot serve HTTP
byte-range requests and so is rejected by the UCSC Browser.

This is the recipe that works, plus the three traps we hit. Read the traps —
two of them are not in the UCSC instructions and cost a day.

## Build the hub with the `trackhub` skill

Max Haeussler maintains a Claude skill that captures the UCSC rules for every
hub file format (bigBed, bigWig, bigMaf, bigChain, 2bit, …). Build and edit the
hub through it — it is the authority, not memory.

- Source: `https://hgwdev.gi.ucsc.edu/~max/anton/skills/trackhub/SKILL.md`
- Installed at: `~/.claude/skills/trackhub/SKILL.md`
- Refresh it before a hub job (Max updates it): re-fetch and overwrite the
  installed copy.

**bigMaf `speciesOrder` — fixed 2026-06-05.** Max found the multiz tracks loaded
but rendered under the wrong sequence names: the assembly names *inside* each
bigMaf were correct, but `speciesOrder` in the per-genome `trackDb.txt` listed
strain aliases (PAM, PvW1, …) the browser could not match. The browser keys on
the first dotted component of each MAF `s` line, and our build wrote query
sources as `GCA_xxxxx.<version>.<chrom>` — so that component is the version-less
accession (`GCA_949152365` from `s GCA_949152365.1.CASCJQ010000001.1`). Fix:
`speciesOrder` now lists the eight version-less `GCA_*` tokens, with
`speciesLabels GCA_900093555="PvP01" …` carrying the strain names for display.
A trackDb-only change — the bigMafs were untouched, so only the eight small
trackDb.txt files were re-uploaded, and hubCheck on the live hub is clean
(0 errors; only optional description-page warnings). The skill's "assembly name"
wording was ambiguous on this point; Max is revising it.

## Prerequisites

- `hubtools` (UCSC kent utility):
  ```
  wget https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/utils/hubtools/hubtools
  chmod a+x hubtools
  ```
- `tuspy` (the TUS client): `pip install tuspy`.
- An apiKey created at genome.ucsc.edu → My Data → Track Hubs → **Track
  Development** (the `#dev` tab), while logged in. The hub lists under the
  **same** account that minted the key.

## `~/.hubtools.conf`

Two lines, mode 600. **The apiKey value must NOT be quoted** (see Trap 2), and
`tusUrl` must point at the `.soe.` host (see Trap 1):

```
apiKey=PASTE_YOUR_KEY_HERE
tusUrl=https://hubspace.soe.ucsc.edu/files
```

```
chmod 600 ~/.hubtools.conf
```

The key is a credential — keep this file 600, never commit it, never paste it
into a tracked file or a ticket.

## Smoke test

```
mkdir /tmp/hubtest && cp one-small.bb /tmp/hubtest/test.bb
cd /tmp/hubtest && hubtools up test
```

Clean exit, no traceback = the path works end to end. The file lands under
My Data → Track Hubs → the hubSpace/upload section, in folder `test`.

## Full upload

Run from the hub root (the directory holding `hub.txt`):

```
cd /path/to/ucsc_hub
hubtools up BRC_Pangenome_Pv_v1
```

`hubtools` walks the whole tree and uploads every non-dotfile, preserving the
relative paths under a top folder named by the argument (`BRC_Pangenome_Pv_v1`).
It records each finished file in `.hubtools.files.json` by mtime, so the
transfer is **resumable** — kill it and re-run, and it skips what is already up.

Observed throughput from our network: ~2.2 MB/s (≈18 Mbit/s); the 4.97 GB
bundle took ~40 min.

Note: `up` uploads *everything* in the tree, so the build helpers that live
beside the data (`*.py`, `*.chain.gz`, `*.as`) ride along. They are not
referenced by `hub.txt` and do no harm — just dead weight. Move them out first
if you want a clean upload.

## Public URL

`hubtools up` does **not** print the public URL and does **not** create a
visible "connected hub" entry — it only stores files. The serving URL is keyed
to your account (a shard prefix plus username). For this hub it is:

```
https://genome.ucsc.edu/hubspace/a6/anekrut/BRC_Pangenome_Pv_v1/hub.txt
```

Load it in the browser with:

```
https://genome.ucsc.edu/cgi-bin/hgTracks?genome=GCA_900093555.2&hubUrl=https://genome.ucsc.edu/hubspace/a6/anekrut/BRC_Pangenome_Pv_v1/hub.txt
```

Your own path is at My Data → Track Hubs (the hubSpace/upload section shows the
loadable URL). That `hub.txt` URL is what you give to collaborators and to
`hubCheck`.

## Validation

```
hubCheck https://hubspace.gi.ucsc.edu/<account-hash>/BRC_Pangenome_Pv_v1/hub.txt
```

The assembled bundle was already validated offline over a byte-range HTTP server
before upload: 0 hard errors, the only messages warnings about optional
per-track description pages.

---

## Traps

### Trap 1 — TLS cert covers only `hubspace.soe.ucsc.edu`

`hubtools` hard-codes the endpoint at line 874:

```python
serverUrl = cfgOption("tusUrl", "https://hubspace.gi.ucsc.edu/files")
```

`hubspace.gi.ucsc.edu` is a CNAME to `hubspace.soe.ucsc.edu` (128.114.119.144),
but the server's certificate carries only the `.soe.` name in its SAN. A request
to the `.gi.` name fails TLS verification:

```
SSL: no alternative certificate subject name matches target host name 'hubspace.gi.ucsc.edu'
```

**Fix:** override the endpoint to the `.soe.` name in `~/.hubtools.conf`:
`tusUrl=https://hubspace.soe.ucsc.edu/files`.

### Trap 2 — a quoted apiKey gives HTTP 500

UCSC's instructions say to write `apiKey="xxxx"` (with quotes). But `hubtools`'
config parser (`parseConf`, line 180) does a plain `split("=", 1)` and does
**not** strip quotes, so it sends the key *including* the quote characters as TUS
metadata. The server's pre-create hook rejects it:

```
TusCommunicationError: ... status 500
```

The unquoted key passes the same hook and returns 200. **Fix:** write the apiKey
with no quotes (`apiKey=E8SS…`). The quote-vs-no-quote 500-vs-200 split is itself
proof the server validates the key.

### Trap 3 — `tusUrl` value must also be unquoted

Same parser, same cause: `tusUrl="https://…"` is read with the quotes attached,
and `requests` then throws `InvalidSchema: No connection adapters were found for
'"https://…"'`. Write `tusUrl=https://…` bare.

---

## History

The endpoint was unreachable on 2026-06-03 (NXDOMAIN on `hubspace.gi.ucsc.edu`)
during a campus DDoS; DNS recovered on 2026-06-04. See `HUBSPACE_UPLOAD_ISSUE.md`
for that record and `EMAIL_TO_MAX.md` for the two papercuts worth fixing upstream
(cert SAN, quote parsing).
