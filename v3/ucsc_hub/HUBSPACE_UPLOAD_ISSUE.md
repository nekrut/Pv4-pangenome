# hubtools upload — `hubspace.gi.ucsc.edu` does not resolve (NXDOMAIN)

> **RESOLVED 2026-06-04.** The NXDOMAIN was transient — campus was under a DDoS
> on 2026-06-03 and `hubspace.gi.ucsc.edu` dropped out of public DNS. It
> recovered the next day (CNAME → `hubspace.soe.ucsc.edu`, 128.114.119.144), and
> the full hub is now staged. Two further traps surfaced once DNS returned (TLS
> cert SAN, quoted apiKey → HTTP 500). The working procedure and all three fixes
> are in **`HUBSPACE_STAGING.md`**; the upstream-feedback note is
> **`EMAIL_TO_MAX.md`**. This file is kept as the record of the original outage.

**Date:** 2026-06-03
**Reporter:** Anton Nekrutenko (anton@nekrut.org)
**Goal:** stage the *P. vivax* 8-strain track hub (~16 GB of bigMaf/bigChain/bigBed) on UCSC hubSpace via `hubtools up`, replacing the current Dropbox hosting.

## Summary

`hubtools up` fails at the TUS upload step because its hard-coded server host
**`hubspace.gi.ucsc.edu` returns NXDOMAIN in public DNS** — from both our local
resolver and Google `8.8.8.8`. Every other UCSC host we use resolves normally.
Config, apiKey, file walk, and `tuspy` all work; the only failure is name
resolution of the upload endpoint.

## Setup (all confirmed working)

- `hubtools` fetched from `https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/utils/hubtools/hubtools`, `chmod a+x`.
- `~/.hubtools.conf` contains `apiKey="…"` (key created at genome.ucsc.edu → My Data → Track Hubs → Track Development; mode 600). `hubtools` reads this file first (line 210), `~/.hg.conf` only as fallback (line 214).
- `tuspy` (tusclient) already installed under python 3.13.
- Quota raised to 1 TB on our account.

## What we ran (Max's recommended smoke test)

```
mkdir /tmp/hubtest
cp one-small.bigChain.bb /tmp/hubtest/      # 52 KB
cd /tmp/hubtest
hubtools up test
```

## Failure

```
tusclient.exceptions.TusCommunicationError:
HTTPSConnectionPool(host='hubspace.gi.ucsc.edu', port=443):
Max retries exceeded with url: /files
(Caused by NameResolutionError("Failed to resolve 'hubspace.gi.ucsc.edu'
([Errno -2] Name or service not known)"))
```

The endpoint is hard-coded in `hubtools` (line 874):

```python
serverUrl = cfgOption("tusUrl", "https://hubspace.gi.ucsc.edu/files")
```

## Diagnosis — it is the host, not our network

DNS lookups for `hubspace.gi.ucsc.edu` (local resolver `127.0.0.53` and public `8.8.8.8`):

```
** server can't find hubspace.gi.ucsc.edu: NXDOMAIN     (127.0.0.53)
** server can't find hubspace.gi.ucsc.edu: NXDOMAIN     (8.8.8.8)
```

`curl https://hubspace.gi.ucsc.edu/files` → `curl: (6) Could not resolve host`.

Other UCSC hosts resolve fine from the same machine at the same time:

| Host                        | Resolves?      | DNS record                          |
| --------------------------- | -------------- | ----------------------------------- |
| `genome.ucsc.edu`           | yes            | 128.114.119.131 / .132 (A)          |
| `hgdownload.soe.ucsc.edu`   | yes            | 128.114.119.163 (A)                 |
| `hgwdev.gi.ucsc.edu`        | yes            | 2607:f5f0:136:1::26 (AAAA only)     |
| `gi.ucsc.edu`               | yes            | 2607:f5f0:136:1::20 (AAAA only)     |
| `hubspace.gi.ucsc.edu`      | **NXDOMAIN**   | none (neither A nor AAAA)           |

The `hgHubConnect` web UI on `genome.ucsc.edu` returns HTTP 200, so the web side
is reachable — only the `hubspace` upload host is unresolvable. Note the
`gi.ucsc.edu` zone publishes **only IPv6 (AAAA)** records, while `hubspace`
publishes nothing.

## Questions for UCSC

1. Is `hubspace.gi.ucsc.edu` the correct public upload endpoint, or has it moved? If moved, what `tusUrl` should we set in `~/.hubtools.conf`?
2. Is the host intentionally internal-only (UCSC network / hgwdev), so external `hubtools up` is not yet supported for our account?
3. The `gi.ucsc.edu` zone is IPv6-only — does the upload endpoint require IPv6 connectivity from the client?

## What we can do once the endpoint resolves

Re-run `hubtools up test` (the 52 KB smoke test) to confirm the file lands under
My Data → Track Hubs → Hub Upload, then `hubtools up <hubName>` over the full
`ucsc_hub/` bundle (resumable — `hubtools` skips already-uploaded files by mtime,
so the 16 GB transfer can be interrupted and resumed).

Our UCSC username: **(to fill in)**.
