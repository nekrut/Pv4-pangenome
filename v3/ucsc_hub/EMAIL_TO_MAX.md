Subject: hubtools up — resolved, plus two papercuts worth fixing

Hi Max,

The hub is up. The whole P. vivax pangenome (BRC_Pangenome_Pv_v1, ~5 GB, 253
files) is staged on hubSpace and replaces our Dropbox hosting — thanks for the
quota bump and the recipe.

The NXDOMAIN I was about to write you about turned out to be the DDoS:
hubspace.gi.ucsc.edu dropped out of public DNS on the 3rd and came back the next
day. Not your bug. But once it resolved, two traps cost me a day, and neither is
in the instructions — flagging in case you want to fix them in hubtools or the
help text:

1. TLS cert SAN. hubtools hard-codes https://hubspace.gi.ucsc.edu/files (line
   874). hubspace.gi is a CNAME to hubspace.soe.ucsc.edu, but the server cert
   carries only the .soe. name, so the .gi. request fails verification:
   "no alternative certificate subject name matches target host name". I worked
   around it with tusUrl=https://hubspace.soe.ucsc.edu/files in ~/.hubtools.conf.
   Either the cert needs hubspace.gi in its SAN, or the default tusUrl should be
   the .soe. name.

2. Quoted apiKey -> HTTP 500. The help text (and the errAbort message at line
   882) tells users to run echo 'apiKey="xxxx"' >> ~/.hubtools.conf. But
   parseConf (line 180) does a plain split("=",1) and does not strip quotes, so
   hubtools sends the key with the quote characters as TUS metadata and the
   pre-create hook 500s. The unquoted key returns 200. So the documented form is
   the one that fails. Either strip surrounding quotes in parseConf, or fix the
   docs to show the key without quotes. Same applies to any quoted tusUrl —
   requests then throws InvalidSchema on the leading quote.

Everything else was smooth: the smoke test, the resumable walk, the mtime cache.
The bundle is hubCheck-clean (0 errors; only warnings about optional per-track
description pages).

Thanks again,
Anton
