# 8-way PGGB graph — pointers to v2 build

The actual graph files live in v2/ to avoid duplication. Recreate via the v2 build pipeline if needed.

| File | Target |
|------|--------|
| pv.og | /media/anton/data/sandbox/Pv4/v2/pggb_out/pggb_input.fa.gz.f705205.c28ecf8.dcad0e6.smooth.fix.og |
| pv.gfa | /media/anton/data/sandbox/Pv4/v2/pggb_out/pggb_input.fa.gz.f705205.c28ecf8.dcad0e6.smooth.fix.gfa |
| pv.gbz | /media/anton/data/sandbox/Pv4/v3/inputs/pggb/pv.gbz |

Recreate command (in v2/):
```bash
pggb -i v2/pggb_in/pv_panSN.fa.gz -n 8 -s 5000 -p 90 -t 32 -o v2/pggb_out/
```
