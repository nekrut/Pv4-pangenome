# Directories still needing sudo for removal

These are intermediate dirs containing root-owned files (created during docker runs).
Need `sudo rm -rf` to clean up before commit.

```bash
sudo rm -rf \
  /media/anton/data/sandbox/Pv4/v3/intermediate_files \
  /media/anton/data/sandbox/Pv4/v3/TOGA2_ref_annotation_17:47_19.05.26_8cdfad94da \
  /media/anton/data/sandbox/Pv4/v3/toga2_run_17:48_19.05.26 \
  /media/anton/data/sandbox/Pv4/v3/work/02b_triage \
  /media/anton/data/sandbox/Pv4/v3/work/02c_toga \
  /media/anton/data/sandbox/Pv4/v3/work/08_genespace
```

All recoverable via:
- `work/02b_triage/`: `scripts/phase_c2_triage.py`
- `work/02c_toga/*-as-ref/`: `scripts/run_phase_c3_toga2.sh`
- `work/08_genespace/`: skipped (substituted with OrthoFinder3 — see `work/08_orthofinder/`)
- `intermediate_files/`, `TOGA2_ref_annotation_*/`, `toga2_run_*/`: TOGA2 internal workdirs

## Additional root-/uid 2000-owned dirs found mid-cleanup

```bash
sudo rm -rf \
  /media/anton/data/sandbox/Pv4/v3/work/00_inventory/busco/busco_downloads
sudo find /media/anton/data/sandbox/Pv4/v3/work/00_inventory/busco -name "*.log" -size +500k | xargs sudo gzip --best
```

- `work/00_inventory/busco/busco_downloads/`: 1.2 GB of BUSCO HMM lineage data (uid 2000 from docker). Re-downloadable automatically on next BUSCO run.
- `work/00_inventory/busco/*_proteins/logs/busco.log`: 4 of 8 logs are root-owned, gzip blocked.
