# Run 2 missing-tree reproduction

The 2022 trees under `flowcorr` were made with the original `TrackAnalyzer`
implementation imported in Git commit `93af001`; its matching header first
appears in commit `ab9e5f1`.  The source snapshot in `src_run2`
contains only the tested plugin files from those commits.  That implementation stores
reconstructed AK8 jets with `pT > 350 GeV` and `|eta| < 2.4`, stores generated
AK8 jets using the configured equivalent cuts, records HLT decisions without
filtering events, and writes the `analyzer/trackTree` schema consumed by the
Run 2 correction macros.

Use `run2_mc_cfg.py` with the historical analyzer in a CMSSW_10_6_30 area.
The era-specific global tags are:

- `2016APV`: `106X_mcRun2_asymptotic_preVFP_v11`
- `2016`: `106X_mcRun2_asymptotic_v17`
- `2017`: `106X_mc2017_realistic_v9`
- `2018`: `106X_upgrade2018_realistic_v16_L1v1`

An isolated CMSSW_10_6_30 test on ten events from the missing 2018
QCD 800--1000 sample produced `analyzer/trackTree` with exactly the same 72
branch names and branch types as a representative 2022 tree.

Following analysis-note Table 4, the required missing-tree reproduction stops
at QCD 2400--3200 and excludes QCD 3200--infinity.  The required missing trees
therefore cover 4,032 official MiniAOD files and 214,867,000 events.  Based on
a representative 2022 output, they would require approximately 1.9 TB of EOS
storage, with sample- and event-content-dependent variation.  Do not submit this production until its
scope, output location, and storage budget have been approved.

## CRAB production

`crab_run2_missing.py` selects one sample at a time from
`missing_samples.json`.  It uses `FileBased` splitting with one MiniAOD file per
job, does not apply a constituent-pT cut, does not apply pileup reweighting, and
writes to:

```
/store/group/phys_heavyions/xiaoyul/Run2_root_out/revisit/trees
```

The five explicit submission commands are:

```bash
RUN2_SAMPLE=2016APV_qcd_800to1000 crab submit -c crab_run2_missing.py
RUN2_SAMPLE=2016_qcd_800to1000 crab submit -c crab_run2_missing.py
RUN2_SAMPLE=2017_qcd_800to1000 crab submit -c crab_run2_missing.py
RUN2_SAMPLE=2018_qcd_800to1000 crab submit -c crab_run2_missing.py
RUN2_SAMPLE=2018_qcd_600to800 crab submit -c crab_run2_missing.py
```

Run those commands only from a built CMSSW_10_6_30 area after checking the
expanded configuration and EOS destination.  Submission is intentionally not
automated as a loop so that the requested sample list is visible in the shell
history and each task can be approved independently.

