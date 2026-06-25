# pO W/Z Combine fits — fit-stage reference

> **The full end-to-end runbook lives in the analysis repo:
> `pO_analysis/README.md`** (skim → MC norm → ABCD QCD → structured Combine
> inputs → **fit (this stage = Module 4)** → fitted yields → observables).
> This file is just the fit-stage quick reference; see Module 4 there for the
> whole procedure, the lxplus split workflow, and how the outputs feed
> `charge_asym.C` / `FBratio.C` / `observables.C`.

Branch: all pO code is on `zheng/po-analysis` (`main` is stock Combine —
`git checkout zheng/po-analysis` first). The fit needs `cmsenv`
(`combine` + `text2workspace.py` on `PATH`).

## Run

```bash
cd HiggsAnalysis-CombinedLimit/test
cmsenv
./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|all] [--dry-run] [--no-postfit]
#   PO_PLOTS=/path/to/pO_analysis/plotting/plots   (else --plots-dir, else autodetect)
```

| arg / option | meaning |
|---|---|
| `mu` / `ele` / `both` | channel(s) (default `both`) |
| `perbin` | 48 per-(charge,y) W regions (lab + FB) |
| `incl`   | `Wp_incl Wm_incl W_incl Z_incl` |
| `combined` | the simultaneous `WZ` fit only |
| `all`    | perbin + incl + combined (default) |
| `--dry-run` | build datacards only (no `cmsenv` needed) |
| `--no-postfit` | skip postfit plots |

Per region: `text2workspace` → `combine -M FitDiagnostics --saveShapes
--saveWithUncertainties` into `pO_fit_out/<chan>/{datacards,fits/<region>,postfit,
summary}`. Per-region failures are logged and skipped (not fatal).

## Pieces (all under `test/`)
- `run_pO_fits.sh` — master driver.
- `my_script/make_pO_datacards.sh` — generates all 53 datacards/channel.
- `my_script/extract_pO_yields.C` — `fit_s` → `<chan>_W_yields.csv`,
  `<chan>_summary.csv`, `<chan>_fitted_yields.root` (single-bin
  `h_mt_W{p,m}_y0..11(_FB)` with Sumw2 = fit error).
- `my_script/make_yields_from_csv.C` — rebuild the `.root` from the CSV if the
  former came out empty (no fit re-run).
- `my_script/draw_postfit_pO.C` — postfit data/MC, same cosmetics as `mtandmet.C`.
- `sync_lxplus.sh` — `upload` inputs+scripts / `download` results, one SSH auth.

## Fit model
- `signal` → POI `r` (per-region W yield). Discriminant = **PF MET shape**.
- EWK `z/ztau/wtau` → one shared `ewk_norm` rateParam (relative MC composition
  LOCKED; overall EWK scale floats).
- ABCD `qcd` → free `qcd_norm`.
- Combined `WZ`: shared `eff_lumi` multiplies W `signal` + Z `zsig`; the Z peak
  pins it (cancels in W/Z ratios). Sanity-check `eff_lumi ≈ 1` afterwards.

All templates are absolutely normalized (`k_s = A·σ·L/N_gen`) — no area norm.

## lxplus
`./sync_lxplus.sh upload` (4 inputs + scripts → lxplus), fit there, then
`./sync_lxplus.sh download [--postfit]`. Defaults:
`FORK_LX=/afs/cern.ch/user/z/zheng/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit`,
`ANA_LX=/afs/cern.ch/user/z/zheng/pO_analysis` (override via env). See
`pO_analysis/README.md` Module 4 for the full split workflow.

## Superseded (kept for reference, not used)
`run_fit.sh`, `my_script/make_combine_input{,_Z}.C`,
`my_script/testdatacard_{inclusive,Zmumu,Zee}.txt`,
`my_script/draw_postfit_{inclusive,Zmumu,Zee}.C`.
