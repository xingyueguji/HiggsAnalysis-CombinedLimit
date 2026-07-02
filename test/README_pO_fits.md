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
./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|all] [--dry-run] [--no-postfit] [--draw-only]
#   PO_PLOTS=/path/to/pO_analysis/plotting/plots   (else --plots-dir, else autodetect)
```

| arg / option | meaning |
|---|---|
| `mu` / `ele` / `both` | channel(s) (default `both`) |
| `perbin` | 48 per-(charge,y) W regions (lab + FB), each fitted **simultaneously with `Z_incl`** (two-channel card) |
| `incl`   | `Wp_incl Wm_incl W_incl Z_incl` (standalone) |
| `combined` | the simultaneous `WZ` (`W_incl`+`Z_incl`) fit only |
| `all`    | perbin + incl + combined (default) |
| `--dry-run` | build datacards only (no `cmsenv` needed) |
| `--no-postfit` | skip postfit plots |
| `--draw-only` | redraw postfit plots from EXISTING fits (no `combine`/`cmsenv`, only `root`) — e.g. after cosmetic changes to `draw_postfit_pO.C`. Respects channel+mode; needs the `fits/` tree from a previous run (not pulled by `sync_lxplus.sh download` — redraw where the fits ran, then `download --postfit`) |

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

## Fit model (two-parameter, 2026-07-01)
- **Two MC scales per fit**: the POI **`r` = all W-related MC** (W `signal` +
  `wtau`, plus the `w`/`wtau` backgrounds under the Z peak in simultaneous
  cards) and **`dy_norm` = all DY-related MC** (`z` + `ztau` + the Z signal in
  simultaneous cards). Relative composition WITHIN each group stays LOCKED by
  the absolute `k_s` templates. Implemented by reusing `r` as a `rateParam` on
  the W backgrounds (signal is index 0 → `r` scales it automatically). Fitted W
  yield = `r`×signal-prefit. Discriminant = **PF MET shape**.
- Standalone `Z_incl` card: roles flip — the POI `r` IS the DY scale
  (`signal`+`ztau`), and the W backgrounds get a free `w_norm`.
- ABCD `qcd` → its own free `qcd_norm`.
- **Simultaneous cards** (`WZ` AND every per-bin W card): two fit channels
  (the W region + `Z_incl`). `r` scales the W-related in both channels;
  the shared `dy_norm` scales the DY-related in both — the high-purity Z peak
  pins it (this **replaces the old shared `eff_lumi`**). Sanity-check
  `dy_norm ≈ 1` afterwards. If `combine_input_Z.root` is missing, per-bin
  cards fall back to standalone W-only and `Z_incl`/`WZ` are skipped.

All templates are absolutely normalized (`k_s = A·σ·L/N_gen`) — no area norm.
Systematics deliberately deferred (statistics-only fits for now).

## Diagnosing fit quality

Each postfit plot (`draw_postfit_pO.C`) now prints, in the top-right info box:
- `#chi^{2}/ndf = … (p=…)` — Poisson (Baker-Cousins) goodness-of-fit of data vs
  the postfit total, robust at low counts (ndf = bins used − floating params).
  χ²/ndf ≈ 1 and a non-tiny p mean good agreement.
- `r = … ± …` — the fitted signal strength for that region.
- `DY norm` / `W norm` / `QCD norm` — the fitted `dy_norm` / `w_norm` /
  `qcd_norm`, whichever float in that fit (two-parameter model diagnostics;
  all should sit near 1 except `qcd_norm`, which is genuinely free).
- a **red** `status N, covQ M` only when the fit did **not** converge cleanly
  (want `status 0`, `covQ 3` = full accurate covariance) — absence of red = OK.

Deeper checks (beyond the plot):
- `pO_fit_out/<chan>/fits/<region>/fit.log` — the combine log (warnings, MINUIT).
- `fitDiagnostics_<region>.root` — `fit_s` (`RooFitResult`: `status()`,
  `covQual()`, `edm()`); nuisance pulls via
  `python diffNuisances.py fitDiagnostics_<region>.root`.
- rigorous GoF p-value with toys:
  `combine -M GoodnessOfFit <ws> --algo saturated -t 200` (+ the same on data).

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
