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
./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|simfit|all] [--dry-run] [--no-postfit] [--draw-only] [--asimov]

# DEFAULT (2026-08-04) = simfit, the GRAND SIMULTANEOUS FIT: one likelihood per
# binning variant (lab, fb) with all 48 W channels ({mu,ele} x {Wp,Wm} x y0..11)
# + BOTH Z peaks. 25 POIs: r_<C>_y<i> (24, mu/e SHARED) + one global r_Z on all
# DY-related MC; qcd_norm free per W channel; w/wtau under Z frozen at MC.
# --asimov adds a prefit-Asimov closure fit (all POIs must return 1; the
# extraction prints PASS/FAIL). Outputs: pO_fit_out<suffix>/simfit/summary/
# {comb_W_yields.csv, comb_summary.csv, comb_fitted_yields.root(+h_cov_yield[_FB])}.
./run_pO_fits.sh --asimov

# W discriminant variants (2026-07-30): --disc met|leppt|leppt_mt40 (default met).
# leppt / leppt_mt40 read combine_input_W_leppt[_mt40].root and write to
# pO_fit_out_leppt[_mt40]/; the Z channel and the fit model are unchanged and
# qcd_norm stays free (the pT variants lack the low-MET QCD anchor, so expect a
# weaker qcd_norm constraint / larger r-qcd correlation).
./run_pO_fits.sh both all --disc leppt_mt40

# *** WARNING: carry the SAME --disc through the WHOLE workflow of a variant ***
# The out-trees have no discriminant marker (the coupling is only the dir
# suffix), so mixed steps corrupt/mislabel silently:
#   - --draw-only MUST repeat the same --disc (it picks the out-tree AND the
#     x-title; forget it and lepton-pT plots get relabeled "PF MET (GeV)").
#     Same rule whenever --out is used.
#   - sync_lxplus.sh download needs NO flag (sweeps all three out-trees).
#   - observables: run the analysis repo's analysis/run_observables.sh <disc>
#     (2026-08-03) -- it feeds the MATCHING tree automatically and writes
#     disc-tagged files + per-disc plot folders. (By hand, remember the histos
#     inside <chan>_fitted_yields.root are named identically across variants
#     (h_yield_*) -- the tree name is the only label.)
# Full workflow + physics notes: the analysis repo's README.md, Module 4.
#   PO_PLOTS=/path/to/pO_analysis/plotting/plots   (else --plots-dir, else autodetect)
```

| arg / option | meaning |
|---|---|
| `mu` / `ele` / `both` | channel(s) for the LEGACY modes (default `both`; `simfit` is always μ+e and ignores it) |
| `simfit` | **DEFAULT.** The grand simultaneous fit (lab + fb workspaces, 25 POIs, μ/e shared) → `pO_fit_out<suffix>/simfit/` |
| `perbin` | legacy: 48 per-(charge,y) W regions (lab + FB), each fitted **simultaneously with `Z_incl`** (two-channel card) |
| `incl`   | legacy: `Wp_incl Wm_incl W_incl Z_incl` (standalone) |
| `combined` | legacy: the simultaneous `WZ` (`W_incl`+`Z_incl`) fit only |
| `all`    | the whole legacy per-flavour pipeline + simfit (when channel = `both`) |
| `--dry-run` | build datacards only (no `cmsenv` needed) |
| `--no-postfit` | skip postfit plots |
| `--asimov` | (simfit) also fit the prefit Asimov dataset per variant — closure: every POI = 1 |
| `--draw-only` | redraw postfit plots from EXISTING fits (no `combine`/`cmsenv`, only `root`) — e.g. after cosmetic changes to `draw_postfit_pO.C`. Respects channel+mode; needs the `fits/` tree from a previous run (not pulled by `sync_lxplus.sh download` — redraw where the fits ran, then `download --postfit`) |

Per region: `text2workspace` → `combine -M FitDiagnostics --saveShapes
--saveWithUncertainties` into `pO_fit_out/<chan>/{datacards,fits/<region>,postfit,
summary}`. Per-region failures are logged and skipped (not fatal).

## Pieces (all under `test/`)
- `run_pO_fits.sh` — master driver (simfit + legacy modes).
- `my_script/make_pO_simfit_cards.sh` — **simfit**: writes the two 50-channel
  grand-fit datacards (`datacard_simfit_{lab,fb}.txt`) + the multiSignalModel
  map files (`t2w_maps_simfit_{lab,fb}.txt` — THE definition of the 25-POI
  model). No `combineCards.py`; `--dry-run` friendly.
- `my_script/extract_pO_simfit.C` — **simfit**: one `fit_s` per variant → all
  POIs + the r-correlation matrix → `comb_W_yields.csv`, `comb_summary.csv`
  (incl. covariance-propagated Wp/Wm/W sums + Asimov closure), and
  `comb_fitted_yields.root` (`h_yield_W{p,m}_y0..11(_FB)`, yields = r×(S_mu+S_ele),
  plus the 24×24 covariance TH2Ds `h_cov_yield[_FB]`, order [Wp_y0..11, Wm_y0..11]).
- `my_script/make_pO_datacards.sh` — legacy: generates all 53 datacards/channel.
- `my_script/extract_pO_yields.C` — legacy: `fit_s` → `<chan>_W_yields.csv`,
  `<chan>_summary.csv`, `<chan>_fitted_yields.root` (single-bin
  `h_mt_W{p,m}_y0..11(_FB)` with Sumw2 = fit error).
- `my_script/make_yields_from_csv.C` — rebuild the `.root` from the CSV if the
  former came out empty (no fit re-run).
- `my_script/draw_postfit_pO.C` — postfit data/MC, same cosmetics as
  `mtandmet.C`; optional trailing args (poi/dy/qcd parameter names + ndf) let
  the simfit pass its per-channel parameter names — defaults = legacy behavior.
- `sync_lxplus.sh` — `upload` inputs+scripts / `download` results (legacy
  per-chan + `simfit/` subtree), one SSH auth.

## Fit model — simfit (2026-08-04, DEFAULT)
- One likelihood per binning variant (lab / fb — the same events rebinned, so
  never combined with each other): 48 W channels + `mu_Z_incl` + `ele_Z_incl`.
- **`r_<C>_y<i>`** (24 POIs): scales `signal`+`wtau` of that (charge, y) bin in
  BOTH flavours' channels — μ/e shared (lepton universality; the relative μ/e
  acceptance×efficiency comes from MC — lepton SFs not applied yet).
- **`r_Z`**: one global scale on ALL DY-related MC (`z`/`ztau` in every W
  channel + `zsig`/`ztau` under both Z peaks). DY rapidity dependence across W
  bins is fixed from MC; the Z peaks pin the normalization.
- `qcd_norm_<channel>`: free per W channel (48). `w`/`wtau` under the Z peaks:
  frozen at absolute MC (0.03–0.06 events — negligible by measurement).
- Implemented via `multiSignalModel` maps (see `t2w_maps_simfit_*.txt`);
  `FitDiagnostics --skipBOnlyFit` (a b-only fit with all POIs at 0 is
  meaningless). Statistical gain vs legacy: the Z data enters ONCE (the legacy
  per-bin cards each re-fit the same Z data, ignoring the induced
  correlations), and the full POI covariance feeds the downstream A/R_FB
  errors.

## Legacy fit model (two-parameter, 2026-07-01)
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
