# pO W/Z Combine fits — fit-stage reference

> **The full end-to-end runbook lives in the analysis repo:
> `pO_analysis/README.md`** (skim → MC norm → ABCD QCD → structured Combine
> inputs → **fit (this stage = Module 4)** → fitted yields → observables).
> This file is just the fit-stage quick reference; see Module 4 there for the
> whole procedure, the lxplus split workflow, and how the outputs feed
> `fiducial_yields.C` (the r's → r × σ_gen; no observable is built on the
> fitted counts) → `charge_asym.C` / `FBratio.C` / `observables.C`.

Branch: all pO code is on `zheng/po-analysis` (`main` is stock Combine —
`git checkout zheng/po-analysis` first). The fit needs `cmsenv`
(`combine` + `text2workspace.py` on `PATH`).

## Run

```bash
cd HiggsAnalysis-CombinedLimit/test
cmsenv
./run_pO_fits.sh [mu|ele|both] [simfit|flavfit|all] [--dry-run] [--no-postfit] [--draw-only] [--asimov] [--no-statonly] [--no-contour] [--extract-only]

# flavfit (2026-09-22) = the PER-FLAVOUR simultaneous fits: the simfit model
# below restricted to ONE lepton flavour -- its 24 W channels + its own Z peak
# (+ its 6 ABCD CRs) per binning variant, 25 POIs (r_<C>_y<i> + r_Z fitted by
# that flavour alone), and otherwise EXACTLY the grand fit's treatment: the
# nuisances that act on that flavour (lumi, its QCD rows, the LHE shapes, muSF
# in the muon fit only), --statonly companion (the stat error), --contour scan,
# --asimov closure, the same extraction. `both flavfit` = one mu-only AND one
# e-only fit; outputs pO_fit_out<suffix>/simfit_{mu,ele}/ with summary files
# simfit_<flav>_{W_yields.csv,summary.csv,fitted_yields.root}. `all` = the
# grand simfit + both flavfits. The analysis repo's run_observables.sh overlays
# mu vs e (plots/flavfit/). The LEGACY per-flavour per-bin modes
# perbin/incl/combined were REMOVED on 2026-09-22 (flavfit replaces them).
#
# DEFAULT (2026-08-04) = simfit, the GRAND SIMULTANEOUS FIT: one likelihood per
# binning variant (lab, fb) with all 48 W channels ({mu,ele} x {Wp,Wm} x y0..11)
# + BOTH Z peaks. 25 POIs: r_<C>_y<i> (24, mu/e SHARED) + one global r_Z on all
# DY-related MC; w/wtau under Z frozen at MC. QCD via QCD_MODE:
#   lnN (default, 2026-08-17): qcd_rate_{mu,ele}_{Wp,Wm} lnN at the ABCD
#     prediction (kappa mu 1.15 / ele 1.20; QCD_LNN_MU/QCD_LNN_ELE);
#   free: the pre-2026-08-17 48 per-channel qcd_norm rateParams;
#   abcd (2026-08-23, --disc leppt_mt40 ONLY): the IN-FIT ABCD -- 12 counting
#     CR channels + free scales qcd_s{B,C,D}_<F>_<C> + the formula rateParam
#     (sB*sC/sD) on the SR qcd_abcd template; EWK subtraction rides the POIs
#     (CRB DY -> r_Z, CRB W -> per-y w_y* mapped to the r's; QCD_WCR=frozen
#     freezes it); reduced residual kappas QCD_ABCD_LNN_MU=1.09/ELE=1.15.
# --asimov adds a prefit-Asimov closure fit (all POIs must return 1 -- and in
# abcd mode the 12 CR scales too, and every shape nuisance at 0; the extraction
# prints PASS/FAIL). Outputs: pO_fit_out<suffix>/simfit/summary/
# {comb_W_yields.csv, comb_summary.csv, comb_fitted_yields.root(+h_cov_yield[_FB])}.
# LHE shape systematics (2026-09-07): the inputs carry <proc>_{nPDF,qcdScale,
# alphaS}Up/Down (+ a <input>_systs.txt sidecar listing them); the cards get one
# `shape` row each on the MC columns + the `lhe` group. LHE_SYST=auto (default,
# the sidecars' common list) | off (cards as before) | nPDF,alphaS (subset).
# Pulls/constraints -> comb_summary.csv <name>_theta rows; stat-only comparison:
# combine ... --freezeNuisanceGroups lhe
# Lepton-SF shape systematics (2026-09-14, muon first): the muon inputs list
# ONE combined nuisance muSF (skim/muon_sf.h: ID, ISO and inclusive trigger SF
# shifts added in quadrature per bin) -> a row with entries on the muon columns
# only (the generator uses the UNION of the four sidecars with per-input
# process lists), group `lepsf`. If the analysis repo ships the three sources
# separately (muID/muIso/muTrig) instead, SF_TRIG_CORR=auto (default) follows
# the sidecar's '#! muTrig corr' directive (coherent for an inclusive trigger
# SF; perbin -> muTrig split into muTrig_y0..11 with `nuisance edit rename`).
# Stat-only: --freezeNuisanceGroups lhe,lepsf
# Stat/syst split of every POI error (2026-09-15): the extractor conditions the
# nominal fit's covariance matrix on the constrained nuisances (Schur
# complement = the frozen-nuisance result, no refit) -> rErr_stat (19th CSV
# col), simfit_<B>_stat summary rows, h_cov_yield[_FB]_stat; --statonly adds
# the frozen-nuisance refit as a cross-check (simfit_<B>_statonlyfit rows);
# --extract-only re-runs just the extraction on an existing fits/ tree.
./run_pO_fits.sh --asimov
./run_pO_fits.sh both flavfit --asimov                                  # mu-only + e-only fits
./run_pO_fits.sh both all --asimov                                      # grand + mu-only + e-only
QCD_MODE=abcd ./run_pO_fits.sh both simfit --disc leppt_mt40 --asimov   # in-fit ABCD
LHE_SYST=off ./run_pO_fits.sh both simfit --disc leppt_mt40             # no theory shape nuisances
SF_TRIG_CORR=perbin QCD_MODE=abcd ./run_pO_fits.sh both simfit --disc leppt_mt40   # robustness: muTrig decorrelated per y bin

# W discriminant variants (2026-07-30): --disc met|leppt|leppt_mt40 (default leppt_mt40 since 2026-08-16).
# leppt / leppt_mt40 read combine_input_W_leppt[_mt40].root and write to
# pO_fit_out_leppt[_mt40]/; the Z channel is unchanged. For the pT variants the
# discriminant lacks the low-MET in-fit QCD anchor, which is exactly what
# QCD_MODE=abcd restores for leppt_mt40 (explicit CR channels).
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
#     inside <fit>_fitted_yields.root are named identically across variants
#     (h_yield_*) -- the tree name is the only label.)
# Full workflow + physics notes: the analysis repo's README.md, Module 4.
#   PO_PLOTS=/path/to/pO_analysis/plotting/plots   (else --plots-dir, else autodetect)
```

| arg / option | meaning |
|---|---|
| `mu` / `ele` / `both` | the flavour(s) of `flavfit` (default `both` = one μ-only AND one e-only fit); `simfit` is always μ+e and ignores it |
| `simfit` | **DEFAULT.** The grand simultaneous fit (lab + fb workspaces, 25 POIs, μ/e shared) → `pO_fit_out<suffix>/simfit/`, summary files `comb_*` |
| `flavfit` | the per-flavour simultaneous fit(s): the same model and passes, one flavour per likelihood → `pO_fit_out<suffix>/simfit_{mu,ele}/`, summary files `simfit_<flav>_*` |
| `all`    | `simfit` (when channel = `both`) + `flavfit` for the channel's flavour(s) |
| `perbin` / `incl` / `combined` | **removed 2026-09-22** (the legacy per-flavour per-bin pipeline) — exits with a pointer to `flavfit` |
| `--dry-run` | build datacards only (no `cmsenv` needed) |
| `--no-postfit` | skip postfit plots |
| `--asimov` | (simfit/flavfit) also fit the prefit Asimov dataset per variant — closure: every POI = 1 |
| `--no-statonly` | (simfit/flavfit) SKIP the frozen-nuisance companion fit. It runs by default (2026-09-15b) and IS the source of the quoted stat error (all constrained nuisances frozen at their post-fit values → `fitDiagnostics_simfit_<B>_statonly.root`); without it the extractor falls back to the nominal fit's covariance conditioned on the constrained nuisances (Schur complement = the Gaussian-exact frozen-nuisance result), which it computes either way as the cross-check → `rErr_stat` (19th CSV col), `simfit_<B>_stat` summary rows, `h_cov_yield[_FB]_stat`, `h_cov_poi[_FB]_stat` |
| `--no-contour` | (simfit/flavfit) SKIP the profiled (σ_W, σ_Z) scan (on by default; additive, own workspace + `contour/contour_<B>/`) |
| `--extract-only` | (simfit/flavfit) re-run only the extraction on an EXISTING `fits/` tree (only `root`, no `cmsenv`) — e.g. on a downloaded lxplus fit after an extractor change. Prefit integrals from the work-dir input copies, else the analysis plots dir; the extractor checks them against the fit's own `shapes_prefit` and WARNs if they are not the inputs that were fitted. `sync_lxplus.sh download` pulls the nominal, `_statonly` and `_asimov` fitDiagnostics so the local re-extraction is complete |
| `--draw-only` | redraw postfit plots from EXISTING fits (no `combine`/`cmsenv`, only `root`) — e.g. after cosmetic changes to `draw_postfit_pO.C`. Respects channel+mode; needs the `fits/` tree from a previous run (not pulled by `sync_lxplus.sh download` — redraw where the fits ran, then `download --postfit`) |

Per fit and binning variant: `text2workspace` (multiSignalModel maps) →
`combine -M FitDiagnostics --saveShapes --saveWithUncertainties
--skipBOnlyFit` into `pO_fit_out<suffix>/<fit>/{datacards,fits/simfit_<B>,
contour/contour_<B>,postfit,summary}` with `<fit>` = `simfit` (grand),
`simfit_mu`, `simfit_ele`. A failed pass is logged and skipped (not fatal).

## Pieces (all under `test/`)
- `run_pO_fits.sh` — master driver (modes `simfit`, `flavfit`, `all`).
- `my_script/make_pO_simfit_cards.sh` — writes the two grand-fit datacards
  (`datacard_simfit_{lab,fb}.txt`, 50 channels / 62 in abcd mode) + the
  multiSignalModel map files (`t2w_maps_simfit_{lab,fb}.txt` — THE definition of
  the 25-POI model; `_sigma` twins for the contour). `SIMFIT_FLAVS=mu|ele`
  (set by `flavfit`) writes the per-flavour card instead: that flavour's
  channels, QCD rows and sidecar-listed shape rows only (the grand card is
  byte-identical with the default `"mu ele"`). No `combineCards.py`;
  `--dry-run` friendly.
- `my_script/extract_pO_simfit.C` — one `fit_s` per variant → all POIs + the
  r-correlation matrix → `<tag>_W_yields.csv`, `<tag>_summary.csv` (incl.
  covariance-propagated Wp/Wm/W sums + Asimov closure), and
  `<tag>_fitted_yields.root` (`h_yield_W{p,m}_y0..11(_FB)`, yields = r×S summed
  over the fit's flavours, plus the 24×24 `h_cov_yield[_FB]` and the 25×25
  `h_cov_poi[_FB]`, each with a `_stat` twin). `<tag>` = `comb` (grand) or
  `simfit_<flav>`; the trailing `flavours` / `outTag` arguments select it.
- `my_script/draw_postfit_pO.C` — postfit data/MC, same cosmetics as
  `mtandmet.C`; optional trailing args (poi/dy/qcd parameter names + ndf) let
  the simfit pass its per-channel parameter names.
- `run_pO_impacts.sh` — impacts + correlation plots; `--fit comb|mu|ele|all`
  (`all` = the grand fit + both flavfits in one call; absent fits are skipped).
- `sync_lxplus.sh` — `upload` inputs+scripts / `download` results (the
  `simfit/` tree and the `simfit_{mu,ele}/` flavfit trees), one SSH auth.
- (Removed 2026-09-22 with the legacy per-bin pipeline: `make_pO_datacards.sh`,
  `extract_pO_yields.C`, `make_yields_from_csv.C` — in git history.)

## Fit model — simfit (2026-08-04, DEFAULT)
- One likelihood per binning variant (lab / fb — the same events rebinned, so
  never combined with each other): 48 W channels + `mu_Z_incl` + `ele_Z_incl`.
- **`r_<C>_y<i>`** (24 POIs): scales `signal`+`wtau` of that (charge, y) bin in
  BOTH flavours' channels — μ/e shared (lepton universality; the relative μ/e
  acceptance×efficiency comes from MC — muon SFs applied in the analysis
  repo's skim since 2026-09-14, electron SFs not yet).
- **`r_Z`**: one global scale on ALL DY-related MC (`z`/`ztau` in every W
  channel + `zsig`/`ztau` under both Z peaks). DY rapidity dependence across W
  bins is fixed from MC; the Z peaks pin the normalization.
- QCD (`QCD_MODE`): **lnN** (default) — `qcd_rate_{mu,ele}_{Wp,Wm}` lnN at the
  ABCD prediction; **free** — `qcd_norm_<channel>` per W channel (48);
  **abcd** (leppt_mt40 only) — in-fit ABCD: CR channels `<F>_<C>_CR{B,C,D}`,
  free scales `qcd_s{B,C,D}_<F>_<C>`, formula rateParam `(@0*@1/@2)` on the SR
  `qcd_abcd` template, reduced residual lnN. `w`/`wtau` under the Z peaks:
  frozen at absolute MC (0.03–0.06 events — negligible by measurement).
- LHE shape systematics (2026-09-07, `LHE_SYST`): `nPDF` (EPPS21 via LHAPDF's
  `PDFSet.uncertainty()`), `qcdScale` (μR/μF envelope), `alphaS` (0.119/0.117
  members) — one `shape` row each, flag `1` on `signal/z/ztau/wtau` (W) and
  `zsig/w/wtau/ztau` (Z), `-` on the data-driven `qcd` and all CR columns;
  templates `<dir>/<proc>_$SYSTEMATIC` (Up/Down) from the analysis repo's
  `combine_input_*.root`, listed in the `<input>_systs.txt` sidecar the
  generator reads. One name = one θ for the whole card (fully correlated). Group
  `lhe` for `--freezeNuisanceGroups`. Sidecar line `lheSysts` → the extractor
  reports pulls + constraints and includes them in the Asimov closure.
- Lepton-SF shape systematics (2026-09-14, muon first): ONE combined nuisance
  `muSF` from the analysis repo's `skim/muon_sf.h` (pp POG ID/ISO SFs and the
  MB-derived inclusive trigger SF, their ±1σ shifts added in quadrature per
  bin), listed by the muon sidecars only → `1` on the muon W/Z MC columns,
  `-` on every electron column (the generator uses the union of the four
  sidecars with per-input process lists); group `lepsf`; same `shape`
  treatment as the theory nuisances. Should the analysis repo ship the three
  sources separately (`muID`, `muIso`, `muTrig`), the muTrig correlation model
  follows the sidecar directive `#! muTrig corr` (`SF_TRIG_CORR=auto`):
  `coherent` for an inclusive trigger SF, `perbin` → split into 12 nuisances
  `muTrig_y<i>` with `nuisance edit rename` (same histograms). Sidecar
  `lheSysts` lists all nuisance names, `sfTrigCorr` the resolved choice.
- Implemented via `multiSignalModel` maps (see `t2w_maps_simfit_*.txt`);
  `FitDiagnostics --skipBOnlyFit` (a b-only fit with all POIs at 0 is
  meaningless). Statistical gain vs legacy: the Z data enters ONCE (the legacy
  per-bin cards each re-fit the same Z data, ignoring the induced
  correlations), and the full POI covariance feeds the downstream A/R_FB
  errors.

## Fit model — flavfit (2026-09-22, the per-flavour fits)
- The simfit model above with ONE lepton flavour per likelihood: that
  flavour's 24 W channels + its Z peak (+ its 6 CR channels in abcd mode).
- Still 25 POIs: `r_<C>_y<i>` now measured by that flavour alone, and its own
  `r_Z` from its own Z peak (so the μ-only and e-only r's are independent
  measurements; the grand fit forces one shared r).
- Every nuisance that acts on the flavour, exactly as in the grand fit:
  `lumi`, its two `qcd_rate_<F>_{Wp,Wm}` rows, the LHE shapes, `muSF` (μ only).
  The other flavour's rows are not written (they would touch nothing).
- Same passes (nominal, `--statonly`, `--contour`, `--asimov`) and extraction.
- The legacy per-flavour per-bin model (two-parameter `r` + `dy_norm` + free
  `qcd_norm`, 48 separate two-channel cards per flavour, no systematics) was
  removed on 2026-09-22.

All templates are absolutely normalized (`k_s = A·σ·L/N_gen`) — no area norm.

## Diagnosing fit quality

Each postfit plot (`draw_postfit_pO.C`) now prints, in the top-right info box:
- `#chi^{2}/ndf = … (p=…)` — Poisson (Baker-Cousins) goodness-of-fit of data vs
  the postfit total, robust at low counts (ndf = bins used − floating params).
  χ²/ndf ≈ 1 and a non-tiny p mean good agreement.
- `r = … ± …` — the fitted signal strength of that channel's bin (`r_<C>_y<i>`).
- `DY norm` / `QCD` — the fitted `r_Z` and the channel's QCD parameter (in
  lnN/abcd mode the shared `qcd_rate_<F>_<C>` pull).
- a **red** `status N, covQ M` only when the fit did **not** converge cleanly
  (want `status 0`, `covQ 3` = full accurate covariance) — absence of red = OK.

Deeper checks (beyond the plot):
- `pO_fit_out<suffix>/<fit>/fits/simfit_<B>/fit.log` (+ `fit_statonly.log`,
  `fit_asimov.log`) — the combine logs; `summary/extract_<fit>.log` — the
  extraction record (stat source, cross-check, closure, WARNs).
- `fitDiagnostics_simfit_<B>.root` — `fit_s` (`RooFitResult`: `status()`,
  `covQual()`, `edm()`); nuisance pulls via
  `python diffNuisances.py fitDiagnostics_simfit_<B>.root`.
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
