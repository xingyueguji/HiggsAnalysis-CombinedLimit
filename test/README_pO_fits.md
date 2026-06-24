# pO W/Z Combine fits — runbook

End-to-end recipe for the rapidity-binned (+ FB) and simultaneous W/Z Combine
fits. The pipeline spans **two repos**:

| repo | role | needs |
|------|------|-------|
| `pO_analysis` (the analysis repo) | makes the skims, ABCD QCD, and the **structured Combine inputs** | ROOT (any recent; tested 6.32) |
| `HiggsAnalysis-CombinedLimit` (this fork, branch `zheng/po-analysis`) | runs the fits, extracts yields, draws postfit plots | **`cmsenv`** (combine + text2workspace.py) |

> Branch check: this fork's `main` is stock Combine. All pO code lives on
> `zheng/po-analysis` — `git checkout zheng/po-analysis` first.

The fit step is the only one that needs `cmsenv`. Everything upstream is plain
ROOT and can be done on the same machine or on lxplus; just make sure the
`combine_input_*.root` files are visible from wherever you run `combine`.

---

## TL;DR

```bash
# --- in pO_analysis (plain ROOT) ---
cd pO_analysis/skim       && ./run_all.sh all && ./run_ngen.sh           # 1. skims + N_gen
cd ../correction          && root -l -q 'qcd_abcd.C+' && root -l -q 'qcd_abcd.C+(true)'   # 2. ABCD QCD (mu, ele)
cd ../plotting            && for a in 'mtandmet.C+(false)' 'mtandmet.C+(true)' \
                                       'dileptonpeak.C+(false)' 'dileptonpeak.C+(true)'; do \
                              root -l -q -b "$a"; done                    # 3. structured Combine inputs

# --- in the Combine fork (cmsenv!) ---
cd HiggsAnalysis-CombinedLimit/test
cmsenv                                                                    # combine on PATH
./run_pO_fits.sh both all                                                 # 4. fit everything

# --- back in pO_analysis (plain ROOT) ---
cd pO_analysis/analysis
root -l -q 'charge_asym.C+("<fork>/test/pO_fit_out/mu/summary/mu_fitted_yields.root")'   # 5. observables
root -l -q 'FBratio.C+("<fork>/test/pO_fit_out/mu/summary/mu_fitted_yields.root")'
```

---

## Step 1 — Skims + N_gen  (`pO_analysis/skim/`)

```bash
cd pO_analysis/skim
./run_all.sh Wmu      # W -> mu nu   -> WToMuNu_pO_PFMet_*_hist.root
./run_all.sh Wel      # W -> e  nu   -> WToElecNu_pO_PFMet_*_hist.root
./run_all.sh Zmm      # Z -> mu mu   -> ZToMuMu_pO2025_*_hist.root
./run_all.sh Zee      # Z -> e e     -> ZToEE_pO2025_*_hist.root
#   (or: ./run_all.sh all)
./run_ngen.sh         # -> skim/rootfile/ngen.root   (MC normalization denominator)
```

`ngen.root` is required: `skim/mc_norm.h::MCScale` reads it for the absolute
scale `k_s = A·σ·L/N_gen`. Without it MCScale returns 1.0 (warns) and the inputs
are NOT absolutely normalized.

## Step 2 — ABCD QCD templates  (`pO_analysis/correction/`)

```bash
cd pO_analysis/correction
root -l -q 'qcd_abcd.C+'        # muon     -> rootfile/qcd_abcd_mu.root
root -l -q 'qcd_abcd.C+(true)'  # electron -> rootfile/qcd_abcd_ele.root
```

These are the data-driven low-MET QCD templates the W fit uses. If they are
missing, Step 3 warns and the W inputs omit `qcd` (the fit then has no QCD
background — wrong for electrons especially).

## Step 3 — Structured Combine inputs  (`pO_analysis/plotting/`)

```bash
cd pO_analysis/plotting
root -l -q -b 'mtandmet.C+(false)'      # muon  W -> plots/combine_input_W.root
root -l -q -b 'mtandmet.C+(true)'       # ele   W -> plots/Elec/combine_input_W.root
root -l -q -b 'dileptonpeak.C+(false)'  # Z mumu  -> plots/combine_input_Z.root
root -l -q -b 'dileptonpeak.C+(true)'   # Z ee    -> plots/Elec/combine_input_Z.root
```

Each `combine_input_W.root` has **one TDirectory per fit region**, each holding
the 6 **absolute** templates (`data_obs/signal/z/ztau/wtau/qcd`):

```
Wp_lab_y0 … Wp_lab_y11   Wm_lab_y0 … Wm_lab_y11     (standard lab-frame bins)
Wp_fb_y0  … Wp_fb_y11    Wm_fb_y0  … Wm_fb_y11      (FB-symmetric bins)
Wp_incl  Wm_incl  W_incl                            (charge-inclusive in y)
```
`combine_input_Z.root` has a single `Z_incl/` dir (`data_obs/signal/w/wtau/ztau`).

Quick check of a file:
```bash
root -l plots/combine_input_W.root
root [1] .ls                       # lists the 51 region dirs
root [2] W_incl->cd(); .ls         # data_obs/signal/z/ztau/wtau/qcd
```

## Step 4 — Run the fits  (this fork, **cmsenv**)

```bash
cd HiggsAnalysis-CombinedLimit/test
cmsenv
./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|all] [options]
```

| arg / option | meaning | default |
|---|---|---|
| `mu` / `ele` / `both` | channel(s) | `both` |
| `perbin` | the 48 per-(charge,y) W regions (lab + FB) | — |
| `incl`   | `Wp_incl Wm_incl W_incl Z_incl` | — |
| `combined` | only the simultaneous `WZ` fit | — |
| `all`    | perbin + incl + combined | `all` |
| `--dry-run` | build datacards + check inputs only (**no cmsenv needed**) | off |
| `--no-postfit` | skip postfit plots (faster) | off |
| `--plots-dir DIR` | analysis plots dir (else `$PO_PLOTS`, else autodetect) | autodetect |
| `--out DIR` | output root | `test/pO_fit_out` |

The driver autodetects the analysis plots dir (local Mac path, then the lxplus
AFS path). If neither applies, point it explicitly:
```bash
PO_PLOTS=/path/to/pO_analysis/plotting/plots ./run_pO_fits.sh both all
#   or
./run_pO_fits.sh both all --plots-dir /path/to/pO_analysis/plotting/plots
```

Per region it runs `text2workspace.py` → `combine -M FitDiagnostics --saveShapes
--saveWithUncertainties` and tolerates per-region failures (logs, then continues).

Examples:
```bash
./run_pO_fits.sh both all --dry-run     # sanity: just make the 53 datacards/chan
./run_pO_fits.sh mu perbin              # only the muon per-bin W fits
./run_pO_fits.sh both combined          # only the simultaneous W+Z fits
```

## Step 5 — Outputs

```
test/pO_fit_out/<chan>/                 # chan = mu | ele
  combine_input_W.root  combine_input_Z.root   # copied inputs (so cards self-resolve)
  datacards/datacard_<region>.txt
  fits/<region>/   fitDiagnostics_<region>.root, workspace.root, fit.log, t2w.log
  postfit/<region>.png/.pdf             # data/MC, same cosmetics as mtandmet.C
  summary/
    <chan>_W_yields.csv      # per-(charge,binning,y): r, rErr, signal_prefit,
                             #   fitted_yield, fitted_yield_err, qcd_norm, ewk_norm
    <chan>_summary.csv       # Wp_incl/Wm_incl/W_incl/Z_incl + WZ combined (r, eff_lumi)
    <chan>_fitted_yields.root  # h_mt_W{p,m}_y{0..11}(_FB) single-bin histos
```

`fitted_yield = r × (prefit signal integral)`, `err = rErr × (same)`.

## Step 6 — Feed the analysis observables  (`pO_analysis/analysis/`)

`<chan>_fitted_yields.root` contains exactly the histogram names
`charge_asym.C` and `FBratio.C` read, with the fit uncertainty in Sumw2 — so the
**raw→fitted swap is just the input path** (no macro edits):

```bash
cd pO_analysis/analysis
F=<fork>/test/pO_fit_out/mu/summary/mu_fitted_yields.root
root -l -q "charge_asym.C+(\"$F\",\"../skim/rootfile/charge_asym_fit_mu.root\")"
root -l -q "FBratio.C+(\"$F\",\"../skim/rootfile/FBratio_fit_mu.root\")"
#   ... and likewise for ele_fitted_yields.root
```
(For cross-sections, read the absolute fitted yields straight from
`<chan>_W_yields.csv` / `<chan>_summary.csv`.)

---

## Fit model (what the datacards encode)

Decided 2026-06-24. Per W fit region (MET discriminant):

- **`signal` → POI `r`** — the per-region W yield strength we extract.
- **EWK `z`,`ztau`,`wtau` → one shared `ewk_norm` rateParam** — relative MC
  composition LOCKED (they're absolute `k_s` templates); only the overall EWK
  normalization floats.
- **ABCD `qcd` → free `qcd_norm` rateParam** — data-driven, least trusted.

Simultaneous **W+Z** (`datacard_WZ.txt`, 2 channels): a shared **`eff_lumi`**
rateParam multiplies BOTH the W `signal` and the Z signal (renamed `zsig`, not a
POI). The high-purity Z peak pins `eff_lumi`; it cancels in W/Z ratios. The
per-bin W fits stay independent — the common scale cancels in charge-asym and F/B
ratios anyway, so the Z control matters mainly for the absolute cross-section.

All templates are **absolutely** normalized (`k_s = A·σ·L/N_gen`) — there is NO
area normalization (the old `make_combine_input*.C` `data_int/mc_sum` rescale is
gone).

## Customization

- **Retune the model** (priors, ranges, freeze QCD, add lnN systematics): edit
  `my_script/make_pO_datacards.sh` (`gen_W_card` / `gen_Z_card` /
  `gen_WZ_combined_card`). No re-skim, no re-make-input.
- **Re-bin / change the QCD split**: those live in `pO_analysis/plotting/mtandmet.C`
  (rapidity edges mirror `skim/skim_common.h::kYEdges` / `kYEdgesFB`); re-run Step 3.
- **Cosmetics** of the postfit plots come from `my_script/plotting_helper.C`
  (synced from the analysis repo's `plotting/plotting_helper.C`) via
  `SaveNicePlot1D_WithBkg`.

## Troubleshooting

- `combine / text2workspace.py not on PATH -- did you cmsenv?` → run `cmsenv` in
  your CMSSW area (or `--dry-run` to only build datacards).
- `analysis plots dir not found` → run Step 3, or pass `--plots-dir` / `$PO_PLOTS`.
- A few **tail FB bins** can have an all-zero `qcd` template (no QCD there); that
  region may warn or fail in `text2workspace` — the driver skips it and continues,
  the yield comes out 0, and `charge_asym`/`FBratio` guard against it.
- Check `eff_lumi ≈ 1` in `<chan>_summary.csv` after the combined fit — it
  validates that the Z control broke the `eff_lumi`↔`r` degeneracy cleanly.

## Superseded (kept for reference, not used by `run_pO_fits.sh`)

`run_fit.sh`, `my_script/make_combine_input{,_Z}.C` (area-normalized + Rayleigh
`pdfbkg`), `my_script/testdatacard_{inclusive,Zmumu,Zee}.txt`,
`my_script/draw_postfit_{inclusive,Zmumu,Zee}.C`.
