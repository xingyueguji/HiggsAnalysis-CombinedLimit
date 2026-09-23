#!/usr/bin/env bash
# =============================================================================
# run_pO_fits.sh -- end-to-end pO W/Z Combine fits, automated.
#
# For each fit it:
#   1. locates the structured Combine inputs from the analysis repo
#        <plots>/combine_input_W<disc>.root , <plots>/combine_input_Z.root      (muon)
#        <plots>/Elec/combine_input_{W<disc>,Z}.root                            (electron)
#      and copies them (+ their *_systs.txt sidecars) into the fit's work dir
#   2. writes the datacard + multiSignalModel maps per binning variant (lab, fb)
#      [my_script/make_pO_simfit_cards.sh]
#   3. runs text2workspace + combine -M FitDiagnostics (+ the --statonly
#      companion, the --contour scan, and the --asimov closure when asked)
#   4. extracts POIs, fitted yields, the stat/syst split and the covariances
#      [my_script/extract_pO_simfit.C]
#   5. draws postfit data/MC plots in the SAME cosmetics as plotting/mtandmet.C
#      [my_script/draw_postfit_pO.C]
#
# Output tree, one work dir per fit (the grand fit and each per-flavour fit):
#   pO_fit_out<suffix>/<fit>/{combine_input_*.root, datacards/, fits/simfit_<B>/,
#                             contour/contour_<B>/, postfit/, summary/}
#   <fit> = simfit (grand, files comb_*) | simfit_mu | simfit_ele (files simfit_<flav>_*)
#
# The LEGACY per-flavour per-bin pipeline (modes perbin/incl/combined:
# 48 separate two-channel cards per flavour, the old r + dy_norm + free
# qcd_norm model with no systematics) was REMOVED on 2026-09-22 (user
# decision) -- flavfit is the per-flavour fit now.
#
# Usage:
#   ./run_pO_fits.sh [mu|ele|both] [simfit|flavfit|all] [options]
#     channel  (default both)   mode (default simfit)
#     all = the grand simfit + BOTH per-flavour fits, in one go
#
#   Mode 'flavfit' (2026-09-22) = the PER-FLAVOUR SIMULTANEOUS FITS: the simfit
#   model below restricted to ONE lepton flavour -- that flavour's 24 W channels
#   + its own Z peak (+ its 6 ABCD control regions) in one likelihood per
#   binning variant, 25 POIs (r_<C>_y<i> + r_Z fitted by that flavour ALONE),
#   and EXACTLY the grand fit's treatment otherwise: the same nuisances where
#   they act on that flavour (lumi, its QCD lnN rows, the LHE shapes, muSF in
#   the muon fit only), the same three passes (nominal + --statonly companion =
#   the stat error + --contour scan), --asimov closure, extraction with the
#   stat/syst split and the 25x25 POI covariance, postfit plots. One fit per
#   flavour of the channel argument (both -> mu AND ele, run one after the
#   other); outputs under <out>/simfit_mu/ and <out>/simfit_ele/, summary files
#   simfit_<flav>_{W_yields.csv,summary.csv,fitted_yields.root}. The point:
#   compare mu with e (analysis/run_observables.sh overlays them) while the
#   electron SFs are not applied -- the grand fit forces one r on both.
#
#   Mode 'simfit' (2026-08-04, the DEFAULT) = the GRAND SIMULTANEOUS FIT: all 48 (flavour,
#   charge, y-bin) W channels + BOTH Z peaks in ONE likelihood per binning
#   variant (lab, fb).  25 POIs: r_<C>_y<i> (24, mu/e SHARED) + one global r_Z
#   scaling all DY-related MC; global lumi lnN on all MC.  QCD: the DEFAULT is
#   DISC-DEPENDENT since 2026-09-15b -- QCD_MODE=abcd (the IN-FIT ABCD) for the
#   primary --disc leppt_mt40, QCD_MODE=lnN (log-normal-constrained at the ABCD
#   prediction per flavour+charge) for met/leppt, where the qcd_abcd template
#   and the CR dirs do not exist.  Env QCD_MODE/QCD_LNN_MU/QCD_LNN_ELE/LUMI_LNN
#   -> make_pO_simfit_cards.sh and an explicit QCD_MODE always wins
#   (QCD_MODE=lnN on leppt_mt40 = the like-for-like comparison; QCD_MODE=free
#   restores the 48 free qcd_norm rateParams).  QCD_MODE=abcd
#   (2026-08-23, --disc leppt_mt40 ONLY) = the IN-FIT ABCD: 12 counting CR
#   channels + free scales qcd_s{B,C,D}_<F>_<C> + the formula rateParam
#   (sB*sC/sD) on the SR qcd_abcd template, so the QCD normalization floats
#   with the CR data and the EWK subtraction rides the POIs; env
#   QCD_ABCD_LNN_MU/QCD_ABCD_LNN_ELE = the reduced residual kappas,
#   QCD_WCR=float|frozen picks the CRB W treatment); w/wtau under the Z
#   peaks frozen at absolute MC.  LHE shape systematics (2026-09-07):
#   nPDF / qcdScale / alphaS `shape` rows from the <proc>_<syst>Up/Down
#   templates the inputs carry (env LHE_SYST=auto|off|list -> the card
#   generator; the inputs' *_systs.txt sidecars travel with the copies).
#   Cross-flavour by construction, so the channel
#   argument is ignored (mode 'all' runs simfit only when channel = both; a
#   single-flavour fit is mode 'flavfit').
#   Outputs under <out>/simfit/ (comb_* files; yields are mu+e combined).
#   options:
#     --disc met|leppt|leppt_mt40
#                       W discriminant (default leppt_mt40 = the PRIMARY one;
#                       the help said "met" until 2026-09-15b, stale since the
#                       2026-08-16 switch -- the code default is DISC above).
#                       met        = PF MET shape (the backup variant)
#                       leppt      = lepton pT, plain W selection
#                       leppt_mt40 = lepton pT with the pT>25 && m_T>40 selection
#                       Reads combine_input_W[_leppt[_mt40]].root and writes to
#                       pO_fit_out[_leppt[_mt40]]/ (unless --out). Z channel and
#                       fit model identical.  NB the pT variants lack the
#                       low-MET QCD anchor, which is exactly why the QCD lnN
#                       (vs the old free qcd_norm) matters most for them.
#     --dry-run         build datacards + check inputs only (no cmsenv needed)
#     --no-postfit      skip the postfit plots (faster)
#     --draw-only       redraw postfit plots from EXISTING fits (no combine run;
#                       use after cosmetic changes to draw_postfit_pO.C)
#     --asimov          (simfit/flavfit) also run a prefit-Asimov closure fit per
#                       variant (-t -1): every fitted POI must come back at 1
#     --no-statonly     (simfit/flavfit) SKIP the frozen-nuisance companion fit.
#                       It is ON by default (2026-09-15b) and is THE source of
#                       the quoted statistical error: all constrained nuisances
#                       frozen at their post-fit values, so the POI errors that
#                       come back are the stat component and downstream
#                       syst = sqrt(total^2 - stat^2). Skipping it makes the
#                       extractor fall back to the conditioned covariance of the
#                       nominal fit (Gaussian-exact, no extra fit) -- which it
#                       computes either way, as the cross-check.
#     --no-contour      (simfit/flavfit) SKIP the profiled (sigma_W, sigma_Z) scan.
#                       ON by default (2026-09-15b). ADDITIVE: same datacard, an
#                       extra map file that promotes the rapidity-inclusive
#                       sigma_W to a POI, its own workspace and its own output
#                       dir contour/contour_<B>/ -- it overwrites nothing from
#                       the nominal fit and the extraction never reads it.
#                       Needs skim/output/gen_xsec_fid.txt; missing -> the pass
#                       is skipped with a WARN, the other two still run.
#     --extract-only    (simfit/flavfit) re-run ONLY the extraction on an EXISTING
#                       fits/ tree (no combine; e.g. on a downloaded fit after an
#                       extractor change). Prefit integrals from the input copies
#                       in the work dir, else from the analysis plots dir -- the
#                       extractor checks them against the fit's own shapes_prefit
#     --plots-dir DIR   analysis plots dir (else $PO_PLOTS, else autodetect)
#     --out DIR         output root (default test/pO_fit_out)
#
# Needs `cmsenv` (combine + text2workspace.py on PATH) for the actual fits.
# bash-3.2 safe (macOS stock bash): no associative arrays.
# =============================================================================
set -uo pipefail   # NOT -e: one failed fit/pass must not abort the other fits

HERE="$(cd "$(dirname "$0")" && pwd)"
MYS="$HERE/my_script"

# A simfit run does THREE fits per binning variant by default (2026-09-15b):
#   (1) the nominal fit          -> fits/simfit_<B>/fitDiagnostics_simfit_<B>.root
#   (2) the --statonly companion -> ..._statonly.root  (THE stat error source)
#   (3) the --contour scan       -> contour/contour_<B>/  (sigma_W as a POI)
# (2) and (3) are additive: separate output files, nothing of (1) is overwritten.
# Turn them off with --no-statonly / --no-contour.
CHAN_ARG="both"; MODE="simfit"; DRYRUN=0; DO_POSTFIT=1; DRAWONLY=0; ASIMOV=0; STATONLY=1; EXTRACTONLY=0
CONTOUR=1; CONTOUR_POINTS="${CONTOUR_POINTS:-2500}"
DISC="leppt_mt40"; OUT_SET=0
OUTROOT="$HERE/pO_fit_out"
PO_PLOTS="${PO_PLOTS:-}"
PO_PLOTS_DEFAULTS="/Users/zhenghuang/pO_analysis/plotting/plots /afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots"

while [ $# -gt 0 ]; do
  case "$1" in
    mu|ele|both)               CHAN_ARG="$1" ;;
    simfit|flavfit|all)        MODE="$1" ;;
    perbin|incl|combined)
      echo "[ERROR] mode '$1' belonged to the legacy per-flavour per-bin pipeline, REMOVED 2026-09-22."
      echo "        The per-flavour fit is now:  ./run_pO_fits.sh [mu|ele|both] flavfit"
      exit 1 ;;
    --dry-run)                 DRYRUN=1 ;;
    --no-postfit)              DO_POSTFIT=0 ;;
    --draw-only)               DRAWONLY=1 ;;
    --asimov)                  ASIMOV=1 ;;
    --statonly)                STATONLY=1 ;;   # (the default since 2026-09-15b; kept as a no-op)
    --no-statonly)             STATONLY=0 ;;   # stat then falls back to the conditioned covariance
    --extract-only)            EXTRACTONLY=1 ;;
    --contour)                 CONTOUR=1 ;;    # (the default since 2026-09-15b; kept as a no-op)
    --no-contour)              CONTOUR=0 ;;    # skip the (sigmaW, r_Z) profiled 2D scan
    --plots-dir)               shift; PO_PLOTS="${1:-}" ;;
    --disc)                    shift; DISC="${1:?--disc needs a value: met|leppt|leppt_mt40}" ;;
    --out)                     shift; OUTROOT="${1:?--out needs a directory}"; OUT_SET=1 ;;
    -h|--help)                 sed -n '3,/bash-3.2 safe/p' "$0"; exit 0 ;;
    *) echo "[ERROR] unknown arg: $1"; exit 1 ;;
  esac
  shift
done

# ---- resolve the discriminant variant ---------------------------------------
case "$DISC" in
  met)        DSUF="";            XT="PF MET (GeV)"        ;;
  leppt)      DSUF="_leppt";      XT="Lepton p_{T} (GeV)"  ;;
  leppt_mt40) DSUF="_leppt_mt40"; XT="Lepton p_{T} (GeV)"  ;;
  *) echo "[ERROR] --disc must be met|leppt|leppt_mt40 (got '$DISC')"; exit 1 ;;
esac
YT="Events / 2.0 GeV"   # MET and lepton-pT templates are both 2 GeV bins
WINNAME="combine_input_W${DSUF}.root"
# unsuffixed out-tree for met (its historical name); suffixed tree for the variants (unless --out given)
if [ "$OUT_SET" -eq 0 ]; then OUTROOT="$HERE/pO_fit_out${DSUF}"; fi

# ---- locate the analysis plots dir (not needed for --draw-only: it reuses ---
# ---- the combine_input_*.root copies already in the work dir; for ----------
# ---- --extract-only it is only the FALLBACK when those copies are absent) ---
if [ "$DRAWONLY" -eq 0 ]; then
  if [ -z "$PO_PLOTS" ]; then
    for d in $PO_PLOTS_DEFAULTS; do
      if [ -f "$d/$WINNAME" ]; then PO_PLOTS="$d"; break; fi
    done
  fi
  if [ -z "$PO_PLOTS" ] || [ ! -d "$PO_PLOTS" ]; then
    if [ "$EXTRACTONLY" -eq 1 ]; then
      echo "[WARN] analysis plots dir not found -- --extract-only will rely on the input copies in the work dir."
      PO_PLOTS=""
    else
      echo "[ERROR] analysis plots dir not found. Set --plots-dir or \$PO_PLOTS to the"
      echo "        dir containing $WINNAME (run plotting/mtandmet.C +"
      echo "        dileptonpeak.C first)."
      exit 2
    fi
  fi
  [ -n "$PO_PLOTS" ] && echo "[run_pO_fits] plots dir : $PO_PLOTS"
fi
echo "[run_pO_fits] channel(s): $CHAN_ARG    mode: $MODE    disc: $DISC    dry-run: $DRYRUN    draw-only: $DRAWONLY    extract-only: $EXTRACTONLY    asimov: $ASIMOV    statonly (stat source): $STATONLY    contour: $CONTOUR"

# ---- cmsenv check -----------------------------------------------------------
HAVE_COMBINE=1
command -v combine          >/dev/null 2>&1 || HAVE_COMBINE=0
command -v text2workspace.py >/dev/null 2>&1 || HAVE_COMBINE=0
if [ "$DRYRUN" -eq 0 ] && [ "$DRAWONLY" -eq 0 ] && [ "$EXTRACTONLY" -eq 0 ] && [ "$HAVE_COMBINE" -eq 0 ]; then
  echo "[ERROR] combine / text2workspace.py not on PATH -- did you cmsenv?"
  echo "        (run with --dry-run to only build datacards.)"
  exit 3
fi
command -v root >/dev/null 2>&1 || { echo "[WARN] root not on PATH: extraction/postfit will be skipped."; }
if { [ "$DRAWONLY" -eq 1 ] || [ "$EXTRACTONLY" -eq 1 ]; } && ! command -v root >/dev/null 2>&1; then
  echo "[ERROR] --draw-only / --extract-only need root on PATH."; exit 3
fi

# =============================================================================
# simfit -- the GRAND SIMULTANEOUS FIT (2026-08-04).  ONE likelihood per
# binning variant (lab / fb): 48 W channels ({mu,ele} x {Wp,Wm} x y0..11) +
# both Z-inclusive peaks.  25 POIs (r_<C>_y<i> shared mu/e + global r_Z),
# QCD lnN per (flavour, charge) + global lumi lnN (2026-08-17; QCD_MODE=free
# restores per-channel free qcd_norm), w/wtau under Z frozen -- see
# my_script/make_pO_simfit_cards.sh for the model definition (card + t2w maps).
#
# The SAME functions run the per-flavour fits (mode flavfit, 2026-09-22):
# run_simfit_set "<flavours>" sets the work dir and the flavour list they all
# read -- "mu ele" -> simfit/ (the grand fit, comb_* files), mu -> simfit_mu/,
# ele -> simfit_ele/ (simfit_<flav>_* files). Nothing else differs, which is
# the point: the per-flavour results get exactly the grand fit's treatment.
# Globals set there: SFLAVS (space list), SFLAVS_CSV, SNAME (work-dir name),
# STAG (summary-file prefix), SWORK/SDCD/SFITS/SPOST/SSUMM/SCONT, and the input
# copies AWMU/AZMU/AWEL/AZEL ("none" for a flavour not in the fit).
# =============================================================================

fit_simfit() {  # $1 = lab | fb : workspace (multiSignalModel) + FitDiagnostics
  B="$1"
  card="$SDCD/datacard_simfit_${B}.txt"; maps="$SDCD/t2w_maps_simfit_${B}.txt"
  if [ ! -f "$card" ] || [ ! -f "$maps" ]; then echo "  [skip] no simfit card/maps for $B"; return; fi
  RD="$SFITS/simfit_$B"; mkdir -p "$RD"
  # the 25-POI model: one --PO map=... per line of the maps file
  PO=(-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose)
  while IFS= read -r m; do [ -n "$m" ] && PO+=(--PO "$m"); done < "$maps"
  (
    cd "$RD" || exit 1
    text2workspace.py "$card" -o workspace.root "${PO[@]}" >t2w.log 2>&1 \
      || { echo "  [FAIL t2w] simfit_$B (see $RD/t2w.log)"; exit 1; }
    # no --rMin/--rMax: there is no POI named 'r'; ranges come from the maps.
    # --skipBOnlyFit: a b-only fit (all 25 POIs at 0) is meaningless here.
    combine -M FitDiagnostics workspace.root \
            --saveShapes --saveWithUncertainties --skipBOnlyFit \
            -n "_simfit_${B}" --cminDefaultMinimizerStrategy 0 >fit.log 2>&1 \
      || { echo "  [FAIL fit] simfit_$B (see $RD/fit.log)"; exit 1; }
    if [ "$STATONLY" -eq 1 ]; then
      # OPTIONAL frozen-nuisance companion (--statonly, 2026-09-14/15): re-fit
      # with every CONSTRAINED nuisance frozen at its POST-FIT value (the
      # Combine breakdown recipe: freeze at the best fit, not at 0), so the
      # POIs land at the same minimum and their errors are the statistical
      # component. Since 2026-09-15 this is a CROSS-CHECK only: the extractor
      # derives the same stat component from the nominal fit's covariance
      # matrix (conditioned on the constrained nuisances -- Gaussian-exact, no
      # refit) and, when this file exists, prints the maximal deviation between
      # the two. The unconstrained rateParams (the in-fit ABCD scales, qcd_norm
      # in free mode) keep floating in both -- they are data-driven statistics.
      SETNUIS=$(root -l -b -q -e 'TFile f("fitDiagnostics_simfit_'"$B"'.root"); auto fr = (RooFitResult*)f.Get("fit_s"); TString s; if (fr) for (auto p : fr->floatParsFinal()) { TString n = p->GetName(); if (n.BeginsWith("r_") || n.BeginsWith("qcd_s") || n.BeginsWith("qcd_norm")) continue; s += TString::Format("%s%s=%.10g", s.Length() ? "," : "", n.Data(), ((RooRealVar*)p)->getVal()); } printf("SETNUIS %s\n", s.Data());' 2>/dev/null | awk '$1=="SETNUIS"{print $2}')
      if [ -n "$SETNUIS" ]; then SETOPT=(--setParameters "$SETNUIS"); echo "  [statonly] simfit_$B: nuisances frozen at their post-fit values ($(echo "$SETNUIS" | tr ',' '\n' | wc -l | tr -d ' ') parameters)"
      else SETOPT=(); echo "  [statonly] WARN simfit_$B: could not read the post-fit nuisance values -> frozen at their nominal values"; fi
      combine -M FitDiagnostics workspace.root \
              --skipBOnlyFit --freezeParameters allConstrainedNuisances "${SETOPT[@]}" \
              -n "_simfit_${B}_statonly" --cminDefaultMinimizerStrategy 0 >fit_statonly.log 2>&1 \
        || { echo "  [FAIL statonly] simfit_$B (see $RD/fit_statonly.log)"; exit 1; }
    fi
    if [ "$ASIMOV" -eq 1 ]; then
      # Prefit S+B Asimov closure: every fitted POI must return 1.
      # NB plain `-t -1` generates the BACKGROUND-ONLY Asimov (src/Combine.cc:844
      # falls back to mc_bonly when neither --expectSignal nor --setParameters is
      # given, and the auto-built b-only model sets ALL POIs to 0 -- killing every
      # POI-scaled template incl. wtau/z/ztau, so all POIs fit to ~0; seen
      # 2026-08-05). With >1 POI combine itself says to use --setParameters, so
      # inject ALL 25 POIs = 1 explicitly (the qcd_norm rateParams already
      # generate at their init value 1).
      SETPARS="r_Z=1"
      for C in Wp Wm; do for iy in $(seq 0 11); do SETPARS="$SETPARS,r_${C}_y${iy}=1"; done; done
      combine -M FitDiagnostics workspace.root \
              --skipBOnlyFit -t -1 --setParameters "$SETPARS" \
              -n "_simfit_${B}_asimov" --cminDefaultMinimizerStrategy 0 >fit_asimov.log 2>&1 \
        || { echo "  [FAIL asimov] simfit_$B (see $RD/fit_asimov.log)"; exit 1; }
    fi
  ) && echo "  [ok] simfit_$B"
}

# =============================================================================
# fit_contour -- the PROFILED (sigma_W, sigma_Z) contour (--contour, 2026-09-15)
#
# Same 62-channel datacard, different map file: t2w_maps_simfit_<B>_sigma.txt
# reparametrizes the 24 per-bin POIs as (sigmaW, 23 shape parameters) so the
# rapidity-inclusive fiducial W cross section is a POI in its own right
# (make_pO_simfit_cards.sh::write_sigma_maps). Profile likelihood is invariant
# under a reparametrization, so the best fit is IDENTICAL to the nominal one --
# what this buys is the exact profiled 2D region, of which the covariance
# ellipse that plotting/xsec_contour.C builds from the nominal fit is the
# Gaussian approximation.
#
# sigma_Z = r_Z * sigma_gen-fid,Z is a one-to-one rescaling of r_Z, so r_Z is
# scanned as-is and the Z axis is converted downstream -- no second POI needed.
#
# Two combine calls: `--algo singles` locates the minimum and gives the 1D
# errors, which set the grid window (+-4 sigma); then `--algo grid` scans.
# =============================================================================
fit_contour() {  # $1 = lab | fb
  B="$1"
  card="$SDCD/datacard_simfit_${B}.txt"; smaps="$SDCD/t2w_maps_simfit_${B}_sigma.txt"
  if [ ! -f "$card" ] || [ ! -f "$smaps" ]; then echo "  [skip] no simfit card / sigma maps for $B"; return; fi
  RD="$SCONT/contour_$B"; mkdir -p "$RD"
  PO=(-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose)
  while IFS= read -r m; do [ -n "$m" ] && PO+=(--PO "$m"); done < "$smaps"
  (
    cd "$RD" || exit 1
    text2workspace.py "$card" -o workspace_sigma.root "${PO[@]}" >t2w.log 2>&1 \
      || { echo "  [FAIL t2w] contour_$B (see $RD/t2w.log)"; exit 1; }
    combine -M MultiDimFit workspace_sigma.root --algo singles \
            -P sigmaW -P r_Z --floatOtherPOIs 1 --saveFitResult \
            -n "_contourfit_${B}" --cminDefaultMinimizerStrategy 0 >fit_singles.log 2>&1 \
      || { echo "  [FAIL singles] contour_$B (see $RD/fit_singles.log)"; exit 1; }
    # window = best fit +- 4 sigma (floored at 0 for sigmaW), from the fit result
    RNG=$(root -l -b -q -e 'TFile f("multidimfit_contourfit_'"$B"'.root"); auto fr = (RooFitResult*)f.Get("fit_mdf"); if (fr) { auto *s = (RooRealVar*)fr->floatParsFinal().find("sigmaW"); auto *z = (RooRealVar*)fr->floatParsFinal().find("r_Z"); if (s && z) printf("RNG sigmaW=%.6g,%.6g:r_Z=%.6g,%.6g BEST %.6g %.6g %.6g %.6g\n", s->getVal()-4*s->getError()>0 ? s->getVal()-4*s->getError() : 0.0, s->getVal()+4*s->getError(), z->getVal()-4*z->getError()>0 ? z->getVal()-4*z->getError() : 0.0, z->getVal()+4*z->getError(), s->getVal(), s->getError(), z->getVal(), z->getError()); }' 2>/dev/null | awk '$1=="RNG"{print $2}')
    BEST=$(root -l -b -q -e 'TFile f("multidimfit_contourfit_'"$B"'.root"); auto fr = (RooFitResult*)f.Get("fit_mdf"); if (fr) { auto *s = (RooRealVar*)fr->floatParsFinal().find("sigmaW"); auto *z = (RooRealVar*)fr->floatParsFinal().find("r_Z"); if (s && z) printf("BEST sigmaW %.5f +/- %.5f | r_Z %.5f +/- %.5f\n", s->getVal(), s->getError(), z->getVal(), z->getError()); }' 2>/dev/null | grep '^BEST' || true)
    [ -n "$BEST" ] && echo "  [contour_$B] $BEST"
    if [ -z "$RNG" ]; then
      echo "  [FAIL range] contour_$B: could not read sigmaW/r_Z from multidimfit_contourfit_${B}.root"; exit 1
    fi
    combine -M MultiDimFit workspace_sigma.root --algo grid --points "$CONTOUR_POINTS" \
            -P sigmaW -P r_Z --floatOtherPOIs 1 --saveNLL \
            --setParameterRanges "$RNG" \
            -n "_contour_${B}" --cminDefaultMinimizerStrategy 0 >fit_grid.log 2>&1 \
      || { echo "  [FAIL grid] contour_$B (see $RD/fit_grid.log)"; exit 1; }
    echo "  [contour_$B] grid: $CONTOUR_POINTS points over $RNG"

    # ---- the STAT-ONLY profiled contour (2026-09-15c, user request) ----------
    # Same scan with every CONSTRAINED nuisance frozen at its POST-FIT value --
    # the identical recipe the --statonly companion fit uses, so this is the
    # profiled twin of the dashed Gaussian stat-only ellipse. The post-fit
    # values come from the NOMINAL fit (the nuisances are named the same in
    # both workspaces; only the POI parametrization differs), and the
    # unconstrained rateParams (qcd_s*) keep floating -- they are statistics.
    # Its own +-4 sigma window: the stat region is ~3x smaller per axis, so
    # re-using the total window would spend most of the grid outside it.
    SND="$SFITS/simfit_$B/fitDiagnostics_simfit_${B}.root"
    if [ -f "$SND" ]; then
      SETNUIS=$(root -l -b -q -e 'TFile f("'"$SND"'"); auto fr = (RooFitResult*)f.Get("fit_s"); TString s; if (fr) for (auto p : fr->floatParsFinal()) { TString n = p->GetName(); if (n.BeginsWith("r_") || n.BeginsWith("a_") || n == "sigmaW" || n.BeginsWith("qcd_s") || n.BeginsWith("qcd_norm")) continue; s += TString::Format("%s%s=%.10g", s.Length() ? "," : "", n.Data(), ((RooRealVar*)p)->getVal()); } printf("SETNUIS %s\n", s.Data());' 2>/dev/null | awk '$1=="SETNUIS"{print $2}')
      if [ -n "$SETNUIS" ]; then SETOPT=(--setParameters "$SETNUIS"); else SETOPT=();
        echo "  [contour_$B] WARN could not read the post-fit nuisance values -> frozen at their nominal values"; fi
      combine -M MultiDimFit workspace_sigma.root --algo singles \
              -P sigmaW -P r_Z --floatOtherPOIs 1 --saveFitResult \
              --freezeParameters allConstrainedNuisances "${SETOPT[@]}" \
              -n "_contourstatfit_${B}" --cminDefaultMinimizerStrategy 0 >fit_singles_stat.log 2>&1 \
        || { echo "  [FAIL singles-stat] contour_$B (see $RD/fit_singles_stat.log)"; exit 1; }
      RNGS=$(root -l -b -q -e 'TFile f("multidimfit_contourstatfit_'"$B"'.root"); auto fr = (RooFitResult*)f.Get("fit_mdf"); if (fr) { auto *s = (RooRealVar*)fr->floatParsFinal().find("sigmaW"); auto *z = (RooRealVar*)fr->floatParsFinal().find("r_Z"); if (s && z) printf("RNG sigmaW=%.6g,%.6g:r_Z=%.6g,%.6g\n", s->getVal()-4*s->getError()>0 ? s->getVal()-4*s->getError() : 0.0, s->getVal()+4*s->getError(), z->getVal()-4*z->getError()>0 ? z->getVal()-4*z->getError() : 0.0, z->getVal()+4*z->getError()); }' 2>/dev/null | awk '$1=="RNG"{print $2}')
      if [ -z "$RNGS" ]; then
        echo "  [contour_$B] WARN no stat-only range -> stat contour skipped"
      else
        combine -M MultiDimFit workspace_sigma.root --algo grid --points "$CONTOUR_POINTS" \
                -P sigmaW -P r_Z --floatOtherPOIs 1 --saveNLL \
                --setParameterRanges "$RNGS" \
                --freezeParameters allConstrainedNuisances "${SETOPT[@]}" \
                -n "_contourstat_${B}" --cminDefaultMinimizerStrategy 0 >fit_grid_stat.log 2>&1 \
          || { echo "  [FAIL grid-stat] contour_$B (see $RD/fit_grid_stat.log)"; exit 1; }
        echo "  [contour_$B] stat-only grid: $CONTOUR_POINTS points over $RNGS"
      fi
    else
      echo "  [contour_$B] WARN $SND absent -> stat-only contour skipped"
    fi
  ) && echo "  [ok] contour_$B"
}

simfit_extract() {  # POIs + yields (summed over the fit's flavours) + covariance -> $SSUMM (root only)
  # lnN kappas + qcd mode the cards were built with (sidecar from
  # make_pO_simfit_cards.sh; missing sidecar / 0 entries -> legacy
  # free-rateParam extraction path; missing qcdMode line -> legacy sidecar,
  # mode inferred from the kappas inside the extractor)
  # lheSysts (2026-09-07): the LHE shape nuisances in the cards ("none"/absent
  # -> none); the extractor reports their pulls + closure.
  # Also called by --extract-only (2026-09-15) on an existing fits/ tree.
  KQM=0; KQE=0; KLU=0; QMODE=""; LHES=""; KF="$SDCD/qcd_lnn_kappas.txt"
  if [ -f "$KF" ]; then
    KQM=$(awk '$1=="kQcdMu"{print $2}' "$KF");  KQM="${KQM:-0}"
    KQE=$(awk '$1=="kQcdEle"{print $2}' "$KF"); KQE="${KQE:-0}"
    KLU=$(awk '$1=="kLumi"{print $2}' "$KF");   KLU="${KLU:-0}"
    QMODE=$(awk '$1=="qcdMode"{print $2}' "$KF"); QMODE="${QMODE:-}"
    LHES=$(awk '$1=="lheSysts"{print $2}' "$KF"); LHES="${LHES:-}"
    [ "$LHES" = "none" ] && LHES=""
    # flavours (2026-09-22): the cards must be the ones of THIS tree (a sidecar
    # without the line predates per-flavour fits = the grand mu,ele card)
    KFL=$(awk '$1=="flavours"{print $2}' "$KF"); KFL="${KFL:-mu,ele}"
    if [ "$KFL" != "$SFLAVS_CSV" ]; then
      echo "[extract] WARN the cards' sidecar says flavours '$KFL' but this is the $SNAME tree ($SFLAVS_CSV)"
      echo "          -- cards from another fit in this work dir? The extraction uses $SFLAVS_CSV."
    fi
  else
    echo "[extract] WARN no sidecar $KF -- legacy (free-rateParam) extraction assumed"
  fi
  if command -v root >/dev/null 2>&1; then
    # tee the FULL extraction output into summary/ (2026-09-15c): the stat-source
    # line, the statonly-vs-conditioned cross-check, the Asimov closure and every
    # WARN are printed, not stored in the CSVs -- and the fit runs on lxplus, so
    # without this they only ever exist in that terminal. summary/ is downloaded
    # wholesale, so the log travels with the results. (Repo convention: if a stage
    # produces numbers anyone quotes, it gets a log.)
    # (log named after the work dir: extract_simfit.log for the grand fit, as
    #  before; extract_simfit_<flav>.log for the per-flavour ones)
    root -b -q "$MYS/extract_pO_simfit.C(\"$SFITS\",\"$AWMU\",\"$AWEL\",\"$SSUMM\",$KQM,$KQE,$KLU,\"$QMODE\",\"$LHES\",\"$SFLAVS_CSV\",\"$STAG\")" 2>&1 \
      | tee "$SSUMM/extract_${SNAME}.log" \
      | grep -E "\[extract-simfit\]|\[asimov\]|WARN|FAIL" || true
    echo "[extract] full log: $SSUMM/extract_${SNAME}.log"
  fi
}

simfit_postfit_all() {  # postfit data/MC per channel of the fit ($SFLAVS), both variants
  # QCD info-box param: lnN mode (2026-08-17 default) -> the shared per
  # (flavour, charge) nuisance qcd_rate_<F>_<C>, whose displayed value is the
  # PULL theta (scale = kappa^theta); free mode -> the per-channel qcd_norm.
  QLNN=0
  if [ -f "$SDCD/qcd_lnn_kappas.txt" ]; then
    KQTEST=$(awk '$1=="kQcdMu"{print $2}' "$SDCD/qcd_lnn_kappas.txt")
    [ -n "$KQTEST" ] && [ "$KQTEST" != "0" ] && QLNN=1
  fi
  for B in lab fb; do
    fd="$SFITS/simfit_$B/fitDiagnostics_simfit_${B}.root"
    [ -f "$fd" ] || continue
    for F in $SFLAVS; do
      if [ "$F" = "ele" ]; then AW="$AWEL"; AZ="$AZEL"; WL="W #rightarrow e #nu"; ZL="Z #rightarrow e e"
      else                      AW="$AWMU"; AZ="$AZMU"; WL="W #rightarrow #mu #nu"; ZL="Z #rightarrow #mu #mu"; fi
      for C in Wp Wm; do
        for iy in $(seq 0 11); do
          R="${C}_${B}_y${iy}"; CH="${F}_${R}"
          QN="qcd_norm_${CH}"; [ "$QLNN" -eq 1 ] && QN="qcd_rate_${F}_${C}"
          # info box: this bin's POI, the global r_Z (shown as DY norm), the QCD
          # param (lnN mode: the shared nuisance -> value shown is the PULL);
          # ndf uses the ~3 params that shape this channel.
          root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"$CH\",\"$AW\",\"$R\",\"$SPOST/$CH\",\"$XT\",\"$YT\",\"$WL\",\"$CH (${SNAME} postfit)\",true,\"r_${C}_y${iy}\",\"r_Z\",\"$QN\",3)" >/dev/null 2>&1
        done
      done
      # the Z peak as seen by this variant's grand fit (r_Z only; W bkg frozen)
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"${F}_Z_incl\",\"$AZ\",\"Z_incl\",\"$SPOST/${F}_Z_incl_${B}\",\"m_{ll} (GeV)\",\"Events / 1.0 GeV\",\"$ZL\",\"Z incl (${SNAME} ${B} postfit)\",false,\"r_Z\",\"none\",\"none\",1)" >/dev/null 2>&1
    done
  done
}

run_simfit_set() {  # $1 = the flavours in ONE likelihood: "mu ele" (grand) | mu | ele (flavfit)
  SFLAVS="$1"
  case "$SFLAVS" in
    "mu ele") SNAME="simfit";     STAG="comb"
              STITLE="simfit: grand simultaneous fit (mu + ele)" ;;
    mu)       SNAME="simfit_mu";  STAG="simfit_mu"
              STITLE="simfit_mu: MUON-ONLY simultaneous fit (flavfit)" ;;
    ele)      SNAME="simfit_ele"; STAG="simfit_ele"
              STITLE="simfit_ele: ELECTRON-ONLY simultaneous fit (flavfit)" ;;
    *)        echo "[ERROR] run_simfit_set: unknown flavour set '$SFLAVS'"; return ;;
  esac
  SFLAVS_CSV=$(echo "$SFLAVS" | tr ' ' ',')
  echo ""
  echo "================ $STITLE ================"
  SWORK="$OUTROOT/$SNAME"; SDCD="$SWORK/datacards"; SFITS="$SWORK/fits"; SPOST="$SWORK/postfit"; SSUMM="$SWORK/summary"
  SCONT="$SWORK/contour"
  # input copies of the flavours IN the fit; "none" for the other one (neither
  # the card generator nor the extractor reads it)
  AWMU="none"; AZMU="none"; AWEL="none"; AZEL="none"
  case " $SFLAVS " in *" mu "*)  AWMU="$SWORK/combine_input_W_mu.root";  AZMU="$SWORK/combine_input_Z_mu.root" ;; esac
  case " $SFLAVS " in *" ele "*) AWEL="$SWORK/combine_input_W_ele.root"; AZEL="$SWORK/combine_input_Z_ele.root" ;; esac

  # ---- draw-only: redraw simfit postfit plots from an EXISTING run -----------
  if [ "$DRAWONLY" -eq 1 ]; then
    miss=0
    for f in "$AWMU" "$AWEL"; do [ "$f" = none ] || [ -f "$f" ] || miss=1; done
    if [ "$miss" -eq 1 ]; then
      echo "[ERROR] $SWORK input copies missing -- no previous $SNAME run (run the full fit first)."
      return
    fi
    nfd=$(find "$SFITS" -name 'fitDiagnostics_simfit_*.root' 2>/dev/null | wc -l | tr -d ' ')
    if [ "$nfd" -eq 0 ]; then
      echo "[ERROR] no fitDiagnostics_simfit_*.root under $SFITS -- --draw-only needs an earlier"
      echo "        simfit run (NB 'sync_lxplus.sh download' does NOT pull fits/; redraw where the fits ran)."
      return
    fi
    mkdir -p "$SPOST"
    echo "[draw-only] redrawing $SNAME postfit plots ..."
    simfit_postfit_all
    echo "[done] $SNAME -> $SPOST"
    return
  fi

  # ---- extract-only: re-run the extraction on an EXISTING fits/ tree ---------
  # (2026-09-15; no combine needed -- e.g. on a downloaded lxplus fit after an
  # extractor change.) The prefit signal integrals come from the input copies
  # in the work dir when present (as on lxplus), else from the analysis plots
  # dir; the extractor compares them with the fitted channels' shapes_prefit
  # and WARNs when they are not the inputs that were fitted.
  if [ "$EXTRACTONLY" -eq 1 ]; then
    if ! ls "$SFITS"/simfit_*/fitDiagnostics_simfit_*.root >/dev/null 2>&1; then
      echo "[ERROR] no fitDiagnostics_simfit_*.root under $SFITS -- --extract-only needs an earlier"
      echo "        run here (or 'sync_lxplus.sh download')."
      return
    fi
    if [ "$AWMU" != none ] && [ ! -f "$AWMU" ]; then
      if [ -n "$PO_PLOTS" ] && [ -f "$PO_PLOTS/$WINNAME" ]; then
        AWMU="$PO_PLOTS/$WINNAME"; echo "[extract-only] no muon W input copy in $SWORK -> $AWMU"
      else echo "[ERROR] --extract-only: no muon W input (neither a work-dir copy nor the plots dir)"; return; fi
    fi
    if [ "$AWEL" != none ] && [ ! -f "$AWEL" ]; then
      if [ -n "$PO_PLOTS" ] && [ -f "$PO_PLOTS/Elec/$WINNAME" ]; then
        AWEL="$PO_PLOTS/Elec/$WINNAME"; echo "[extract-only] no electron W input copy in $SWORK -> $AWEL"
      else echo "[ERROR] --extract-only: no electron W input (neither a work-dir copy nor the plots dir)"; return; fi
    fi
    mkdir -p "$SSUMM"
    echo "[extract-only] re-extracting from $SFITS ..."
    simfit_extract
    echo "[done] extract-only -> $SSUMM/${STAG}_W_yields.csv (+ ${STAG}_summary.csv, ${STAG}_fitted_yields.root)"
    return
  fi

  WMU_SRC="$PO_PLOTS/$WINNAME";      ZMU_SRC="$PO_PLOTS/combine_input_Z.root"
  WEL_SRC="$PO_PLOTS/Elec/$WINNAME"; ZEL_SRC="$PO_PLOTS/Elec/combine_input_Z.root"
  miss=0
  for pair in "$AWMU|$WMU_SRC" "$AZMU|$ZMU_SRC" "$AWEL|$WEL_SRC" "$AZEL|$ZEL_SRC"; do
    [ "${pair%%|*}" = none ] && continue          # flavour not in this fit
    [ -f "${pair#*|}" ] || { echo "[ERROR] $SNAME input missing: ${pair#*|}"; miss=1; }
  done
  if [ "$miss" -eq 1 ]; then
    echo "[ERROR] $SNAME needs the W AND Z inputs of every flavour it fits (no Z fallback: r_Z is pinned by the peak) -- skipped."
    return
  fi

  mkdir -p "$SWORK" "$SDCD" "$SFITS" "$SPOST" "$SSUMM"
  # --contour (2026-09-15): ask the card generator for the SECOND map file that
  # promotes the rapidity-inclusive sigma_W to a POI. It needs the gen fiducial
  # cross sections, which live in the ANALYSIS repo next to the plots dir
  # (skim/output/gen_xsec_fid.txt, written by skim/gen_xsec.C and uploaded by
  # sync_lxplus.sh). The card itself is unchanged either way.
  if [ "$CONTOUR" -eq 1 ]; then
    export SIGMA_POI=1
    export GEN_XSEC_FID="${GEN_XSEC_FID:-$PO_PLOTS/../../skim/output/gen_xsec_fid.txt}"
    if [ ! -f "$GEN_XSEC_FID" ]; then
      # NOT fatal: the contour is one of three passes and the nominal fit must
      # not die because an optional input is missing.
      echo "[WARN] --contour needs the gen fiducial sidecar, not found: $GEN_XSEC_FID"
      echo "       run skim/gen_xsec.C (analysis repo), then sync_lxplus.sh upload; or set GEN_XSEC_FID."
      echo "       -> contour pass SKIPPED; the nominal and stat-only fits run as usual."
      CONTOUR=0; unset SIGMA_POI
    else
      mkdir -p "$SCONT"
      echo "[contour] sigma_W promoted to a POI; gen sidecar: $GEN_XSEC_FID"
    fi
  fi
  # LHE shape-systematics sidecars (2026-09-07, <input minus .root>_systs.txt):
  # travel with the inputs so the card generator finds them next to the copies
  # (absent -> no shape rows; a stale copy is removed so it cannot lie).
  # (only the flavours in the fit: a per-flavour work dir holds one flavour's
  #  inputs, so nothing of the other can leak into its cards)
  for pair in "$WMU_SRC|$AWMU" "$ZMU_SRC|$AZMU" "$WEL_SRC|$AWEL" "$ZEL_SRC|$AZEL"; do
    psrc="${pair%%|*}"; pdst="${pair##*|}"
    [ "$pdst" = none ] && continue
    if [ -f "${psrc%.root}_systs.txt" ]; then cp -f "${psrc%.root}_systs.txt" "${pdst%.root}_systs.txt"
    else rm -f "${pdst%.root}_systs.txt"; fi
    cp -f "$psrc" "$pdst"
  done

  # absolute-path datacards so combine resolves shapes from any CWD.
  # The generator can refuse (e.g. QCD_MODE=abcd with a non-leppt_mt40 disc);
  # the script runs without -e, so check explicitly or the fit stage would run
  # on stale/absent cards. SIMFIT_FLAVS restricts the card to this fit's
  # flavours (the grand card is byte-identical to the pre-2026-09-22 one).
  if ! SIMFIT_FLAVS="$SFLAVS" /bin/bash "$MYS/make_pO_simfit_cards.sh" "$AWMU" "$AZMU" "$AWEL" "$AZEL" "$SDCD" "$DISC"; then
    echo "[ERROR] $SNAME datacard generation failed -- $SNAME skipped."
    return
  fi

  if [ "$DRYRUN" -eq 1 ]; then
    echo "[dry-run] $SNAME datacards + t2w maps in $SDCD ; skipping fits."
    return
  fi

  for B in lab fb; do fit_simfit "$B"; done
  if [ "$CONTOUR" -eq 1 ]; then for B in lab fb; do fit_contour "$B"; done; fi

  # ---- extract POIs + yields (summed over the fit's flavours) + covariance ----
  simfit_extract

  # ---- postfit plots (per channel of the fit: 2 variants x 25 per flavour) ----
  if [ "$DO_POSTFIT" -eq 1 ] && command -v root >/dev/null 2>&1; then
    echo "[postfit] drawing $SNAME postfit plots (2 variants x 25 channels per flavour: $SFLAVS) ..."
    simfit_postfit_all
  fi

  echo "[done] $SNAME -> $SWORK"
  echo "       yields CSV     : $SSUMM/${STAG}_W_yields.csv"
  echo "       POI summary    : $SSUMM/${STAG}_summary.csv"
  echo "       analysis input : $SSUMM/${STAG}_fitted_yields.root  (+ h_cov_yield[_FB], h_cov_poi[_FB])"
  echo "       -> the observables chain for this discriminant, locally:"
  echo "          (analysis repo) analysis/run_observables.sh $DISC"
}

case "$CHAN_ARG" in
  both) CHANS="mu ele" ;;
  *)    CHANS="$CHAN_ARG" ;;
esac
case "$MODE" in
  simfit)    # the grand fit, mu + e in one likelihood
    [ "$CHAN_ARG" != "both" ] && echo "[note] simfit is the mu+e GRAND fit -- channel '$CHAN_ARG' ignored (a $CHAN_ARG-only fit: ./run_pO_fits.sh $CHAN_ARG flavfit)."
    run_simfit_set "mu ele" ;;
  flavfit)   # one simultaneous fit per flavour of the channel argument
    for c in $CHANS; do run_simfit_set "$c"; done ;;
  all)       # every fit: the grand fit + the per-flavour fits
    if [ "$CHAN_ARG" = "both" ]; then run_simfit_set "mu ele"
    else echo "[note] mode 'all' with channel '$CHAN_ARG': the grand simfit needs both flavours -> skipped ($CHAN_ARG flavfit still runs)."; fi
    for c in $CHANS; do run_simfit_set "$c"; done ;;
esac
echo ""
echo "[run_pO_fits] all done. Output under: $OUTROOT"
