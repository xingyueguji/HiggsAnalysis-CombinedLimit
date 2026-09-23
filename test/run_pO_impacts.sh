#!/usr/bin/env bash
# =============================================================================
# run_pO_impacts.sh -- nuisance IMPACT plots + fit COVARIANCE/correlation plots
# for the grand simultaneous fit (simfit) or, with --fit mu|ele, for one of the
# per-flavour fits (flavfit).  Run AFTER ./run_pO_fits.sh: it needs each
# variant's workspace.root (+ fitDiagnostics for the covariance) under
# pO_fit_out<suffix>/<fit dir>/fits/simfit_{lab,fb}/.
#
# Impacts (combineTool.py -M Impacts + plotImpacts.py, both SHIPPED with
# Combine v10 in <fork>/scripts/ -- no CombineHarvester checkout needed):
#   1. initial MultiDimFit (all 25 POIs profiled together)
#   2. per-nuisance fits: each constrained nuisance fixed at its +-1sigma
#      postfit values, everything else re-profiled (2026-08-17 lnN model:
#      qcd_rate_{mu,ele}_{Wp,Wm} x4 + lumi = 5 nuisances -> ~11 fits/variant)
#   3. impacts_simfit_<B>.json -> one PDF per requested POI (plotImpacts --POI)
#
# Covariance (my_script/plot_pO_cov.C, needs only root -- can also be re-run
# locally on a downloaded tree, where fitDiagnostics is absent: it then skips
# the parameter matrix and still draws the yield correlation from summary/):
#   - full floating-parameter correlation matrix from fit_s (POIs + nuisances)
#   - the 24x24 fitted-YIELD correlation from summary/<prefix>_fitted_yields.root
#     (comb_ for the grand fit, simfit_<flav>_ for the per-flavour ones)
#     (h_cov_yield / h_cov_yield_FB, the matrices the downstream error
#      propagation actually uses)
#
# Usage (impacts need cmsenv; --cov-only needs only root):
#   ./run_pO_impacts.sh [--disc met|leppt|leppt_mt40] [--variant lab|fb|both]
#                       [--fit comb|mu|ele|all]       (default: comb = the grand fit)
#                       [--pois "r_Wp_y3 r_Z ..."]    (default: ALL 25 POIs)
#                       [--impacts-only | --cov-only] [--dry-run]
# --fit (2026-09-22) picks the fit: comb = pO_fit_out<suffix>/simfit/, mu / ele
# = the per-flavour fits of run_pO_fits.sh mode flavfit (simfit_mu/, simfit_ele/;
# same layout, same 25 POIs, summary file simfit_<flav>_fitted_yields.root).
# --fit all = comb + mu + ele in one go (= the fit driver's mode all); a list
# ("mu,ele") picks several. A fit whose work dir does not exist (e.g. flavfit
# never ran) is skipped with a note, so `all` also works on a grand-fit-only tree.
# Outputs:  pO_fit_out<suffix>/<fit dir>/impacts/   (json + per-POI PDFs)
#           pO_fit_out<suffix>/<fit dir>/cov/       (correlation-matrix png+pdf)
# Both directories are pulled back to the Mac by  ./sync_lxplus.sh download.
# bash-3.2 safe (macOS stock bash): no associative arrays.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
MYS="$HERE/my_script"

DISC="leppt_mt40"; VARIANTS="lab fb"; POIS="all"; DO_IMP=1; DO_COV=1; DRY=0; FIT="comb"
while [ $# -gt 0 ]; do
  case "$1" in
    --disc)    shift; DISC="${1:-}";;
    --fit)     shift; FIT="${1:-}";;
    --variant) shift
               case "${1:-both}" in
                 both) VARIANTS="lab fb";;
                 lab|fb) VARIANTS="$1";;
                 *) echo "[ERROR] --variant must be lab|fb|both"; exit 1;;
               esac;;
    --pois)    shift; POIS="${1:-all}";;
    --impacts-only) DO_COV=0;;
    --cov-only)     DO_IMP=0;;
    --dry-run) DRY=1;;
    # the WHOLE header, however long it grows (a fixed range truncates it)
    -h|--help) sed -n '3,/^# ===/p' "$0" | sed '$d'; exit 0;;
    *) echo "[ERROR] unknown option: $1"; exit 1;;
  esac
  shift
done

case "$DISC" in
  met)        SFX="";;
  leppt)      SFX="_leppt";;
  leppt_mt40) SFX="_leppt_mt40";;
  *) echo "[ERROR] unknown --disc '$DISC' (met|leppt|leppt_mt40)"; exit 1;;
esac
# the fit(s): all = the grand fit + both per-flavour fits; a comma- or
# space-separated list picks several (validated up front, before any work)
case "$FIT" in
  all) FITS="comb mu ele";;
  *)   FITS="$(echo "$FIT" | tr ',' ' ')";;
esac
[ -n "$(echo $FITS)" ] || { echo "[ERROR] --fit needs comb|mu|ele|all"; exit 1; }
for F in $FITS; do
  case "$F" in
    comb|mu|ele) ;;
    *) echo "[ERROR] unknown --fit '$F' (comb|mu|ele|all, or a list like mu,ele)"; exit 1;;
  esac
done

# default POI list = the 25 POIs of the simfit model
if [ "$POIS" = "all" ]; then
  POIS="r_Z"
  for C in Wp Wm; do
    for iy in 0 1 2 3 4 5 6 7 8 9 10 11; do POIS="$POIS r_${C}_y${iy}"; done
  done
fi

run() { echo "+ $*"; if [ "$DRY" -eq 0 ]; then "$@"; fi; }

err=0; DONE_FITS=""
# one fit: its work dir + the prefix of its summary files (run_pO_fits.sh)
do_fit() {
  local F="$1" MODE
  case "$F" in
    comb) FDIR="simfit";     FPRE="comb";       MODE="simfit";;
    mu)   FDIR="simfit_mu";  FPRE="simfit_mu";  MODE="mu flavfit";;
    ele)  FDIR="simfit_ele"; FPRE="simfit_ele"; MODE="ele flavfit";;
  esac
  SWORK="$HERE/pO_fit_out${SFX}/$FDIR"
  IMP="$SWORK/impacts"; COV="$SWORK/cov"
  # an absent fit is a note, not an error: --fit all on a tree where flavfit
  # never ran must still do the grand fit
  if [ ! -d "$SWORK" ]; then
    echo "[skip] fit $F: no $SWORK (run ./run_pO_fits.sh $MODE --disc $DISC first)"
    return
  fi

  for B in $VARIANTS; do
    RD="$SWORK/fits/simfit_$B"
    WS="$RD/workspace.root"
    FD="$RD/fitDiagnostics_simfit_${B}.root"

    # ---- impacts ------------------------------------------------------------
    if [ "$DO_IMP" -eq 1 ]; then
      if [ ! -f "$WS" ]; then
        echo "[skip] impacts $FDIR/$B: no $WS (run ./run_pO_fits.sh $MODE --disc $DISC first)"
      elif ! command -v combineTool.py >/dev/null 2>&1; then
        echo "[ERROR] combineTool.py not on PATH -- impacts need cmsenv"; err=1
      else
        echo "== impacts $FDIR/simfit_$B ($DISC) =="
        WD="$IMP/wd_$B"; mkdir -p "$WD"     # contains the many intermediate roots
        JSON="$IMP/impacts_simfit_${B}.json"
        (
          cd "$WD" || exit 1
          run combineTool.py -M Impacts -d "$WS" -m 125 \
              --robustFit 1 --cminDefaultMinimizerStrategy 0 --doInitialFit \
            && run combineTool.py -M Impacts -d "$WS" -m 125 \
                --robustFit 1 --cminDefaultMinimizerStrategy 0 --doFits \
            && run combineTool.py -M Impacts -d "$WS" -m 125 -o "$JSON"
        ) || { echo "[FAIL impacts] $FDIR/simfit_$B (see $WD)"; err=1; continue; }
        if [ "$DRY" -eq 0 ] && [ ! -f "$JSON" ]; then
          echo "[FAIL impacts] $FDIR/simfit_$B: no $JSON produced"; err=1; continue
        fi
        for P in $POIS; do
          run plotImpacts.py -i "$JSON" -o "$IMP/impacts_simfit_${B}_${P}" --POI "$P" \
            || { echo "[warn] plotImpacts failed for POI $P ($FDIR/$B)"; err=1; }
        done
        echo "  -> $IMP/impacts_simfit_${B}_<poi>.pdf"
      fi
    fi

    # ---- covariance / correlation plots ------------------------------------
    if [ "$DO_COV" -eq 1 ]; then
      if ! command -v root >/dev/null 2>&1; then
        echo "[ERROR] root not on PATH -- covariance plots need it"; err=1
      else
        echo "== covariance $FDIR/simfit_$B ($DISC) =="
        mkdir -p "$COV"
        SUMR="$SWORK/summary/${FPRE}_fitted_yields.root"
        run root -b -q "$MYS/plot_pO_cov.C(\"$FD\",\"$SUMR\",\"$B\",\"$COV\")" \
          || { echo "[FAIL cov] $FDIR/simfit_$B"; err=1; }
      fi
    fi
  done
  DONE_FITS="$DONE_FITS $FDIR"
}

for F in $FITS; do do_fit "$F"; done

echo ""
if [ "$err" -ne 0 ]; then echo "[run_pO_impacts] FINISHED WITH ERRORS (see above)"; exit 1; fi
if [ -z "$DONE_FITS" ]; then
  echo "[run_pO_impacts] nothing done: no work dir for --fit $FIT under pO_fit_out${SFX}/"; exit 1
fi
for d in $DONE_FITS; do echo "[run_pO_impacts] done -> pO_fit_out${SFX}/$d/{impacts,cov}"; done
echo "  pull them to the Mac with:  ./sync_lxplus.sh download"
