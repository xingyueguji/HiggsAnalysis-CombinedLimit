#!/usr/bin/env bash
# =============================================================================
# run_pO_impacts.sh -- nuisance IMPACT plots + fit COVARIANCE/correlation plots
# for the grand simultaneous fit (simfit).  Run AFTER ./run_pO_fits.sh: it needs
# each variant's workspace.root (+ fitDiagnostics for the covariance) under
# pO_fit_out<suffix>/simfit/fits/simfit_{lab,fb}/.
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
#   - the 24x24 fitted-YIELD correlation from summary/comb_fitted_yields.root
#     (h_cov_yield / h_cov_yield_FB, the matrices the downstream error
#      propagation actually uses)
#
# Usage (impacts need cmsenv; --cov-only needs only root):
#   ./run_pO_impacts.sh [--disc met|leppt|leppt_mt40] [--variant lab|fb|both]
#                       [--pois "r_Wp_y3 r_Z ..."]    (default: ALL 25 POIs)
#                       [--impacts-only | --cov-only] [--dry-run]
# Outputs:  pO_fit_out<suffix>/simfit/impacts/   (json + per-POI PDFs)
#           pO_fit_out<suffix>/simfit/cov/       (correlation-matrix png+pdf)
# Both directories are pulled back to the Mac by  ./sync_lxplus.sh download.
# bash-3.2 safe (macOS stock bash): no associative arrays.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
MYS="$HERE/my_script"

DISC="leppt_mt40"; VARIANTS="lab fb"; POIS="all"; DO_IMP=1; DO_COV=1; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --disc)    shift; DISC="${1:-}";;
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
    -h|--help) sed -n '2,31p' "$0"; exit 0;;
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
SWORK="$HERE/pO_fit_out${SFX}/simfit"
IMP="$SWORK/impacts"; COV="$SWORK/cov"

# default POI list = the 25 POIs of the simfit model
if [ "$POIS" = "all" ]; then
  POIS="r_Z"
  for C in Wp Wm; do
    for iy in 0 1 2 3 4 5 6 7 8 9 10 11; do POIS="$POIS r_${C}_y${iy}"; done
  done
fi

run() { echo "+ $*"; if [ "$DRY" -eq 0 ]; then "$@"; fi; }

err=0
for B in $VARIANTS; do
  RD="$SWORK/fits/simfit_$B"
  WS="$RD/workspace.root"
  FD="$RD/fitDiagnostics_simfit_${B}.root"

  # ---- impacts --------------------------------------------------------------
  if [ "$DO_IMP" -eq 1 ]; then
    if [ ! -f "$WS" ]; then
      echo "[skip] impacts $B: no $WS (run ./run_pO_fits.sh --disc $DISC first)"
    elif ! command -v combineTool.py >/dev/null 2>&1; then
      echo "[ERROR] combineTool.py not on PATH -- impacts need cmsenv"; err=1
    else
      echo "== impacts simfit_$B ($DISC) =="
      WD="$IMP/wd_$B"; mkdir -p "$WD"     # contains the many intermediate roots
      JSON="$IMP/impacts_simfit_${B}.json"
      (
        cd "$WD" || exit 1
        run combineTool.py -M Impacts -d "$WS" -m 125 \
            --robustFit 1 --cminDefaultMinimizerStrategy 0 --doInitialFit \
          && run combineTool.py -M Impacts -d "$WS" -m 125 \
              --robustFit 1 --cminDefaultMinimizerStrategy 0 --doFits \
          && run combineTool.py -M Impacts -d "$WS" -m 125 -o "$JSON"
      ) || { echo "[FAIL impacts] simfit_$B (see $WD)"; err=1; continue; }
      if [ "$DRY" -eq 0 ] && [ ! -f "$JSON" ]; then
        echo "[FAIL impacts] simfit_$B: no $JSON produced"; err=1; continue
      fi
      for P in $POIS; do
        run plotImpacts.py -i "$JSON" -o "$IMP/impacts_simfit_${B}_${P}" --POI "$P" \
          || { echo "[warn] plotImpacts failed for POI $P ($B)"; err=1; }
      done
      echo "  -> $IMP/impacts_simfit_${B}_<poi>.pdf"
    fi
  fi

  # ---- covariance / correlation plots --------------------------------------
  if [ "$DO_COV" -eq 1 ]; then
    if ! command -v root >/dev/null 2>&1; then
      echo "[ERROR] root not on PATH -- covariance plots need it"; err=1
    else
      echo "== covariance simfit_$B ($DISC) =="
      mkdir -p "$COV"
      SUMR="$SWORK/summary/comb_fitted_yields.root"
      run root -b -q "$MYS/plot_pO_cov.C(\"$FD\",\"$SUMR\",\"$B\",\"$COV\")" \
        || { echo "[FAIL cov] simfit_$B"; err=1; }
    fi
  fi
done

echo ""
if [ "$err" -ne 0 ]; then echo "[run_pO_impacts] FINISHED WITH ERRORS (see above)"; exit 1; fi
echo "[run_pO_impacts] done -> $IMP , $COV"
echo "  pull them to the Mac with:  ./sync_lxplus.sh download"
