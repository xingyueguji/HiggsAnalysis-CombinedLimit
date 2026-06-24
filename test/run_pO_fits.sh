#!/usr/bin/env bash
# =============================================================================
# run_pO_fits.sh -- end-to-end pO W/Z Combine fits, automated.
#
# For each lepton channel (muon / electron) it:
#   1. locates the structured Combine inputs from the analysis repo
#        <plots>/combine_input_W.root , <plots>/combine_input_Z.root      (muon)
#        <plots>/Elec/combine_input_{W,Z}.root                            (electron)
#   2. generates ALL datacards (per-(charge,y) lab + FB, per-charge incl,
#      W_incl, Z_incl, and the simultaneous W+Z card)        [make_pO_datacards.sh]
#   3. runs text2workspace + combine -M FitDiagnostics per region, each in its
#      own clean output subdir
#   4. extracts fitted signal yields -> CSV + h_mt_W{p,m}_y{..}(_FB) histograms
#      that analysis/charge_asym.C and analysis/FBratio.C read directly
#      [extract_pO_yields.C]
#   5. draws postfit data/MC plots in the SAME cosmetics as plotting/mtandmet.C
#      [draw_postfit_pO.C]
#
# Output tree (clean, one dir per region):
#   pO_fit_out/<chan>/{combine_input_*.root, datacards/, fits/<region>/,
#                      postfit/, summary/<chan>_W_yields.csv,
#                      <chan>_summary.csv, <chan>_fitted_yields.root}
#
# Usage:
#   ./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|all] [options]
#     channel  (default both)   mode (default all)
#   options:
#     --dry-run         build datacards + check inputs only (no cmsenv needed)
#     --no-postfit      skip the postfit plots (faster)
#     --plots-dir DIR   analysis plots dir (else $PO_PLOTS, else autodetect)
#     --out DIR         output root (default test/pO_fit_out)
#
# Needs `cmsenv` (combine + text2workspace.py on PATH) for the actual fits.
# bash-3.2 safe (macOS stock bash): no associative arrays.
# =============================================================================
set -uo pipefail   # NOT -e: per-bin fit failures must not abort the whole loop

HERE="$(cd "$(dirname "$0")" && pwd)"
MYS="$HERE/my_script"

CHAN_ARG="both"; MODE="all"; DRYRUN=0; DO_POSTFIT=1
OUTROOT="$HERE/pO_fit_out"
PO_PLOTS="${PO_PLOTS:-}"
PO_PLOTS_DEFAULTS="/Users/zhenghuang/pO_analysis/plotting/plots /afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots"

while [ $# -gt 0 ]; do
  case "$1" in
    mu|ele|both)               CHAN_ARG="$1" ;;
    perbin|incl|combined|all)  MODE="$1" ;;
    --dry-run)                 DRYRUN=1 ;;
    --no-postfit)              DO_POSTFIT=0 ;;
    --plots-dir)               shift; PO_PLOTS="${1:-}" ;;
    --out)                     shift; OUTROOT="${1:-$OUTROOT}" ;;
    -h|--help)                 sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "[ERROR] unknown arg: $1"; exit 1 ;;
  esac
  shift
done

# ---- locate the analysis plots dir -----------------------------------------
if [ -z "$PO_PLOTS" ]; then
  for d in $PO_PLOTS_DEFAULTS; do
    if [ -f "$d/combine_input_W.root" ]; then PO_PLOTS="$d"; break; fi
  done
fi
if [ -z "$PO_PLOTS" ] || [ ! -d "$PO_PLOTS" ]; then
  echo "[ERROR] analysis plots dir not found. Set --plots-dir or \$PO_PLOTS to the"
  echo "        dir containing combine_input_W.root (run plotting/mtandmet.C +"
  echo "        dileptonpeak.C first)."
  exit 2
fi
echo "[run_pO_fits] plots dir : $PO_PLOTS"
echo "[run_pO_fits] channel(s): $CHAN_ARG    mode: $MODE    dry-run: $DRYRUN"

# ---- cmsenv check -----------------------------------------------------------
HAVE_COMBINE=1
command -v combine          >/dev/null 2>&1 || HAVE_COMBINE=0
command -v text2workspace.py >/dev/null 2>&1 || HAVE_COMBINE=0
if [ "$DRYRUN" -eq 0 ] && [ "$HAVE_COMBINE" -eq 0 ]; then
  echo "[ERROR] combine / text2workspace.py not on PATH -- did you cmsenv?"
  echo "        (run with --dry-run to only build datacards.)"
  exit 3
fi
command -v root >/dev/null 2>&1 || { echo "[WARN] root not on PATH: extraction/postfit will be skipped."; }

# ---- region list for a mode -------------------------------------------------
build_regions() {  # echoes space-separated region labels (excludes the WZ combo)
  m="$1"; out=""
  if [ "$m" = "perbin" ] || [ "$m" = "all" ]; then
    for C in Wp Wm; do for B in lab fb; do for iy in $(seq 0 11); do out="$out ${C}_${B}_y${iy}"; done; done; done
  fi
  if [ "$m" = "incl" ] || [ "$m" = "all" ]; then
    out="$out Wp_incl Wm_incl W_incl Z_incl"
  fi
  echo "$out"
}
want_combined() { [ "$1" = "combined" ] || [ "$1" = "all" ]; }

# ---- one fit region: datacard -> workspace -> FitDiagnostics ----------------
fit_region() {  # $1 = region label, $2 = fits base dir, $3 = datacards dir
  R="$1"; FITS="$2"; DCD="$3"
  card="$DCD/datacard_${R}.txt"
  if [ ! -f "$card" ]; then echo "  [skip] no datacard for $R"; return; fi
  RD="$FITS/$R"; mkdir -p "$RD"
  (
    cd "$RD" || exit 1
    text2workspace.py "$card" -o workspace.root >t2w.log 2>&1 || { echo "  [FAIL t2w] $R (see $RD/t2w.log)"; exit 1; }
    combine -M FitDiagnostics workspace.root \
            --saveShapes --saveWithUncertainties \
            -n "_${R}" --rMin 0 --rMax 20 \
            --cminDefaultMinimizerStrategy 0 >fit.log 2>&1 \
      || { echo "  [FAIL fit] $R (see $RD/fit.log)"; exit 1; }
  ) && echo "  [ok] $R"
}

# ---- postfit plot for a region ---------------------------------------------
postfit_region() {  # $1=region(label/fitChannel) $2=fits $3=post $4=absW $5=absZ $6=chlabel $7=zlabel
  R="$1"; FITS="$2"; POST="$3"; AW="$4"; AZ="$5"; WL="$6"; ZL="$7"
  fd="$FITS/$R/fitDiagnostics_${R}.root"
  [ -f "$fd" ] || return
  case "$R" in
    Z_incl)
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"Z_incl\",\"$AZ\",\"Z_incl\",\"$POST/$R\",\"m_{ll} (GeV)\",\"Events / 1.0 GeV\",\"$ZL\",\"$R (postfit)\",false)" >/dev/null 2>&1 ;;
    WZ)
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"Wincl\",\"$AW\",\"W_incl\",\"$POST/WZ_Wincl\",\"PF MET (GeV)\",\"Events / 2.0 GeV\",\"$WL\",\"W+Z combined: W (postfit)\",true)"  >/dev/null 2>&1
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"Zincl\",\"$AZ\",\"Z_incl\",\"$POST/WZ_Zincl\",\"m_{ll} (GeV)\",\"Events / 1.0 GeV\",\"$ZL\",\"W+Z combined: Z (postfit)\",false)" >/dev/null 2>&1 ;;
    *)
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"$R\",\"$AW\",\"$R\",\"$POST/$R\",\"PF MET (GeV)\",\"Events / 2.0 GeV\",\"$WL\",\"$R (postfit)\",true)" >/dev/null 2>&1 ;;
  esac
}

# ---- one channel ------------------------------------------------------------
run_channel() {
  chan="$1"
  if [ "$chan" = "ele" ]; then
    WIN_SRC="$PO_PLOTS/Elec/combine_input_W.root"; ZIN_SRC="$PO_PLOTS/Elec/combine_input_Z.root"
    WL="W #rightarrow e #nu"; ZL="Z #rightarrow e e"
  else
    WIN_SRC="$PO_PLOTS/combine_input_W.root"; ZIN_SRC="$PO_PLOTS/combine_input_Z.root"
    WL="W #rightarrow #mu #nu"; ZL="Z #rightarrow #mu #mu"
  fi
  echo ""
  echo "================ channel: $chan ================"
  if [ ! -f "$WIN_SRC" ]; then echo "[ERROR] missing $WIN_SRC"; return; fi
  if [ ! -f "$ZIN_SRC" ]; then echo "[WARN] missing $ZIN_SRC (Z + combined fits will be skipped)"; fi

  WORK="$OUTROOT/$chan"; DCD="$WORK/datacards"; FITS="$WORK/fits"; POST="$WORK/postfit"; SUMM="$WORK/summary"
  mkdir -p "$WORK" "$DCD" "$FITS" "$POST" "$SUMM"
  cp -f "$WIN_SRC" "$WORK/combine_input_W.root"
  [ -f "$ZIN_SRC" ] && cp -f "$ZIN_SRC" "$WORK/combine_input_Z.root"
  ABS_W="$WORK/combine_input_W.root"; ABS_Z="$WORK/combine_input_Z.root"

  # absolute-path datacards so combine resolves shapes from any CWD
  /bin/bash "$MYS/make_pO_datacards.sh" "$ABS_W" "$ABS_Z" "$DCD"

  if [ "$DRYRUN" -eq 1 ]; then
    echo "[dry-run] datacards in $DCD ; skipping fits."
    return
  fi

  REGIONS="$(build_regions "$MODE")"
  echo "[fit] regions: $(echo $REGIONS | wc -w) ; combined: $(want_combined "$MODE" && echo yes || echo no)"
  for R in $REGIONS; do fit_region "$R" "$FITS" "$DCD"; done
  if want_combined "$MODE" && [ -f "$ABS_Z" ]; then fit_region "WZ" "$FITS" "$DCD"; fi

  # ---- extract machine-readable yields + analysis histos ----
  if command -v root >/dev/null 2>&1; then
    root -b -q "$MYS/extract_pO_yields.C(\"$chan\",\"$FITS\",\"$ABS_W\",\"$ABS_Z\",\"$SUMM\")" 2>&1 | grep -E "\[extract\]|WARN" || true
  fi

  # ---- postfit plots ----
  if [ "$DO_POSTFIT" -eq 1 ] && command -v root >/dev/null 2>&1; then
    echo "[postfit] drawing ..."
    for R in $REGIONS; do postfit_region "$R" "$FITS" "$POST" "$ABS_W" "$ABS_Z" "$WL" "$ZL"; done
    if want_combined "$MODE" && [ -f "$ABS_Z" ]; then postfit_region "WZ" "$FITS" "$POST" "$ABS_W" "$ABS_Z" "$WL" "$ZL"; fi
  fi

  echo "[done] $chan -> $WORK"
  echo "       yields CSV : $SUMM/${chan}_W_yields.csv"
  echo "       summary    : $SUMM/${chan}_summary.csv"
  echo "       analysis in: $SUMM/${chan}_fitted_yields.root"
  echo "       -> feed analysis macros, e.g.:"
  echo "          charge_asym(\"$SUMM/${chan}_fitted_yields.root\")"
  echo "          FBratio(\"$SUMM/${chan}_fitted_yields.root\")"
}

case "$CHAN_ARG" in
  both) CHANS="mu ele" ;;
  *)    CHANS="$CHAN_ARG" ;;
esac
for c in $CHANS; do run_channel "$c"; done
echo ""
echo "[run_pO_fits] all done. Output under: $OUTROOT"
