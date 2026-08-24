#!/usr/bin/env bash
# =============================================================================
# run_pO_fits.sh -- end-to-end pO W/Z Combine fits, automated.
#
# For each lepton channel (muon / electron) it:
#   1. locates the structured Combine inputs from the analysis repo
#        <plots>/combine_input_W.root , <plots>/combine_input_Z.root      (muon)
#        <plots>/Elec/combine_input_{W,Z}.root                            (electron)
#   2. generates ALL datacards (per-(charge,y) lab + FB -- each a TWO-channel
#      card fitted simultaneously with Z_incl -- per-charge incl, W_incl,
#      Z_incl, and the simultaneous W+Z card).  Two-parameter model: POI 'r'
#      scales all W-related MC, 'dy_norm' all DY-related MC (shared with the
#      Z peak in the simultaneous cards), 'qcd_norm' the data-driven QCD.
#      [make_pO_datacards.sh]
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
#   ./run_pO_fits.sh [mu|ele|both] [perbin|incl|combined|simfit|all] [options]
#     channel  (default both)   mode (default simfit)
#
#   Mode 'simfit' (2026-08-04, the DEFAULT) = the GRAND SIMULTANEOUS FIT: all 48 (flavour,
#   charge, y-bin) W channels + BOTH Z peaks in ONE likelihood per binning
#   variant (lab, fb).  25 POIs: r_<C>_y<i> (24, mu/e SHARED) + one global r_Z
#   scaling all DY-related MC; QCD lnN-constrained at the ABCD prediction per
#   (flavour, charge) + global lumi lnN on all MC (2026-08-17 default; env
#   QCD_MODE/QCD_LNN_MU/QCD_LNN_ELE/LUMI_LNN -> make_pO_simfit_cards.sh;
#   QCD_MODE=free restores the 48 free qcd_norm rateParams; QCD_MODE=abcd
#   (2026-08-23, --disc leppt_mt40 ONLY) = the IN-FIT ABCD: 12 counting CR
#   channels + free scales qcd_s{B,C,D}_<F>_<C> + the formula rateParam
#   (sB*sC/sD) on the SR qcd_abcd template, so the QCD normalization floats
#   with the CR data and the EWK subtraction rides the POIs; env
#   QCD_ABCD_LNN_MU/QCD_ABCD_LNN_ELE = the reduced residual kappas,
#   QCD_WCR=float|frozen picks the CRB W treatment); w/wtau under the Z
#   peaks frozen at absolute MC.  Cross-flavour by construction, so the channel
#   argument is ignored (mode 'all' runs simfit only when channel = both).
#   Outputs under <out>/simfit/ (comb_* files; yields are mu+e combined).
#   The legacy per-bin pipeline (perbin/incl/combined) is UNCHANGED and stays
#   runnable for comparison (it refits the same Z data in every per-bin card);
#   mode 'all' = the full legacy per-flavour pipeline PLUS simfit.
#   options:
#     --disc met|leppt|leppt_mt40
#                       W discriminant (default met = PF MET shape).
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
#     --asimov          (simfit only) also run a prefit-Asimov closure fit per
#                       variant (-t -1): every fitted POI must come back at 1
#     --plots-dir DIR   analysis plots dir (else $PO_PLOTS, else autodetect)
#     --out DIR         output root (default test/pO_fit_out)
#
# Needs `cmsenv` (combine + text2workspace.py on PATH) for the actual fits.
# bash-3.2 safe (macOS stock bash): no associative arrays.
# =============================================================================
set -uo pipefail   # NOT -e: per-bin fit failures must not abort the whole loop

HERE="$(cd "$(dirname "$0")" && pwd)"
MYS="$HERE/my_script"

CHAN_ARG="both"; MODE="simfit"; DRYRUN=0; DO_POSTFIT=1; DRAWONLY=0; ASIMOV=0
DISC="leppt_mt40"; OUT_SET=0
OUTROOT="$HERE/pO_fit_out"
PO_PLOTS="${PO_PLOTS:-}"
PO_PLOTS_DEFAULTS="/Users/zhenghuang/pO_analysis/plotting/plots /afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots"

while [ $# -gt 0 ]; do
  case "$1" in
    mu|ele|both)                      CHAN_ARG="$1" ;;
    perbin|incl|combined|simfit|all)  MODE="$1" ;;
    --dry-run)                 DRYRUN=1 ;;
    --no-postfit)              DO_POSTFIT=0 ;;
    --draw-only)               DRAWONLY=1 ;;
    --asimov)                  ASIMOV=1 ;;
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
# legacy output path for met; suffixed tree for the variants (unless --out given)
if [ "$OUT_SET" -eq 0 ]; then OUTROOT="$HERE/pO_fit_out${DSUF}"; fi

# ---- locate the analysis plots dir (not needed for --draw-only: it reuses ---
# ---- the combine_input_*.root copies already in the work dir) ---------------
if [ "$DRAWONLY" -eq 0 ]; then
  if [ -z "$PO_PLOTS" ]; then
    for d in $PO_PLOTS_DEFAULTS; do
      if [ -f "$d/$WINNAME" ]; then PO_PLOTS="$d"; break; fi
    done
  fi
  if [ -z "$PO_PLOTS" ] || [ ! -d "$PO_PLOTS" ]; then
    echo "[ERROR] analysis plots dir not found. Set --plots-dir or \$PO_PLOTS to the"
    echo "        dir containing $WINNAME (run plotting/mtandmet.C +"
    echo "        dileptonpeak.C first)."
    exit 2
  fi
  echo "[run_pO_fits] plots dir : $PO_PLOTS"
fi
echo "[run_pO_fits] channel(s): $CHAN_ARG    mode: $MODE    disc: $DISC    dry-run: $DRYRUN    draw-only: $DRAWONLY    asimov: $ASIMOV"

# ---- cmsenv check -----------------------------------------------------------
HAVE_COMBINE=1
command -v combine          >/dev/null 2>&1 || HAVE_COMBINE=0
command -v text2workspace.py >/dev/null 2>&1 || HAVE_COMBINE=0
if [ "$DRYRUN" -eq 0 ] && [ "$DRAWONLY" -eq 0 ] && [ "$HAVE_COMBINE" -eq 0 ]; then
  echo "[ERROR] combine / text2workspace.py not on PATH -- did you cmsenv?"
  echo "        (run with --dry-run to only build datacards.)"
  exit 3
fi
command -v root >/dev/null 2>&1 || { echo "[WARN] root not on PATH: extraction/postfit will be skipped."; }
if [ "$DRAWONLY" -eq 1 ] && ! command -v root >/dev/null 2>&1; then
  echo "[ERROR] --draw-only needs root on PATH."; exit 3
fi

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
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"Wincl\",\"$AW\",\"W_incl\",\"$POST/WZ_Wincl\",\"$XT\",\"$YT\",\"$WL\",\"W+Z fit (postfit)\",true)"  >/dev/null 2>&1
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"Zincl\",\"$AZ\",\"Z_incl\",\"$POST/WZ_Zincl\",\"m_{ll} (GeV)\",\"Events / 1.0 GeV\",\"$ZL\",\"W+Z fit (postfit)\",false)" >/dev/null 2>&1 ;;
    *)
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"$R\",\"$AW\",\"$R\",\"$POST/$R\",\"$XT\",\"$YT\",\"$WL\",\"$R (postfit)\",true)" >/dev/null 2>&1 ;;
  esac
}

# ---- one channel ------------------------------------------------------------
run_channel() {
  chan="$1"
  if [ "$chan" = "ele" ]; then
    WIN_SRC="$PO_PLOTS/Elec/$WINNAME"; ZIN_SRC="$PO_PLOTS/Elec/combine_input_Z.root"
    WL="W #rightarrow e #nu"; ZL="Z #rightarrow e e"
  else
    WIN_SRC="$PO_PLOTS/$WINNAME"; ZIN_SRC="$PO_PLOTS/combine_input_Z.root"
    WL="W #rightarrow #mu #nu"; ZL="Z #rightarrow #mu #mu"
  fi
  echo ""
  echo "================ channel: $chan ================"
  WORK="$OUTROOT/$chan"; DCD="$WORK/datacards"; FITS="$WORK/fits"; POST="$WORK/postfit"; SUMM="$WORK/summary"
  ABS_W="$WORK/combine_input_W.root"; ABS_Z="$WORK/combine_input_Z.root"

  # ---- draw-only: redraw postfit plots from an EXISTING fit run --------------
  # Reuses the work-dir input copies (physical axes) + fits/<region>/ from the
  # previous run; regenerates nothing else.
  if [ "$DRAWONLY" -eq 1 ]; then
    if [ ! -f "$ABS_W" ]; then
      echo "[ERROR] $ABS_W missing -- no previous fit run for '$chan' (run the full pipeline first)."
      return
    fi
    nfd=$(find "$FITS" -name 'fitDiagnostics_*.root' 2>/dev/null | wc -l | tr -d ' ')
    if [ "$nfd" -eq 0 ]; then
      echo "[ERROR] no fitDiagnostics_*.root under $FITS -- --draw-only needs the fits/"
      echo "        tree from an earlier run (NB 'sync_lxplus.sh download' does NOT pull fits/;"
      echo "        redraw where the fits ran, then download --postfit)."
      return
    fi
    mkdir -p "$POST"
    REGIONS="$(build_regions "$MODE")"
    echo "[draw-only] $nfd fit result(s) under fits/; redrawing $(echo $REGIONS | wc -w | tr -d ' ') region(s) + combined: $(want_combined "$MODE" && echo yes || echo no) ..."
    for R in $REGIONS; do postfit_region "$R" "$FITS" "$POST" "$ABS_W" "$ABS_Z" "$WL" "$ZL"; done
    if want_combined "$MODE" && [ -f "$ABS_Z" ]; then postfit_region "WZ" "$FITS" "$POST" "$ABS_W" "$ABS_Z" "$WL" "$ZL"; fi
    echo "[done] $chan -> $POST"
    return
  fi

  if [ ! -f "$WIN_SRC" ]; then echo "[ERROR] missing $WIN_SRC"; return; fi
  if [ ! -f "$ZIN_SRC" ]; then echo "[WARN] missing $ZIN_SRC (Z + combined fits will be skipped)"; fi

  mkdir -p "$WORK" "$DCD" "$FITS" "$POST" "$SUMM"
  cp -f "$WIN_SRC" "$WORK/combine_input_W.root"
  [ -f "$ZIN_SRC" ] && cp -f "$ZIN_SRC" "$WORK/combine_input_Z.root"

  # absolute-path datacards so combine resolves shapes from any CWD
  /bin/bash "$MYS/make_pO_datacards.sh" "$ABS_W" "$ABS_Z" "$DCD" "$DISC"

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

# =============================================================================
# simfit -- the GRAND SIMULTANEOUS FIT (2026-08-04).  ONE likelihood per
# binning variant (lab / fb): 48 W channels ({mu,ele} x {Wp,Wm} x y0..11) +
# both Z-inclusive peaks.  25 POIs (r_<C>_y<i> shared mu/e + global r_Z),
# QCD lnN per (flavour, charge) + global lumi lnN (2026-08-17; QCD_MODE=free
# restores per-channel free qcd_norm), w/wtau under Z frozen -- see
# my_script/make_pO_simfit_cards.sh for the model definition (card + t2w maps).
# Cross-flavour by construction, so it lives OUTSIDE run_channel().
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

simfit_postfit_all() {  # postfit data/MC per channel of the grand fit, both variants
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
    for F in mu ele; do
      if [ "$F" = "ele" ]; then AW="$AWEL"; AZ="$AZEL"; WL="W #rightarrow e #nu"; ZL="Z #rightarrow e e"
      else                      AW="$AWMU"; AZ="$AZMU"; WL="W #rightarrow #mu #nu"; ZL="Z #rightarrow #mu #mu"; fi
      for C in Wp Wm; do
        for iy in $(seq 0 11); do
          R="${C}_${B}_y${iy}"; CH="${F}_${R}"
          QN="qcd_norm_${CH}"; [ "$QLNN" -eq 1 ] && QN="qcd_rate_${F}_${C}"
          # info box: this bin's POI, the global r_Z (shown as DY norm), the QCD
          # param (lnN mode: the shared nuisance -> value shown is the PULL);
          # ndf uses the ~3 params that shape this channel.
          root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"$CH\",\"$AW\",\"$R\",\"$SPOST/$CH\",\"$XT\",\"$YT\",\"$WL\",\"$CH (simfit postfit)\",true,\"r_${C}_y${iy}\",\"r_Z\",\"$QN\",3)" >/dev/null 2>&1
        done
      done
      # the Z peak as seen by this variant's grand fit (r_Z only; W bkg frozen)
      root -b -q "$MYS/draw_postfit_pO.C(\"$fd\",\"${F}_Z_incl\",\"$AZ\",\"Z_incl\",\"$SPOST/${F}_Z_incl_${B}\",\"m_{ll} (GeV)\",\"Events / 1.0 GeV\",\"$ZL\",\"Z incl (simfit ${B} postfit)\",false,\"r_Z\",\"none\",\"none\",1)" >/dev/null 2>&1
    done
  done
}

run_simfit() {
  echo ""
  echo "================ simfit: grand simultaneous fit (mu + ele) ================"
  SWORK="$OUTROOT/simfit"; SDCD="$SWORK/datacards"; SFITS="$SWORK/fits"; SPOST="$SWORK/postfit"; SSUMM="$SWORK/summary"
  AWMU="$SWORK/combine_input_W_mu.root";  AZMU="$SWORK/combine_input_Z_mu.root"
  AWEL="$SWORK/combine_input_W_ele.root"; AZEL="$SWORK/combine_input_Z_ele.root"

  # ---- draw-only: redraw simfit postfit plots from an EXISTING run -----------
  if [ "$DRAWONLY" -eq 1 ]; then
    if [ ! -f "$AWMU" ] || [ ! -f "$AWEL" ]; then
      echo "[ERROR] $SWORK input copies missing -- no previous simfit run (run the full simfit first)."
      return
    fi
    nfd=$(find "$SFITS" -name 'fitDiagnostics_simfit_*.root' 2>/dev/null | wc -l | tr -d ' ')
    if [ "$nfd" -eq 0 ]; then
      echo "[ERROR] no fitDiagnostics_simfit_*.root under $SFITS -- --draw-only needs an earlier"
      echo "        simfit run (NB 'sync_lxplus.sh download' does NOT pull fits/; redraw where the fits ran)."
      return
    fi
    mkdir -p "$SPOST"
    echo "[draw-only] redrawing simfit postfit plots ..."
    simfit_postfit_all
    echo "[done] simfit -> $SPOST"
    return
  fi

  WMU_SRC="$PO_PLOTS/$WINNAME";      ZMU_SRC="$PO_PLOTS/combine_input_Z.root"
  WEL_SRC="$PO_PLOTS/Elec/$WINNAME"; ZEL_SRC="$PO_PLOTS/Elec/combine_input_Z.root"
  miss=0
  for f in "$WMU_SRC" "$ZMU_SRC" "$WEL_SRC" "$ZEL_SRC"; do
    [ -f "$f" ] || { echo "[ERROR] simfit input missing: $f"; miss=1; }
  done
  if [ "$miss" -eq 1 ]; then
    echo "[ERROR] simfit needs BOTH flavours' W AND Z inputs (no Z fallback: r_Z is pinned by the peaks) -- skipped."
    return
  fi

  mkdir -p "$SWORK" "$SDCD" "$SFITS" "$SPOST" "$SSUMM"
  cp -f "$WMU_SRC" "$AWMU"; cp -f "$ZMU_SRC" "$AZMU"
  cp -f "$WEL_SRC" "$AWEL"; cp -f "$ZEL_SRC" "$AZEL"

  # absolute-path datacards so combine resolves shapes from any CWD.
  # The generator can refuse (e.g. QCD_MODE=abcd with a non-leppt_mt40 disc);
  # the script runs without -e, so check explicitly or the fit stage would run
  # on stale/absent cards.
  if ! /bin/bash "$MYS/make_pO_simfit_cards.sh" "$AWMU" "$AZMU" "$AWEL" "$AZEL" "$SDCD" "$DISC"; then
    echo "[ERROR] simfit datacard generation failed -- simfit skipped."
    return
  fi

  if [ "$DRYRUN" -eq 1 ]; then
    echo "[dry-run] simfit datacards + t2w maps in $SDCD ; skipping fits."
    return
  fi

  for B in lab fb; do fit_simfit "$B"; done

  # ---- extract POIs + mu+e-combined yields + covariance ----
  # lnN kappas + qcd mode the cards were built with (sidecar from
  # make_pO_simfit_cards.sh; missing sidecar / 0 entries -> legacy
  # free-rateParam extraction path; missing qcdMode line -> legacy sidecar,
  # mode inferred from the kappas inside the extractor)
  KQM=0; KQE=0; KLU=0; QMODE=""; KF="$SDCD/qcd_lnn_kappas.txt"
  if [ -f "$KF" ]; then
    KQM=$(awk '$1=="kQcdMu"{print $2}' "$KF");  KQM="${KQM:-0}"
    KQE=$(awk '$1=="kQcdEle"{print $2}' "$KF"); KQE="${KQE:-0}"
    KLU=$(awk '$1=="kLumi"{print $2}' "$KF");   KLU="${KLU:-0}"
    QMODE=$(awk '$1=="qcdMode"{print $2}' "$KF"); QMODE="${QMODE:-}"
  fi
  if command -v root >/dev/null 2>&1; then
    root -b -q "$MYS/extract_pO_simfit.C(\"$SFITS\",\"$AWMU\",\"$AWEL\",\"$SSUMM\",$KQM,$KQE,$KLU,\"$QMODE\")" 2>&1 \
      | grep -E "\[extract-simfit\]|\[asimov\]|WARN|FAIL" || true
  fi

  # ---- postfit plots (per channel of the grand fit: 2 variants x 50) ----
  if [ "$DO_POSTFIT" -eq 1 ] && command -v root >/dev/null 2>&1; then
    echo "[postfit] drawing simfit postfit plots (2 variants x 50 channels) ..."
    simfit_postfit_all
  fi

  echo "[done] simfit -> $SWORK"
  echo "       combined yields CSV : $SSUMM/comb_W_yields.csv"
  echo "       POI summary         : $SSUMM/comb_summary.csv"
  echo "       analysis input      : $SSUMM/comb_fitted_yields.root  (+ h_cov_yield[_FB])"
  echo "       -> feed analysis macros, e.g.:"
  echo "          charge_asym(\"$SSUMM/comb_fitted_yields.root\")"
  echo "          FBratio(\"$SSUMM/comb_fitted_yields.root\")"
}

case "$CHAN_ARG" in
  both) CHANS="mu ele" ;;
  *)    CHANS="$CHAN_ARG" ;;
esac
if [ "$MODE" != "simfit" ]; then
  for c in $CHANS; do run_channel "$c"; done
fi
if [ "$MODE" = "simfit" ]; then
  [ "$CHAN_ARG" != "both" ] && echo "[note] simfit always uses BOTH flavours; channel arg '$CHAN_ARG' ignored."
  run_simfit
elif [ "$MODE" = "all" ]; then
  if [ "$CHAN_ARG" = "both" ]; then run_simfit
  else echo "[note] mode 'all' with channel '$CHAN_ARG': simfit needs both flavours -> skipped."; fi
fi
echo ""
echo "[run_pO_fits] all done. Output under: $OUTROOT"
