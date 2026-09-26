#!/usr/bin/env bash
# =============================================================================
# run_pO_idiso_sf.sh -- the ELECTRON ID+ISO SCALE-FACTOR FIT (2026-09-24).
#
# A study SEPARATE from the nominal fit stream: its own inputs, cards, output
# tree and extraction. Nothing here reads or writes pO_fit_out*/ (which
# analysis/run_observables.sh and plotting/postfit_incl.C consume), and
# run_pO_fits.sh knows nothing about it.
#
# What it measures (the analysis repo's correction/idiso_sf_common.h has the
# full definition): the efficiency of eleMVAIdWP90 && eleMVAIsoWP90 on
# W -> e nu events, eps = N_W(pass) / N_W(pass + fail), with N_W the POST-FIT W
# yields of an electron-only W+Z fit (pass and fail W channels per charge and
# coarse bin + the Z_PP peak, which pins the DY), and SF = eps_data / eps_MC, to
# be compared with the EGM 2025Prompt wp90iso SF.
# + the INCLUSIVE Z TAG-AND-PROBE (2026-09-25), its OWN small fit (tag "ztnp",
# run once): Z_PP (both legs pass, r_ZPP) and Z_PF (the probe fails, r_ZPF) ->
# eps_Z = 2 N_PP / (2 N_PP + N_PF). Not in the W likelihood: a shared DY scale
# let the W channels pull it (SF_Z moved 0.96 -> 1.01 across the W variants).
# The W channels are fitted in the lepton pT at m_T > 40 with the NOMINAL IN-FIT
# ABCD (QCD_MODE=abcd of run_pO_fits.sh; user 2026-09-25: "always the in-fit
# ABCD"): per channel the counting CRs CRB/CRC/CRD in the likelihood, free
# scales sB/sC/sD, the SR QCD scaled by (@0*@1/@2), and the residual lnN on the
# SR QCD only (env IDISO_QCD_ABCD_LNN, default 1.15 = the nominal electron value;
# 'none' drops it).
#
# Per input variant <tag> = <scheme>_<disc>[_sblo|_sbhi]:
#   1. copies <inputs>/combine_input_idiso_<tag>.root (+ _meta.txt) into the work dir
#   2. writes the card + multiSignalModel maps   [my_script/make_pO_idiso_sf_cards.sh]
#   3. text2workspace + combine -M FitDiagnostics (+ the --asimov closure)
#   4. extracts r, eps_data, eps_MC, SF          [my_script/extract_pO_idiso_sf.C]
#   5. postfit plots per channel                 [my_script/draw_postfit_pO.C]
# Output tree: <out>/<tag>/{combine_input_idiso_<tag>.root, datacards/, fits/,
#                            summary/, postfit/}      (default <out> = test/pO_idiso_sf_out)
#
# Usage (under cmsenv; --dry-run needs neither combine nor cmsenv):
#   ./run_pO_idiso_sf.sh [--scheme incl|abseta|pt|ztnp|all] [--disc leppt_mt40]
#                        [--sbwin full|lo|hi|all] [--asimov] [--dry-run]
#                        [--extract-only] [--no-postfit] [--inputs DIR] [--out DIR]
#   defaults: every scheme, every sideband window (9 W fits) + the Z
#             tag-and-probe fit
#   scheme  incl = one bin; abseta = |eta_SC| 0-0.8, 0.8-1.4442, 1.566-2.0, 2.0-2.4;
#           pt = 25-35, 35-50, >50 GeV  (EGM's own edges, merged);
#           ztnp = the Z tag-and-probe (m_ee; --disc/--sbwin do not apply)
#   disc    leppt_mt40 only: the lepton pT at m_T > 40 (the nominal in-fit-ABCD
#           discriminant; MET and m_T cannot be fitted -- m_T is the ABCD axis)
#   sbwin   the relIso sideband of regions C/D (and of the SR QCD shape): full
#           0.3-1.0 (nominal), lo 0.3-0.6, hi 0.6-1.0 -- the window variation
#           measures how the result depends on the choice of sideband
#   --inputs  dir holding combine_input_idiso_*.root; default: env IDISO_INPUTS,
#             else the first existing of the analysis repo's
#             correction/rootfile/idiso_sf_ele (local, lxplus)
#   --extract-only  redo only the extraction on an existing fits/ tree (root only)
# bash-3.2 safe.
# =============================================================================
set -uo pipefail   # NOT -e: one failed fit must not abort the others

HERE="$(cd "$(dirname "$0")" && pwd)"
MYS="$HERE/my_script"

SCHEMES="incl abseta pt ztnp"; DISCS="leppt_mt40"; WINS="full lo hi"
ASIMOV=0; DRYRUN=0; EXTRACTONLY=0; DO_POSTFIT=1
OUTROOT="$HERE/pO_idiso_sf_out"
INDIR="${IDISO_INPUTS:-}"
INDIR_DEFAULTS="/Users/zhenghuang/pO_analysis/correction/rootfile/idiso_sf_ele /afs/cern.ch/user/z/zheng/pO_analysis/correction/rootfile/idiso_sf_ele"

while [ $# -gt 0 ]; do
  case "$1" in
    --scheme) shift; case "${1:-}" in all) SCHEMES="incl abseta pt ztnp" ;; incl|abseta|pt|ztnp) SCHEMES="$1" ;;
                                      *) echo "[ERROR] --scheme incl|abseta|pt|ztnp|all"; exit 1 ;; esac ;;
    --disc)   shift; case "${1:-}" in leppt_mt40|all) DISCS="leppt_mt40" ;;
                                      *) echo "[ERROR] --disc leppt_mt40 (met/mt: m_T is the in-fit ABCD axis)"; exit 1 ;; esac ;;
    --sbwin)  shift; case "${1:-}" in all) WINS="full lo hi" ;; full|lo|hi) WINS="$1" ;;
                                      *) echo "[ERROR] --sbwin full|lo|hi|all"; exit 1 ;; esac ;;
    --asimov)       ASIMOV=1 ;;
    --dry-run)      DRYRUN=1 ;;
    --extract-only) EXTRACTONLY=1 ;;
    --no-postfit)   DO_POSTFIT=0 ;;
    --inputs) shift; INDIR="${1:?--inputs needs a directory}" ;;
    --out)    shift; OUTROOT="${1:?--out needs a directory}" ;;
    -h|--help) sed -n '3,/^# ===/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "[ERROR] unknown arg: $1"; exit 1 ;;
  esac
  shift
done

if [ -z "$INDIR" ]; then
  for d in $INDIR_DEFAULTS; do
    if [ -d "$d" ] && ls "$d"/combine_input_idiso_*.root >/dev/null 2>&1; then INDIR="$d"; break; fi
  done
fi
if [ "$EXTRACTONLY" -eq 0 ] && { [ -z "$INDIR" ] || [ ! -d "$INDIR" ]; }; then
  echo "[ERROR] no inputs dir with combine_input_idiso_*.root -- run the analysis repo's"
  echo "        correction/run_idiso_sf.sh (skim + inputs), sync them, or pass --inputs DIR."
  exit 2
fi
echo "[run_pO_idiso_sf] inputs: ${INDIR:-<work-dir copies>}   out: $OUTROOT"
echo "[run_pO_idiso_sf] schemes: $SCHEMES   discs: $DISCS   sideband windows: $WINS   asimov: $ASIMOV   dry-run: $DRYRUN"
echo "[run_pO_idiso_sf] QCD: in-fit ABCD, residual lnN IDISO_QCD_ABCD_LNN=${IDISO_QCD_ABCD_LNN:-1.15}"
KAPX="${IDISO_QCD_ABCD_LNN:-1.15}"; [ "$KAPX" = "none" ] && KAPX="1.0"   # the extractor's kappa (1 = no residual)

HAVE_COMBINE=1
command -v combine          >/dev/null 2>&1 || HAVE_COMBINE=0
command -v text2workspace.py >/dev/null 2>&1 || HAVE_COMBINE=0
if [ "$DRYRUN" -eq 0 ] && [ "$EXTRACTONLY" -eq 0 ] && [ "$HAVE_COMBINE" -eq 0 ]; then
  echo "[ERROR] combine / text2workspace.py not on PATH -- did you cmsenv?  (--dry-run only builds the cards)"
  exit 3
fi
HAVE_ROOT=1; command -v root >/dev/null 2>&1 || HAVE_ROOT=0
if [ "$EXTRACTONLY" -eq 1 ] && [ "$HAVE_ROOT" -eq 0 ]; then echo "[ERROR] --extract-only needs root"; exit 3; fi

# the fit tags: <scheme>_<disc>[_sblo|_sbhi] per W variant, and the one "ztnp"
TAGS=""
for S in $SCHEMES; do
  if [ "$S" = "ztnp" ]; then TAGS="$TAGS ztnp"; continue; fi
  for D in $DISCS; do
    for WN in $WINS; do
      case "$WN" in full) WS="" ;; lo) WS="_sblo" ;; hi) WS="_sbhi" ;; esac
      TAGS="$TAGS ${S}_${D}${WS}"
    done
  done
done

nfail=0
for TAG in $TAGS; do
      WORK="$OUTROOT/$TAG"; DCD="$WORK/datacards"; FITS="$WORK/fits"; SUMM="$WORK/summary"; POST="$WORK/postfit"
      IN="$WORK/combine_input_idiso_${TAG}.root"; META="${IN%.root}_meta.txt"
      echo ""
      echo "================ ID+iso SF fit: $TAG ================"
      mkdir -p "$DCD" "$FITS" "$SUMM" "$POST"

      if [ "$EXTRACTONLY" -eq 0 ]; then
        SRC="$INDIR/combine_input_idiso_${TAG}.root"
        if [ ! -f "$SRC" ] || [ ! -f "${SRC%.root}_meta.txt" ]; then
          echo "  [skip] no input $SRC (+ _meta.txt)"; nfail=$((nfail+1)); continue
        fi
        cp -f "$SRC" "$IN"; cp -f "${SRC%.root}_meta.txt" "$META"
        /bin/bash "$MYS/make_pO_idiso_sf_cards.sh" "$IN" "$DCD" "$TAG" || { echo "  [FAIL cards] $TAG"; nfail=$((nfail+1)); continue; }
        [ "$DRYRUN" -eq 1 ] && continue

        CARD="$DCD/datacard_idiso_${TAG}.txt"; MAPS="$DCD/t2w_maps_idiso_${TAG}.txt"
        PO=(-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose)
        while IFS= read -r m; do [ -n "$m" ] && PO+=(--PO "$m"); done < "$MAPS"
        (
          cd "$FITS" || exit 1
          text2workspace.py "$CARD" -o workspace.root "${PO[@]}" >t2w.log 2>&1 \
            || { echo "  [FAIL t2w] $TAG (see $FITS/t2w.log)"; exit 1; }
          # --skipBOnlyFit: every POI at 0 would switch the W off -- meaningless here
          combine -M FitDiagnostics workspace.root --saveShapes --saveWithUncertainties --skipBOnlyFit \
                  -n "_idiso_${TAG}" --cminDefaultMinimizerStrategy 0 >fit.log 2>&1 \
            || { echo "  [FAIL fit] $TAG (see $FITS/fit.log)"; exit 1; }
          if [ "$ASIMOV" -eq 1 ]; then
            # prefit S+B Asimov closure: every POI must return 1, i.e. SF = 1 exactly.
            # All POIs are set explicitly -- a plain -t -1 generates the B-ONLY Asimov
            # (all POIs at 0), see the note in run_pO_fits.sh.
            SETPARS=$(sed -n 's/^map=.*:\(r_[A-Za-z0-9_]*\)\[.*/\1=1/p' "$MAPS" | sort -u | paste -sd, -)
            combine -M FitDiagnostics workspace.root --skipBOnlyFit -t -1 --setParameters "$SETPARS" \
                    -n "_idiso_${TAG}_asimov" --cminDefaultMinimizerStrategy 0 >fit_asimov.log 2>&1 \
              || { echo "  [FAIL asimov] $TAG (see $FITS/fit_asimov.log)"; exit 1; }
          fi
        ) && echo "  [ok] fit $TAG" || { nfail=$((nfail+1)); continue; }
      fi

      FD="$FITS/fitDiagnostics_idiso_${TAG}.root"; FDA="$FITS/fitDiagnostics_idiso_${TAG}_asimov.root"
      [ -f "$FDA" ] || FDA="none"
      if [ "$HAVE_ROOT" -eq 1 ] && [ -f "$FD" ]; then
        # the full console output IS the record (fit quality, closure, every WARN); summary/ is downloaded
        root -b -q "$MYS/extract_pO_idiso_sf.C(\"$FD\",\"$FDA\",\"$IN\",\"$META\",\"$SUMM\",\"$TAG\",$KAPX)" 2>&1 \
          | tee "$SUMM/extract_idiso_sf_${TAG}.log" | grep -E "\[idiso-sf\]|\[asimov\]|\[expected\]|WARN|FAIL" || true
      elif [ ! -f "$FD" ]; then
        echo "  [skip extract] no $FD"
      fi

      if [ "$DO_POSTFIT" -eq 1 ] && [ "$HAVE_ROOT" -eq 1 ] && [ -f "$FD" ] && [ "$EXTRACTONLY" -eq 0 ]; then
        XT=$(sed -n 's/^xtitle //p' "$META"); YT=$(sed -n 's/^ytitle //p' "$META")
        NB=$(awk '$1=="nbins"{print $2}' "$META")
        CARD="$DCD/datacard_idiso_${TAG}.txt"
        if [ "$TAG" = "ztnp" ]; then
          ZCATS="PP PF"
        else
          ZCATS="PP"
          for C in Wp Wm; do
            k=0
            while [ "$k" -lt "$NB" ]; do
              for CAT in pass fail; do
                R="${C}_${CAT}_k${k}"; CH="ele_${R}"
                # the QCD label of the info box: the free rateParam only exists in free mode
                QN="none"; grep -q "^qcd_norm_${CH} " "$CARD" 2>/dev/null && QN="qcd_norm_${CH}"
                root -b -q "$MYS/draw_postfit_pO.C(\"$FD\",\"$CH\",\"$IN\",\"$R\",\"$POST/$CH\",\"$XT\",\"$YT\",\"W #rightarrow e#nu, ID+iso SF fit\",\"$CH ($TAG postfit)\",true,\"r_${R}\",\"r_Z\",\"$QN\",3)" >/dev/null 2>&1
              done
              k=$((k+1))
            done
          done
        fi
        for ZC in $ZCATS; do
          # the DY scale: r_Z (the W fits' anchor) or the tag-and-probe's own r_ZPP / r_ZPF
          if [ "$TAG" = "ztnp" ]; then ZP="r_Z${ZC}"; else ZP="r_Z"; fi
          [ "$ZC" = "PP" ] && ZT="Z #rightarrow ee TnP, both pass" || ZT="Z #rightarrow ee TnP, probe fails"
          root -b -q "$MYS/draw_postfit_pO.C(\"$FD\",\"ele_Z_${ZC}\",\"$IN\",\"Z_${ZC}\",\"$POST/ele_Z_${ZC}\",\"m_{ee} (GeV)\",\"Events / 1.0 GeV\",\"$ZT\",\"Z ${ZC} ($TAG postfit)\",false,\"$ZP\",\"none\",\"bkg_norm_ele_Z_${ZC}\",2)" >/dev/null 2>&1
        done
        echo "  [ok] postfit plots in $POST"
      fi
done

echo ""
echo "[run_pO_idiso_sf] done, failures=$nfail   (results: $OUTROOT/<tag>/summary/idiso_sf_<tag>.csv)"
[ "$nfail" -eq 0 ]
