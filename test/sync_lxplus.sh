#!/usr/bin/env bash
# =============================================================================
# sync_lxplus.sh -- push the pO Combine inputs + fit scripts to lxplus, and pull
# the fit results back, for the split workflow (build inputs locally -> fit on
# lxplus -> observables locally).  See README_pO_fits.md "Step 4b".
#
# Usage:
#   ./sync_lxplus.sh upload           # inputs + scripts (everything needed to fit)
#   ./sync_lxplus.sh upload-inputs    # only the structured Combine input files
#                                     # (4 required MET/Z + up to 4 optional
#                                     #  lepton-pT variant W files, if built)
#   ./sync_lxplus.sh upload-scripts   # only the pipeline scripts
#   ./sync_lxplus.sh download         # pull summary/ (fitted yields + CSVs), both chans,
#                                     # from ALL discriminant out-trees present
#                                     # (pO_fit_out, pO_fit_out_leppt, pO_fit_out_leppt_mt40)
#   ./sync_lxplus.sh download --postfit   # also pull the postfit plots (skipped if absent)
# Options (any command):
#   --chan mu|ele   restrict to one channel (default: both)
#   --dry-run       show what rsync would do, transfer nothing (-n)
#
# ONE password prompt per run: a single shared SSH connection (ControlMaster) is
# opened up front and reused by every rsync/ssh.  Missing remote dirs are skipped,
# not errored.  Paths default to the user's lxplus layout; override via env, e.g.
#   LX=me@lxplus.cern.ch FORK_LX=/afs/.../HiggsAnalysis-CombinedLimit ./sync_lxplus.sh upload
#
# Tip: to avoid passwords entirely on lxplus, get a Kerberos ticket first:
#   kinit zheng@CERN.CH      (then ssh uses GSSAPI; ControlMaster still 1 connection)
# bash-3.2 safe.
# =============================================================================
set -uo pipefail

# ---- remote (lxplus) paths --------------------------------------------------
LX="${LX:-zheng@lxplus.cern.ch}"
ANA_LX="${ANA_LX:-/afs/cern.ch/user/z/zheng/pO_analysis}"
FORK_LX="${FORK_LX:-/afs/cern.ch/user/z/zheng/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit}"

# ---- local paths (this script lives in <fork>/test) -------------------------
HERE="$(cd "$(dirname "$0")" && pwd)"            # <fork>/test
FORK_LOCAL="$(cd "$HERE/.." && pwd)"             # <fork>
ANA_LOCAL="${ANA_LOCAL:-/Users/zhenghuang/pO_analysis}"

# ---- args -------------------------------------------------------------------
CMD="${1:-}"; [ $# -gt 0 ] && shift
CHANS="mu ele"; DRY=""; POSTFIT=0
while [ $# -gt 0 ]; do
  case "$1" in
    --chan) shift; CHANS="${1:-}";;
    --postfit) POSTFIT=1;;
    --dry-run|-n) DRY="-n";;
    -h|--help) sed -n '2,29p' "$0"; exit 0;;
    *) echo "[ERROR] unknown option: $1"; exit 1;;
  esac
  shift
done

# ---- one shared SSH connection (ControlMaster) ------------------------------
CM="${TMPDIR:-/tmp}/cm_pO_$$"                     # short control-socket path
SSH_OPTS="-o ControlMaster=auto -o ControlPath=$CM -o ControlPersist=300 -o ConnectTimeout=20"
RSH="ssh $SSH_OPTS"
close_master() { ssh $SSH_OPTS -O exit "$LX" 2>/dev/null || true; rm -f "$CM" 2>/dev/null || true; }
open_master() {
  echo "[sync] opening one shared SSH connection to $LX (single auth for the whole run) ..."
  ssh $SSH_OPTS "$LX" true || { echo "[ERROR] cannot connect to $LX"; exit 1; }
  trap close_master EXIT
}
rexists() { ssh $SSH_OPTS "$LX" "test -e '$1'" >/dev/null 2>&1; }   # 0 if remote path exists
rmkdir()  { ssh $SSH_OPTS "$LX" "mkdir -p '$1'"; }
run()     { echo "+ rsync $*"; rsync $DRY -avz -e "$RSH" "$@"; }

err=0

upload_inputs() {
  echo "== upload structured Combine inputs -> $LX:$ANA_LX/plotting/ =="
  local missing=0
  # required: the PF-MET W inputs + the Z inputs
  for f in plots/combine_input_W.root plots/combine_input_Z.root \
           plots/Elec/combine_input_W.root plots/Elec/combine_input_Z.root; do
    [ -f "$ANA_LOCAL/plotting/$f" ] || { echo "[MISS] $ANA_LOCAL/plotting/$f"; missing=1; }
  done
  if [ "$missing" -eq 1 ]; then
    echo "[ERROR] missing input(s) -- run Steps 1-3 (skim, run_ngen, qcd_abcd, mtandmet/dileptonpeak) first."
    err=1; return
  fi
  # optional: the lepton-pT discriminant variants (2026-07-30). Uploaded when
  # present; a missing variant is only a note, not an error.
  SEND="plots/combine_input_W.root plots/combine_input_Z.root plots/Elec/combine_input_W.root plots/Elec/combine_input_Z.root"
  for f in plots/combine_input_W_leppt.root plots/combine_input_W_leppt_mt40.root \
           plots/Elec/combine_input_W_leppt.root plots/Elec/combine_input_W_leppt_mt40.root; do
    if [ -f "$ANA_LOCAL/plotting/$f" ]; then SEND="$SEND $f"
    else echo "[note] optional variant not built (skipped): $f"; fi
  done
  rmkdir "$ANA_LX/plotting"
  ( cd "$ANA_LOCAL/plotting" && run --relative $SEND "$LX:$ANA_LX/plotting/" ) || err=1
}

upload_scripts() {
  echo "== upload pipeline scripts -> $LX:$FORK_LX/test/ =="
  rmkdir "$FORK_LX/test/my_script"
  run "$FORK_LOCAL/test/run_pO_fits.sh" \
      "$FORK_LOCAL/test/run_pO_impacts.sh" \
      "$LX:$FORK_LX/test/" || err=1
  run "$FORK_LOCAL/test/my_script/make_pO_datacards.sh" \
      "$FORK_LOCAL/test/my_script/make_pO_simfit_cards.sh" \
      "$FORK_LOCAL/test/my_script/extract_pO_yields.C" \
      "$FORK_LOCAL/test/my_script/extract_pO_simfit.C" \
      "$FORK_LOCAL/test/my_script/make_yields_from_csv.C" \
      "$FORK_LOCAL/test/my_script/draw_postfit_pO.C" \
      "$FORK_LOCAL/test/my_script/plot_pO_cov.C" \
      "$FORK_LOCAL/test/my_script/plotting_helper.C" \
      "$FORK_LOCAL/test/my_script/CMS_lumi.C" \
      "$FORK_LOCAL/test/my_script/CMS_lumi.h" \
      "$LX:$FORK_LX/test/my_script/" || err=1
}

download_results() {
  # one out-tree per discriminant: "" (PF MET), _leppt, _leppt_mt40
  local got=0
  for sfx in "" "_leppt" "_leppt_mt40"; do
    local tree="pO_fit_out${sfx}"
    echo "== download fit results <- $LX:$FORK_LX/test/$tree/ =="
    for c in $CHANS; do
      local rsum="$FORK_LX/test/$tree/$c/summary"
      if rexists "$rsum"; then
        mkdir -p "$FORK_LOCAL/test/$tree/$c/summary"
        run "$LX:$rsum/" "$FORK_LOCAL/test/$tree/$c/summary/" && got=1 || err=1
      else
        echo "[skip] no remote $tree/$c/summary (fit not run for that channel/discriminant?)"
      fi
      if [ "$POSTFIT" -eq 1 ]; then
        local rpost="$FORK_LX/test/$tree/$c/postfit"
        if rexists "$rpost"; then
          mkdir -p "$FORK_LOCAL/test/$tree/$c/postfit"
          run "$LX:$rpost/" "$FORK_LOCAL/test/$tree/$c/postfit/" || err=1
        else
          echo "[skip] no remote $tree/$c/postfit (ran with --no-postfit, or plots not made)"
        fi
      fi
    done
    # the grand simultaneous fit (simfit, 2026-08-04): flavourless, so outside
    # the per-channel loop; skipped silently when not run for this disc
    local rsimsum="$FORK_LX/test/$tree/simfit/summary"
    if rexists "$rsimsum"; then
      mkdir -p "$FORK_LOCAL/test/$tree/simfit/summary"
      run "$LX:$rsimsum/" "$FORK_LOCAL/test/$tree/simfit/summary/" && got=1 || err=1
    else
      echo "[skip] no remote $tree/simfit/summary (simfit not run for that discriminant?)"
    fi
    if [ "$POSTFIT" -eq 1 ]; then
      local rsimpost="$FORK_LX/test/$tree/simfit/postfit"
      if rexists "$rsimpost"; then
        mkdir -p "$FORK_LOCAL/test/$tree/simfit/postfit"
        run "$LX:$rsimpost/" "$FORK_LOCAL/test/$tree/simfit/postfit/" || err=1
      fi
    fi
    # impacts + covariance plots (run_pO_impacts.sh, 2026-08-17): pulled whenever
    # present -- the json + per-POI PDFs and correlation matrices, but NOT the
    # wd_* intermediate fit files (many higgsCombine*.root, useless locally)
    local d
    for d in impacts cov; do
      local rdir="$FORK_LX/test/$tree/simfit/$d"
      if rexists "$rdir"; then
        mkdir -p "$FORK_LOCAL/test/$tree/simfit/$d"
        run --exclude 'wd_*' "$LX:$rdir/" "$FORK_LOCAL/test/$tree/simfit/$d/" || err=1
      fi
    done
  done
  [ "$got" -eq 0 ] && echo "[warn] nothing downloaded -- did the fit run on lxplus yet?"
}

case "$CMD" in
  upload)         open_master; upload_inputs; upload_scripts ;;
  upload-inputs)  open_master; upload_inputs ;;
  upload-scripts) open_master; upload_scripts ;;
  download)       open_master; download_results ;;
  *) echo "usage: $0 {upload|upload-inputs|upload-scripts|download} [--chan mu|ele] [--postfit] [--dry-run]"; exit 1 ;;
esac

echo ""
if [ "$err" -ne 0 ]; then
  echo "[sync_lxplus] FINISHED WITH ERRORS (see above)."; exit 1
fi
case "$CMD" in
  upload|upload-inputs|upload-scripts)
    echo "[sync_lxplus] upload done. Next, on lxplus:"
    echo "    ssh $LX"
    echo "    cd <CMSSW>/src && cmsenv && cd $FORK_LX/test"
    echo "    # DEFAULT = the grand simultaneous fit (simfit, mu+e in one likelihood):"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh --asimov"
    echo "    # legacy per-flavour pipeline + simfit together:"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both all"
    echo "    # PRIMARY discriminant (2026-08-17): lepton pT with mT>40"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both simfit --asimov --disc leppt_mt40"
    echo "    # backup: PF MET"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both simfit --asimov"
    echo "    # then impacts + covariance plots (pulled by 'download'):"
    echo "    ./run_pO_impacts.sh --disc leppt_mt40" ;;
  download)
    echo "[sync_lxplus] download done. Next, locally (for the lepton-pT variants,"
    echo "swap pO_fit_out for pO_fit_out_leppt or pO_fit_out_leppt_mt40):"
    echo "    cd $ANA_LOCAL/analysis"
    for c in $CHANS; do
      echo "    root -l -q 'charge_asym.C+(\"$FORK_LOCAL/test/pO_fit_out/$c/summary/${c}_fitted_yields.root\")'"
      echo "    root -l -q 'FBratio.C+(\"$FORK_LOCAL/test/pO_fit_out/$c/summary/${c}_fitted_yields.root\")'"
    done ;;
esac
