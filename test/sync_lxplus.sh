#!/usr/bin/env bash
# =============================================================================
# sync_lxplus.sh -- push the pO Combine inputs + fit scripts to lxplus, and pull
# the fit results back, for the split workflow (build inputs locally -> fit on
# lxplus -> observables locally).  See README_pO_fits.md "Step 4b".
#
# Usage:
#   ./sync_lxplus.sh upload           # inputs + scripts (everything needed to fit)
#   ./sync_lxplus.sh upload-inputs    # only the 4 structured Combine input files
#   ./sync_lxplus.sh upload-scripts   # only the pipeline scripts
#   ./sync_lxplus.sh download         # pull summary/ (fitted yields + CSVs), both chans
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
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
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
  for f in plots/combine_input_W.root plots/combine_input_Z.root \
           plots/Elec/combine_input_W.root plots/Elec/combine_input_Z.root; do
    [ -f "$ANA_LOCAL/plotting/$f" ] || { echo "[MISS] $ANA_LOCAL/plotting/$f"; missing=1; }
  done
  if [ "$missing" -eq 1 ]; then
    echo "[ERROR] missing input(s) -- run Steps 1-3 (skim, run_ngen, qcd_abcd, mtandmet/dileptonpeak) first."
    err=1; return
  fi
  rmkdir "$ANA_LX/plotting"
  ( cd "$ANA_LOCAL/plotting" && run --relative \
      plots/combine_input_W.root      plots/combine_input_Z.root \
      plots/Elec/combine_input_W.root plots/Elec/combine_input_Z.root \
      "$LX:$ANA_LX/plotting/" ) || err=1
}

upload_scripts() {
  echo "== upload pipeline scripts -> $LX:$FORK_LX/test/ =="
  rmkdir "$FORK_LX/test/my_script"
  run "$FORK_LOCAL/test/run_pO_fits.sh" "$LX:$FORK_LX/test/" || err=1
  run "$FORK_LOCAL/test/my_script/make_pO_datacards.sh" \
      "$FORK_LOCAL/test/my_script/extract_pO_yields.C" \
      "$FORK_LOCAL/test/my_script/draw_postfit_pO.C" \
      "$FORK_LOCAL/test/my_script/plotting_helper.C" \
      "$FORK_LOCAL/test/my_script/CMS_lumi.C" \
      "$FORK_LOCAL/test/my_script/CMS_lumi.h" \
      "$LX:$FORK_LX/test/my_script/" || err=1
}

download_results() {
  echo "== download fit results <- $LX:$FORK_LX/test/pO_fit_out/ =="
  local got=0
  for c in $CHANS; do
    local rsum="$FORK_LX/test/pO_fit_out/$c/summary"
    if rexists "$rsum"; then
      mkdir -p "$FORK_LOCAL/test/pO_fit_out/$c/summary"
      run "$LX:$rsum/" "$FORK_LOCAL/test/pO_fit_out/$c/summary/" && got=1 || err=1
    else
      echo "[skip] no remote summary/ for '$c' (fit not run for that channel?)"
    fi
    if [ "$POSTFIT" -eq 1 ]; then
      local rpost="$FORK_LX/test/pO_fit_out/$c/postfit"
      if rexists "$rpost"; then
        mkdir -p "$FORK_LOCAL/test/pO_fit_out/$c/postfit"
        run "$LX:$rpost/" "$FORK_LOCAL/test/pO_fit_out/$c/postfit/" || err=1
      else
        echo "[skip] no remote postfit/ for '$c' (ran with --no-postfit, or plots not made)"
      fi
    fi
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
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both all" ;;
  download)
    echo "[sync_lxplus] download done. Next, locally:"
    echo "    cd $ANA_LOCAL/analysis"
    for c in $CHANS; do
      echo "    root -l -q 'charge_asym.C+(\"$FORK_LOCAL/test/pO_fit_out/$c/summary/${c}_fitted_yields.root\")'"
      echo "    root -l -q 'FBratio.C+(\"$FORK_LOCAL/test/pO_fit_out/$c/summary/${c}_fitted_yields.root\")'"
    done ;;
esac
