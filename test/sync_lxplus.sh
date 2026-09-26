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
#   ./sync_lxplus.sh download         # pull, for every out-tree present
#                                     # (pO_fit_out[_leppt[_mt40]]) and for each
#                                     # fit in it -- simfit/ (grand) and, when
#                                     # run, simfit_mu/ + simfit_ele/ (flavfit,
#                                     # 2026-09-22; same layout, same pulls):
#                                     #   summary/   fitted yields + CSVs
#                                     #   datacards/ cards + t2w maps (incl. the
#                                     #              _sigma contour maps) + the
#                                     #              kappa sidecar, AS FITTED
#                                     #   impacts/ cov/ contour/   (minus wd_* and
#                                     #              workspace_sigma.root)
#                                     #   fits/fitDiagnostics_simfit_<B>{,_statonly,_asimov}.root
#                                     #              + the per-variant *.log
#                                     # (summary/ carries extract_simfit.log --
#                                     #  the stat-source line, the statonly-vs-
#                                     #  conditioned cross-check and the Asimov
#                                     #  closure, which live nowhere else)
#   ./sync_lxplus.sh download --postfit   # also pull the postfit plots (skipped if absent)
#   ./sync_lxplus.sh upload-idiso     # the ELECTRON ID+ISO SF study (2026-09-25), a
#                                     # SEPARATE stream: its inputs (the analysis
#                                     # repo's correction/rootfile/idiso_sf_ele/)
#                                     # + run_pO_idiso_sf.sh and its two scripts.
#                                     # upload / download never touch it.
#   ./sync_lxplus.sh download-idiso   # pull test/pO_idiso_sf_out/ (summary CSVs +
#                                     # extraction logs, cards as fitted, combine
#                                     # logs, fitDiagnostics; --postfit adds the
#                                     # postfit plots; never the workspaces)
# Options (any command):
#   --chan mu|ele   restrict to one channel (default: both; also selects which
#                   flavfit tree, simfit_mu/ or simfit_ele/, is downloaded)
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
    # print the WHOLE header block, however long it grows: from the line after
    # the opening '# ===' rule to the closing one. A fixed line range silently
    # truncates the help every time the header is extended -- as it just did.
    -h|--help) sed -n '3,/^# ===/p' "$0" | sed '$d'; exit 0;;
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
  # the LHE shape-systematics sidecars (2026-09-07, <input>_systs.txt): the card
  # generator reads them next to the inputs, so they must travel together
  local sc
  for f in $SEND; do
    sc="${f%.root}_systs.txt"
    [ -f "$ANA_LOCAL/plotting/$sc" ] && SEND="$SEND $sc"
  done
  rmkdir "$ANA_LX/plotting"
  ( cd "$ANA_LOCAL/plotting" && run --relative $SEND "$LX:$ANA_LX/plotting/" ) || err=1
  # the gen FIDUCIAL cross sections (2026-09-15, skim/gen_xsec.C): needed by
  # run_pO_fits.sh --contour, which bakes them into the reparametrized map file
  # that promotes the rapidity-inclusive sigma_W to a POI. Lives outside
  # plotting/, hence its own transfer; absent is only a note (--contour then
  # errors with the path it looked for).
  if [ -f "$ANA_LOCAL/skim/output/gen_xsec_fid.txt" ]; then
    rmkdir "$ANA_LX/skim/output"
    run "$ANA_LOCAL/skim/output/gen_xsec_fid.txt" "$LX:$ANA_LX/skim/output/" || err=1
  else
    echo "[note] no skim/output/gen_xsec_fid.txt (run skim/gen_xsec.C) -- --contour will not work remotely"
  fi
}

upload_scripts() {
  echo "== upload pipeline scripts -> $LX:$FORK_LX/test/ =="
  rmkdir "$FORK_LX/test/my_script"
  run "$FORK_LOCAL/test/run_pO_fits.sh" \
      "$FORK_LOCAL/test/run_pO_impacts.sh" \
      "$LX:$FORK_LX/test/" || err=1
  # (the legacy per-bin scripts make_pO_datacards.sh / extract_pO_yields.C /
  #  make_yields_from_csv.C were removed with that pipeline on 2026-09-22)
  run "$FORK_LOCAL/test/my_script/make_pO_simfit_cards.sh" \
      "$FORK_LOCAL/test/my_script/extract_pO_simfit.C" \
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
    # (the legacy per-channel trees <tree>/{mu,ele}/ are no longer pulled: that
    #  pipeline was removed on 2026-09-22)
    # the grand simultaneous fit (simfit, 2026-08-04), skipped silently when not
    # run for this disc, + the per-flavour simultaneous fits (flavfit,
    # 2026-09-22): simfit_mu/ and simfit_ele/ have EXACTLY the grand fit's
    # layout, so everything below is pulled for each of them too (restricted
    # by --chan)
    local sim flavtrees=""
    for c in $CHANS; do flavtrees="$flavtrees simfit_$c"; done
    for sim in simfit $flavtrees; do
    local rsimsum="$FORK_LX/test/$tree/$sim/summary"
    if rexists "$rsimsum"; then
      mkdir -p "$FORK_LOCAL/test/$tree/$sim/summary"
      run "$LX:$rsimsum/" "$FORK_LOCAL/test/$tree/$sim/summary/" && got=1 || err=1
    else
      echo "[skip] no remote $tree/$sim/summary ($sim not run for that discriminant?)"
    fi
    if [ "$POSTFIT" -eq 1 ]; then
      local rsimpost="$FORK_LX/test/$tree/$sim/postfit"
      if rexists "$rsimpost"; then
        mkdir -p "$FORK_LOCAL/test/$tree/$sim/postfit"
        run "$LX:$rsimpost/" "$FORK_LOCAL/test/$tree/$sim/postfit/" || err=1
      fi
    fi
    # datacards (2026-08-24): the cards + t2w maps + qcd_lnn_kappas.txt sidecar
    # ACTUALLY FITTED on lxplus. Local --dry-run regenerates the local copies
    # with local env defaults (e.g. QCD_MODE=lnN instead of abcd), so without
    # this pull the local cards/sidecar silently disagree with the downloaded
    # summary. Plus impacts + covariance plots (run_pO_impacts.sh, 2026-08-17):
    # pulled whenever present -- the json + per-POI PDFs and correlation
    # matrices, but NOT the wd_* intermediate fit files (many
    # higgsCombine*.root, useless locally)
    # 'contour' (--contour, 2026-09-15) holds the profiled (sigmaW, r_Z) scan;
    # only the higgsCombine*/multidimfit* results are wanted, not the
    # reparametrized workspaces (workspace_sigma.root is ~100 MB and useless
    # locally -- plotting/xsec_contour.C reads only the scan tree).
    local d
    for d in datacards impacts cov contour; do
      local rdir="$FORK_LX/test/$tree/$sim/$d"
      if rexists "$rdir"; then
        mkdir -p "$FORK_LOCAL/test/$tree/$sim/$d"
        run --exclude 'wd_*' --exclude 'workspace_sigma.root' \
            "$LX:$rdir/" "$FORK_LOCAL/test/$tree/$sim/$d/" || err=1
      fi
    done
    # the FitDiagnostics results themselves (2026-09-07): with LHE shape
    # nuisances in the fit the postfit shapes are no longer prefit x scale, so
    # plotting/postfit_incl.C reads shapes_fit_s from these files. Only the two
    # fitDiagnostics_simfit_<B>.root (not the workspaces/logs/higgsCombine*).
    # (+ the _statonly companion -- since 2026-09-15b THE source of the quoted
    # statistical error, run by default -- and the _asimov closure fit, so a
    # local `run_pO_fits.sh --extract-only` on the downloaded tree reproduces
    # the complete summary incl. the stat/syst split and the closure rows)
    # NB `kind` -- NOT `sfx`, which is the OUTER discriminant loop's variable.
    # Re-using it here works only because bash pre-expands a for-loop's word
    # list, and that is far too subtle to rely on.
    # the per-variant combine logs (t2w.log / fit.log / fit_statonly.log /
    # fit_asimov.log): small text, and the only place a failed or badly
    # converged fit explains itself -- the .root files do not.
    local Blog
    for Blog in lab fb; do
      if rexists "$FORK_LX/test/$tree/$sim/fits/simfit_$Blog"; then
        mkdir -p "$FORK_LOCAL/test/$tree/$sim/fits/simfit_$Blog"
        run --include '*.log' --exclude '*' \
            "$LX:$FORK_LX/test/$tree/$sim/fits/simfit_$Blog/" \
            "$FORK_LOCAL/test/$tree/$sim/fits/simfit_$Blog/" || err=1
      fi
    done
    # (inside every fit tree the per-variant names stay simfit_<B> /
    #  fitDiagnostics_simfit_<B>*.root -- the tree dir is the fit's identity)
    local B rfd kind
    for B in lab fb; do
      for kind in "" "_statonly" "_asimov"; do
        rfd="$FORK_LX/test/$tree/$sim/fits/simfit_$B/fitDiagnostics_simfit_${B}${kind}.root"
        if rexists "$rfd"; then
          mkdir -p "$FORK_LOCAL/test/$tree/$sim/fits/simfit_$B"
          run "$LX:$rfd" "$FORK_LOCAL/test/$tree/$sim/fits/simfit_$B/" || err=1
        fi
      done
    done
    done   # sim (the grand fit + the flavfit trees)
  done
  [ "$got" -eq 0 ] && echo "[warn] nothing downloaded -- did the fit run on lxplus yet?"
}

# ---- the ELECTRON ID+ISO SF study (2026-09-25): a SEPARATE stream -------------
# its own inputs, scripts and out-tree (test/pO_idiso_sf_out/); nothing above
# reads or writes any of it, and these two functions touch nothing else.
upload_idiso() {
  local d="$ANA_LOCAL/correction/rootfile/idiso_sf_ele"
  echo "== upload the ID+iso SF inputs -> $LX:$ANA_LX/correction/rootfile/idiso_sf_ele/ =="
  if ! ls "$d"/combine_input_idiso_*.root >/dev/null 2>&1; then
    echo "[ERROR] no $d/combine_input_idiso_*.root -- run correction/run_idiso_sf.sh (skim + inputs) first."
    err=1; return
  fi
  rmkdir "$ANA_LX/correction/rootfile/idiso_sf_ele"
  # the inputs travel with their _meta.txt sidecars (the card generator reads nbins there)
  ( cd "$d" && run combine_input_idiso_*.root combine_input_idiso_*_meta.txt \
                   "$LX:$ANA_LX/correction/rootfile/idiso_sf_ele/" ) || err=1
  echo "== upload the ID+iso SF scripts -> $LX:$FORK_LX/test/ =="
  rmkdir "$FORK_LX/test/my_script"
  run "$FORK_LOCAL/test/run_pO_idiso_sf.sh" "$LX:$FORK_LX/test/" || err=1
  # + the shared postfit plotter and its helpers (unchanged files are skipped by rsync)
  run "$FORK_LOCAL/test/my_script/make_pO_idiso_sf_cards.sh" \
      "$FORK_LOCAL/test/my_script/extract_pO_idiso_sf.C" \
      "$FORK_LOCAL/test/my_script/draw_postfit_pO.C" \
      "$FORK_LOCAL/test/my_script/plotting_helper.C" \
      "$FORK_LOCAL/test/my_script/CMS_lumi.C" \
      "$FORK_LOCAL/test/my_script/CMS_lumi.h" \
      "$LX:$FORK_LX/test/my_script/" || err=1
}

download_idiso() {
  local r="$FORK_LX/test/pO_idiso_sf_out"
  echo "== download the ID+iso SF results <- $LX:$r/ =="
  if ! rexists "$r"; then echo "[warn] no remote $r -- did run_pO_idiso_sf.sh run yet?"; return; fi
  mkdir -p "$FORK_LOCAL/test/pO_idiso_sf_out"
  # everything per tag (summary/, datacards/, the combine logs, fitDiagnostics_idiso_*.root)
  # EXCEPT the workspaces, combine's higgsCombine* files and the input copies (the
  # originals are local). Plain --exclude patterns: macOS's rsync has no '***'.
  if [ "$POSTFIT" -eq 1 ]; then
    run --exclude 'workspace.root' --exclude 'higgsCombine*' --exclude 'combine_input_idiso_*.root' \
        "$LX:$r/" "$FORK_LOCAL/test/pO_idiso_sf_out/" || err=1
  else
    run --exclude 'workspace.root' --exclude 'higgsCombine*' --exclude 'combine_input_idiso_*.root' \
        --exclude 'postfit' "$LX:$r/" "$FORK_LOCAL/test/pO_idiso_sf_out/" || err=1
  fi
}

case "$CMD" in
  upload)         open_master; upload_inputs; upload_scripts ;;
  upload-inputs)  open_master; upload_inputs ;;
  upload-scripts) open_master; upload_scripts ;;
  download)       open_master; download_results ;;
  upload-idiso)   open_master; upload_idiso ;;
  download-idiso) open_master; download_idiso ;;
  *) echo "usage: $0 {upload|upload-inputs|upload-scripts|download|upload-idiso|download-idiso} [--chan mu|ele] [--postfit] [--dry-run]"; exit 1 ;;
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
    echo "    # DEFAULT: simfit (mu+e in one likelihood), --disc leppt_mt40,"
    echo "    # QCD_MODE=abcd, and THREE fits per binning variant --"
    echo "    # nominal + --statonly (the stat error) + --contour (sigma_W as a POI):"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both simfit --asimov"
    echo "    # the mu-only and e-only fits (flavfit: the same model, one flavour each):"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both flavfit --asimov"
    echo "    # or ALL THREE (grand + mu-only + e-only) in one go:"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both all --asimov"
    echo "    # backup discriminant, PF MET (defaults to QCD_MODE=lnN):"
    echo "    PO_PLOTS=$ANA_LX/plotting/plots ./run_pO_fits.sh both simfit --asimov --disc met"
    echo "    # then impacts + covariance plots (pulled by 'download'; --fit mu|ele for flavfit):"
    echo "    ./run_pO_impacts.sh --disc leppt_mt40" ;;
  download)
    echo "[sync_lxplus] download done. Next, locally -- ONE command runs the whole"
    echo "observables chain (charge asymmetry, F/B, sigma = r x sigma_gen, and the"
    echo "(sigma_W, sigma_Z) contour; + the mu-vs-e overlays when the flavfit trees"
    echo "exist) for that discriminant:"
    echo "    cd $ANA_LOCAL/analysis"
    echo "    ./run_observables.sh leppt_mt40      # or: met | all"
    echo "Check in its output that the contour came from the scan, not the ellipse:"
    echo "    [profiled] using .../higgsCombine_contour_lab.MultiDimFit.mH*.root"
    echo "    [profiled] scan min vs covariance best fit: d(sigma_W) = +0.0000 nb" ;;
  upload-idiso)
    echo "[sync_lxplus] ID+iso SF upload done. Next, on lxplus:"
    echo "    ssh $LX"
    echo "    cd <CMSSW>/src && cmsenv && cd $FORK_LX/test"
    echo "    # the 9 W variants (3 binnings x 3 sideband windows; lepton pT at m_T > 40, the nominal"
    echo "    # in-fit ABCD) + the Z tag-and-probe fit 'ztnp', each with its Asimov closure:"
    echo "    ./run_pO_idiso_sf.sh --asimov" ;;
  download-idiso)
    echo "[sync_lxplus] ID+iso SF download done: $FORK_LOCAL/test/pO_idiso_sf_out/<tag>/summary/"
    echo "    idiso_sf_<tag>.csv (per bin; the Z tag-and-probe under <tag> = ztnp),"
    echo "    extract_idiso_sf_<tag>.log (fit quality, Asimov closure, the in-fit ABCD multipliers)" ;;
esac
