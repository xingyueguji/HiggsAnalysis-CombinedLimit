#!/usr/bin/env bash
# =============================================================================
# make_pO_simfit_cards.sh -- generate the GRAND SIMULTANEOUS FIT datacards,
# one per binning variant (lab / fb): ALL (flavour, charge, y-bin) W channels
# plus BOTH Z-inclusive peaks in ONE likelihood (2026-08-04 scheme).
#
# Channels per card (50):
#   {mu,ele} x {Wp,Wm} x y0..y11   ->  48 W channels, named <F>_<C>_<B>_y<i>
#   mu_Z_incl, ele_Z_incl          ->   2 Z channels
#
# Fit model: 2N+1 = 25 POIs (N = 12 rapidity bins), applied at text2workspace
# time via multiSignalModel -- the map file written NEXT TO each card
# (t2w_maps_simfit_<B>.txt, one 'map=' per line) is the single definition:
#   r_<C>_y<i>  scales the W-related MC (signal + wtau) of rapidity bin i,
#               charge C, in BOTH flavours' channels  -> 24 POIs, mu/e SHARED
#   r_Z         scales ALL DY-related MC everywhere (z + ztau in every W
#               channel, zsig + ztau under both Z peaks) -> the "+1".  The DY
#               rapidity dependence across W bins is taken from MC (trusted);
#               only this one global normalization floats, pinned by the peaks.
#   QCD (2026-08-17, default): log-normal-constrained at the ABCD prediction --
#   ONE lnN nuisance per (flavour, charge), qcd_rate_{mu,ele}_{Wp,Wm},
#   correlated across the 12 y bins of that flavour+charge (T and the template
#   are measured once, inclusively in y; uncorrelated across charge/flavour so
#   the charge asymmetry gets no unearned cancellation).  QCD_MODE=free restores
#   the pre-2026-08-17 model: 48 free qcd_norm_<channel> rateParams.
#   QCD_MODE=abcd (2026-08-23, leppt_mt40 ONLY): the IN-FIT ABCD -- 12 extra
#   counting CR channels <F>_<C>_{CRB,CRC,CRD} ((relIso x m_T) plane: CRB
#   iso-pass m_T<30, CRC anti-iso m_T>40 = the SR template's source region,
#   CRD anti-iso m_T<30; 1-bin templates written by plotting/mtandmet.C).
#   Free scales qcd_s{B,C,D}_<F>_<C> float the three QCD counts and the SR qcd
#   (shapes path qcd_abcd, total = B0*C40/D0) is scaled by the FORMULA
#   rateParam (@0*@1/@2) -- Combine's documented ABCD pattern
#   (docs/part2/settinguptheanalysis.md) -- so the normalization floats with
#   the CR data and the EWK subtraction rides the POIs: CRB z/ztau are swept
#   into r_Z by the catch-all map, the CRB W content is 12 per-y processes
#   w_y0..11 mapped to r_<C>_y<i> (QCD_WCR=float, default) or one frozen wfix
#   (QCD_WCR=frozen).  The qcd_rate_* lnN rows stay, with the REDUCED kappas
#   (window (+) FF-shift only -- stat is profiled in-fit, plane transport does
#   not apply; derived in correction/logs/qcd_abcd_*.log, 2026-08-23).
#   lumi lnN (2026-08-17): one global nuisance on EVERY MC template (signal, z,
#   ztau, wtau, zsig, and the CR wfix/w_y*/ewk -- NOT the data-driven qcd),
#   +-3% on L = 46.5 nb^-1.
#   w, wtau under the Z peaks: FROZEN at absolute MC (0.03-0.06 events vs
#   ~250-370 signal; decision 2026-08-04 -- no scaling parameter; the lumi lnN
#   does ride on them, consistently with "everything MC scales with L").
#
# The 50-channel card (62 in abcd mode) is written DIRECTLY (no
# combineCards.py), so --dry-run works without cmsenv.  Process indices are
# consistent across channels:
#   signal = 0 (W signal)   zsig = -1 (DY signal under the Z peaks)
#   z = 1   ztau = 2   wtau = 3   qcd = 4   w = 5
#   wfix = 6   ewk = 7   w_y0..w_y11 = 8..19        (abcd-mode CRs only)
# Note lab and fb are the SAME events rebinned -> they live in SEPARATE cards
# (separate likelihoods), never combined with each other.  In abcd mode the
# same CR channels appear in both cards, BUT the CRB w_y* shapes come from
# w_lab_y*/w_fb_y* respectively -- the per-y split must match that card's POIs.
#
# The legacy per-bin pipeline (make_pO_datacards.sh) is untouched and stays
# runnable for comparison.
#
# Usage: make_pO_simfit_cards.sh <muW> <muZ> <eleW> <eleZ> <outdir> [disc]
# bash-3.2 safe (macOS stock bash): no associative arrays.
# =============================================================================
set -euo pipefail

WMU="${1:?need muon W input file}"
ZMU="${2:?need muon Z input file}"
WEL="${3:?need electron W input file}"
ZEL="${4:?need electron Z input file}"
OUTDIR="${5:?need output dir}"
DISC="${6:-met}"
case "$DISC" in
  met)        DISCLABEL="PF MET" ;;
  leppt)      DISCLABEL="lepton pT" ;;
  leppt_mt40) DISCLABEL="lepton pT (pT>25 && mT>40 selection)" ;;
  *)          DISCLABEL="$DISC" ;;
esac
mkdir -p "$OUTDIR"

# ---- systematics config (env-overridable) -----------------------------------
# QCD lnN kappas derived 2026-08-17 from the qcd_abcd.C printout (per component,
# combined in quadrature; same value for both charges of a flavour):
#   mu : stat 9%  (Terr 8.3%, anti-iso mt40 count 3.9%)  (+) T-plane transport
#        10% (m_T-plane vs MET-plane T: 9.7%/10.8%)      (+) anti-iso tilt  5%
#        = 14.4%  -> 1.15
#   ele: stat 4%  (Terr 3.2%, anti-iso mt40 count 2.4%)  (+) transport 10%
#        (9.1% e+, ~0% e- taken as accidental)           (+) anti-iso tilt 15%
#        = 18.5%  -> 1.20
# Override for robustness scans, e.g.  QCD_LNN_MU=1.3 QCD_LNN_ELE=1.3 ./run_pO_fits.sh ...
QCD_MODE="${QCD_MODE:-lnN}"        # lnN | free | abcd  (free = pre-2026-08-17 model;
                                   #  abcd = in-fit ABCD, leppt_mt40 only, 2026-08-23)
QCD_LNN_MU="${QCD_LNN_MU:-1.15}"
QCD_LNN_ELE="${QCD_LNN_ELE:-1.20}"
# abcd-mode REDUCED kappas (residual on the SR qcd only: anti-iso window (+)
# the fake-factor total shift; stat is profiled by the CR channels and the
# plane-transport row does not apply). From correction/logs/qcd_abcd_*.log
# (2026-08-23): mu window 9.5% (+) FF 5.5% -> 1.09 (tilt- and FF-based agree);
# ele window 7.0% (+) FF 13.2% -> 1.15 (the FF-based value; the <pT>-tilt-based
# one is 1.11 -- the FF shift measures the same iso-pT correlation directly, so
# the larger is used).
QCD_ABCD_LNN_MU="${QCD_ABCD_LNN_MU:-1.09}"
QCD_ABCD_LNN_ELE="${QCD_ABCD_LNN_ELE:-1.15}"
QCD_WCR="${QCD_WCR:-float}"        # abcd-mode CRB W content: float (per-y w_y*
                                   #  mapped to the r POIs) | frozen (wfix at
                                   #  absolute MC; add ~2%/1% residual to kappa)
LUMI_LNN="${LUMI_LNN:-1.03}"       # +-3% on kLumi_invnb = 46.5 nb^-1 (all MC)

if [ "$QCD_MODE" = "abcd" ] && [ "$DISC" != "leppt_mt40" ]; then
  echo "[make_pO_simfit_cards] ERROR: QCD_MODE=abcd requires the leppt_mt40 discriminant" >&2
  echo "  (the qcd_abcd template + CR dirs exist only in combine_input_W_leppt_mt40.root;" >&2
  echo "   the met fit's QCD is already data-constrained in-fit by its low-MET region)" >&2
  exit 2
fi

YBINS="0 1 2 3 4 5 6 7 8 9 10 11"

gen_simfit_card() {  # $1 = lab | fb
  B="$1"
  CARD="$OUTDIR/datacard_simfit_${B}.txt"
  MAPS="$OUTDIR/t2w_maps_simfit_${B}.txt"

  SHAPES=""; RP=""
  BINL="bin        "; OBSL="observation"
  MB="bin     "; MP="process "; MI="process "; MR="rate    "
  # systematics rows (one entry per process column, aligned with MP)
  SQMWP="qcd_rate_mu_Wp   lnN"; SQMWM="qcd_rate_mu_Wm   lnN"
  SQEWP="qcd_rate_ele_Wp  lnN"; SQEWM="qcd_rate_ele_Wm  lnN"
  SLUMI="lumi             lnN"
  NCH=0

  # ---- 48 W channels ----------------------------------------------------------
  # abcd mode fits the A0-normalized template (7th input-file object) so the
  # formula rateParam (sB*sC/sD, init 1) needs no baked constants.
  QPATH="qcd"; [ "$QCD_MODE" = "abcd" ] && QPATH="qcd_abcd"
  for F in mu ele; do
    if [ "$F" = "mu" ]; then WF="$WMU"; else WF="$WEL"; fi
    for C in Wp Wm; do
      for iy in $YBINS; do
        R="${C}_${B}_y${iy}"; CH="${F}_${R}"
        SHAPES="${SHAPES}
shapes data_obs ${CH} ${WF} ${R}/data_obs
shapes signal   ${CH} ${WF} ${R}/signal
shapes z        ${CH} ${WF} ${R}/z
shapes ztau     ${CH} ${WF} ${R}/ztau
shapes wtau     ${CH} ${WF} ${R}/wtau
shapes qcd      ${CH} ${WF} ${R}/${QPATH}"
        BINL="$BINL $CH"; OBSL="$OBSL -1"
        MB="$MB $CH $CH $CH $CH $CH"
        MP="$MP signal z ztau wtau qcd"
        MI="$MI 0 1 2 3 4"
        MR="$MR -1 -1 -1 -1 -1"
        if [ "$QCD_MODE" = "free" ]; then
          RP="${RP}
qcd_norm_${CH} rateParam ${CH} qcd 1 [0,10]"
          SQMWP="$SQMWP - - - - -"; SQMWM="$SQMWM - - - - -"
          SQEWP="$SQEWP - - - - -"; SQEWM="$SQEWM - - - - -"
        else
          # lnN mode: the full ABCD kappa; abcd mode: the REDUCED residual
          # (the same 4 row names, on the SR qcd columns only -- never on the
          # CR qcd, whose yields are measurements the scales float).
          if [ "$QCD_MODE" = "abcd" ]; then
            K="$QCD_ABCD_LNN_MU"; [ "$F" = "ele" ] && K="$QCD_ABCD_LNN_ELE"
          else
            K="$QCD_LNN_MU"; [ "$F" = "ele" ] && K="$QCD_LNN_ELE"
          fi
          H=" - - - - ${K}"; N=" - - - - -"
          case "${F}_${C}" in
            mu_Wp)  SQMWP="$SQMWP$H"; SQMWM="$SQMWM$N"; SQEWP="$SQEWP$N"; SQEWM="$SQEWM$N" ;;
            mu_Wm)  SQMWP="$SQMWP$N"; SQMWM="$SQMWM$H"; SQEWP="$SQEWP$N"; SQEWM="$SQEWM$N" ;;
            ele_Wp) SQMWP="$SQMWP$N"; SQMWM="$SQMWM$N"; SQEWP="$SQEWP$H"; SQEWM="$SQEWM$N" ;;
            ele_Wm) SQMWP="$SQMWP$N"; SQMWM="$SQMWM$N"; SQEWP="$SQEWP$N"; SQEWM="$SQEWM$H" ;;
          esac
        fi
        SLUMI="$SLUMI ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} -"
        NCH=$((NCH+1))
      done
    done
  done

  # ---- 2 Z channels (w/wtau deliberately have NO scaling parameter: frozen) ---
  for F in mu ele; do
    if [ "$F" = "mu" ]; then ZF="$ZMU"; else ZF="$ZEL"; fi
    CH="${F}_Z_incl"
    SHAPES="${SHAPES}
shapes data_obs ${CH} ${ZF} Z_incl/data_obs
shapes zsig     ${CH} ${ZF} Z_incl/signal
shapes w        ${CH} ${ZF} Z_incl/w
shapes wtau     ${CH} ${ZF} Z_incl/wtau
shapes ztau     ${CH} ${ZF} Z_incl/ztau"
    BINL="$BINL $CH"; OBSL="$OBSL -1"
    MB="$MB $CH $CH $CH $CH"
    MP="$MP zsig w wtau ztau"
    MI="$MI -1 5 3 2"
    MR="$MR -1 -1 -1 -1"
    SQMWP="$SQMWP - - - -"; SQMWM="$SQMWM - - - -"
    SQEWP="$SQEWP - - - -"; SQEWM="$SQEWM - - - -"
    SLUMI="$SLUMI ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN}"
    NCH=$((NCH+1))
  done

  # ---- 12 ABCD control-region channels (QCD_MODE=abcd only) -------------------
  # Counting channels (1-bin templates from plotting/mtandmet.C). The free
  # scales qcd_s{B,C,D} float the QCD counts; the SR qcd is scaled by the
  # formula rateParam (@0*@1/@2) -- the ABCD relation inside the likelihood.
  # ALIGNMENT RULE: every new (channel, process) column appends exactly one
  # entry to MB/MP/MI/MR AND to all five systematics rows in the same block.
  if [ "$QCD_MODE" = "abcd" ]; then
    for F in mu ele; do
      if [ "$F" = "mu" ]; then WF="$WMU"; else WF="$WEL"; fi
      for C in Wp Wm; do
        # --- CRB: iso-pass, m_T<30 (EWK ~10-20%: z/ztau ride r_Z; W floats or freezes) ---
        CH="${F}_${C}_CRB"
        SHAPES="${SHAPES}
shapes data_obs ${CH} ${WF} ${C}_CRB/data_obs
shapes qcd      ${CH} ${WF} ${C}_CRB/qcd
shapes z        ${CH} ${WF} ${C}_CRB/z
shapes ztau     ${CH} ${WF} ${C}_CRB/ztau"
        BINL="$BINL $CH"; OBSL="$OBSL -1"
        if [ "$QCD_WCR" = "frozen" ]; then
          SHAPES="${SHAPES}
shapes wfix     ${CH} ${WF} ${C}_CRB/wfix"
          MB="$MB $CH $CH $CH $CH"
          MP="$MP qcd z ztau wfix"
          MI="$MI 4 1 2 6"
          MR="$MR -1 -1 -1 -1"
          SQMWP="$SQMWP - - - -"; SQMWM="$SQMWM - - - -"
          SQEWP="$SQEWP - - - -"; SQEWM="$SQEWM - - - -"
          SLUMI="$SLUMI - ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN}"
        else
          MB="$MB $CH $CH $CH"
          MP="$MP qcd z ztau"
          MI="$MI 4 1 2"
          MR="$MR -1 -1 -1"
          SQMWP="$SQMWP - - -"; SQMWM="$SQMWM - - -"
          SQEWP="$SQEWP - - -"; SQEWM="$SQEWM - - -"
          SLUMI="$SLUMI - ${LUMI_LNN} ${LUMI_LNN}"
          # per-y W content: card process w_y<i> <- histogram w_<B>_y<i>
          # (lab and fb cards MUST wire their own split -- see the header note)
          PIDX=8
          for iy in $YBINS; do
            SHAPES="${SHAPES}
shapes w_y${iy}   ${CH} ${WF} ${C}_CRB/w_${B}_y${iy}"
            MB="$MB $CH"; MP="$MP w_y${iy}"; MI="$MI $PIDX"; MR="$MR -1"
            SQMWP="$SQMWP -"; SQMWM="$SQMWM -"; SQEWP="$SQEWP -"; SQEWM="$SQEWM -"
            SLUMI="$SLUMI ${LUMI_LNN}"
            PIDX=$((PIDX+1))
          done
        fi
        NCH=$((NCH+1))
        # --- CRC (anti-iso, m_T>40) and CRD (anti-iso, m_T<30): one frozen ewk ---
        for RG in CRC CRD; do
          CH="${F}_${C}_${RG}"
          SHAPES="${SHAPES}
shapes data_obs ${CH} ${WF} ${C}_${RG}/data_obs
shapes qcd      ${CH} ${WF} ${C}_${RG}/qcd
shapes ewk      ${CH} ${WF} ${C}_${RG}/ewk"
          BINL="$BINL $CH"; OBSL="$OBSL -1"
          MB="$MB $CH $CH"
          MP="$MP qcd ewk"
          MI="$MI 4 7"
          MR="$MR -1 -1"
          SQMWP="$SQMWP - -"; SQMWM="$SQMWM - -"
          SQEWP="$SQEWP - -"; SQEWM="$SQEWM - -"
          SLUMI="$SLUMI - ${LUMI_LNN}"
          NCH=$((NCH+1))
        done
        # --- the three free scales + the functional SR multiplier -------------
        # The formula's ONLY args are the scales (init 1 each), and the SR
        # template total is B0*C40/D0 exactly, so prefit == the ABCD prediction
        # and Asimov closure demands all three scales = 1.
        RP="${RP}
qcd_sB_${F}_${C} rateParam ${F}_${C}_CRB qcd 1 [0,10]
qcd_sC_${F}_${C} rateParam ${F}_${C}_CRC qcd 1 [0,10]
qcd_sD_${F}_${C} rateParam ${F}_${C}_CRD qcd 1 [0,10]
qcd_abcd_${F}_${C} rateParam ${F}_${C}_${B}_y* qcd (@0*@1/@2) qcd_sB_${F}_${C},qcd_sC_${F}_${C},qcd_sD_${F}_${C}"
      done
    done
  fi

  # ---- systematics block: lumi always; the 4 QCD lnN rows only in lnN mode ----
  SYST="
${SLUMI}"
  if [ "$QCD_MODE" != "free" ]; then
    SYST="${SYST}
${SQMWP}
${SQMWM}
${SQEWP}
${SQEWM}"
  fi
  if [ "$QCD_MODE" = "free" ]; then
    QDESC="qcd_norm free rateParam per W channel"
  elif [ "$QCD_MODE" = "abcd" ]; then
    QDESC="QCD in-fit ABCD: 12 CR channels + (sB*sC/sD) formula rateParam on the SR qcd_abcd; residual lnN mu ${QCD_ABCD_LNN_MU} / ele ${QCD_ABCD_LNN_ELE}; CRB W-part ${QCD_WCR}"
  else
    QDESC="QCD lnN (ABCD) per flavour x charge: mu ${QCD_LNN_MU}, ele ${QCD_LNN_ELE}"
  fi

  cat > "$CARD" <<EOF
# Auto-generated by make_pO_simfit_cards.sh -- GRAND SIMULTANEOUS FIT, ${B} binning,
# ${DISCLABEL} W discriminant.  48 W channels ({mu,ele} x {Wp,Wm} x y0..11) + both
# Z-inclusive peaks in ONE likelihood.  The 25-POI model (r_<C>_y<i> shared mu/e,
# + global r_Z on all DY-related MC) is applied at text2workspace time via
# multiSignalModel with the maps in t2w_maps_simfit_${B}.txt.
# ${QDESC}; lumi lnN ${LUMI_LNN} on all MC (not qcd);
# w/wtau under the Z peaks FROZEN at absolute MC (negligible: <0.1 evt).
imax ${NCH}
jmax *
kmax *
------------${SHAPES}
------------
${BINL}
${OBSL}
------------
${MB}
${MP}
${MI}
${MR}
------------${SYST}${RP}
EOF

  # ---- the model: one 'map=' per line, consumed as --PO map=... at t2w time ---
  # re.match runs each regex against "bin/process"; the trailing /... makes y1
  # safe against y10/y11, and the \$ anchors keep 'z' from matching 'ztau'.
  : > "$MAPS"
  for C in Wp Wm; do
    for iy in $YBINS; do
      echo "map=(mu|ele)_${C}_${B}_y${iy}/(signal|wtau)\$:r_${C}_y${iy}[1,0,10]" >> "$MAPS"
    done
  done
  if [ "$QCD_MODE" = "abcd" ] && [ "$QCD_WCR" != "frozen" ]; then
    # CRB W content floats with the SAME POIs (the in-fit EWK subtraction).
    # The process name ENDS in the y index, so the \$ anchor alone keeps
    # w_y1 from matching w_y10 (no trailing-/ trick needed here).
    for C in Wp Wm; do
      for iy in $YBINS; do
        echo "map=(mu|ele)_${C}_CRB/w_y${iy}\$:r_${C}_y${iy}[1,0,10]" >> "$MAPS"
      done
    done
  fi
  echo "map=.*/(z|ztau|zsig)\$:r_Z[1,0,10]" >> "$MAPS"

  echo "[make_pO_simfit_cards] wrote $(basename "$CARD") (${NCH} channels) + $(basename "$MAPS") ($(wc -l < "$MAPS" | tr -d ' ') maps)"
}

gen_simfit_card lab
gen_simfit_card fb

# Sidecar recording the lnN kappas THESE cards were built with -- read back by
# run_pO_fits.sh at extraction time (single source even if env changes between
# card generation and extraction).  kQcd* = 0 means "free rateParam mode":
# extract_pO_simfit.C then reads the per-channel qcd_norm params (legacy path).
# qcdMode (2026-08-23) makes the mode explicit; sidecars WITHOUT the line are
# legacy (mode inferred: kQcd > 0 -> lnN, = 0 -> free). In abcd mode the kQcd*
# entries hold the REDUCED kappas.
KQM="$QCD_LNN_MU"; KQE="$QCD_LNN_ELE"
if [ "$QCD_MODE" = "free" ]; then KQM=0; KQE=0; fi
if [ "$QCD_MODE" = "abcd" ]; then KQM="$QCD_ABCD_LNN_MU"; KQE="$QCD_ABCD_LNN_ELE"; fi
cat > "$OUTDIR/qcd_lnn_kappas.txt" <<EOF
kQcdMu $KQM
kQcdEle $KQE
kLumi $LUMI_LNN
qcdMode $QCD_MODE
EOF

echo "[make_pO_simfit_cards] done -> ${OUTDIR} (W discriminant: ${DISCLABEL})"
echo "[make_pO_simfit_cards] QCD mode: ${QCD_MODE} (mu ${KQM} / ele ${KQE}); lumi lnN ${LUMI_LNN}"
