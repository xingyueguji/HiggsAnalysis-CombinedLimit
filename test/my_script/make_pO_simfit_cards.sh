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
#   lumi lnN (2026-08-17): one global nuisance on EVERY MC template (signal, z,
#   ztau, wtau, zsig -- NOT the data-driven qcd), +-3% on L = 46.5 nb^-1.
#   w, wtau under the Z peaks: FROZEN at absolute MC (0.03-0.06 events vs
#   ~250-370 signal; decision 2026-08-04 -- no scaling parameter; the lumi lnN
#   does ride on them, consistently with "everything MC scales with L").
#
# The 50-channel card is written DIRECTLY (no combineCards.py), so --dry-run
# works without cmsenv.  Process indices are consistent across channels:
#   signal = 0 (W signal)   zsig = -1 (DY signal under the Z peaks)
#   z = 1   ztau = 2   wtau = 3   qcd = 4   w = 5
# Note lab and fb are the SAME events rebinned -> they live in SEPARATE cards
# (separate likelihoods), never combined with each other.
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
QCD_MODE="${QCD_MODE:-lnN}"        # lnN | free  (free = pre-2026-08-17 model)
QCD_LNN_MU="${QCD_LNN_MU:-1.15}"
QCD_LNN_ELE="${QCD_LNN_ELE:-1.20}"
LUMI_LNN="${LUMI_LNN:-1.03}"       # +-3% on kLumi_invnb = 46.5 nb^-1 (all MC)

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
shapes qcd      ${CH} ${WF} ${R}/qcd"
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
          K="$QCD_LNN_MU"; [ "$F" = "ele" ] && K="$QCD_LNN_ELE"
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
  echo "map=.*/(z|ztau|zsig)\$:r_Z[1,0,10]" >> "$MAPS"

  echo "[make_pO_simfit_cards] wrote $(basename "$CARD") (${NCH} channels) + $(basename "$MAPS") ($(wc -l < "$MAPS" | tr -d ' ') maps)"
}

gen_simfit_card lab
gen_simfit_card fb

# Sidecar recording the lnN kappas THESE cards were built with -- read back by
# run_pO_fits.sh at extraction time (single source even if env changes between
# card generation and extraction).  kQcd* = 0 means "free rateParam mode":
# extract_pO_simfit.C then reads the per-channel qcd_norm params (legacy path).
KQM="$QCD_LNN_MU"; KQE="$QCD_LNN_ELE"
if [ "$QCD_MODE" = "free" ]; then KQM=0; KQE=0; fi
cat > "$OUTDIR/qcd_lnn_kappas.txt" <<EOF
kQcdMu $KQM
kQcdEle $KQE
kLumi $LUMI_LNN
EOF

echo "[make_pO_simfit_cards] done -> ${OUTDIR} (W discriminant: ${DISCLABEL})"
echo "[make_pO_simfit_cards] QCD mode: ${QCD_MODE} (mu ${KQM} / ele ${KQE}); lumi lnN ${LUMI_LNN}"
