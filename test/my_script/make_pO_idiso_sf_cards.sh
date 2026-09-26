#!/usr/bin/env bash
# =============================================================================
# make_pO_idiso_sf_cards.sh -- datacard + multiSignalModel maps of the
# ELECTRON ID+ISO SF FIT (2026-09-24). A study SEPARATE from the nominal W/Z
# fits: nothing here is read by run_pO_fits.sh, and this script reads nothing
# of theirs.
#
# What is measured (definition: the analysis repo's correction/idiso_sf_common.h)
# -- the efficiency of eleMVAIdWP90 && eleMVAIsoWP90 on W -> e nu events,
#   eps = N_W(pass) / N_W(pass + fail),
# with N_W the POST-FIT W yields (the electron QCD is far too large to count).
# pass and fail are DISJOINT event samples, so they are channels of ONE
# likelihood and total = pass + fail needs no channel of its own.
#
# TWO KINDS OF CARD, chosen by the input's scheme (its _meta.txt):
#
# (1) a W input (scheme incl|abseta|pt; the lepton pT at m_T > 40): per charge
# C = Wp|Wm, category CAT = pass|fail and coarse bin k -- R = <C>_<CAT>_k<k>:
#   ele_R        the SR (m_T > 40, pT)       processes signal z ztau wtau qcd
#   ele_R_CRB    1-bin, the category, m_T<30  processes qcd signal wtau z ztau
#   ele_R_CRC    1-bin, sideband, m_T > 40    processes qcd ewk
#   ele_R_CRD    1-bin, sideband, m_T < 30    processes qcd ewk
# plus the DY anchor
#   ele_Z_PP (both Z TnP legs pass)           processes zsig w wtau ztau bkg
# The model (map file, one 'map=' per line, --PO map=... at t2w time):
#   r_R              scales the W MC (signal + wtau) of the SR AND of its CRB --
#                    independent per channel, so the post-fit W counts N = r x S
#                    of pass and fail are what the extraction turns into
#                    eps_data (extract_pO_idiso_sf.C); in CRB it is the in-fit
#                    EWK subtraction of the ABCD
#   r_Z              scales every DY template (z, ztau of the SR and CRB, zsig
#                    and ztau of Z_PP) -- pinned by the Z_PP peak
# QCD = THE NOMINAL IN-FIT ABCD (QCD_MODE=abcd of make_pO_simfit_cards.sh; user
# 2026-09-25: "always the in-fit ABCD, so we don't assume the prefit
# normalization"): the free scales qcd_s{B,C,D}_ele_R ([0,10]) float the three
# CR QCD counts, and the SR qcd (total B0 C0 / D0) is scaled by the formula
# rateParam qcd_abcd_ele_R = (@0*@1/@2) of them; CRC/CRD ewk frozen. The residual
# lnN qcd_rate_ele_<C>_<CAT> (kappa IDISO_QCD_ABCD_LNN, default 1.15 = the
# nominal electron value; 'none' drops it) sits on the SR qcd only, correlated
# across the coarse bins -- never on a CR qcd, whose yields are measurements.
#
# (2) the Z tag-and-probe input (scheme ztnp, 2026-09-25) -- its OWN fit:
#   ele_Z_PP (both legs pass)              processes zsig w wtau ztau bkg
#   ele_Z_PF (tag passes, probe fails)     processes zsig w wtau ztau bkg
#   r_ZPP, r_ZPF     scale the DY (zsig, ztau) of each category
#   -> eps_Z = 2 N_PP / (2 N_PP + N_PF) (a PP event holds two passing probes)
# It is not in the W likelihood on purpose: a shared r_Z let the W pass
# channels (~100 DY events each at low MET, degenerate there with the QCD
# normalization) pull the Z efficiency -- SF_Z moved 0.96 -> 1.01 across the W
# variants, whose Z channels are identical.
#
# Both: bkg_norm_<CH> = free rateParam of the flat Z background, init 1,
# [0,100]; w, wtau under the Z peaks frozen at MC (as in the nominal fit).
# No lumi / theory rows: a common normalization cancels in eps; this is a
# cross-check, not a cross section.
#
# Usage:  make_pO_idiso_sf_cards.sh <combine_input.root> <outdir> <tag>
#   <input minus .root>_meta.txt (written by correction/idiso_sf_inputs.C)
#   provides nbins -- the only thing the card needs to know about the scheme.
# Writes <outdir>/datacard_idiso_<tag>.txt and <outdir>/t2w_maps_idiso_<tag>.txt.
# Process indices follow the nominal generator's table:
#   signal = 0, zsig = -1, z = 1, ztau = 2, wtau = 3, qcd = 4, w = 5
# bash-3.2 safe.
# =============================================================================
set -euo pipefail

if [ $# -lt 3 ]; then
  echo "usage: $0 <combine_input.root> <outdir> <tag>" >&2
  exit 2
fi
IN="$1"; OUTDIR="$2"; TAG="$3"
META="${IN%.root}_meta.txt"
[ -f "$IN" ]   || { echo "[make_pO_idiso_sf_cards] ERROR no input $IN" >&2; exit 2; }
[ -f "$META" ] || { echo "[make_pO_idiso_sf_cards] ERROR no sidecar $META" >&2; exit 2; }
NB=$(awk '$1=="nbins"{print $2}' "$META")
case "$NB" in ''|*[!0-9]*) echo "[make_pO_idiso_sf_cards] ERROR bad nbins '$NB' in $META" >&2; exit 2 ;; esac
SCHEME=$(awk '$1=="scheme"{print $2}' "$META")
DISC=$(awk '$1=="disc"{print $2}' "$META")
SBWIN=$(awk '$1=="sbwin"{print $2}' "$META")
KAPPA="${IDISO_QCD_ABCD_LNN:-1.15}"   # the residual lnN on the SR qcd ('none' = no row)
mkdir -p "$OUTDIR"
CARD="$OUTDIR/datacard_idiso_${TAG}.txt"
MAPS="$OUTDIR/t2w_maps_idiso_${TAG}.txt"

SHAPES=""; RP=""
BINL="bin        "; OBSL="observation"
MB="bin     "; MP="process "; MI="process "; MR="rate    "
NCH=0
if [ "$SCHEME" = "ztnp" ]; then ZTNP=1; WCHG=""; ZCATS="PP PF"; else ZTNP=0; WCHG="Wp Wm"; ZCATS="PP"; fi
for C in $WCHG; do
  k=0
  while [ "$k" -lt "$NB" ]; do
    for CAT in pass fail; do
      R="${C}_${CAT}_k${k}"; CH="ele_${R}"
      SHAPES="${SHAPES}
shapes data_obs ${CH} ${IN} ${R}/data_obs
shapes signal   ${CH} ${IN} ${R}/signal
shapes z        ${CH} ${IN} ${R}/z
shapes ztau     ${CH} ${IN} ${R}/ztau
shapes wtau     ${CH} ${IN} ${R}/wtau
shapes qcd      ${CH} ${IN} ${R}/qcd"
      BINL="$BINL $CH"; OBSL="$OBSL -1"
      MB="$MB $CH $CH $CH $CH $CH"
      MP="$MP signal z ztau wtau qcd"
      MI="$MI 0 1 2 3 4"
      MR="$MR -1 -1 -1 -1 -1"
      NCH=$((NCH+1))
      # --- the in-fit ABCD control regions of this channel (1-bin counting) ---
      SHAPES="${SHAPES}
shapes data_obs ${CH}_CRB ${IN} ${R}_CRB/data_obs
shapes qcd      ${CH}_CRB ${IN} ${R}_CRB/qcd
shapes signal   ${CH}_CRB ${IN} ${R}_CRB/signal
shapes wtau     ${CH}_CRB ${IN} ${R}_CRB/wtau
shapes z        ${CH}_CRB ${IN} ${R}_CRB/z
shapes ztau     ${CH}_CRB ${IN} ${R}_CRB/ztau"
      BINL="$BINL ${CH}_CRB"; OBSL="$OBSL -1"
      MB="$MB ${CH}_CRB ${CH}_CRB ${CH}_CRB ${CH}_CRB ${CH}_CRB"
      MP="$MP qcd signal wtau z ztau"
      MI="$MI 4 0 3 1 2"
      MR="$MR -1 -1 -1 -1 -1"
      NCH=$((NCH+1))
      for RG in CRC CRD; do
        SHAPES="${SHAPES}
shapes data_obs ${CH}_${RG} ${IN} ${R}_${RG}/data_obs
shapes qcd      ${CH}_${RG} ${IN} ${R}_${RG}/qcd
shapes ewk      ${CH}_${RG} ${IN} ${R}_${RG}/ewk"
        BINL="$BINL ${CH}_${RG}"; OBSL="$OBSL -1"
        MB="$MB ${CH}_${RG} ${CH}_${RG}"
        MP="$MP qcd ewk"
        MI="$MI 4 7"
        MR="$MR -1 -1"
        NCH=$((NCH+1))
      done
      # the three free scales + the SR multiplier (@0*@1/@2): prefit == the ABCD
      # prediction, and the Asimov closure demands all three scales = 1
      RP="${RP}
qcd_sB_ele_${R} rateParam ${CH}_CRB qcd 1 [0,10]
qcd_sC_ele_${R} rateParam ${CH}_CRC qcd 1 [0,10]
qcd_sD_ele_${R} rateParam ${CH}_CRD qcd 1 [0,10]
qcd_abcd_ele_${R} rateParam ${CH} qcd (@0*@1/@2) qcd_sB_ele_${R},qcd_sC_ele_${R},qcd_sD_ele_${R}"
    done
    k=$((k+1))
  done
done
NZCOL=0
for ZC in $ZCATS; do
  R="Z_${ZC}"; CH="ele_${R}"; NZCOL=$((NZCOL+5))
  SHAPES="${SHAPES}
shapes data_obs ${CH} ${IN} ${R}/data_obs
shapes zsig     ${CH} ${IN} ${R}/signal
shapes w        ${CH} ${IN} ${R}/w
shapes wtau     ${CH} ${IN} ${R}/wtau
shapes ztau     ${CH} ${IN} ${R}/ztau
shapes bkg      ${CH} ${IN} ${R}/bkg"
  BINL="$BINL $CH"; OBSL="$OBSL -1"
  MB="$MB $CH $CH $CH $CH $CH"
  MP="$MP zsig w wtau ztau bkg"
  MI="$MI -1 5 3 2 6"
  MR="$MR -1 -1 -1 -1 -1"
  RP="${RP}
bkg_norm_${CH} rateParam ${CH} bkg 1 [0,100]"
  NCH=$((NCH+1))
done

# the residual lnN: one row per (charge, category) -- kappa on the SR qcd column
# of every coarse bin of that (charge, category), '-' on every other column. The
# column order is the loop order above: per W channel 5 (SR) + 5 (CRB) + 2 (CRC)
# + 2 (CRD), then 5 per Z channel.
lnn_row() {
  local row="qcd_rate_ele_$1_$2 lnN" c k cat i=0
  for c in Wp Wm; do
    k=0
    while [ "$k" -lt "$NB" ]; do
      for cat in pass fail; do
        if [ "$c" = "$1" ] && [ "$cat" = "$2" ]; then row="$row - - - - $KAPPA"; else row="$row - - - - -"; fi
        row="$row - - - - - - - - -"   # CRB (5), CRC (2), CRD (2): never on a CR
      done
      k=$((k+1))
    done
  done
  while [ "$i" -lt "$NZCOL" ]; do row="$row -"; i=$((i+1)); done
  echo "$row"
}
SYST=""
if [ "$KAPPA" != "none" ] && [ "$ZTNP" -eq 0 ]; then
  for C in Wp Wm; do
    for CAT in pass fail; do
      SYST="${SYST}
$(lnn_row "$C" "$CAT")"
    done
  done
fi

cat > "$CARD" <<EOF
# Auto-generated by make_pO_idiso_sf_cards.sh -- ELECTRON ID+ISO SF FIT, tag ${TAG}
# (scheme ${SCHEME}, discriminant ${DISC}, QCD sideband ${SBWIN}). SEPARATE from the
# nominal fits. W input: per charge and coarse bin a pass and a fail W channel
# (eleMVAIdWP90 && eleMVAIsoWP90 passed / failed by the leading electron; the
# lepton pT at m_T > 40) with its in-fit ABCD control regions CRB/CRC/CRD + the
# Z_PP DY anchor; ztnp input: the Z tag-and-probe (Z_PP both pass, Z_PF fails).
# Model: t2w_maps_idiso_${TAG}.txt. QCD: in-fit ABCD, residual lnN ${KAPPA}
# (W cards only); bkg_norm_* free.
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

# re.match runs each regex against "bin/process"; the trailing / keeps k1 from
# matching k10 (should the schemes ever grow), the \$ keeps 'z' off 'ztau'.
: > "$MAPS"
if [ "$ZTNP" -eq 1 ]; then
  echo "map=ele_Z_PP/(zsig|ztau)\$:r_ZPP[1,0,10]" >> "$MAPS"
  echo "map=ele_Z_PF/(zsig|ztau)\$:r_ZPF[1,0,10]" >> "$MAPS"
else
  for C in Wp Wm; do
    k=0
    while [ "$k" -lt "$NB" ]; do
      for CAT in pass fail; do
        # the SR and its CRB (the in-fit EWK subtraction rides the same r)
        echo "map=ele_${C}_${CAT}_k${k}(_CRB)?/(signal|wtau)\$:r_${C}_${CAT}_k${k}[1,0,10]" >> "$MAPS"
      done
      k=$((k+1))
    done
  done
  echo "map=.*/(z|ztau|zsig)\$:r_Z[1,0,10]" >> "$MAPS"
fi

if [ "$ZTNP" -eq 1 ]; then WHAT="Z tag-and-probe"; else WHAT="QCD in-fit ABCD, residual lnN ${KAPPA}"; fi
echo "[make_pO_idiso_sf_cards] wrote $(basename "$CARD") (${NCH} channels, ${WHAT}) + $(basename "$MAPS") ($(wc -l < "$MAPS" | tr -d ' ') maps)"
