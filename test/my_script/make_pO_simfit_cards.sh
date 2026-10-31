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
# PER-FLAVOUR cards (2026-09-22, env SIMFIT_FLAVS = mu | ele; the default
# "mu ele" is the grand fit above): the SAME model restricted to ONE lepton
# flavour -- its 24 W channels + its own Z peak (+ its 18 CR channels in abcdy
# mode), still 25 POIs (r_<C>_y<i> + r_Z, now fitted by that flavour alone),
# and every nuisance that acts on it: lumi, that flavour's two QCD rows, and
# the shape rows ITS sidecars list (so muSF only in the muon card). Rows that
# would be all '-' (the other flavour's QCD lnN) are not written. Driven by
# run_pO_fits.sh mode `flavfit` (work dirs simfit_mu/, simfit_ele/); the
# other flavour's input arguments are then ignored (pass "none"). With both
# flavours the cards and maps are byte-identical to the pre-2026-09-22 ones.
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
#   QCD -- TWO models since 2026-10-26 (user decision), one per discriminant:
#   met (PF MET, no m_T cut), QCD_MODE=lnN (default) | free: log-normal-
#   constrained at the ABCD prediction -- ONE lnN nuisance per (flavour,
#   charge), qcd_rate_{mu,ele}_{Wp,Wm}, correlated across the 12 y bins of that
#   flavour+charge (T and the template are measured once, inclusively in y;
#   uncorrelated across charge/flavour so the charge asymmetry gets no unearned
#   cancellation; the low-MET bins of the fit pin the level in-fit).
#   QCD_MODE=free = 48 free qcd_norm_<channel> rateParams (pre-2026-08-17).
#   leppt_mt40 (lepton pT, pT>25 && m_T>40), QCD_MODE=abcdy (2026-10-19; its
#   only model since 2026-10-26): the IN-FIT ABCD done PER RAPIDITY GROUP.
#   The groups come from the input's group sidecar <W input minus .root>_abcdy.txt
#   (plotting/mtandmet.C, lines '<lab|fb> <group> <bin indices>'; today
#   fwd = y10,y11 | bwd = y0,y1 | cen = y2..y9, correction/abcd_ygroups.h) --
#   abcdy refuses to run without it, and the muon and electron sidecars must
#   agree. Per (flavour, charge, group g): counting CR channels
#   <F>_<C>_<g>_CR{B,C,D} ((relIso x m_T) plane: CRB iso-pass m_T<30, CRC
#   anti-iso m_T>40 = the SR template's source region, CRD anti-iso m_T<30;
#   1-bin templates from the input dirs <C>_<B>_<g>_CR*, the group's own
#   control regions, FB-binned planes in the fb card), free scales
#   qcd_s{B,C,D}_<F>_<C>_<g> floating the three QCD counts, and the FORMULA
#   rateParam qcd_abcd_<F>_<C>_<g> = (@0*@1/@2) -- Combine's documented ABCD
#   pattern (docs/part2/settinguptheanalysis.md) -- written once PER SR
#   CHANNEL of the group under the same name with the identical formula
#   (Combine builds ONE function from the repeated lines:
#   ModelTools.doRateParams skips a name that already exists). The SR qcd
#   column reads qcd_abcdy (each group carrying its own B_g*C40_g/D_g prefit
#   total), so the normalization floats with the CR data and the EWK
#   subtraction rides the POIs: CRB z/ztau are swept into r_Z by the catch-all
#   map, the CRB W parts w_y<i> of a group map to r_<C>_y<i> (QCD_WCR=float,
#   default; maps and sigma maps name the group channel) or one frozen wfix
#   (QCD_WCR=frozen). The residual lnN on the SR qcd is, since 2026-10-30
#   (user decision, QCD_RATE_SPLIT=group, the default), ONE PER (flavour,
#   charge, group), qcd_rate_<F>_<C>_<g>, with a kappa per (flavour, group)
#   from correction/logs/qcd_abcd_*.log "PER-GROUP RESIDUAL lnN SIZING" (see
#   the QCD_ABCD_LNN_<FLAV>_<GROUP> defaults); until then (and with
#   QCD_RATE_SPLIT=charge) one qcd_rate_<F>_<C> per (flavour, charge), shared
#   by the groups, with the REDUCED kappas QCD_ABCD_LNN_MU/ELE (window only for
#   e, window (+) FF shift for mu -- stat is profiled in-fit). The kappa sidecar
#   gains a qcdGroups line and, when split, a qcdRateGroups line (both read
#   back by run_pO_fits.sh for the extractor).
#   (Removed 2026-10-26: QCD_MODE=abcd = the same in-fit ABCD with ONE set of
#   CRs per (flavour, charge), 2026-08-23..10-26; QCD_MODE=mtsplit = the no-m_T
#   lepton-pT fit, disc leppt, 2026-10-21..26; disc met_mt40, 2026-10-20..26.)
#   lumi lnN (2026-08-17): one global nuisance on EVERY MC template (signal, z,
#   ztau, wtau, zsig, and the CR wfix/w_y*/ewk -- NOT the data-driven qcd),
#   +-3% on L = 46.5 nb^-1.
#   w, wtau under the Z peaks: FROZEN at absolute MC (0.03-0.06 events vs
#   ~250-370 signal; decision 2026-08-04 -- no scaling parameter; the lumi lnN
#   does ride on them, consistently with "everything MC scales with L").
#   ttbar (2026-10-16): the inclusive POWHEG ttbar background, its own column
#   `ttbar` (process id 20) in every W channel, both Z peaks and CRB (abcdy
#   mode; in CRC/CRD it sits inside the frozen `ewk`). FROZEN at absolute MC:
#   no map matches it (scale 1 -- it must ride neither r_<C>_y<i> nor r_Z),
#   only the lumi lnN and its own normalization lnN `ttbar_xsec`, kappa
#   TTBAR_LNN (default 1.15 = the sample's own inclusive scale envelope
#   +13.3/-12.9% (+) EPPS21 ~3% at 68% CL, skim/output/lhe_weights_TT.txt),
#   plus the shape rows its sidecars list (LHE theory -- area-normalized,
#   shape only, so no double counting with ttbar_xsec -- and the lepton SFs).
#   LHE SHAPE systematics (2026-09-07): the Combine inputs carry, for every MC
#   process of every SR region, <process>_<syst>Up/Down built from the
#   generator's per-event LHE weights (plotting/mtandmet.C + dileptonpeak.C
#   from the skim twins; skim/lhe_updown.py):
#     nPDF      EPPS21nlo_CT18Anlo_O16, LHAPDF PDFSet.uncertainty() per bin
#     qcdScale  muR/muF envelope (per-bin max/min over the 9-point grid)
#     alphaS    the alpha_s 0.119 / 0.117 member templates
#   The sidecar next to each input (<input minus .root>_systs.txt) lists what
#   was written; this script uses the UNION of the four sidecars and emits ONE
#   `shape` row each (entries 1 on the MC columns of the inputs whose sidecar
#   lists it, - on the data-driven qcd, on every CR column and on the columns
#   of inputs that do not list it) plus the fifth `shapes` token
#   <dir>/<proc>_$SYSTEMATIC on the MC processes. One row name = one nuisance
#   for the whole card -> fully correlated across channels, flavours, charges
#   and processes (the correlation that cancels in the charge asymmetry). A
#   `lhe` group holds the theory ones for --freezeNuisanceGroups.
#   LEPTON-SF SHAPE systematics (2026-09-14, muon first): the muon inputs carry
#   <process>_muSFUp/Down (skim/muon_sf.h) = ONE combined muon-SF nuisance --
#   the pp-2025 POG ID and ISO SFs and our MB-derived inclusive trigger SF,
#   each varied +-1 sigma, added in quadrature bin by bin (three independent
#   sources -> the exact 1 sigma of the product; user decision: one Up and one
#   Down for the fit). Listed by the muon sidecars only -> `1` on the muon
#   channels' MC columns, `-` on every electron column; group `lepsf`. The
#   analysis repo can alternatively list the three sources separately
#   (muID/muIso/muTrig); everything below then applies per source, and the
#   muTrig correlation model follows the trigger-SF binning chosen
#   in the analysis repo (skim/muon_sf.h; recorded in the muon sidecars as the
#   directive '#! muTrig corr <coherent|perbin>', read when SF_TRIG_CORR=auto,
#   the default): the SF is INCLUSIVE in rapidity since 2026-09-14 (the per-y
#   values are consistent with a flat SF, chi2/ndf 11.5/11) -> ONE coherent
#   nuisance. With a per-y SF (pure-statistics per-bin errors) the directive
#   says perbin and the single template variation is split into 12 independent
#   nuisances muTrig_y<i> with `nuisance edit rename` (one per rapidity bin,
#   both charges; Combine maps each new name back to the SAME
#   <proc>_muTrigUp/Down histograms, no duplicated templates; the Z-peak piece
#   stays under muTrig). SF_TRIG_CORR=perbin|coherent overrides the directive.
#   ELECTRON SFs (2026-09-29): the electron inputs carry <process>_eleSFUp/Down
#   (skim/electron_sf.h: EGM Reco + wp90iso ID+iso + our MB trigger SF, in
#   quadrature) -- the muSF recipe on the electron columns; group `lepsf`.
#   ELECTRON FAKE-FACTOR QCD (2026-09-29): the electron SR qcd template is
#   fake-factor weighted (correction/qcd_abcd.C buildFFTemplate) and its one
#   shape parameter -- the slope of F(pT) -- is a Gaussian-constrained shape
#   nuisance PER CHARGE, qcdFF_ele_Wp / qcdFF_ele_Wm, listed in the electron
#   leppt_mt40 sidecar as 'qcdFF_ele_Wp qcd_abcdy@Wp': entries `1` on the SR
#   qcd column (histogram qcd_abcdy) of that charge's electron channels only
#   (the '@<charge>' token), group `qcdff`. The qcd `shapes` line carries the
#   <dir>/qcd_abcdy_$SYSTEMATIC token for it. A shape row that ends up with no
#   entry at all is dropped with a note.
#     LHE_SYST=auto  (default) the sidecars' union
#     LHE_SYST=off   no shape rows (cards byte-identical to the pre-09-07 ones)
#     LHE_SYST=a,b   an explicit subset (each must be listed by some sidecar)
#
# The 50-channel card (86 in abcdy mode) is written DIRECTLY (no
# combineCards.py), so --dry-run works without cmsenv.  Process indices are
# consistent across channels:
#   signal = 0 (W signal)   zsig = -1 (DY signal under the Z peaks)
#   z = 1   ztau = 2   wtau = 3   qcd = 4   w = 5
#   wfix = 6   ewk = 7   w_y0..w_y11 = 8..19        (abcdy-mode CRs only; w_y<i> = 8+i)
#   ttbar = 20                                      (2026-10-16, every channel but CRC/CRD)
# Note lab and fb are the SAME events rebinned -> they live in SEPARATE cards
# (separate likelihoods), never combined with each other.  In abcdy mode the
# lab card reads the lab group CRs and the fb card the FB-binned ones, and the
# CRB w_y* shapes come from w_lab_y*/w_fb_y* respectively -- the per-y split
# must match that card's POIs.
#
# (The legacy per-bin pipeline, make_pO_datacards.sh, was removed 2026-09-22:
# the per-flavour fit is this generator with SIMFIT_FLAVS, see above.)
#
# Usage: make_pO_simfit_cards.sh <muW> <muZ> <eleW> <eleZ> <outdir> [disc]
#        (SIMFIT_FLAVS=mu|ele: the absent flavour's two inputs are not read)
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
  leppt_mt40) DISCLABEL="lepton pT (pT>25 && mT>40 selection)" ;;
  *) echo "[make_pO_simfit_cards] ERROR disc '$DISC' (leppt_mt40 | met; leppt and met_mt40 were removed 2026-10-26)" >&2; exit 2 ;;
esac
mkdir -p "$OUTDIR"

# ---- systematics config (env-overridable) -----------------------------------
# QCD lnN kappas of the met fit (QCD_MODE=lnN), since 2026-10-26 (user
# decision) = the Oct-22 correction/logs/qcd_abcd_{mu,ele}.log "systematic
# budget on the multijet normalisation", line "lnN kappa EXCLUDING the
# transport row" (max of the two charges; components in quadrature):
#   mu : stat on T 6.9/6.8% (+) sideband stat 3.8/3.8% (+) anti-iso window
#        12.8/9.6% (+) anti-iso tilt 0.6/2.7%  (mu+/mu-)            -> 1.15
#   ele (since 2026-10-28, user decision: the eleMVAIsoWP90-only selection and
#        its pass/fail ABCD plane): stat on T 2.3/2.2% (+) sideband stat
#        0.4/0.4% (+) anti-iso window 15.9/15.7% (+) anti-iso tilt 14.5/14.2%
#        (e+/e-)                                                     -> 1.22
# The plane-transport row (m_T vs MET plane T: 22.8/18.8% mu) is left out
# (user decision; with it mu 1.27). For the muon the log's r-scan has the two
# planes agree at r = 1.21, i.e. near the fitted r, so that row is the prefit W
# normalization and would double-count the signal strength. For the electron
# it is 1.8/1.9% since 2026-10-28 and moves nothing (1.22 with or without).
# History: ele 1.28 2026-10-26..28 (Oct-22, MVAId WP90 + MVAIso WP90 + relIso
# < 0.3: window 26.2/17.6% (+) tilt 7.5/8.0%, transport 14.6/11.6% left out);
# before 2026-10-26 mu 1.15 / e 1.20, derived 2026-08-17 on the pre-WP90
# electron selection.
# Override for robustness scans, e.g.  QCD_LNN_MU=1.3 QCD_LNN_ELE=1.3 ./run_pO_fits.sh ...
# QCD_MODE is DISC-DEPENDENT (2026-10-26, user decision: one multijet model per
# fit): abcdy for leppt_mt40 (its only model -- the qcd_abcdy template and the
# group CR dirs exist only in combine_input_W_leppt_mt40.root), lnN (default)
# or free for met. A mode that does not belong to the discriminant is an error.
if [ -z "${QCD_MODE:-}" ]; then
  if [ "$DISC" = "leppt_mt40" ]; then QCD_MODE=abcdy; else QCD_MODE=lnN; fi
  echo "[make_pO_simfit_cards] QCD_MODE not set -> default '$QCD_MODE' for disc '$DISC'"
fi
QCD_LNN_MU="${QCD_LNN_MU:-1.15}"
QCD_LNN_ELE="${QCD_LNN_ELE:-1.22}"
# abcdy-mode REDUCED kappas (residual on the SR qcd only: anti-iso window (+)
# the fake-factor total shift; stat is profiled by the CR channels and the
# plane-transport row does not apply). Since 2026-10-26 (user decision) from
# the Oct-22 correction/logs/qcd_abcd_*.log (max of the two charges):
#   mu : "in-fit ABCD reduced lnN kappa" = anti-iso window 10.1%/3.0% (+) tilt
#        0.6%/2.7% (mu+/mu-) -> 1.10 (the FF-shift variant gives 1.10 too)
#   ele: the SR template IS fake-factor weighted and its slope uncertainty is
#        the qcdFF_ele_<C> shape nuisance, so neither the FF-shift nor the tilt
#        row remains -- "residual lnN kappa of the FF model" = the anti-iso
#        window variation of the FF-weighted total, 4.5% (e+) / 9.9% (e-)
#        -> 1.099, rounded up to 1.10
# Before 2026-10-26: mu 1.09 (2026-08-23 log, window 9.5% (+) FF 5.5%), ele 1.05
# (2026-09-29 log, FF window 4.5%/4.1%; before 2026-09-29 window (+) FF shift
# -> 1.15).
QCD_ABCD_LNN_MU="${QCD_ABCD_LNN_MU:-1.10}"
QCD_ABCD_LNN_ELE="${QCD_ABCD_LNN_ELE:-1.10}"
# abcdy residual lnN PER RAPIDITY GROUP (2026-10-30, user decision): one
# qcd_rate_<F>_<C>_<g> per (flavour, charge, group) instead of one per
# (flavour, charge) shared by the groups -- the groups' CR -> SR transfer is
# allowed to fail independently (the e-only fit of 2026-10-28 needed the edges'
# multijet x1.1-1.5 and the centre's x0.6-0.8 of the template, and the shared
# row could only move them together). kappa per (flavour, group) from
# correction/logs/qcd_abcd_<lep>.log "PER-GROUP RESIDUAL lnN SIZING" =
# 1 + sqrt(window^2 + (the m_T 30-40 closure beyond its stat error)^2), max of
# the two charges, rounded up; a group with no QCD_ABCD_LNN_<FLAV>_<GROUP>
# falls back to QCD_ABCD_LNN_<FLAV>. QCD_RATE_SPLIT=charge restores the shared
# rows (= the cards of the 2026-10-28 fit).
QCD_RATE_SPLIT="${QCD_RATE_SPLIT:-group}"
QCD_ABCD_LNN_MU_FWD="${QCD_ABCD_LNN_MU_FWD:-1.44}"    # window 18.7%, closure 1.77 +- 0.63 (mu-) -> 43% beyond stat
QCD_ABCD_LNN_MU_BWD="${QCD_ABCD_LNN_MU_BWD:-1.08}"    # window 7.6%, closures within stat
QCD_ABCD_LNN_MU_CEN="${QCD_ABCD_LNN_MU_CEN:-1.16}"    # window 14.3%, closure 0.77 +- 0.17 (mu-)
QCD_ABCD_LNN_ELE_FWD="${QCD_ABCD_LNN_ELE_FWD:-1.25}"  # window 9.4%, closure 1.26 +- 0.13 (e-)
QCD_ABCD_LNN_ELE_BWD="${QCD_ABCD_LNN_ELE_BWD:-1.29}"  # window 10.6%, closure 1.30 +- 0.14 (e-)
QCD_ABCD_LNN_ELE_CEN="${QCD_ABCD_LNN_ELE_CEN:-1.19}"  # window 1.4%, closure 0.80 +- 0.08 (e+)
case "$QCD_RATE_SPLIT" in
  group|charge) ;;
  *) echo "[make_pO_simfit_cards] ERROR QCD_RATE_SPLIT='$QCD_RATE_SPLIT' (group | charge)" >&2; exit 2 ;;
esac
qrate_kappa() {  # $1 = mu|ele, $2 = group -> the per-group kappa (fallback: the flavour's)
  local FU GU v
  FU=$(echo "$1" | tr '[:lower:]' '[:upper:]'); GU=$(echo "$2" | tr '[:lower:]' '[:upper:]')
  eval "v=\${QCD_ABCD_LNN_${FU}_${GU}:-}"
  if [ -z "$v" ]; then
    if [ "$1" = "ele" ]; then v="$QCD_ABCD_LNN_ELE"; else v="$QCD_ABCD_LNN_MU"; fi
  fi
  echo "$v"
}
QCD_WCR="${QCD_WCR:-float}"        # abcdy-mode CRB W content: float (per-y w_y*
                                   #  mapped to the r POIs) | frozen (wfix at
                                   #  absolute MC; add ~2%/1% residual to kappa)
LUMI_LNN="${LUMI_LNN:-1.03}"       # +-3% on kLumi_invnb = 46.5 nb^-1 (all MC)
TTBAR_LNN="${TTBAR_LNN:-1.15}"     # ttbar_xsec: normalization of the frozen ttbar MC
                                   #  (scale +13.3/-12.9% (+) EPPS21 ~3%, 68% CL; see header)
LHE_SYST="${LHE_SYST:-auto}"       # auto | off | comma list (see header)
SF_TRIG_CORR="${SF_TRIG_CORR:-auto}"    # auto (the muon W sidecar's '#! muTrig corr' directive, coherent if
                                        #  absent) | perbin | coherent -- see header
SIGMA_POI="${SIGMA_POI:-0}"        # 1 -> ALSO write t2w_maps_simfit_<B>_sigma.txt, the
                                   #  reparametrization that makes the rapidity-inclusive
                                   #  sigma_W a POI (same card; see write_sigma_maps)
GEN_XSEC_FID="${GEN_XSEC_FID:-}"   # gen fiducial sigma sidecar (skim/output/gen_xsec_fid[_mt40].txt,
                                   #  the set of the discriminant -- run_pO_fits.sh picks it);
                                   #  required by SIGMA_POI=1, ignored otherwise
SIMFIT_FLAVS="${SIMFIT_FLAVS:-mu ele}"  # the lepton flavours in the card (2026-09-22):
                                        #  "mu ele" = the grand fit, mu | ele = a per-flavour fit

# ---- flavours in this card: FLAVS in the canonical order mu, ele (so the
# grand card's column order never depends on how the list was spelled) ------
for f in $SIMFIT_FLAVS; do
  case "$f" in
    mu|ele) ;;
    *) echo "[make_pO_simfit_cards] ERROR SIMFIT_FLAVS='$SIMFIT_FLAVS' (allowed: mu, ele)" >&2; exit 2 ;;
  esac
done
FLAVS=""; USE_MU=0; USE_ELE=0
case " $SIMFIT_FLAVS " in *" mu "*)  FLAVS="$FLAVS mu";  USE_MU=1 ;; esac
case " $SIMFIT_FLAVS " in *" ele "*) FLAVS="$FLAVS ele"; USE_ELE=1 ;; esac
FLAVS="${FLAVS# }"
if [ -z "$FLAVS" ]; then echo "[make_pO_simfit_cards] ERROR SIMFIT_FLAVS is empty" >&2; exit 2; fi
# the flavour alternation of every map regex: "(mu|ele)" keeps the grand
# fit's maps byte-identical; a per-flavour card matches its own channels only
if [ "$USE_MU" = 1 ] && [ "$USE_ELE" = 1 ]; then FLAVRE="(mu|ele)"; else FLAVRE="$FLAVS"; fi
# the inputs of the flavours in use must exist (the card only records paths,
# so a missing file would otherwise surface much later, inside text2workspace)
for pair in "$USE_MU|$WMU" "$USE_MU|$ZMU" "$USE_ELE|$WEL" "$USE_ELE|$ZEL"; do
  [ "${pair%%|*}" = 1 ] || continue
  [ -f "${pair#*|}" ] || { echo "[make_pO_simfit_cards] ERROR input not found: ${pair#*|}" >&2; exit 2; }
done

case "$DISC:$QCD_MODE" in
  leppt_mt40:abcdy|met:lnN|met:free) ;;
  *) echo "[make_pO_simfit_cards] ERROR QCD_MODE='$QCD_MODE' for disc '$DISC' (leppt_mt40: abcdy | met: lnN, free;" >&2
     echo "  QCD_MODE=abcd and mtsplit were removed 2026-10-26)" >&2
     exit 2 ;;
esac

# ---- abcdy (2026-10-19): the rapidity groups, from the W inputs' group sidecar
# <W input minus .root>_abcdy.txt (plotting/mtandmet.C writes it ONLY when every
# abcdy object is in the input): '<lab|fb> <group> <bin indices>' lines. The
# sidecars of the flavours in the card must agree (one set of map regexes);
# every binning must place each of the 12 bins in exactly one group.
YSC=""
if [ "$QCD_MODE" = "abcdy" ]; then
  GSUF="_${QCD_MODE}.txt"
  for pair in "mu|$USE_MU|$WMU" "ele|$USE_ELE|$WEL"; do
    fl="${pair%%|*}"; rest="${pair#*|}"
    [ "${rest%%|*}" = 1 ] || continue
    sc="${rest#*|}"; sc="${sc%.root}${GSUF}"
    if [ ! -f "$sc" ]; then
      echo "[make_pO_simfit_cards] ERROR QCD_MODE=$QCD_MODE: group sidecar not found: $sc" >&2
      echo "  (regenerate the inputs with plotting/run_combine_inputs.sh after correction/run_qcd_abcd.sh)" >&2
      exit 2
    fi
    if [ -z "$YSC" ]; then YSC="$sc"
    elif [ "$(grep -v '^#' "$YSC")" != "$(grep -v '^#' "$sc")" ]; then
      echo "[make_pO_simfit_cards] ERROR QCD_MODE=$QCD_MODE: the group sidecars disagree: $YSC vs $sc" >&2; exit 2
    fi
  done
  for Bchk in lab fb; do
    chk=$(awk -v B="$Bchk" '!/^#/ && NF > 0 && $1 == B {
            ng++; if ($2 !~ /^[A-Za-z0-9]+$/ || seen[$2]++) bad = 1
            for (i = 3; i <= NF; i++) { if ($i !~ /^[0-9]+$/ || $i + 0 > 11) bad = 1; c[$i + 0]++ } }
          END { for (y = 0; y < 12; y++) if (c[y] != 1) bad = 1; print (ng > 0 && !bad) ? "ok" : "bad" }' "$YSC")
    if [ "$chk" != "ok" ]; then
      echo "[make_pO_simfit_cards] ERROR QCD_MODE=$QCD_MODE: $YSC: the $Bchk groups do not place each of the 12 bins in exactly one group" >&2
      exit 2
    fi
  done
  echo "[make_pO_simfit_cards] $QCD_MODE rapidity groups from $(basename "$YSC"):"
  grep -v '^#' "$YSC" | sed 's/^/    /'
fi
ygrp_names() {  # $1 = lab|fb -> the group names, sidecar order
  awk -v B="$1" '!/^#/ && NF > 0 && $1 == B {print $2}' "$YSC"
}
ygrp_bins() {   # $1 = lab|fb, $2 = group -> its bin indices, space separated
  awk -v B="$1" -v G="$2" '!/^#/ && $1 == B && $2 == G {for (i = 3; i <= NF; i++) printf "%s%s", (i > 3 ? " " : ""), $i; print ""}' "$YSC"
}
ygrp_of() {     # $1 = lab|fb, $2 = bin index -> its group
  awk -v B="$1" -v Y="$2" '!/^#/ && $1 == B {for (i = 3; i <= NF; i++) if ($i == Y) {print $2; exit}}' "$YSC"
}

# ---- shape systematics: the UNION of the four sidecars, per-INPUT process lists ----
# Sidecar format (plotting/mtandmet.C, dileptonpeak.C): '# comment' lines and
# '<syst> <hist-process> ...' lines, e.g. 'nPDF signal z ztau wtau' (W inputs)
# or 'nPDF signal ztau w wtau' (Z inputs; hist 'signal' = card process 'zsig').
# Since 2026-09-14 a systematic may be listed by ONE flavour's inputs only (the
# muon SFs muID / muIso / muTrig live in the muon W and Z sidecars): its row
# gets `1` on the columns of the inputs that list it and `-` everywhere else,
# so a card never references a shape that is missing -- the per-input lists
# below are the single source of that. The theory (LHE) systematics are still
# expected in all four inputs (a WARN flags anything else).
sidecar_systs() {  # $1 = input .root -> the syst names listed next to it
  local sc="${1%.root}_systs.txt"
  [ -f "$sc" ] || return 0
  awk '!/^#/ && NF>0 {print $1}' "$sc"
}
sidecar_procs() {  # $1 = input .root, $2 = syst -> its histogram-process list
  local sc="${1%.root}_systs.txt"
  [ -f "$sc" ] || return 0
  awk -v s="$2" '!/^#/ && $1==s {for (i=2;i<=NF;i++) printf "%s ", $i; print ""}' "$sc"
}
LHE_THEORY="nPDF qcdScale alphaS"   # the LHE (theory) families; every other listed name is a lepton-SF one,
                                    # except the qcdFF_* fake-factor slopes (group qcdff, 2026-09-29)
LHENAMES=""   # space-separated shape systematics in use ("" = none), sidecar order
# Only the sidecars of the flavours IN this card count (2026-09-22): a
# per-flavour card must not pick up the other flavour's systematics -- muSF in
# an electron card would be a row of '-' only. Historical order muW eleW muZ
# eleZ, so the grand card's row order is unchanged.
used_sidecar_systs() {
  if [ "$USE_MU" = 1 ];  then sidecar_systs "$WMU"; fi
  if [ "$USE_ELE" = 1 ]; then sidecar_systs "$WEL"; fi
  if [ "$USE_MU" = 1 ];  then sidecar_systs "$ZMU"; fi
  if [ "$USE_ELE" = 1 ]; then sidecar_systs "$ZEL"; fi
}
USED_LABELS=""; ALLONES=""   # e.g. "muW eleW muZ eleZ" / "1111" for the grand card
for pair in "$USE_MU|muW" "$USE_ELE|eleW" "$USE_MU|muZ" "$USE_ELE|eleZ"; do
  if [ "${pair%%|*}" = 1 ]; then USED_LABELS="$USED_LABELS ${pair#*|}"; ALLONES="${ALLONES}1"; fi
done
USED_LABELS="${USED_LABELS# }"
if [ "$LHE_SYST" != "off" ]; then
  ALLS=$( used_sidecar_systs | awk 'NF && !seen[$0]++' )
  for s in $ALLS; do
    case ",$LHE_SYST," in *,auto,*|*",$s,"*) ;; *) continue ;; esac
    LHENAMES="$LHENAMES $s"
    have=""
    for pair in "$USE_MU|$WMU" "$USE_ELE|$WEL" "$USE_MU|$ZMU" "$USE_ELE|$ZEL"; do
      [ "${pair%%|*}" = 1 ] || continue
      case " $(sidecar_systs "${pair#*|}" | tr '\n' ' ') " in *" $s "*) have="${have}1" ;; *) have="${have}0" ;; esac
    done
    case " $LHE_THEORY " in
      *" $s "*) [ "$have" = "$ALLONES" ] || echo "[make_pO_simfit_cards] WARN theory syst '$s' listed by inputs ($USED_LABELS) = $have only -> entries only where listed (inputs out of step?)" >&2 ;;
      *) case "$s" in
           qcdFF_*) echo "[make_pO_simfit_cards] fake-factor QCD syst '$s' listed by inputs ($USED_LABELS) = $have -> the qcd column(s) its sidecar names" ;;
           *) echo "[make_pO_simfit_cards] lepton-SF syst '$s' listed by inputs ($USED_LABELS) = $have -> entries on those columns only" ;;
         esac ;;
    esac
  done
  LHENAMES="${LHENAMES# }"
  if [ "$LHE_SYST" != "auto" ]; then
    for s in $(echo "$LHE_SYST" | tr ',' ' '); do
      case " $LHENAMES " in *" $s "*) ;; *)
        echo "[make_pO_simfit_cards] ERROR LHE_SYST='$s' is not listed in any input's sidecar" >&2; exit 2 ;;
      esac
    done
  fi
fi
NLHE=$(echo "$LHENAMES" | wc -w | tr -d ' ')
# per-syst process lists (histogram names) of EACH input, indexed like LHENAMES
LHEP_W_MU=(); LHEP_W_ELE=(); LHEP_Z_MU=(); LHEP_Z_ELE=()
i=0
for s in $LHENAMES; do
  LHEP_W_MU[$i]="$(sidecar_procs "$WMU" "$s")"
  LHEP_W_ELE[$i]="$(sidecar_procs "$WEL" "$s")"
  LHEP_Z_MU[$i]="$(sidecar_procs "$ZMU" "$s")"
  LHEP_Z_ELE[$i]="$(sidecar_procs "$ZEL" "$s")"
  i=$((i+1))
done
# muTrig correlation model (2026-09-14): decided in the analysis repo together
# with the trigger-SF binning (skim/muon_sf.h kTrigBinning -> kMuTrigCorr) and
# recorded in the muon sidecars as '#! muTrig corr <coherent|perbin>' (a '#!'
# line is a directive, not a systematic -- sidecar_systs() skips '#' lines).
# auto follows it (coherent when absent = the inclusive SF); an explicit env
# value overrides for robustness checks.
SF_TRIG_SRC="env"
if [ "$SF_TRIG_CORR" = "auto" ] && [ "$USE_MU" = 0 ]; then
  SF_TRIG_CORR="coherent"; SF_TRIG_SRC="auto: no muon channels in this card"
elif [ "$SF_TRIG_CORR" = "auto" ]; then
  SF_TRIG_CORR=$(awk '$1=="#!" && $2=="muTrig" && $3=="corr" {print $4; exit}' "${WMU%.root}_systs.txt" 2>/dev/null || true)
  if [ -n "$SF_TRIG_CORR" ]; then SF_TRIG_SRC="auto: muon W sidecar directive"
  else SF_TRIG_CORR="coherent"; SF_TRIG_SRC="auto: no directive in the muon W sidecar -> coherent (only matters for a separate muTrig nuisance)"; fi
fi
case "$SF_TRIG_CORR" in
  perbin|coherent) ;;
  *) echo "[make_pO_simfit_cards] ERROR SF_TRIG_CORR='$SF_TRIG_CORR' (auto | perbin | coherent)" >&2; exit 2 ;;
esac
# theory vs lepton-SF vs fake-factor-QCD vs MET-resolution split (the
# --freezeNuisanceGroups groups; metRes* = the MC MET smearing of the W inputs,
# group metres)
THNAMES=""; SFNAMES=""; QFNAMES=""; MRNAMES=""
for s in $LHENAMES; do
  case " $LHE_THEORY " in
    *" $s "*) THNAMES="$THNAMES $s" ;;
    *) case "$s" in
         qcdFF_*)  QFNAMES="$QFNAMES $s" ;;
         metRes*)  MRNAMES="$MRNAMES $s" ;;
         *)        SFNAMES="$SFNAMES $s" ;;
       esac ;;
  esac
done
THNAMES="${THNAMES# }"; SFNAMES="${SFNAMES# }"; QFNAMES="${QFNAMES# }"; MRNAMES="${MRNAMES# }"
# fifth `shapes` token (Combine appends Up/Down to $SYSTEMATIC); empty when no systs
SY=""; [ "$NLHE" -gt 0 ] && SY='_$SYSTEMATIC'

# Append one entry per shape row for a block of process columns:
#   $1 = W | Z | CR (which input's list applies; CR = never), $2 = mu | ele
#   (which flavour's input), $3.. = the HISTOGRAM names of the columns in MP
#   order (zsig -> signal; the SR qcd column -> $QPATH, i.e. qcd_abcdy in abcdy
#   mode -- the electron sidecar lists the qcdFF rows as 'qcdFF_ele_Wp qcd_abcdy@Wp').
# A sidecar process token <hist>@<charge> (2026-09-29, the per-charge
# fake-factor slopes) matches only while CURC -- set by the calling block to
# the channel's charge, empty outside the W channels -- equals <charge>.
# ALIGNMENT RULE: call it in every block that appends to MB/MP/MI/MR, with
# exactly as many names as columns.
CURC=""
lhe_append() {
  local kind="$1" flav="$2"; shift 2
  local i=0 s p e plist
  for s in $LHENAMES; do
    e=""
    plist=""
    if [ "$kind" = "W" ]; then
      if [ "$flav" = "mu" ]; then plist="${LHEP_W_MU[$i]}"; else plist="${LHEP_W_ELE[$i]}"; fi
    elif [ "$kind" = "Z" ]; then
      if [ "$flav" = "mu" ]; then plist="${LHEP_Z_MU[$i]}"; else plist="${LHEP_Z_ELE[$i]}"; fi
    fi
    for p in "$@"; do
      case " $plist " in
        *" $p "*) e="$e 1" ;;
        *" $p@$CURC "*) if [ -n "$CURC" ]; then e="$e 1"; else e="$e -"; fi ;;
        *) e="$e -" ;;
      esac
    done
    SLHE[$i]="${SLHE[$i]}$e"
    i=$((i+1))
  done
}

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
  STTX="ttbar_xsec       lnN"   # ttbar normalization (2026-10-16), the ttbar columns only
  # LHE shape rows, one per systematic (entries appended by lhe_append)
  SLHE=()
  i=0; for s in $LHENAMES; do SLHE[$i]="$(printf '%-17s%s' "$s" 'shape')"; i=$((i+1)); done
  NCH=0

  # ---- 48 W channels ----------------------------------------------------------
  # abcdy mode fits the in-fit template qcd_abcdy (each rapidity group at its
  # own B_g*C40_g/D_g) so the group's formula rateParam (sB*sC/sD, init 1)
  # needs no baked constants; lnN/free (met) the plain qcd.
  QPATH="qcd"; [ "$QCD_MODE" = "abcdy" ] && QPATH="qcd_abcdy"
  for F in $FLAVS; do
    if [ "$F" = "mu" ]; then WF="$WMU"; else WF="$WEL"; fi
    for C in Wp Wm; do
      for iy in $YBINS; do
        R="${C}_${B}_y${iy}"
        CH="${F}_${R}"
        SHAPES="${SHAPES}
shapes data_obs ${CH} ${WF} ${R}/data_obs
shapes signal   ${CH} ${WF} ${R}/signal${SY:+ ${R}/signal$SY}
shapes z        ${CH} ${WF} ${R}/z${SY:+ ${R}/z$SY}
shapes ztau     ${CH} ${WF} ${R}/ztau${SY:+ ${R}/ztau$SY}
shapes wtau     ${CH} ${WF} ${R}/wtau${SY:+ ${R}/wtau$SY}
shapes ttbar    ${CH} ${WF} ${R}/ttbar${SY:+ ${R}/ttbar$SY}
shapes qcd      ${CH} ${WF} ${R}/${QPATH}${SY:+ ${R}/${QPATH}$SY}"
        BINL="$BINL $CH"; OBSL="$OBSL -1"
        MB="$MB $CH $CH $CH $CH $CH $CH"
        MP="$MP signal z ztau wtau ttbar qcd"
        MI="$MI 0 1 2 3 20 4"
        MR="$MR -1 -1 -1 -1 -1 -1"
        # the qcd column is named by its HISTOGRAM ($QPATH: qcd_abcdy in abcdy
        # mode) so a sidecar row on the data-driven template (qcdFF_*) finds it;
        # CURC lets a '<hist>@<charge>' token restrict it to this charge
        CURC="$C"
        lhe_append W "$F" signal z ztau wtau ttbar "$QPATH"
        CURC=""
        if [ "$QCD_MODE" = "free" ]; then
          RP="${RP}
qcd_norm_${CH} rateParam ${CH} qcd 1 [0,10]"
          SQMWP="$SQMWP - - - - - -"; SQMWM="$SQMWM - - - - - -"
          SQEWP="$SQEWP - - - - - -"; SQEWM="$SQEWM - - - - - -"
        else
          # lnN mode: the full ABCD kappa; abcdy mode: the REDUCED residual, ONE
          # row per (flavour, charge) shared by its rapidity groups (the same 4
          # row names, on the SR qcd columns only -- never on the CR qcd, whose
          # yields are measurements the scales float)
          if [ "$QCD_MODE" = "abcdy" ]; then
            K="$QCD_ABCD_LNN_MU"; [ "$F" = "ele" ] && K="$QCD_ABCD_LNN_ELE"
          else
            K="$QCD_LNN_MU"; [ "$F" = "ele" ] && K="$QCD_LNN_ELE"
          fi
          H=" - - - - - ${K}"; N=" - - - - - -"
          case "${F}_${C}" in
            mu_Wp)  SQMWP="$SQMWP$H"; SQMWM="$SQMWM$N"; SQEWP="$SQEWP$N"; SQEWM="$SQEWM$N" ;;
            mu_Wm)  SQMWP="$SQMWP$N"; SQMWM="$SQMWM$H"; SQEWP="$SQEWP$N"; SQEWM="$SQEWM$N" ;;
            ele_Wp) SQMWP="$SQMWP$N"; SQMWM="$SQMWM$N"; SQEWP="$SQEWP$H"; SQEWM="$SQEWM$N" ;;
            ele_Wm) SQMWP="$SQMWP$N"; SQMWM="$SQMWM$N"; SQEWP="$SQEWP$N"; SQEWM="$SQEWM$H" ;;
          esac
        fi
        SLUMI="$SLUMI ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} -"
        STTX="$STTX - - - - ${TTBAR_LNN} -"
        NCH=$((NCH+1))
      done
    done
  done

  # ---- 2 Z channels (w/wtau deliberately have NO scaling parameter: frozen) ---
  # (one in a per-flavour card: that flavour's own peak pins its r_Z)
  for F in $FLAVS; do
    if [ "$F" = "mu" ]; then ZF="$ZMU"; else ZF="$ZEL"; fi
    CH="${F}_Z_incl"
    SHAPES="${SHAPES}
shapes data_obs ${CH} ${ZF} Z_incl/data_obs
shapes zsig     ${CH} ${ZF} Z_incl/signal${SY:+ Z_incl/signal$SY}
shapes w        ${CH} ${ZF} Z_incl/w${SY:+ Z_incl/w$SY}
shapes wtau     ${CH} ${ZF} Z_incl/wtau${SY:+ Z_incl/wtau$SY}
shapes ztau     ${CH} ${ZF} Z_incl/ztau${SY:+ Z_incl/ztau$SY}
shapes ttbar    ${CH} ${ZF} Z_incl/ttbar${SY:+ Z_incl/ttbar$SY}"
    BINL="$BINL $CH"; OBSL="$OBSL -1"
    MB="$MB $CH $CH $CH $CH $CH"
    MP="$MP zsig w wtau ztau ttbar"
    MI="$MI -1 5 3 2 20"
    MR="$MR -1 -1 -1 -1 -1"
    SQMWP="$SQMWP - - - - -"; SQMWM="$SQMWM - - - - -"
    SQEWP="$SQEWP - - - - -"; SQEWM="$SQEWM - - - - -"
    SLUMI="$SLUMI ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN}"
    STTX="$STTX - - - - ${TTBAR_LNN}"
    lhe_append Z "$F" signal w wtau ztau ttbar
    NCH=$((NCH+1))
  done

  # ---- ABCD control-region channels (QCD_MODE=abcdy: 36) ----------------------
  # Counting channels (1-bin templates from plotting/mtandmet.C). The free
  # scales qcd_s{B,C,D} float the QCD counts; the SR qcd is scaled by the
  # formula rateParam (@0*@1/@2) -- the ABCD relation inside the likelihood.
  # ALIGNMENT RULE: every new (channel, process) column appends exactly one
  # entry to MB/MP/MI/MR AND to all six systematics rows in the same block
  # (lumi, ttbar_xsec, the 4 QCD lnN) AND (via lhe_append CR ..., always '-')
  # to every LHE shape row.
  # One set per rapidity group of the sidecar: input dirs <C>_<B>_<g>_CR*,
  # channels <F>_<C>_<g>_CR*, scales qcd_s*_<F>_<C>_<g>, the CRB W parts of the
  # group's bins only, and the formula repeated on each SR channel of the group
  # under ONE name. (Until 2026-10-26 QCD_MODE=abcd used the same block with
  # ONE rapidity-inclusive pseudo-group per (flavour, charge).)
  if [ "$QCD_MODE" = "abcdy" ]; then
    CRGROUPS=$(ygrp_names "$B")
    for F in $FLAVS; do
      if [ "$F" = "mu" ]; then WF="$WMU"; else WF="$WEL"; fi
      for C in Wp Wm; do
       for G in $CRGROUPS; do
        GS="_${G}"; DP="${C}_${B}_${G}"; GBINS=$(ygrp_bins "$B" "$G")
        # --- CRB: iso-pass, m_T<30 (EWK ~10-20%: z/ztau ride r_Z; W floats or freezes) ---
        CH="${F}_${C}${GS}_CRB"
        SHAPES="${SHAPES}
shapes data_obs ${CH} ${WF} ${DP}_CRB/data_obs
shapes qcd      ${CH} ${WF} ${DP}_CRB/qcd
shapes z        ${CH} ${WF} ${DP}_CRB/z
shapes ztau     ${CH} ${WF} ${DP}_CRB/ztau
shapes ttbar    ${CH} ${WF} ${DP}_CRB/ttbar"
        BINL="$BINL $CH"; OBSL="$OBSL -1"
        # ttbar in CRB (2026-10-16): frozen like the SR ttbar -- lumi + ttbar_xsec
        if [ "$QCD_WCR" = "frozen" ]; then
          SHAPES="${SHAPES}
shapes wfix     ${CH} ${WF} ${DP}_CRB/wfix"
          MB="$MB $CH $CH $CH $CH $CH"
          MP="$MP qcd z ztau ttbar wfix"
          MI="$MI 4 1 2 20 6"
          MR="$MR -1 -1 -1 -1 -1"
          SQMWP="$SQMWP - - - - -"; SQMWM="$SQMWM - - - - -"
          SQEWP="$SQEWP - - - - -"; SQEWM="$SQEWM - - - - -"
          SLUMI="$SLUMI - ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN}"
          STTX="$STTX - - - ${TTBAR_LNN} -"
          lhe_append CR "$F" qcd z ztau ttbar wfix
        else
          MB="$MB $CH $CH $CH $CH"
          MP="$MP qcd z ztau ttbar"
          MI="$MI 4 1 2 20"
          MR="$MR -1 -1 -1 -1"
          SQMWP="$SQMWP - - - -"; SQMWM="$SQMWM - - - -"
          SQEWP="$SQEWP - - - -"; SQEWM="$SQEWM - - - -"
          SLUMI="$SLUMI - ${LUMI_LNN} ${LUMI_LNN} ${LUMI_LNN}"
          STTX="$STTX - - - ${TTBAR_LNN}"
          lhe_append CR "$F" qcd z ztau ttbar
          # per-y W content: card process w_y<i> <- histogram w_<B>_y<i>
          # (lab and fb cards MUST wire their own split -- see the header note)
          # process index 8+i for w_y<i> (only the group's bins)
          for iy in $GBINS; do
            PIDX=$((8+iy))
            SHAPES="${SHAPES}
shapes w_y${iy}   ${CH} ${WF} ${DP}_CRB/w_${B}_y${iy}"
            MB="$MB $CH"; MP="$MP w_y${iy}"; MI="$MI $PIDX"; MR="$MR -1"
            SQMWP="$SQMWP -"; SQMWM="$SQMWM -"; SQEWP="$SQEWP -"; SQEWM="$SQEWM -"
            SLUMI="$SLUMI ${LUMI_LNN}"
            STTX="$STTX -"
            lhe_append CR "$F" w_y${iy}
          done
        fi
        NCH=$((NCH+1))
        # --- CRC (anti-iso, m_T>40) and CRD (anti-iso, m_T<30): one frozen ewk ---
        # (ewk = every MC sample, ttbar included since 2026-10-16: lumi only --
        #  the ttbar part is ~0.1-1% of ewk, so no ttbar_xsec entry here)
        for RG in CRC CRD; do
          CH="${F}_${C}${GS}_${RG}"
          SHAPES="${SHAPES}
shapes data_obs ${CH} ${WF} ${DP}_${RG}/data_obs
shapes qcd      ${CH} ${WF} ${DP}_${RG}/qcd
shapes ewk      ${CH} ${WF} ${DP}_${RG}/ewk"
          BINL="$BINL $CH"; OBSL="$OBSL -1"
          MB="$MB $CH $CH"
          MP="$MP qcd ewk"
          MI="$MI 4 7"
          MR="$MR -1 -1"
          SQMWP="$SQMWP - -"; SQMWM="$SQMWM - -"
          SQEWP="$SQEWP - -"; SQEWM="$SQEWM - -"
          SLUMI="$SLUMI - ${LUMI_LNN}"
          STTX="$STTX - -"
          lhe_append CR "$F" qcd ewk
          NCH=$((NCH+1))
        done
        # --- the three free scales + the functional SR multiplier -------------
        # The formula's ONLY args are the scales (init 1 each), and the group's
        # SR template total is B0_g*C40_g/D0_g exactly [x A0ff/A0 for e], so
        # prefit == the ABCD prediction and Asimov closure demands all three
        # scales = 1.
        RP="${RP}
qcd_sB_${F}_${C}${GS} rateParam ${F}_${C}${GS}_CRB qcd 1 [0,10]
qcd_sC_${F}_${C}${GS} rateParam ${F}_${C}${GS}_CRC qcd 1 [0,10]
qcd_sD_${F}_${C}${GS} rateParam ${F}_${C}${GS}_CRD qcd 1 [0,10]"
        FARGS="qcd_sB_${F}_${C}${GS},qcd_sC_${F}_${C}${GS},qcd_sD_${F}_${C}${GS}"
        # one line per SR channel of the group, SAME name + formula: Combine
        # creates the function once and every listed channel shares it
        for iy in $GBINS; do
          RP="${RP}
qcd_abcd_${F}_${C}${GS} rateParam ${F}_${C}_${B}_y${iy} qcd (@0*@1/@2) ${FARGS}"
        done
       done
      done
    done
  fi

  # ---- abcdy with QCD_RATE_SPLIT=group (2026-10-30): the residual lnN rows per
  # (flavour, charge, rapidity group), built from the FINAL column list (MB/MP),
  # so the blocks above need no extra bookkeeping: entry kappa_g on the SR qcd
  # column of every bin of the group, '-' on every other column (other bins,
  # other charges/flavours, every CR qcd -- whose yields are measurements).
  QGROWS=""; QRKTOK=""
  if [ "$QCD_MODE" = "abcdy" ] && [ "$QCD_RATE_SPLIT" = "group" ]; then
    read -r -a MBA <<< "$MB"; read -r -a MPA <<< "$MP"
    for F in $FLAVS; do
      for G in $(ygrp_names "$B"); do
        QRKTOK="${QRKTOK:+$QRKTOK,}${F}_${G}=$(qrate_kappa "$F" "$G")"
      done
      for C in Wp Wm; do
        for G in $(ygrp_names "$B"); do
          K=$(qrate_kappa "$F" "$G")
          GB=" $(ygrp_bins "$B" "$G") "
          ROW="$(printf '%-24s%s' "qcd_rate_${F}_${C}_${G}" 'lnN')"
          NHIT=0
          j=1
          while [ "$j" -lt "${#MBA[@]}" ]; do
            E="-"
            if [ "${MPA[$j]}" = "qcd" ]; then
              case "${MBA[$j]}" in
                "${F}_${C}_${B}_y"*)
                  iy="${MBA[$j]##*_y}"
                  case "$GB" in *" $iy "*) E="$K"; NHIT=$((NHIT+1)) ;; esac ;;
              esac
            fi
            ROW="$ROW $E"
            j=$((j+1))
          done
          if [ "$NHIT" -eq 0 ]; then
            echo "[make_pO_simfit_cards] ERROR qcd_rate_${F}_${C}_${G} has no SR qcd column ($B card)" >&2
            exit 2
          fi
          QGROWS="${QGROWS}
${ROW}"
        done
      done
    done
  fi

  # ---- systematics block: lumi + ttbar_xsec always; the QCD lnN rows unless free ----
  SYST="
${SLUMI}
${STTX}"
  if [ -n "$QGROWS" ]; then
    SYST="${SYST}${QGROWS}"
  elif [ "$QCD_MODE" != "free" ]; then
    # only the rows of the flavours IN the card: in a per-flavour card the
    # other flavour's two rows would be all '-' (a nuisance touching nothing)
    if [ "$USE_MU" = 1 ]; then
      SYST="${SYST}
${SQMWP}
${SQMWM}"
    fi
    if [ "$USE_ELE" = 1 ]; then
      SYST="${SYST}
${SQEWP}
${SQEWM}"
    fi
  fi
  # shape rows (2026-09-07 LHE, 2026-09-14 lepton SFs, 2026-09-29 the
  # fake-factor QCD slopes) + the groups for --freezeNuisanceGroups: `lhe` =
  # the theory ones, `lepsf` = the lepton-SF ones (after the per-bin split
  # below), `qcdff` = the qcdFF_* slopes. A row with NO entry at all is
  # dropped -- a nuisance touching nothing.
  KEPT=""
  i=0
  for s in $LHENAMES; do
    case " ${SLHE[$i]} " in
      *" 1 "*) SYST="${SYST}
${SLHE[$i]}"; KEPT="$KEPT $s" ;;
      *) echo "[make_pO_simfit_cards] ${B}: shape syst '$s' has no entry in this card (QCD_MODE=${QCD_MODE}) -> row dropped" ;;
    esac
    i=$((i+1))
  done
  keep_only() { local out="" s; for s in $1; do case " $KEPT " in *" $s "*) out="$out $s" ;; esac; done; echo "${out# }"; }
  THK=$(keep_only "$THNAMES"); SFK=$(keep_only "$SFNAMES"); QFK=$(keep_only "$QFNAMES"); MRK=$(keep_only "$MRNAMES")
  # muTrig per-bin decorrelation (2026-09-14, SF_TRIG_CORR=perbin): the MB
  # trigger SF's per-y errors are pure statistics, so the ONE coherent template
  # variation is split into 12 independent nuisances muTrig_y<i> -- one per
  # rapidity bin, both charges, every muon W column of that bin -- with
  # Combine's `nuisance edit rename` (process '*', channel regex, full match),
  # which maps each new name back to the SAME <proc>_muTrigUp/Down histograms
  # (DatacardParser.systematicsShapeMap; no duplicated templates). The Z-peak
  # piece stays under `muTrig` (a two-leg factor of 1 +- 1e-4). The edit lines
  # go at the END of the card (after the rateParams), as Combine documents.
  SFEDIT=""; SFALL="$SFK"
  if [ "$SF_TRIG_CORR" = "perbin" ]; then
    case " $KEPT " in *" muTrig "*)
      for iy in $YBINS; do
        SFEDIT="${SFEDIT}
nuisance edit rename * mu_W[pm]_${B}_y${iy} muTrig muTrig_y${iy}"
        SFALL="$SFALL muTrig_y${iy}"
      done ;;
    esac
  fi
  NUISLIST="$THK $SFALL $QFK $MRK"   # every shape nuisance of the card, post-edit names (for the sidecar)
  LHEDESC="no shape systematics (LHE_SYST=${LHE_SYST})"
  if [ "$NLHE" -gt 0 ]; then
    [ -n "$THK" ] && SYST="${SYST}
lhe group = ${THK}"
    [ -n "$SFALL" ] && SYST="${SYST}
lepsf group = ${SFALL}"
    [ -n "$QFK" ] && SYST="${SYST}
qcdff group = ${QFK}"
    [ -n "$MRK" ] && SYST="${SYST}
metres group = ${MRK}"
    LHEDESC="shape systematics -- theory (group lhe): ${THK:-none}; lepton SFs (group lepsf): ${SFALL:-none}; <proc>_<syst>Up/Down from the input files' sidecars, entries only on the columns of the inputs listing them"
    [ -n "$QFK" ] && LHEDESC="${LHEDESC}; fake-factor QCD slopes (group qcdff): ${QFK} on the SR ${QPATH} of their charge"
    [ -n "$MRK" ] && LHEDESC="${LHEDESC}; MC MET resolution (group metres): ${MRK}"
    [ -n "$SFEDIT" ] && LHEDESC="${LHEDESC}; muTrig decorrelated per rapidity bin via nuisance edit rename (SF_TRIG_CORR=perbin)"
  fi
  # the kappas / CR count of the flavours in the card (grand card: the
  # historical text, byte for byte)
  KDA=""; KDL=""
  if [ "$USE_MU" = 1 ];  then KDA="mu ${QCD_ABCD_LNN_MU}"; KDL="mu ${QCD_LNN_MU}"; fi
  if [ "$USE_ELE" = 1 ]; then KDA="${KDA:+$KDA / }ele ${QCD_ABCD_LNN_ELE}"; KDL="${KDL:+$KDL, }ele ${QCD_LNN_ELE}"; fi
  NCR=$((6 * (USE_MU + USE_ELE)))
  if [ "$QCD_MODE" = "free" ]; then
    QDESC="qcd_norm free rateParam per W channel"
  elif [ "$QCD_MODE" = "abcdy" ]; then
    NGRP=$(ygrp_names "$B" | wc -l | tr -d ' ')
    QDESC="QCD in-fit ABCD PER RAPIDITY GROUP (abcdy; $(ygrp_names "$B" | tr '\n' ' ' | sed 's/ $//')): $((NCR * NGRP)) CR channels + one (sB*sC/sD) formula per (flavour, charge, group) on the SR qcd_abcdy; residual lnN ${KDA} shared by the groups; CRB W-part ${QCD_WCR}"
  else
    QDESC="QCD lnN (ABCD) per flavour x charge: ${KDL}"
  fi
  WCHT="48 W channels ({mu,ele} x {Wp,Wm} x y0..11)"
  if [ "$USE_MU" = 1 ] && [ "$USE_ELE" = 1 ]; then
    CARDHEAD="# Auto-generated by make_pO_simfit_cards.sh -- GRAND SIMULTANEOUS FIT, ${B} binning,
# ${DISCLABEL} W discriminant.  ${WCHT} + both
# Z-inclusive peaks in ONE likelihood.  The 25-POI model (r_<C>_y<i> shared mu/e,
# + global r_Z on all DY-related MC) is applied at text2workspace time via
# multiSignalModel with the maps in t2w_maps_simfit_${B}.txt.
# ${QDESC}; lumi lnN ${LUMI_LNN} on all MC (not qcd);
# w/wtau under the Z peaks FROZEN at absolute MC (negligible: <0.1 evt);
# ttbar FROZEN at absolute MC everywhere, ttbar_xsec lnN ${TTBAR_LNN}."
  else
    CARDHEAD="# Auto-generated by make_pO_simfit_cards.sh -- PER-FLAVOUR SIMULTANEOUS FIT (${FLAVS} only),
# ${B} binning, ${DISCLABEL} W discriminant.  24 W channels (${FLAVS} x {Wp,Wm} x
# y0..11) + the ${FLAVS} Z-inclusive peak in ONE likelihood.  The 25-POI model
# (r_<C>_y<i> of this flavour alone, + r_Z on all DY-related MC) is applied at
# text2workspace time via multiSignalModel with the maps in t2w_maps_simfit_${B}.txt.
# ${QDESC}; lumi lnN ${LUMI_LNN} on all MC (not qcd);
# w/wtau under the Z peak FROZEN at absolute MC (negligible: <0.1 evt);
# ttbar FROZEN at absolute MC everywhere, ttbar_xsec lnN ${TTBAR_LNN}."
  fi

  cat > "$CARD" <<EOF
${CARDHEAD}
# ${LHEDESC}
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
------------${SYST}${RP}${SFEDIT}
EOF

  # ---- the model: one 'map=' per line, consumed as --PO map=... at t2w time ---
  # re.match runs each regex against "bin/process"; the trailing /... makes y1
  # safe against y10/y11, and the \$ anchors keep 'z' from matching 'ztau'.
  # FLAVRE = "(mu|ele)" in the grand card (unchanged), "mu" / "ele" in a
  # per-flavour one, where the other flavour's channels do not exist.
  : > "$MAPS"
  for C in Wp Wm; do
    for iy in $YBINS; do
      echo "map=${FLAVRE}_${C}_${B}_y${iy}/(signal|wtau)\$:r_${C}_y${iy}[1,0,10]" >> "$MAPS"
    done
  done
  if [ "$QCD_MODE" = "abcdy" ] && [ "$QCD_WCR" != "frozen" ]; then
    # CRB W content floats with the SAME POIs (the in-fit EWK subtraction), on
    # the GROUP CRB channel that holds w_y<i>. The process name ENDS in the y
    # index, so the \$ anchor alone keeps w_y1 from matching w_y10.
    for C in Wp Wm; do
      for iy in $YBINS; do
        echo "map=${FLAVRE}_${C}_$(ygrp_of "$B" "$iy")_CRB/w_y${iy}\$:r_${C}_y${iy}[1,0,10]" >> "$MAPS"
      done
    done
  fi
  echo "map=.*/(z|ztau|zsig)\$:r_Z[1,0,10]" >> "$MAPS"

  echo "[make_pO_simfit_cards] wrote $(basename "$CARD") (${NCH} channels) + $(basename "$MAPS") ($(wc -l < "$MAPS" | tr -d ' ') maps)"

  if [ "$SIGMA_POI" = "1" ]; then
    write_sigma_maps "$B"
  else
    # never leave a stale sigma map next to freshly generated cards: it would
    # look like a valid model definition for a card it was not built from
    rm -f "$OUTDIR/t2w_maps_simfit_${B}_sigma.txt"
  fi
  return 0
}

# ---------------------------------------------------------------------------
# SIGMA_POI=1 (2026-09-15): a SECOND map file per variant that reparametrizes
# the SAME datacard so the rapidity-inclusive W cross section is itself a POI.
#
#   sigma_W = Sum_i r_i sigma_gen-fid,i   is a fixed linear combination of the
#   24 per-bin POIs, so promoting it is a pure change of variables:
#
#       POIs   sigmaW, a_<C>_y<i>  (23 free; ONE bin is the anchor, a == 1)
#       funcs  Wnorm = G_anchor + Sum_j a_j G_j
#              r_i   = sigmaW * a_i / Wnorm      (anchor: r = sigmaW / Wnorm)
#
#   => Sum_i r_i G_i == sigmaW identically (verified to 6e-16 in RooFit), and
#   a_i == 1 with sigmaW = Sum_i G_i reproduces r_i == 1, so the Asimov closure
#   point is unchanged. The profile likelihood is invariant under a bijective
#   reparametrization, so the best fit and the PROFILED interval of sigma_W are
#   the exact ones -- which is the whole point: the Hesse/covariance ellipse
#   that plotting/xsec_contour.C draws from the nominal fit is the GAUSSIAN
#   approximation to this, and this is what validates it (and what a genuine
#   `MultiDimFit --algo grid -P sigmaW -P r_Z` contour needs).
#
# Chosen over the "redefine one bin as sigma_total - sum(others)" subtraction
# because the subtraction sends that bin's yield NEGATIVE once sigmaW is scanned
# ~1.4 sigma down (the largest bin carries only ~5% of the total), which is well
# inside a 95% CL contour; the ratio form keeps every r_i > 0 by construction.
#
# The per-bin gen fiducial cross sections G_i come from the analysis repo's
# skim/output/gen_xsec_fid[_mt40].txt (written by skim/gen_xsec.C; the _mt40
# set = m_T > 40, for leppt_mt40, since 2026-10-13) -- they must be
# baked in because sigma_W is that specific combination. NOTE they are the
# NOMINAL-PDF values, the same approximation plotting/xsec_fiducial.C makes.
#
# The DATACARD is untouched: POIs live in the map, so one card serves both
# parametrizations and nothing about the nominal fit changes.
# ---------------------------------------------------------------------------
write_sigma_maps() {
  B="$1"
  SMAPS="$OUTDIR/t2w_maps_simfit_${B}_sigma.txt"
  if [ ! -f "$GEN_XSEC_FID" ]; then
    echo "[make_pO_simfit_cards][ERROR] SIGMA_POI=1 but gen sidecar not found: $GEN_XSEC_FID"
    echo "                              (run skim/gen_xsec.C, then sync_lxplus.sh upload)"
    return 1
  fi
  # abcdy (2026-10-19): "<bin>:<group> ..." of this binning, so the CRB pattern
  # names the group channel holding w_y<bin>; empty in every other mode
  local GM="" iyg
  if [ "$QCD_MODE" = "abcdy" ]; then
    for iyg in $YBINS; do GM="$GM $iyg:$(ygrp_of "$B" "$iyg")"; done
  fi
  awk -v B="$B" -v QMODE="$QCD_MODE" -v QWCR="$QCD_WCR" -v FRE="$FLAVRE" -v GM="$GM" '
    BEGIN { ng = split(GM, gma, " "); for (k = 1; k <= ng; k++) { split(gma[k], kv, ":"); grp[kv[1]] = kv[2] } }
    /^#/ { next }
    $1 == B { key = $2 "_y" $3; g[key] = $4 + 0; order[++n] = key; tot += $4 }
    END {
      if (n != 24) { printf "[ERROR] %d (not 24) %s rows in the gen sidecar\n", n, B > "/dev/stderr"; exit 2 }
      anchor = order[n]
      printf "map=__poi_decl__/__poi_decl__:sigmaW[%.6f,0,%.6f]\n", tot, 4 * tot
      terms = ""; args = ""; k = 0
      for (i = 1; i <= n; i++) {
        if (order[i] == anchor) continue
        printf "map=__poi_decl__/__poi_decl__:a_%s[1,0,10]\n", order[i]
        terms = terms sprintf("%s@%d*%.6f", (k ? "+" : ""), k, g[order[i]])
        args  = args  sprintf("%sa_%s", (k ? "," : ""), order[i])
        k++
      }
      # Wnorm must be declared BEFORE the r_* expressions that use it:
      # MultiSignalModel runs factory statements in map-file order.
      printf "map=__fn_decl__/__fn_decl__:Wnorm=expr;;Wnorm(\"%.6f+%s\",%s)\n", g[anchor], terms, args
      for (i = 1; i <= n; i++) {
        key = order[i]; split(key, p, "_y"); C = p[1]; iy = p[2]
        # ONE line per POI carrying BOTH patterns (multiSignalModel splits the
        # map list on commas): re-listing a factory-defined name on a second
        # plain map= line would push it into the POI set as a bare variable and
        # collide with the RooFormulaVar.
        pat = sprintf("%s_%s_%s_y%s/(signal|wtau)$", FRE, C, B, iy)
        if (QMODE == "abcdy" && QWCR != "frozen")
          pat = pat sprintf(",%s_%s_%s_CRB/w_y%s$", FRE, C, grp[iy], iy)
        if (key == anchor)
          printf "map=%s:r_%s=expr;;r_%s(\"@0/@1\",sigmaW,Wnorm)\n", pat, key, key
        else
          printf "map=%s:r_%s=expr;;r_%s(\"@0*@1/@2\",sigmaW,a_%s,Wnorm)\n", pat, key, key, key
      }
      print "map=.*/(z|ztau|zsig)$:r_Z[1,0,10]"
      printf "[make_pO_simfit_cards] %s: sigmaW anchored on %s, Sum_i sigma_gen,i = %.4f nb\n",
             B, anchor, tot > "/dev/stderr"
    }' "$GEN_XSEC_FID" > "$SMAPS" || { echo "[make_pO_simfit_cards][ERROR] sigma-map generation failed for $B"; rm -f "$SMAPS"; return 1; }
  echo "[make_pO_simfit_cards] wrote $(basename "$SMAPS") ($(wc -l < "$SMAPS" | tr -d ' ') maps, sigmaW + 23 a + r_Z = 25 POIs)"
}

gen_simfit_card lab
gen_simfit_card fb

# Sidecar recording the lnN kappas THESE cards were built with -- read back by
# run_pO_fits.sh at extraction time (single source even if env changes between
# card generation and extraction).  kQcd* = 0 means "free rateParam mode":
# extract_pO_simfit.C then reads the per-channel qcd_norm params (legacy path).
# qcdMode (2026-08-23) makes the mode explicit; sidecars WITHOUT the line are
# legacy (mode inferred: kQcd > 0 -> lnN, = 0 -> free). In abcdy mode the kQcd*
# entries hold the REDUCED kappas.
KQM="$QCD_LNN_MU"; KQE="$QCD_LNN_ELE"
if [ "$QCD_MODE" = "free" ]; then KQM=0; KQE=0; fi
if [ "$QCD_MODE" = "abcdy" ]; then KQM="$QCD_ABCD_LNN_MU"; KQE="$QCD_ABCD_LNN_ELE"; fi
# lheSysts (2026-09-07): ALL shape nuisances in the cards, comma-joined, by
# their POST-EDIT names (theory + lepton SFs incl. the per-bin muTrig_y<i> +
# since 2026-09-29 the qcdFF_* fake-factor slopes; "none" when off/absent)
# -- the extractor reports their pulls and includes them in the Asimov
# closure; postfit_incl.C switches to shapes_fit_s on it.
# The key keeps its 2026-09-07 name for the readers; sfTrigCorr records the
# muTrig correlation choice.
LHELIST=$(echo "$NUISLIST" | tr -s ' ' ',' | sed 's/^,//; s/,$//'); LHELIST="${LHELIST:-none}"
# flavours (2026-09-22): the lepton flavours of THESE cards, comma-joined so the
# file stays one key + one token per line (postfit_incl.C reads it as pairs).
# The extractor is told the same list by run_pO_fits.sh and checks it; a
# sidecar without the line predates per-flavour fits = mu,ele.
cat > "$OUTDIR/qcd_lnn_kappas.txt" <<EOF
kQcdMu $KQM
kQcdEle $KQE
kLumi $LUMI_LNN
kTtbar $TTBAR_LNN
qcdMode $QCD_MODE
lheSysts $LHELIST
sfTrigCorr $SF_TRIG_CORR
flavours $(echo "$FLAVS" | tr ' ' ',')
EOF
# qcdGroups (abcdy only, 2026-10-19): the rapidity groups THESE cards were built
# with, ONE token (postfit_incl.C reads the file as key/value pairs):
#   lab:fwd=10,11/bwd=0,1/cen=2,...,9;fb:fwd=...
# run_pO_fits.sh passes it to the extractor (per-group scales + multipliers).
# Absent in lnN/free mode (the met fit).
if [ "$QCD_MODE" = "abcdy" ]; then
  QGTOK=""
  for Bq in lab fb; do
    part=""
    for g in $(ygrp_names "$Bq"); do part="${part:+$part/}$g=$(ygrp_bins "$Bq" "$g" | tr ' ' ',')"; done
    QGTOK="${QGTOK:+$QGTOK;}$Bq:$part"
  done
  echo "qcdGroups $QGTOK" >> "$OUTDIR/qcd_lnn_kappas.txt"
fi
# qcdRateGroups (abcdy + QCD_RATE_SPLIT=group, 2026-10-30): the per-group
# residual kappas THESE cards carry, ONE token "<flav>_<group>=<kappa>,...";
# run_pO_fits.sh passes it to the extractor (it then reads the per-group
# qcd_rate_<F>_<C>_<g> pulls). Absent = one shared row per (flavour, charge).
if [ -n "$QRKTOK" ]; then
  echo "qcdRateGroups $QRKTOK" >> "$OUTDIR/qcd_lnn_kappas.txt"
fi

echo "[make_pO_simfit_cards] done -> ${OUTDIR} (flavours: ${FLAVS}; W discriminant: ${DISCLABEL})"
KMSG=""
if [ "$USE_MU" = 1 ];  then KMSG="mu ${KQM}"; fi
if [ "$USE_ELE" = 1 ]; then KMSG="${KMSG:+$KMSG / }ele ${KQE}"; fi
[ -n "$QRKTOK" ] && KMSG="residual lnN per rapidity group: ${QRKTOK}"
echo "[make_pO_simfit_cards] QCD mode: ${QCD_MODE} (${KMSG}); lumi lnN ${LUMI_LNN}; ttbar_xsec lnN ${TTBAR_LNN}"
echo "[make_pO_simfit_cards] shape systematics: theory [${THK:-none}] lepton SF [${SFALL:-none}] fake-factor QCD [${QFK:-none}] (LHE_SYST=${LHE_SYST}, SF_TRIG_CORR=${SF_TRIG_CORR} from ${SF_TRIG_SRC})"
