#!/usr/bin/env bash
#
# End-to-end Combine pipeline for one channel.
#
# Usage:
#   ./run_fit.sh W       # W -> mu nu
#   ./run_fit.sh Zmu     # Z -> mu mu
#   ./run_fit.sh Ze      # Z -> e e

set -euo pipefail

CH="${1:-}"
if [[ -z "$CH" ]]; then
    echo "usage: $0 <W|Zmu|Ze>"
    exit 1
fi

case "$CH" in
    W)
        INPUT_MACRO='./my_script/make_combine_input.C'
        DATACARD='./my_script/testdatacard_inclusive.txt'
        WORKSPACE='./my_script/W_inclusive_workspace.root'
        FIT_NAME='_w_inclusivefit'
        PLOT_MACRO='./my_script/draw_postfit_inclusive.C'
        ;;
    Zmu)
        INPUT_MACRO='./my_script/make_combine_input_Z.C(0)'
        DATACARD='./my_script/testdatacard_Zmumu.txt'
        WORKSPACE='./my_script/Z_mumu_workspace.root'
        FIT_NAME='_z_mumufit'
        PLOT_MACRO='./my_script/draw_postfit_Zmumu.C'
        ;;
    Ze)
        INPUT_MACRO='./my_script/make_combine_input_Z.C(1)'
        DATACARD='./my_script/testdatacard_Zee.txt'
        WORKSPACE='./my_script/Z_ee_workspace.root'
        FIT_NAME='_z_eefit'
        PLOT_MACRO='./my_script/draw_postfit_Zee.C'
        ;;
    *)
        echo "unknown channel: $CH (use W, Zmu, or Ze)"
        exit 1
        ;;
esac

# Preflight
command -v text2workspace.py >/dev/null || { echo "ERROR: did you cmsenv?"; exit 1; }
command -v combine           >/dev/null || { echo "ERROR: did you cmsenv?"; exit 1; }

echo ""
echo "==== [$CH] 1/4  build combine input ===="
root -b -q "$INPUT_MACRO"

echo ""
echo "==== [$CH] 2/4  text2workspace ===="
text2workspace.py "$DATACARD" -o "$WORKSPACE"

echo ""
echo "==== [$CH] 3/4  FitDiagnostics ===="
combine -M FitDiagnostics "$WORKSPACE" \
        --saveShapes --saveWithUncertainties \
        -n "$FIT_NAME"

echo ""
echo "==== [$CH] 4/4  draw postfit ===="
mkdir -p ./plots_postfit
root -b -q "$PLOT_MACRO"

echo ""
echo "Done [$CH]:"
echo "  fit  = fitDiagnostics${FIT_NAME}.root"
echo "  plot = ./plots_postfit/${CH}_postfit.{png,pdf}"