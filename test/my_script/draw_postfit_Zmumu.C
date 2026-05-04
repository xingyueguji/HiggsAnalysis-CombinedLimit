#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "plotting_helper.C"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

// ----------------------------------------------------------------
//  Helpers — identical to draw_postfit_inclusive.C
//  (GetFitPar, RemapToReferenceBinning, MakePullHist,
//   SaveNicePlot1D_WithBkg_Postfit_Pull)
//  Copy them verbatim from your W file.
// ----------------------------------------------------------------
static double GetFitPar(RooFitResult* fr, const char* name, double& err) { /* ... */ }
static TH1*   RemapToReferenceBinning(TH1* hIn, TH1* hRef, const char* newName) { /* ... */ }
static TH1D*  MakePullHist(TH1* hData, TH1* hPost, const char* name = "h_pull") { /* ... */ }
static void   SaveNicePlot1D_WithBkg_Postfit_Pull(/* same signature, same body */) { /* ... */ }

// ----------------------------------------------------------------
//  Main
// ----------------------------------------------------------------
void draw_postfit_Zmumu() {
    const std::string fitFile  = "fitDiagnostics_z_mumufit.root";
    const std::string dataFile = "combine_input_Zmumu.root";
    const std::string outDir   = "./plots_postfit";

    gSystem->mkdir(outDir.c_str(), kTRUE);

    TFile* fFit  = TFile::Open(fitFile.c_str(),  "READ");
    if (!fFit  || fFit ->IsZombie()) { std::cerr << "[ERROR] cannot open " << fitFile  << "\n"; return; }
    TFile* fData = TFile::Open(dataFile.c_str(), "READ");
    if (!fData || fData->IsZombie()) { std::cerr << "[ERROR] cannot open " << dataFile << "\n"; return; }

    // -------- data --------
    TH1D* hData = dynamic_cast<TH1D*>(fData->Get("data_obs"));
    if (!hData) { std::cerr << "[ERROR] no data_obs in " << dataFile << "\n"; return; }
    hData = dynamic_cast<TH1D*>(hData->Clone("data_obs_plot"));
    hData->SetDirectory(nullptr);

    // -------- postfit shapes (s+b fit) --------
    TH1* hSignal_raw = dynamic_cast<TH1*>(fFit->Get("shapes_fit_s/ch1/signal"));
    TH1* hW_raw      = dynamic_cast<TH1*>(fFit->Get("shapes_fit_s/ch1/w"));
    TH1* hWtau_raw   = dynamic_cast<TH1*>(fFit->Get("shapes_fit_s/ch1/wtau"));
    TH1* hZtau_raw   = dynamic_cast<TH1*>(fFit->Get("shapes_fit_s/ch1/ztau"));
    TH1* hTotal_raw  = dynamic_cast<TH1*>(fFit->Get("shapes_fit_s/ch1/total"));
    if (!hSignal_raw || !hW_raw || !hWtau_raw || !hZtau_raw || !hTotal_raw) {
        std::cerr << "[ERROR] Missing one or more postfit shapes in " << fitFile << "\n"; return;
    }

    TH1* hSignal = RemapToReferenceBinning(hSignal_raw, hData, "signal_postfit");
    TH1* hW      = RemapToReferenceBinning(hW_raw,      hData, "w_postfit");
    TH1* hWtau   = RemapToReferenceBinning(hWtau_raw,   hData, "wtau_postfit");
    TH1* hZtau   = RemapToReferenceBinning(hZtau_raw,   hData, "ztau_postfit");
    TH1* hTotal  = RemapToReferenceBinning(hTotal_raw,  hData, "total_postfit");
    for (TH1* h : {hSignal, hW, hWtau, hZtau, hTotal}) if (h) h->SetDirectory(nullptr);

    // -------- fit parameters (rateParam version) --------
    RooFitResult* fit_s = dynamic_cast<RooFitResult*>(fFit->Get("fit_s"));
    if (!fit_s) { std::cerr << "[ERROR] Cannot find fit_s\n"; return; }

    double er=0, ew=0, ewt=0, ezt=0;
    double r         = GetFitPar(fit_s, "r",         er);
    double w_rate    = GetFitPar(fit_s, "w_rate",    ew);
    double wtau_rate = GetFitPar(fit_s, "wtau_rate", ewt);
    double ztau_rate = GetFitPar(fit_s, "ztau_rate", ezt);

    // -------- style --------
    PlotStyle ps;
    ps.drawOpt   = "E";
    ps.showStats = false;
    ps.logy      = true;
    ps.boxX1 = 0.2; ps.boxX2 = 0.5; ps.boxY1 = 0.49; ps.boxY2 = 0.78;

    PlotTuner commonTuner = [&](TCanvas* c, TH1* h) {
        (void)c;
        if (!h) return;
        if (ps.logy) h->SetMinimum(1.0);
        double ymax = hData->GetMaximum();
        if (hTotal) ymax = std::max(ymax, hTotal->GetMaximum());
        h->SetMaximum(10.25 * ymax);
    };

    std::vector<TH1*>         bkgs  = { hSignal, hW, hWtau, hZtau };
    std::vector<std::string>  names = { "Signal", "W^{+}/W^{-}", "W^{+}/W^{-} #tau", "DY #tau" };

    std::vector<std::string> box = {
        Form("Data: %.0f",                hData ->Integral(1, hData ->GetNbinsX())),
        Form("Postfit total: %.1f",       hTotal->Integral(1, hTotal->GetNbinsX())),
        Form("r = %.3f #pm %.3f",         r, er),
        Form("w rate = %.3f #pm %.3f",    w_rate,    ew),
        Form("wtau rate = %.3f #pm %.3f", wtau_rate, ewt),
        Form("ztau rate = %.3f #pm %.3f", ztau_rate, ezt),
    };

    SaveNicePlot1D_WithBkg_Postfit_Pull(
        hData, bkgs, names, hTotal,
        outDir + "/Zmumu_postfit",
        "m_{#mu#mu} (GeV)",
        "Events / 1.0 GeV",          // adjust if your binning differs
        "",
        "Z #rightarrow #mu #mu",
        "inclusive postfit",
        box, ps, commonTuner);

    fFit ->Close();
    fData->Close();
}