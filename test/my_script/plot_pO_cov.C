// =============================================================================
// plot_pO_cov.C -- correlation-matrix plots for the grand simultaneous fit
// (simfit).  Two matrices per binning variant:
//
//   (1) corr_params_simfit_<B>: the FULL floating-parameter correlation matrix
//       from fit_s in fitDiagnostics_simfit_<B>.root -- POIs + constrained
//       nuisances (2026-08-17 lnN model: 24 r + r_Z + 4 qcd_rate + lumi = 30;
//       in QCD_MODE=free cards the 48 qcd_norm rateParams appear instead).
//       Skipped with a WARN when the fitDiagnostics file is absent (e.g. when
//       re-running locally on a downloaded tree, where fits/ is not synced).
//
//   (2) corr_yield_<B>: the 24x24 fitted-YIELD correlation, converted from
//       h_cov_yield (lab) / h_cov_yield_FB (fb) in summary/comb_fitted_yields
//       .root -- THE matrix the downstream error propagation (charge_asym,
//       FBratio, xsec sums) actually uses.  This one works locally too.
//
// Driven by test/run_pO_impacts.sh; standalone:
//   root -b -q 'plot_pO_cov.C("<fitDiagnostics>","<comb_fitted_yields.root>","lab","<outDir>")'
// =============================================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include <cmath>
#include <iostream>

namespace {

// one canvas cosmetic for every matrix: diverging palette, z locked to [-1,1]
void DrawCorr(TH2 *h, const TString &title, const TString &outStem)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(87); // kTemperatureMap: blue -1 .. white 0 .. red +1
  h->SetTitle(title);
  h->GetZaxis()->SetRangeUser(-1.0, 1.0);
  const double lbl = (h->GetNbinsX() > 30) ? 0.012 : 0.02;
  h->GetXaxis()->SetLabelSize(lbl);
  h->GetYaxis()->SetLabelSize(lbl);
  h->GetXaxis()->LabelsOption("v");
  TCanvas c("c_corr", "c_corr", 1000, 900);
  c.SetLeftMargin(0.16);
  c.SetBottomMargin(0.16);
  c.SetRightMargin(0.13);
  c.SetTopMargin(0.08);
  h->Draw("COLZ");
  c.SaveAs(outStem + ".png");
  c.SaveAs(outStem + ".pdf");
}

} // namespace

void plot_pO_cov(const char *fitDiag,    // fits/simfit_<B>/fitDiagnostics_simfit_<B>.root
                 const char *yieldsRoot, // summary/comb_fitted_yields.root
                 const char *variant,    // "lab" | "fb"
                 const char *outDir)
{
  gSystem->mkdir(outDir, kTRUE);
  const TString B(variant);
  bool didAny = false;

  // ---- (1) full parameter correlation matrix from fit_s ---------------------
  if (!gSystem->AccessPathName(fitDiag)) {
    TFile *fd = TFile::Open(fitDiag, "READ");
    RooFitResult *fr = (fd && !fd->IsZombie()) ? (RooFitResult *)fd->Get("fit_s") : nullptr;
    if (fr) {
      TH2 *hp = fr->correlationHist(Form("corr_params_%s", variant));
      hp->SetDirectory(nullptr);
      DrawCorr(hp, Form("Floating-parameter correlation (simfit %s);;", variant),
               TString(outDir) + "/corr_params_simfit_" + B);
      std::cout << "[plot_pO_cov] wrote corr_params_simfit_" << B << " ("
                << hp->GetNbinsX() << " parameters)\n";
      didAny = true;
    } else {
      std::cerr << "[plot_pO_cov] WARN no fit_s in " << fitDiag << "\n";
    }
    if (fd) { fd->Close(); delete fd; }
  } else {
    std::cerr << "[plot_pO_cov] WARN fitDiagnostics absent (" << fitDiag
              << ") -- parameter matrix skipped (fits/ not synced locally?)\n";
  }

  // ---- (2) fitted-yield correlation from h_cov_yield[_FB] -------------------
  const char *covName = (B == "fb") ? "h_cov_yield_FB" : "h_cov_yield";
  TFile *fy = TFile::Open(yieldsRoot, "READ");
  TH2D *hc = (fy && !fy->IsZombie()) ? (TH2D *)fy->Get(covName) : nullptr;
  if (hc) {
    // clone-then-convert: never mutate the file-owned histogram
    TH2D *hr = (TH2D *)hc->Clone(Form("corr_yield_%s", variant));
    hr->SetDirectory(nullptr);
    const int n = hc->GetNbinsX();
    for (int i = 1; i <= n; ++i)
      for (int j = 1; j <= n; ++j) {
        const double vii = hc->GetBinContent(i, i), vjj = hc->GetBinContent(j, j);
        const double d = (vii > 0 && vjj > 0) ? hc->GetBinContent(i, j) / std::sqrt(vii * vjj) : 0.0;
        hr->SetBinContent(i, j, d);
      }
    DrawCorr(hr, Form("Fitted-yield correlation, %s (r covariance incl. lnN nuisances);;", covName),
             TString(outDir) + "/corr_yield_" + B);
    std::cout << "[plot_pO_cov] wrote corr_yield_" << B << " (" << n << " bins)\n";
    didAny = true;
  } else {
    std::cerr << "[plot_pO_cov] WARN no " << covName << " in " << yieldsRoot << "\n";
  }
  if (fy) { fy->Close(); delete fy; }

  if (!didAny) std::cerr << "[plot_pO_cov] ERROR: nothing plotted\n";
}
