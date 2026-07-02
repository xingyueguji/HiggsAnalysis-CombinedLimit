// =============================================================================
// draw_postfit_pO.C  --  generic postfit data/MC plot for ONE pO fit region,
// styled IDENTICALLY to the analysis-repo plots (plotting/mtandmet.C): it routes
// through the same SaveNicePlot1D_WithBkg (translucent alpha=0.65 fills,
// color-matched outlines width 3, CMS_lumi) from the synced plotting_helper.C.
//
// It reads the POSTFIT shapes (shapes_fit_s/<fitChannel>/<process>) from the
// FitDiagnostics output and remaps them onto the physical x-axis (MET or mass)
// taken from the input file's data_obs, then stacks signal + backgrounds with
// data overlaid.  The stack sum IS the postfit total.
//
// Usage (under cmsenv):
//   root -b -q 'draw_postfit_pO.C("<fitDiag.root>","<fitChannel>","<input.root>",
//                                 "<inputRegion>","<outNoExt>","PF MET (GeV)",
//                                 "Events / 2.0 GeV","W #rightarrow #mu #nu",
//                                 "Wp lab y3 (postfit)", true)'
// =============================================================================
#include "plotting_helper.C" // synced from analysis repo (alpha=0.65 + width-3 stack)
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>

// Poisson (Baker-Cousins) chi2 between data and the postfit total: well-defined
// at low/zero counts (unlike Pearson). Empty (no data, no model) bins are
// skipped; nbinsUsed feeds ndf = nbinsUsed - n_floating_params.
static double BakerCousinsChi2(TH1 *d, TH1 *t, int &nbinsUsed) {
  nbinsUsed = 0;
  double chi2 = 0.0;
  const int n = std::min(d->GetNbinsX(), t->GetNbinsX());
  for (int i = 1; i <= n; ++i) {
    double di = d->GetBinContent(i);
    double ti = t->GetBinContent(i);
    if (ti <= 1e-9 && di <= 0.0) continue; // no information in this bin
    if (ti < 1e-9) ti = 1e-9;              // floor to keep the log finite
    double term = ti - di;
    if (di > 0.0) term += di * std::log(di / ti);
    chi2 += 2.0 * term;
    ++nbinsUsed;
  }
  return chi2;
}

static TH1D *RemapToRef(TH1 *ref, TH1 *src, const char *nm) {
  if (!ref || !src) return 0;
  TH1D *h = (TH1D *)ref->Clone(nm);
  h->SetDirectory(0);
  h->Reset();
  int n = std::min(src->GetNbinsX(), h->GetNbinsX());
  if (src->GetNbinsX() != h->GetNbinsX())
    std::cerr << "[WARN] " << nm << ": postfit nbins " << src->GetNbinsX()
              << " != ref " << h->GetNbinsX() << " (copying " << n << ")\n";
  for (int b = 1; b <= n; ++b) {
    h->SetBinContent(b, src->GetBinContent(b));
    h->SetBinError(b, src->GetBinError(b));
  }
  return h;
}

void draw_postfit_pO(const char *fitDiagFile,
                     const char *fitChannel,   // folder under shapes_fit_s
                     const char *inputFile,    // structured combine input
                     const char *inputRegion,  // TDirectory holding data_obs (physical axis)
                     const char *outNoExt,
                     const char *xTitle,
                     const char *yTitle,
                     const char *subTitle1,
                     const char *subTitle2,
                     bool isW = true)
{
  TFile *ff = TFile::Open(fitDiagFile, "READ");
  if (!ff || ff->IsZombie()) { std::cerr << "[ERROR] cannot open " << fitDiagFile << "\n"; return; }
  TFile *fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) { std::cerr << "[ERROR] cannot open " << inputFile << "\n"; return; }

  TH1 *ref = (TH1 *)fin->Get(TString::Format("%s/data_obs", inputRegion));
  if (!ref) { std::cerr << "[ERROR] no " << inputRegion << "/data_obs in " << inputFile << "\n"; return; }
  TH1D *hData = (TH1D *)ref->Clone("postfit_data");
  hData->SetDirectory(0);

  // process list + legend labels (stack order matches mtandmet.C)
  std::vector<std::string> procs, labels;
  if (isW) {
    procs  = {"signal", "z", "ztau", "wtau", "qcd"};
    labels = {"W signal", "DY", "DY #tau", "W #tau", "QCD (ABCD)"};
  } else {
    procs  = {"signal", "w", "wtau", "ztau"};
    labels = {"Z signal", "W^{+}/W^{-}", "W^{+}/W^{-} #tau", "DY #tau"};
  }

  std::vector<TH1 *> bkgs;
  std::vector<std::string> names;
  for (size_t i = 0; i < procs.size(); ++i) {
    TH1 *src = (TH1 *)ff->Get(TString::Format("shapes_fit_s/%s/%s", fitChannel, procs[i].c_str()));
    if (!src && i == 0) // signal slot: the simultaneous W+Z fit renames the Z signal 'zsig'
      src = (TH1 *)ff->Get(TString::Format("shapes_fit_s/%s/zsig", fitChannel));
    TH1D *h = RemapToRef(ref, src, TString::Format("postfit_%s", procs[i].c_str()));
    if (h) { bkgs.push_back(h); names.push_back(labels[i]); }
  }
  TH1 *tot = (TH1 *)ff->Get(TString::Format("shapes_fit_s/%s/total", fitChannel));
  TH1D *hTot = RemapToRef(ref, tot, "postfit_total");

  // ---- fit-quality diagnostics shown on the plot ----
  std::vector<std::string> box = {
      Form("Data: %.0f", hData->Integral()),
      Form("Postfit total: %.0f", hTot ? hTot->Integral() : 0.0)};

  RooFitResult *fr = (RooFitResult *)ff->Get("fit_s");
  if (hTot) {
    int nUsed = 0;
    const double chi2 = BakerCousinsChi2(hData, hTot, nUsed); // Poisson GoF
    const int nfloat = fr ? fr->floatParsFinal().getSize() : 0;
    int ndf = nUsed - nfloat;
    if (ndf < 1) ndf = (nUsed > 0 ? nUsed : 1);
    box.push_back(Form("#chi^{2}/ndf = %.2f, p = %.2f", chi2 / ndf, TMath::Prob(chi2, ndf)));
  }
  if (fr) {
    RooRealVar *rv = (RooRealVar *)fr->floatParsFinal().find("r");
    const bool bad = (fr->status() != 0 || fr->covQual() < 3); // not converged / bad covariance
    if (rv)
      box.push_back(bad
        ? Form("r = %.3f #pm %.3f  #color[2]{(status %d, covQ %d)}",
               rv->getVal(), rv->getError(), fr->status(), fr->covQual())
        : Form("r = %.3f #pm %.3f", rv->getVal(), rv->getError()));
    else if (bad)
      box.push_back(Form("#color[2]{fit status %d, covQ %d}", fr->status(), fr->covQual()));
    // Two-parameter model diagnostics: whichever of these float in this fit
    // (dy_norm in W + simultaneous cards, w_norm in the standalone Z card,
    // qcd_norm in W cards).  The info box auto-sizes to its line count.
    // Display labels avoid '_' (TLatex would render it as a subscript).
    const char *pars[3]  = {"dy_norm", "w_norm", "qcd_norm"};
    const char *plabs[3] = {"DY norm", "W norm", "QCD norm"};
    for (int ip = 0; ip < 3; ++ip) {
      RooRealVar *x = (RooRealVar *)fr->floatParsFinal().find(pars[ip]);
      if (x) box.push_back(Form("%s = %.3f #pm %.3f", plabs[ip], x->getVal(), x->getError()));
    }
  }

  PlotStyle ps;
  ps.drawOpt = "hist";
  ps.showStats = false;
  ps.logy = isW;            // log-y for the MET tails; linear for the Z peak
  ps.normBkgToData = false; // ABSOLUTE postfit yields -- never area-normalize
  ps.headerX = 0.56;        // channel header shifted left so it fits in-frame
  ps.boxTextSize = 0.028;   // smaller fit-result / info text
  ps.boxX1 = 0.56; ps.boxX2 = 0.93; // info box upper-right, contained in the frame
  ps.boxY1 = 0.56; ps.boxY2 = 0.76; // sits below the header
  ps.legX1 = 0.70; ps.legY1 = 0.15; // legend -> lower-right (away from the box)
  ps.legX2 = 0.93; ps.legY2 = 0.48;
  PlotTuner tuner = [&](TCanvas *c, TH1 *h) {
    (void)c; if (!h) return;
    if (ps.logy) h->SetMinimum(1.0);
    h->SetMaximum((ps.logy ? 10.25 : 1.4) * h->GetMaximum());
  };

  SaveNicePlot1D_WithBkg(hData, bkgs, names, outNoExt, xTitle, yTitle,
                         "", subTitle1, subTitle2, box, ps, tuner);

  ff->Close(); fin->Close();
  std::cout << "[postfit] " << outNoExt << ".png\n";
}
