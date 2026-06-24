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
#include <algorithm>
#include <vector>
#include <string>

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
    TH1D *h = RemapToRef(ref, src, TString::Format("postfit_%s", procs[i].c_str()));
    if (h) { bkgs.push_back(h); names.push_back(labels[i]); }
  }
  TH1 *tot = (TH1 *)ff->Get(TString::Format("shapes_fit_s/%s/total", fitChannel));
  TH1D *hTot = RemapToRef(ref, tot, "postfit_total");

  std::vector<std::string> box = {
      Form("Data: %.0f", hData->Integral()),
      Form("Postfit total: %.0f", hTot ? hTot->Integral() : 0.0)};

  PlotStyle ps;
  ps.drawOpt = "hist";
  ps.showStats = false;
  ps.logy = isW;            // log-y for the MET tails; linear for the Z peak
  ps.boxY1 = 0.62; ps.boxY2 = 0.82;
  ps.normBkgToData = false; // ABSOLUTE postfit yields -- never area-normalize
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
