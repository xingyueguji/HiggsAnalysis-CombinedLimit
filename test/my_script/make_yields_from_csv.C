// =============================================================================
// make_yields_from_csv.C -- rebuild <chan>_fitted_yields.root from the
// <chan>_W_yields.csv that extract_pO_yields.C already wrote. Use this when you
// have the CSV but the .root came out empty (no need to re-run the fits), or to
// regenerate the analysis-consumer histos locally.
//
// Produces the single-bin h_mt_W{p,m}_y{0..11}(_FB) histograms (full integral =
// fitted yield, Sumw2 error = fit uncertainty) that analysis/charge_asym.C and
// analysis/FBratio.C read directly.
//
// Usage:
//   root -b -q 'make_yields_from_csv.C("mu_W_yields.csv","mu_fitted_yields.root")'
// =============================================================================
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

void make_yields_from_csv(const char *csvPath, const char *outRoot) {
  std::ifstream in(csvPath);
  if (!in) { std::cerr << "[ERROR] cannot open " << csvPath << "\n"; return; }
  TFile *fy = TFile::Open(outRoot, "RECREATE");
  if (!fy || fy->IsZombie()) { std::cerr << "[ERROR] cannot create " << outRoot << "\n"; return; }

  std::string line;
  std::getline(in, line); // header: region,charge,binning,ybin,r,rErr,signal_prefit,fitted_yield,fitted_yield_err,...
  int n = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::vector<std::string> c;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) c.push_back(tok);
    if (c.size() < 9) continue;
    const std::string charge = c[1], binning = c[2];
    int iy = std::atoi(c[3].c_str());
    double y = std::atof(c[7].c_str());   // fitted_yield
    double e = std::atof(c[8].c_str());   // fitted_yield_err
    TString nm = TString::Format("h_mt_%s_y%d%s", charge.c_str(), iy, (binning == "fb" ? "_FB" : ""));
    fy->cd();
    TH1D *h = new TH1D(nm, nm, 1, 0.0, 1.0);
    h->Sumw2();
    h->SetBinContent(1, y);
    h->SetBinError(1, e);
    h->SetDirectory(fy);
    fy->WriteTObject(h, nm, "Overwrite");
    ++n;
  }
  fy->Close(); delete fy;
  std::cout << "[csv->yields] wrote " << n << " histos to " << outRoot << "\n";
}
