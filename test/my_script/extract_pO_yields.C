// =============================================================================
// extract_pO_yields.C  --  turn the per-region FitDiagnostics outputs into
//   (a) a machine-readable CSV of fitted signal yields per (charge, rapidity bin)
//   (b) single-bin histograms named EXACTLY as analysis/charge_asym.C and
//       analysis/FBratio.C expect, so those macros consume FITTED yields with no
//       code change (they take the input file as their first argument):
//         lab bins -> h_mt_Wp_y{0..11}    / h_mt_Wm_y{0..11}
//         FB  bins -> h_mt_Wp_y{0..11}_FB / h_mt_Wm_y{0..11}_FB
//       Each is a 1-bin TH1D whose full-range integral == fitted signal yield and
//       whose Sumw2 error == the fit uncertainty (charge_asym/FBratio read the
//       error via TH1::IntegralAndError, so Sumw2 carries the right sigma).
//
// Fitted signal yield = r * (prefit signal integral), error = rErr * (same),
// because the POI 'r' scales the signal normalization linearly.  For the
// combined W+Z fit the W signal also carries the shared eff_lumi (r*eff_lumi).
//
// Usage (run under cmsenv, after the fits):
//   root -b -q 'extract_pO_yields.C("mu","<fitsDir>","<W.root>","<Z.root>","<outDir>")'
// =============================================================================
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include <fstream>
#include <iostream>

namespace {
struct FitVals { bool ok; double r,rE,qn,qnE,en,enE,el,elE; };

FitVals readFit(const TString &path) {
  FitVals v; v.ok=false; v.r=v.rE=v.qn=v.qnE=v.en=v.enE=v.el=v.elE=0;
  TFile *f = TFile::Open(path);
  if (!f || f->IsZombie()) { std::cerr << "[WARN] cannot open fit file: " << path << "\n"; return v; }
  RooFitResult *fr = (RooFitResult *)f->Get("fit_s");
  if (fr) {
    const RooArgList &ps = fr->floatParsFinal();
    auto g = [&](const char *n, double &val, double &err) {
      RooRealVar *x = (RooRealVar *)ps.find(n);
      if (x) { val = x->getVal(); err = x->getError(); }
    };
    v.ok = (ps.find("r") != 0);
    g("r", v.r, v.rE); g("qcd_norm", v.qn, v.qnE);
    g("ewk_norm", v.en, v.enE); g("eff_lumi", v.el, v.elE);
  } else std::cerr << "[WARN] no 'fit_s' in " << path << " (fit failed?)\n";
  f->Close(); delete f;
  return v;
}
} // namespace

void extract_pO_yields(const char *chan,        // "mu" or "ele" (label only)
                       const char *fitsDir,     // <workdir>/fits
                       const char *wInputFile,  // structured combine_input_W.root
                       const char *zInputFile,  // structured combine_input_Z.root
                       const char *outDir)      // <workdir>/summary
{
  gSystem->mkdir(outDir, kTRUE);

  TFile *wIn = TFile::Open(wInputFile, "READ");
  TFile *zIn = TFile::Open(zInputFile, "READ");
  auto sigPrefit = [&](TFile *fin, const TString &region) -> double {
    if (!fin || fin->IsZombie()) return -1;
    TH1 *h = (TH1 *)fin->Get(region + "/signal");
    return h ? h->Integral() : -1; // bins 1..N, matches Combine rate convention
  };

  TString yieldsRoot = TString::Format("%s/%s_fitted_yields.root", outDir, chan);
  TFile *fy = TFile::Open(yieldsRoot, "RECREATE");
  auto makeYieldHist = [&](const TString &name, double y, double e) {
    TH1D *h = new TH1D(name, name, 1, 0.0, 1.0);
    h->Sumw2();
    h->SetBinContent(1, y);
    h->SetBinError(1, e);
    h->SetDirectory(fy);
    h->Write(name, TObject::kOverwrite);
  };

  std::ofstream csv(TString::Format("%s/%s_W_yields.csv", outDir, chan).Data());
  csv << "region,charge,binning,ybin,r,rErr,signal_prefit,fitted_yield,fitted_yield_err,"
         "qcd_norm,qcd_normErr,ewk_norm,ewk_normErr\n";

  const char *charges[2]  = {"Wp", "Wm"};
  const char *binnings[2] = {"lab", "fb"};
  for (int ic = 0; ic < 2; ++ic)
    for (int ib = 0; ib < 2; ++ib)
      for (int iy = 0; iy < 12; ++iy) {
        TString R   = TString::Format("%s_%s_y%d", charges[ic], binnings[ib], iy);
        TString fit = TString::Format("%s/%s/fitDiagnostics_%s.root", fitsDir, R.Data(), R.Data());
        FitVals v   = readFit(fit);
        double Isig = sigPrefit(wIn, R);
        double y = (v.ok && Isig > 0) ? v.r * Isig : 0.0;
        double e = (v.ok && Isig > 0) ? v.rE * Isig : 0.0;
        csv << R << "," << charges[ic] << "," << binnings[ib] << "," << iy << ","
            << v.r << "," << v.rE << "," << Isig << "," << y << "," << e << ","
            << v.qn << "," << v.qnE << "," << v.en << "," << v.enE << "\n";
        TString hname = TString::Format("h_mt_%s_y%d%s", charges[ic], iy, (ib == 1 ? "_FB" : ""));
        makeYieldHist(hname, y, e);
      }
  csv.close();
  fy->Close(); delete fy;
  std::cout << "[extract] wrote " << yieldsRoot << "\n";
  std::cout << "[extract] wrote " << outDir << "/" << chan << "_W_yields.csv\n";

  // ---- inclusive + Z + combined summary -------------------------------------
  std::ofstream scsv(TString::Format("%s/%s_summary.csv", outDir, chan).Data());
  scsv << "fit,param,value,error,signal_prefit,fitted_yield,fitted_yield_err\n";
  auto dumpW = [&](const TString &R) {
    TString fit = TString::Format("%s/%s/fitDiagnostics_%s.root", fitsDir, R.Data(), R.Data());
    FitVals v = readFit(fit);
    double Isig = sigPrefit(wIn, R);
    double y = (v.ok && Isig > 0) ? v.r * Isig : 0, e = (v.ok && Isig > 0) ? v.rE * Isig : 0;
    scsv << R << ",r," << v.r << "," << v.rE << "," << Isig << "," << y << "," << e << "\n";
  };
  dumpW("Wp_incl"); dumpW("Wm_incl"); dumpW("W_incl");
  { // Z standalone (POI r scales Z signal)
    TString R = "Z_incl";
    FitVals v = readFit(TString::Format("%s/%s/fitDiagnostics_%s.root", fitsDir, R.Data(), R.Data()));
    double Isig = sigPrefit(zIn, R);
    double y = (v.ok && Isig > 0) ? v.r * Isig : 0, e = (v.ok && Isig > 0) ? v.rE * Isig : 0;
    scsv << "Z_incl,r," << v.r << "," << v.rE << "," << Isig << "," << y << "," << e << "\n";
  }
  { // combined W+Z: r (W strength) and the shared eff_lumi
    FitVals v = readFit(TString::Format("%s/WZ/fitDiagnostics_WZ.root", fitsDir));
    double Isig = sigPrefit(wIn, "W_incl");
    double y = (Isig > 0) ? v.r * v.el * Isig : 0; // W signal incl shared eff_lumi
    scsv << "WZ_combined,r," << v.r << "," << v.rE << "," << Isig << "," << y << ",0\n";
    scsv << "WZ_combined,eff_lumi," << v.el << "," << v.elE << ",,,\n";
  }
  scsv.close();
  std::cout << "[extract] wrote " << outDir << "/" << chan << "_summary.csv\n";

  if (wIn) { wIn->Close(); delete wIn; }
  if (zIn) { zIn->Close(); delete zIn; }
}
