// =============================================================================
// extract_pO_simfit.C -- extraction for the GRAND SIMULTANEOUS FIT (simfit,
// 2026-08-04): ONE FitDiagnostics result per binning variant (lab / fb) holds
// ALL 25 POIs -- r_<C>_y<i> (24, mu/e shared) + the global DY scale r_Z --
// plus the QCD normalization parameters and (2026-08-17) the lumi nuisance.
//
// QCD modes (2026-08-17): the trailing kappa arguments select how the QCD
// normalization is read back, matching how make_pO_simfit_cards.sh built the
// cards (run_pO_fits.sh reads them from the qcd_lnn_kappas.txt sidecar):
//   kQcdMu/kQcdEle > 1 -> lnN mode: the fit has 4 shared nuisances
//     qcd_rate_{mu,ele}_{Wp,Wm} (theta, N(0,1)-constrained); the CSV qcd
//     columns then carry the MULTIPLIER kappa^theta (err = mult*ln(kappa)*dtheta),
//     identical semantics to the legacy columns (the factor on the qcd template).
//   kQcdMu/kQcdEle = 0 -> legacy free-rateParam mode: per-channel qcd_norm_* read.
//   qcdMode = "abcd" (2026-08-23, 8th argument; from the sidecar's qcdMode
//     line) -> the IN-FIT ABCD: the SR qcd is scaled by the FORMULA rateParam
//     (sB*sC/sD) x the residual kappa^theta. A formula rateParam is a
//     RooFormulaVar, ABSENT from floatParsFinal, so the multiplier is
//     evaluated here from its floating constituents qcd_s{B,C,D}_<F>_<C>
//     (+ theta), with the error propagated through the full correlation
//     matrix. The CSV qcd columns keep their meaning: the factor on the qcd
//     template written in the input file (which in abcd mode is qcd_abcd,
//     total B0*C40/D0). kQcdMu/kQcdEle then hold the REDUCED kappas.
//     An empty/absent qcdMode falls back to the kappa-based inference
//     (legacy sidecars): kQcd > 1 -> lnN, = 0 -> free.
//   kLumi > 1 -> the global 'lumi' lnN exists; its multiplier (+err) is appended
//     as two extra CSV columns (lumi,lumiErr; 1,0 when off) -- postfit_incl.C
//     uses it to keep the "prefit x fitted scale" reconstruction exact.
//   The trailing qcd_model CSV column (18th, 2026-08-23) records which of the
//     three modes produced the row, so readers (postfit_incl.C) know whether
//     the multiplier applies to the `qcd` or the `qcd_abcd` template.
//   lheSysts (2026-09-07, 9th argument; the sidecar's lheSysts line, comma
//     list, ""/"none" = none): the LHE shape nuisances (nPDF, qcdScale,
//     alphaS) -- their pulls AND post-fit constraints go to comb_summary.csv
//     (<name>_theta rows: value = pull, error = constraint; error < 1 means the
//     data constrained that shape) and they join the Asimov closure (theta = 0).
// STAT COMPONENT -- the 19th CSV column rErr_stat, the comb_summary.csv rows
//   simfit_<B>_stat and the matrices h_cov_yield[_FB]_stat +
//   h_cov_poi[_FB]_stat. rErr / h_cov_yield / h_cov_poi stay the TOTAL
//   (profiled) errors, so downstream syst = sqrt(total^2 - stat^2).
//   SOURCE, user decision 2026-09-15b: the --statonly COMPANION FIT
//   (fitDiagnostics_simfit_<B>_statonly.root -- every constrained nuisance
//   frozen at its POST-FIT value, the Combine breakdown recipe run
//   numerically), now produced by default. When it is absent the extraction
//   falls back to the CONDITIONED covariance of fit_s (ComputeStatCov below:
//   the Schur complement = the Gaussian-exact equivalent of that refit,
//   without the refit), and the conditioning is always computed anyway as the
//   CROSS-CHECK -- the two agreeing says the likelihood is parabolic in the
//   POI directions. One line reports which source was used; everything reads
//   through the statErr/statCov lambdas, so the choice is made in one place.
//   NB the conditioning guarantees stat <= total (a PSD term is subtracted); a
//   separate refit does NOT, so POIs with stat > total are counted and warned
//   about rather than silently floored by a downstream max(0, .).
// PREFIT-S CHECK (2026-09-15): the prefit signal integrals S are read from the
//   two W input files and, when the nominal fit carries shapes_prefit
//   (--saveShapes), compared per region with the sum of the fitted channels'
//   prefit signal shapes -- a WARN means the input files are NOT the ones that
//   were fitted (e.g. --extract-only on regenerated inputs).
// Nuisance thetas (pulls) are also dumped to comb_summary.csv (<name>_theta
// rows) -- THE check that the ABCD prediction and its assigned uncertainty are
// consistent with the data (|pull| ~> 1 means kappa too small or template biased).
// A final sweep prints any floating parameter of fit_s not reported above, so
// a new nuisance can never be silently invisible.
// PER-FLAVOUR FITS (2026-09-22, run_pO_fits.sh mode flavfit): the trailing
//   `flavours` argument ("mu,ele" = the grand fit, default; "mu" or "ele")
//   names the lepton flavours the fit contains, and `outTag` the prefix of the
//   three output files ("comb" for the grand fit, as always; "simfit_mu" /
//   "simfit_ele" for the per-flavour ones). Only the W input files of those
//   flavours are opened, S sums over them alone (so the yields are that
//   flavour's), the prefit-S check compares with that many fitted channels,
//   and the per-flavour QCD reporting (pulls, CR scales, the CSV qcd columns)
//   covers them alone -- the CSV keeps its 19-column layout, with the absent
//   flavour's two qcd columns written as 0,0 ("not in this fit").
//
// It writes into <outDir> (both variants into the same files; the file names
// below are for outTag = "comb" -- per-flavour fits use their own prefix):
//   (a) comb_W_yields.csv -- one row per (charge, binning, y bin): r, rErr,
//       the mu+e SUMMED prefit signal integral S, and the COMBINED fitted
//       yield Y = r * S.  With the mu/e-shared r the per-flavour yields are
//       100% correlated, so the combined yield is the primary result.  The
//       first 9 columns keep the legacy <chan>_W_yields.csv layout
//       (c[7]/c[8] = yield/err are load-bearing for downstream CSV readers);
//       the per-flavour qcd norms and r_Z follow.
//   (b) comb_fitted_yields.root -- single-bin h_yield_W{p,m}_y{0..11}[_FB]
//       (+ deprecated h_mt_* aliases): EXACTLY the names analysis/charge_asym.C
//       and analysis/FBratio.C read (pass this file as their first argument).
//       PLUS the 24x24 fitted-yield covariance matrices h_cov_yield (lab) /
//       h_cov_yield_FB (fb) built from the fit correlation matrix, fixed order
//       [Wp_y0..y11, Wm_y0..y11] (axis bin labels set accordingly), so
//       downstream error propagation can include the r-correlation cross terms
//       (bins are correlated through the shared r_Z and, weakly, the QCD).
//   (c) comb_summary.csv -- every POI (value, error, prefit S, yield),
//       fit status/covQual, covariance-propagated inclusive sums (Wp/Wm/W),
//       and -- when the --asimov closure fit exists -- the Asimov-fitted POIs
//       (closure: every POI must come back at 1).
//
// Usage (run under cmsenv, after run_pO_fits.sh simfit):
//   root -b -q 'extract_pO_simfit.C("<fitsDir>","<muW.root>","<eleW.root>","<outDir>")'
// (run_pO_fits.sh passes every argument, incl. flavours/outTag for flavfit;
//  the W file of a flavour not in the fit may be "none")
// =============================================================================
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

namespace {

struct Par { bool ok; double v, e; };

Par GetPar(const RooFitResult *fr, const TString &name) {
  Par p; p.ok = false; p.v = p.e = 0.0;
  if (!fr) return p;
  RooRealVar *x = (RooRealVar *)fr->floatParsFinal().find(name);
  if (x) { p.ok = true; p.v = x->getVal(); p.e = x->getError(); }
  return p;
}

// lnN nuisance theta -> multiplicative scale kappa^theta with propagated error.
// !ok when the nuisance is absent from fit_s (or kappa is not a valid lnN).
Par LnNScale(const RooFitResult *fr, double kappa, const TString &name) {
  Par p; p.ok = false; p.v = 1.0; p.e = 0.0;
  if (kappa <= 1.0) return p;
  Par th = GetPar(fr, name);
  if (!th.ok) return p;
  p.ok = true;
  p.v = std::pow(kappa, th.v);
  p.e = p.v * std::log(kappa) * th.e;
  return p;
}

// abcd mode: the SR qcd multiplier M = kappa^theta * sB*sC/sD. The formula
// rateParam itself is a RooFormulaVar (not in floatParsFinal), so M is
// evaluated from its floating constituents with the error propagated through
// the full correlation matrix:  dM/dp_i / M = (ln kappa, 1/sB, 1/sC, -1/sD).
// Parameters with zero error (or absent) drop out of the propagation.
Par AbcdScale(const RooFitResult *fr, double kappa, const char *flav, const char *chg) {
  Par p; p.ok = false; p.v = 1.0; p.e = 0.0;
  if (!fr) return p;
  const TString nB = TString::Format("qcd_sB_%s_%s", flav, chg);
  const TString nC = TString::Format("qcd_sC_%s_%s", flav, chg);
  const TString nD = TString::Format("qcd_sD_%s_%s", flav, chg);
  const Par sB = GetPar(fr, nB), sC = GetPar(fr, nC), sD = GetPar(fr, nD);
  if (!sB.ok || !sC.ok || !sD.ok || sB.v <= 0.0 || sC.v <= 0.0 || sD.v <= 0.0)
    return p;
  const bool hasK = (kappa > 1.0);
  const TString nT = TString::Format("qcd_rate_%s_%s", flav, chg);
  Par th; th.ok = false; th.v = 0.0; th.e = 0.0;
  if (hasK) th = GetPar(fr, nT); // absent -> theta 0 with no error contribution
  p.ok = true;
  p.v = (hasK && th.ok ? std::pow(kappa, th.v) : 1.0) * sB.v * sC.v / sD.v;
  const TString nm[4] = {nT, nB, nC, nD};
  const double  gr[4] = {(hasK && th.ok) ? std::log(kappa) : 0.0,
                         1.0 / sB.v, 1.0 / sC.v, -1.0 / sD.v};
  const double  er[4] = {th.e, sB.e, sC.e, sD.e};
  double var = 0.0;
  for (int a = 0; a < 4; ++a) {
    if (er[a] <= 0.0) continue;
    for (int b = 0; b < 4; ++b) {
      if (er[b] <= 0.0) continue;
      const double rho = (a == b) ? 1.0 : fr->correlation(nm[a].Data(), nm[b].Data());
      var += gr[a] * gr[b] * rho * er[a] * er[b];
    }
  }
  p.e = p.v * std::sqrt(var > 0.0 ? var : 0.0);
  return p;
}

// nullptr when the file does not exist / has no fit_s; caller owns *f.
RooFitResult *OpenFitS(const TString &path, TFile *&f) {
  f = nullptr;
  if (gSystem->AccessPathName(path)) return nullptr; // kTRUE = NOT accessible
  f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) { if (f) { f->Close(); delete f; f = nullptr; } return nullptr; }
  return (RooFitResult *)f->Get("fit_s");
}

// ---- STATISTICAL component of the post-fit covariance, from fit_s ALONE ------
// (2026-09-15.) The error Minuit reports for a POI is the TOTAL (profiled) one:
// every nuisance is free to move while the POI is displaced. Its statistical
// part is the error the same POI would have with every CONSTRAINED nuisance
// held fixed at its post-fit value -- and that number is already contained in
// the nominal fit's covariance matrix. With V = the full Hesse covariance over
// all floating parameters, split into the STATISTICAL block k (the POIs and the
// unconstrained rateParams: r_*, qcd_s*, qcd_norm* -- data-driven quantities,
// e.g. the in-fit ABCD scales) and the nuisance block n (everything else), the
// covariance of k with n FIXED is the Schur complement
//     V_stat = V_kk - V_kn V_nn^{-1} V_nk      ( = (H_kk)^{-1},  H = V^{-1} ),
// i.e. exactly what a refit with `--freezeParameters allConstrainedNuisances`
// at the post-fit values returns in the Gaussian (Hesse) approximation -- the
// Combine "breakdown" recipe -- without the refit. Checked on the 2026-09-14
// fit: the two formulas agree to 1e-16 and sqrt(V_ii) reproduces every quoted
// error to 2e-4 (the errors ARE the Hesse diagonal). The subtracted term is
// positive semi-definite, so stat <= total per parameter by construction and
// syst = sqrt(total^2 - stat^2) is always real. The classification is BY NAME
// (the same rule as run_pO_fits.sh's frozen-refit one-liner) and is printed, so
// a new unconstrained parameter under another name is visible as
// "conditioned out" instead of silently shrinking the stat errors.
bool IsStatParam(const TString &n) {
  return n.BeginsWith("r_") || n.BeginsWith("qcd_s") || n.BeginsWith("qcd_norm");
}

struct StatCov {
  bool ok = false;
  std::vector<TString> kept, nuis; // statistical parameters / conditioned-out nuisances
  std::map<TString, int> idx;      // kept name -> row of V
  TMatrixD V;                      // kept x kept statistical covariance
  double err(const TString &n) const {
    std::map<TString, int>::const_iterator it = idx.find(n);
    return (it == idx.end()) ? -1.0 : std::sqrt(std::max(0.0, V(it->second, it->second)));
  }
  double cov(const TString &a, const TString &b) const {
    std::map<TString, int>::const_iterator ia = idx.find(a), ib = idx.find(b);
    return (ia == idx.end() || ib == idx.end()) ? 0.0 : V(ia->second, ib->second);
  }
};

StatCov ComputeStatCov(const RooFitResult *fr) {
  StatCov sc;
  if (!fr) return sc;
  const RooArgList &pars = fr->floatParsFinal();
  const int n = pars.getSize();
  const TMatrixDSym &Vall = fr->covarianceMatrix(); // ordered like floatParsFinal
  if (n == 0 || Vall.GetNrows() != n) return sc;
  std::vector<int> ik, in;
  for (int i = 0; i < n; ++i) {
    const TString nm = pars.at(i)->GetName();
    if (IsStatParam(nm)) { ik.push_back(i); sc.kept.push_back(nm); }
    else                 { in.push_back(i); sc.nuis.push_back(nm); }
  }
  const int nk = (int)ik.size(), nn = (int)in.size();
  if (nk == 0) return sc;
  TMatrixD Vkk(nk, nk);
  for (int a = 0; a < nk; ++a)
    for (int b = 0; b < nk; ++b) Vkk(a, b) = Vall(ik[a], ik[b]);
  if (nn > 0) {
    TMatrixD Vkn(nk, nn), Vnn(nn, nn);
    for (int a = 0; a < nk; ++a)
      for (int c = 0; c < nn; ++c) Vkn(a, c) = Vall(ik[a], in[c]);
    for (int c = 0; c < nn; ++c)
      for (int d = 0; d < nn; ++d) Vnn(c, d) = Vall(in[c], in[d]);
    double det = 0.0;
    Vnn.Invert(&det);
    if (det == 0.0 || !std::isfinite(det)) return sc; // singular nuisance block: give up
    const TMatrixD Vnk(TMatrixD::kTransposed, Vkn);
    Vkk -= Vkn * Vnn * Vnk;
  }
  sc.V.ResizeTo(nk, nk);
  sc.V = Vkk;
  for (int a = 0; a < nk; ++a) sc.idx[sc.kept[a]] = a;
  sc.ok = true;
  return sc;
}

} // namespace

void extract_pO_simfit(const char *fitsDir,  // <workdir>/fits (simfit_lab/, simfit_fb/)
                       const char *muWFile,  // muon structured combine_input_W.root
                       const char *eleWFile, // electron structured combine_input_W.root
                       const char *outDir,   // <workdir>/summary
                       double kQcdMu  = 0.0, // lnN kappas the cards were built with
                       double kQcdEle = 0.0, // (0 = free-rateParam legacy mode;
                       double kLumi   = 0.0, //  see qcd_lnn_kappas.txt sidecar)
                       const char *qcdMode = "", // "abcd"|"lnN"|"free"|"" (sidecar
                                                 //  qcdMode line; "" = infer from kappas)
                       const char *lheSysts = "", // comma list of the LHE shape
                                                  //  nuisances in the cards ("" = none)
                       const char *flavours = "mu,ele", // lepton flavours IN the fit
                                                        //  (2026-09-22: "mu" | "ele" = flavfit)
                       const char *outTag = "comb")     // output-file prefix: comb | simfit_mu | simfit_ele
{
  gSystem->mkdir(outDir, kTRUE);

  // ---- which lepton flavours this fit contains ([0] = mu, [1] = ele) --------
  bool useFl[2] = {false, false};
  {
    TString s(flavours);
    TObjArray *toks = s.Tokenize(",");
    for (int i = 0; i < toks->GetEntries(); ++i) {
      TString t = ((TObjString *)toks->At(i))->GetString();
      t.ReplaceAll(" ", "");
      if (t == "mu") useFl[0] = true;
      else if (t == "ele") useFl[1] = true;
      else if (t != "") std::cerr << "[extract-simfit] WARN unknown flavour '" << t << "' ignored\n";
    }
    delete toks;
  }
  const int nFl = (useFl[0] ? 1 : 0) + (useFl[1] ? 1 : 0);
  if (nFl == 0) {
    std::cerr << "[extract-simfit] WARN no valid flavour in '" << flavours << "' -- aborting extraction\n";
    return;
  }
  std::cout << "[extract-simfit] fit '" << outTag << "': flavours "
            << (useFl[0] ? "mu " : "") << (useFl[1] ? "ele" : "") << "\n";

  const bool isAbcd = (TString(qcdMode) == "abcd");
  std::vector<TString> lheNames; // the LHE shape nuisances (sidecar lheSysts line)
  {
    TString s(lheSysts);
    if (s != "" && s != "none") {
      TObjArray *toks = s.Tokenize(",");
      for (int i = 0; i < toks->GetEntries(); ++i) {
        TString t = ((TObjString *)toks->At(i))->GetString();
        t.ReplaceAll(" ", "");
        if (t != "") lheNames.push_back(t);
      }
      delete toks;
    }
  }
  // the qcd_model tag stamped on every CSV row (readers pick qcd vs qcd_abcd)
  const char *qModel = isAbcd ? "abcd" : ((kQcdMu > 1.0 || kQcdEle > 1.0) ? "lnN" : "free");

  // the W input files of the flavours IN the fit (a flavour not in it is never
  // opened -- run_pO_fits.sh passes "none" for it)
  TFile *wmu = useFl[0] ? TFile::Open(muWFile, "READ") : nullptr;
  TFile *wel = useFl[1] ? TFile::Open(eleWFile, "READ") : nullptr;
  if ((useFl[0] && (!wmu || wmu->IsZombie())) || (useFl[1] && (!wel || wel->IsZombie()))) {
    std::cerr << "[extract-simfit] WARN cannot open W input file(s) -- aborting extraction\n";
    return;
  }
  TFile *wIn[2] = {wmu, wel};
  auto sigPrefit = [&](TFile *fin, const TString &region) -> double {
    TH1 *h = (TH1 *)fin->Get(region + "/signal");
    return h ? h->Integral() : -1.0; // bins 1..N, matches Combine rate convention
  };

  TString yieldsRoot = TString::Format("%s/%s_fitted_yields.root", outDir, outTag);
  TFile *fy = TFile::Open(yieldsRoot, "RECREATE");
  auto makeYieldHist = [&](const TString &name, double y, double e) {
    fy->cd(); // OpenFitS opens/closes files in the loop: make fy current first
    TH1D *h = new TH1D(name, name, 1, 0.0, 1.0);
    h->Sumw2();
    h->SetBinContent(1, y);
    h->SetBinError(1, e);
    h->SetDirectory(fy);
    fy->WriteTObject(h, name, "Overwrite");
  };

  // 19th column rErr_stat (2026-09-14/15): the STATISTICAL component of the POI
  // error, from the nominal fit's covariance conditioned on the constrained
  // nuisances (ComputeStatCov); -1 only when fit_s carries no usable covariance.
  // rErr stays the TOTAL (profiled) error; syst = sqrt(rErr^2 - rErr_stat^2).
  std::ofstream csv(TString::Format("%s/%s_W_yields.csv", outDir, outTag).Data());
  csv << "region,charge,binning,ybin,r,rErr,signal_prefit,fitted_yield,fitted_yield_err,"
         "qcd_norm_mu,qcd_norm_muErr,qcd_norm_ele,qcd_norm_eleErr,r_Z,r_ZErr,lumi,lumiErr,"
         "qcd_model,rErr_stat\n";

  std::ofstream scsv(TString::Format("%s/%s_summary.csv", outDir, outTag).Data());
  scsv << "fit,param,value,error,signal_prefit,fitted_yield,fitted_yield_err\n";

  const char *charges[2]  = {"Wp", "Wm"};
  const char *binnings[2] = {"lab", "fb"};
  const int NB = 12, NPOI = 24;

  for (int ib = 0; ib < 2; ++ib) {
    const char *B = binnings[ib];
    TString fitName = TString::Format("simfit_%s", B);
    TFile *fd = nullptr;
    RooFitResult *fr = OpenFitS(TString::Format("%s/%s/fitDiagnostics_%s.root",
                                                fitsDir, fitName.Data(), fitName.Data()), fd);
    if (!fr) {
      std::cerr << "[extract-simfit] WARN no fit_s for " << fitName << " -- variant skipped\n";
      if (fd) { fd->Close(); delete fd; }
      continue;
    }
    Par rz = GetPar(fr, "r_Z");
    std::cout << "[extract-simfit] " << fitName << ": status " << fr->status()
              << ", covQual " << fr->covQual()
              << ", r_Z = " << Form("%.4f +/- %.4f", rz.v, rz.e) << "\n";
    // the STATISTICAL component of every POI error, from fit_s alone (2026-09-15):
    // the full post-fit covariance conditioned on the constrained nuisances
    const StatCov sc = ComputeStatCov(fr);
    if (sc.ok) {
      std::cout << "[extract-simfit] " << fitName << ": stat component from the fit covariance: "
                << sc.kept.size() << " statistical parameters kept (POIs + unconstrained rateParams), "
                << sc.nuis.size() << " constrained nuisances conditioned out:";
      for (size_t i = 0; i < sc.nuis.size(); ++i) std::cout << " " << sc.nuis[i];
      std::cout << "\n";
    } else
      std::cout << "[extract-simfit] WARN " << fitName
                << ": no usable covariance matrix in fit_s -> conditioned stat component unavailable\n";
    // the frozen-nuisance companion fit (run_pO_fits.sh --statonly, ON by
    // default since 2026-09-15b): every CONSTRAINED nuisance frozen at its
    // POST-FIT value, so the POI errors that come back ARE the statistical
    // component -- the Combine breakdown recipe, run numerically.
    TFile *fds = nullptr;
    RooFitResult *frs = OpenFitS(TString::Format("%s/%s/fitDiagnostics_%s_statonly.root",
                                                 fitsDir, fitName.Data(), fitName.Data()), fds);
    if (frs) {
      Par rzs = GetPar(frs, "r_Z");
      std::cout << "[extract-simfit] " << fitName << "_statonlyfit: status " << frs->status()
                << ", covQual " << frs->covQual()
                << ", r_Z = " << Form("%.4f +/- %.4f vs %.4f conditioned", rzs.v, rzs.e, sc.err("r_Z")) << "\n";
      if (frs->status() != 0 || frs->covQual() < 3)
        std::cout << "[extract-simfit] WARN " << fitName << "_statonlyfit quality flags -- inspect fit_statonly.log\n";
      scsv << fitName << "_statonlyfit,r_Z," << rzs.v << "," << rzs.e << ",,,\n";
      scsv << fitName << "_statonlyfit,fit_status," << frs->status() << "," << frs->covQual() << ",,,\n";
    }

    // ---- WHICH source the STATISTICAL component comes from -------------------
    // USER DECISION 2026-09-15b: the --statonly companion FIT is the primary
    // source; syst = sqrt(total^2 - stat^2) downstream. The conditioned
    // covariance (ComputeStatCov, the Schur complement of fit_s) is the exact
    // Gaussian answer and needs no refit, so it stays as the FALLBACK when the
    // companion is absent AND as the cross-check printed below -- the two
    // agreeing is the statement that the likelihood is parabolic in the POI
    // directions. Everything downstream (rErr_stat, h_cov_yield[_FB]_stat,
    // h_cov_poi[_FB]_stat, the simfit_<B>_stat rows and the inclusive sums)
    // reads ONLY through statErr/statCov, so the choice is made in one place.
    //
    // NB unlike the conditioning, a separate refit does NOT guarantee
    // stat <= total (different minimizer path, non-parabolic directions), so
    // every consumer of syst = sqrt(total^2 - stat^2) must clamp -- and the
    // violations are counted and reported here rather than silently floored.
    const bool statFromFit = (frs != nullptr);
    auto statAvail = [&]() { return statFromFit || sc.ok; };
    auto statErr = [&](const TString &n) -> double {
      if (statFromFit) { Par p = GetPar(frs, n); return p.ok ? p.e : -1.0; }
      return sc.ok ? sc.err(n) : -1.0;
    };
    auto statCov = [&](const TString &a, const TString &b) -> double {
      if (statFromFit) {
        Par pa = GetPar(frs, a), pb = GetPar(frs, b);
        if (!pa.ok || !pb.ok) return 0.0;
        return ((a == b) ? 1.0 : frs->correlation(a.Data(), b.Data())) * pa.e * pb.e;
      }
      return sc.ok ? sc.cov(a, b) : 0.0;
    };
    std::cout << "[extract-simfit] " << fitName << ": STAT component from "
              << (statFromFit ? "the --statonly companion FIT (frozen nuisances)"
                              : (sc.ok ? "the CONDITIONED covariance (no companion fit found)"
                                       : "NOTHING -- no companion fit and no usable covariance"))
              << "\n";
    if (statAvail()) {
      const double ezs = statErr("r_Z");
      std::cout << "[extract-simfit] " << fitName << ": r_Z = "
                << Form("%.4f +/- %.4f total, +/- %.4f stat -> syst %.4f", rz.v, rz.e, ezs,
                        std::sqrt(std::max(0.0, rz.e * rz.e - ezs * ezs)))
                << "\n";
      scsv << fitName << "_stat,r_Z," << rz.v << "," << ezs << ",,,\n";
    }
    if (fr->status() != 0 || fr->covQual() < 3)
      std::cout << "[extract-simfit] WARN " << fitName
                << " quality flags (status!=0 or covQual<3) -- inspect fit.log\n";

    // ---- constrained nuisances (2026-08-17 lnN model): pulls + multipliers --
    // The pull table is the ABCD-consistency check: |theta| ~> 1 means the
    // assigned kappa is too small or the template normalization is biased.
    Par lum; lum.ok = false; lum.v = 1.0; lum.e = 0.0;
    if (kLumi > 1.0) {
      lum = LnNScale(fr, kLumi, "lumi");
      Par lth = GetPar(fr, "lumi");
      if (lum.ok) {
        std::cout << "[extract-simfit] " << fitName << ": lumi pull = "
                  << Form("%.3f +/- %.3f (scale %.4f +/- %.4f)", lth.v, lth.e, lum.v, lum.e) << "\n";
        scsv << fitName << ",lumi_theta," << lth.v << "," << lth.e << ",,,\n";
      } else {
        std::cerr << "[extract-simfit] WARN 'lumi' nuisance not in fit_s (kLumi="
                  << kLumi << " given)\n";
      }
    }
    if (kQcdMu > 1.0 || kQcdEle > 1.0) {
      for (int ifl = 0; ifl < 2; ++ifl) {
        if (!useFl[ifl]) continue; // flavour not in this fit
        const char *flav = (ifl == 0) ? "mu" : "ele";
        const double kap = (ifl == 0) ? kQcdMu : kQcdEle;
        if (kap <= 1.0) continue;
        for (int ic2 = 0; ic2 < 2; ++ic2) {
          const TString nn = TString::Format("qcd_rate_%s_%s", flav, charges[ic2]);
          Par th = GetPar(fr, nn);
          if (!th.ok) { std::cerr << "[extract-simfit] WARN nuisance " << nn << " not in fit_s\n"; continue; }
          std::cout << "[extract-simfit] " << fitName << ": " << nn << " pull = "
                    << Form("%.3f +/- %.3f (scale %.4f)", th.v, th.e, std::pow(kap, th.v)) << "\n";
          scsv << fitName << "," << nn << "_theta," << th.v << "," << th.e << ",,,\n";
        }
      }
    }

    // ---- abcd mode: the CR scales + the assembled SR multiplier -------------
    // sB/sC/sD are the fitted QCD yields of the three CRs relative to the
    // prefit EWK-subtracted counts; the SR multiplier sB*sC/sD (x kappa^theta)
    // is what the CSV qcd columns carry. Postfit QCD yield = multiplier x
    // Sum_y Int(qcd_abcd) = multiplier x B0*C40/D0.
    if (isAbcd) {
      for (int ifl = 0; ifl < 2; ++ifl) {
        if (!useFl[ifl]) continue; // flavour not in this fit
        const char *flav = (ifl == 0) ? "mu" : "ele";
        const double kap = (ifl == 0) ? kQcdMu : kQcdEle;
        for (int ic2 = 0; ic2 < 2; ++ic2) {
          const char *scl[3] = {"B", "C", "D"};
          double sv[3] = {0.0, 0.0, 0.0};
          bool allok = true;
          for (int is2 = 0; is2 < 3; ++is2) {
            const TString nn = TString::Format("qcd_s%s_%s_%s", scl[is2], flav, charges[ic2]);
            Par s = GetPar(fr, nn);
            if (!s.ok) {
              std::cerr << "[extract-simfit] WARN CR scale " << nn << " not in fit_s\n";
              allok = false;
              continue;
            }
            sv[is2] = s.v;
            scsv << fitName << "," << nn << "," << s.v << "," << s.e << ",,,\n";
          }
          Par m = AbcdScale(fr, kap, flav, charges[ic2]);
          if (allok && m.ok) {
            std::cout << "[extract-simfit] " << fitName << ": qcd ABCD " << flav << "_"
                      << charges[ic2] << ": sB=" << Form("%.3f", sv[0])
                      << " sC=" << Form("%.3f", sv[1]) << " sD=" << Form("%.3f", sv[2])
                      << " -> SR multiplier = " << Form("%.4f +/- %.4f", m.v, m.e) << "\n";
            scsv << fitName << ",qcd_abcd_mult_" << flav << "_" << charges[ic2] << ","
                 << m.v << "," << m.e << ",,,\n";
          }
        }
      }
    }

    // ---- LHE shape nuisances (2026-09-07): pull + constraint ---------------
    // theta ~ N(0,1) a priori: value = pull (data preferred a shifted template),
    // error = post-fit constraint (< 1 only if the data measure that shape --
    // for a theory nuisance that is worth knowing, it means the data are
    // tuning e.g. muR/muF).
    for (size_t il = 0; il < lheNames.size(); ++il) {
      Par th = GetPar(fr, lheNames[il]);
      if (!th.ok) {
        std::cerr << "[extract-simfit] WARN LHE nuisance " << lheNames[il] << " not in fit_s\n";
        continue;
      }
      std::cout << "[extract-simfit] " << fitName << ": " << lheNames[il] << " pull = "
                << Form("%.3f +/- %.3f", th.v, th.e)
                << (th.e < 0.9 ? "  (constrained by the data)" : "") << "\n";
      scsv << fitName << "," << lheNames[il] << "_theta," << th.v << "," << th.e << ",,,\n";
    }

    // ---- safety net: every floating parameter not reported above -----------
    {
      const RooArgList &fp = fr->floatParsFinal();
      for (int i = 0; i < fp.getSize(); ++i) {
        const RooRealVar *x = (const RooRealVar *)fp.at(i);
        const TString n = x->GetName();
        bool known = n.BeginsWith("r_") || n == "lumi" || n.BeginsWith("qcd_rate_") ||
                     n.BeginsWith("qcd_s") || n.BeginsWith("qcd_norm_");
        for (size_t il = 0; il < lheNames.size() && !known; ++il) known = (n == lheNames[il]);
        if (!known)
          std::cout << "[extract-simfit] WARN " << fitName << ": floating parameter " << n << " = "
                    << Form("%.4f +/- %.4f", x->getVal(), x->getError())
                    << " is not reported in the CSVs (add it to the sidecar / this extractor)\n";
      }
    }

    // ---- POIs + prefit integrals, fixed order [Wp_y0..11, Wm_y0..11] --------
    std::vector<TString> pois, regs;
    std::vector<double> rv(NPOI, 0), re(NPOI, 0), S(NPOI, 0), res(NPOI, -1.0);
    std::vector<double> rvf(NPOI, 0), ref(NPOI, -1.0); // the companion fit's r / rErr
    std::vector<double> rec(NPOI, -1.0);               // the CONDITIONED stat error (cross-check)
    std::vector<bool> ok(NPOI, false), okf(NPOI, false);
    double rmin = 1e30, rmax = -1e30, maxDr = 0.0, maxDe = 0.0, maxDS = 0.0;
    int nStatGtTot = 0; double worstStatGtTot = 0.0; TString worstStatPoi;
    const char *flavs[2] = {"mu", "ele"};
    for (int ic = 0; ic < 2; ++ic)
      for (int iy = 0; iy < NB; ++iy) {
        const int k = ic * NB + iy;
        pois.push_back(TString::Format("r_%s_y%d", charges[ic], iy));
        regs.push_back(TString::Format("%s_%s_y%d", charges[ic], B, iy));
        Par p = GetPar(fr, pois[k]);
        if (p.ok && statAvail()) res[k] = statErr(pois[k]); // THE stat component (source chosen above)
        if (sc.ok && p.ok) rec[k] = sc.err(pois[k]);        // conditioned value, for the cross-check
        if (frs) { // the companion fit's own r (same minimum expected)
          Par ps = GetPar(frs, pois[k]);
          if (ps.ok) {
            okf[k] = true; ref[k] = ps.e; rvf[k] = ps.v;
            if (p.ok) maxDr = std::max(maxDr, std::fabs(ps.v - p.v));
            if (rec[k] > 0) maxDe = std::max(maxDe, std::fabs(ps.e / rec[k] - 1.0));
          }
        }
        // syst = sqrt(total^2 - stat^2) needs stat <= total; the conditioning
        // guarantees it, a separate refit does not -- count and report instead
        // of letting a downstream max(0, .) hide it
        if (p.ok && res[k] > 0 && res[k] > p.e) {
          ++nStatGtTot;
          const double rel = res[k] / p.e - 1.0;
          if (rel > worstStatGtTot) { worstStatGtTot = rel; worstStatPoi = pois[k]; }
        }
        // S = the prefit signal summed over the flavours IN the fit (both for the
        // grand fit, that flavour's alone for a per-flavour one)
        S[k] = 0;
        for (int jf = 0; jf < 2; ++jf) {
          if (!useFl[jf]) continue;
          const double s = sigPrefit(wIn[jf], regs[k]);
          if (s < 0) std::cerr << "[extract-simfit] WARN missing prefit " << flavs[jf] << " signal for " << regs[k] << "\n";
          else S[k] += s;
        }
        ok[k] = p.ok && S[k] > 0;
        // prefit-S check against the fitted channels' saved prefit signal shapes
        // (shapes_prefit/<F>_<C>_<B>_y<i>/signal, present with --saveShapes): the
        // input files must be the ones that were fitted
        {
          double sfit = 0; int nf = 0;
          for (int jf = 0; jf < 2; ++jf) {
            if (!useFl[jf]) continue;
            TH1 *hp = (TH1 *)fd->Get(TString::Format("shapes_prefit/%s_%s_%s_y%d/signal", flavs[jf], charges[ic], B, iy));
            if (hp) { sfit += hp->Integral(); ++nf; }
          }
          if (nf == nFl && S[k] > 0) maxDS = std::max(maxDS, std::fabs(sfit / S[k] - 1.0));
        }
        rv[k] = p.v; re[k] = p.e;
        if (!p.ok) {
          std::cerr << "[extract-simfit] WARN POI " << pois[k] << " not in fit_s\n";
        } else {
          if (p.v < rmin) rmin = p.v;
          if (p.v > rmax) rmax = p.v;
          if (p.v <= 1e-3 || p.v >= 9.99)
            std::cerr << "[extract-simfit] WARN " << pois[k] << " = " << p.v
                      << " at its range boundary -- error unreliable\n";
        }

        const double y = ok[k] ? rv[k] * S[k] : 0.0;
        const double e = ok[k] ? re[k] * S[k] : 0.0;
        const TString suffix = TString::Format("%s_y%d%s", charges[ic], iy, (ib == 1 ? "_FB" : ""));
        makeYieldHist("h_yield_" + suffix, y, e);
        makeYieldHist("h_mt_" + suffix, y, e); // deprecated alias

        // qcd factor on this channel's template: abcd mode -> the evaluated
        // formula (sB*sC/sD) x kappa^theta (on the qcd_abcd template); lnN mode
        // -> kappa^theta of the shared (flavour, charge) nuisance; legacy free
        // mode -> the per-channel rateParam
        Par qmu = isAbcd ? AbcdScale(fr, kQcdMu, "mu", charges[ic])
                  : (kQcdMu > 1.0)
                      ? LnNScale(fr, kQcdMu, TString::Format("qcd_rate_mu_%s", charges[ic]))
                      : GetPar(fr, TString::Format("qcd_norm_mu_%s", regs[k].Data()));
        Par qel = isAbcd ? AbcdScale(fr, kQcdEle, "ele", charges[ic])
                  : (kQcdEle > 1.0)
                      ? LnNScale(fr, kQcdEle, TString::Format("qcd_rate_ele_%s", charges[ic]))
                      : GetPar(fr, TString::Format("qcd_norm_ele_%s", regs[k].Data()));
        // a flavour not in this fit has no QCD parameter: 0,0 = "not in this fit"
        // (the helpers would return a plausible-looking 1 for an absent one)
        if (!useFl[0]) { qmu.v = 0.0; qmu.e = 0.0; }
        if (!useFl[1]) { qel.v = 0.0; qel.e = 0.0; }
        csv << regs[k] << "," << charges[ic] << "," << B << "," << iy << ","
            << rv[k] << "," << re[k] << "," << S[k] << "," << y << "," << e << ","
            << qmu.v << "," << qmu.e << "," << qel.v << "," << qel.e << ","
            << rz.v << "," << rz.e << "," << lum.v << "," << lum.e << ","
            << qModel << "," << res[k] << "\n";
        scsv << fitName << "," << pois[k] << "," << rv[k] << "," << re[k] << ","
             << S[k] << "," << y << "," << e << "\n";
        if (ok[k] && res[k] >= 0) // stat component: same r, stat error (source above)
          scsv << fitName << "_stat," << pois[k] << "," << rv[k] << "," << res[k] << ","
               << S[k] << "," << y << "," << res[k] * S[k] << "\n";
        if (okf[k])
          scsv << fitName << "_statonlyfit," << pois[k] << "," << rvf[k] << "," << ref[k] << ","
               << S[k] << "," << rvf[k] * S[k] << "," << ref[k] * S[k] << "\n";
      }
    std::cout << "[extract-simfit] " << fitName << ": 24 x r in ["
              << Form("%.3f, %.3f", rmin, rmax) << "]\n";
    if (maxDS > 1e-5)
      std::cout << "[extract-simfit] WARN " << fitName << ": prefit signal integrals of the input files differ from the"
                << " fitted channels' shapes_prefit by up to " << Form("%.2e", maxDS)
                << " (relative) -- are these the inputs that were fitted?\n";
    if (statAvail()) {
      // typical stat/total error ratio of the 24 POIs (the rest is the profiled
      // systematics: QCD lnN, lumi, the theory shapes, the lepton SFs)
      double sumRatio = 0; int nRatio = 0;
      for (int k = 0; k < NPOI; ++k) if (ok[k] && res[k] > 0 && re[k] > 0) { sumRatio += res[k] / re[k]; ++nRatio; }
      std::cout << "[extract-simfit] " << fitName << ": mean rErr_stat / rErr over the 24 POIs = "
                << Form("%.3f", nRatio ? sumRatio / nRatio : 0.0) << "\n";
    }
    if (nStatGtTot > 0)
      std::cout << "[extract-simfit] WARN " << fitName << ": " << nStatGtTot
                << " POI(s) with stat > total (worst " << worstStatPoi << " by "
                << Form("%+.2f%%", 100 * worstStatGtTot)
                << ") -- syst = sqrt(total^2 - stat^2) is floored at 0 there."
                << " The conditioning cannot do this; a separate refit can, so it points at a"
                << " non-parabolic direction or a minimizer difference. Inspect fit_statonly.log.\n";
    if (frs) {
      // the companion must sit at the same minimum (nuisances frozen at their
      // post-fit values) and, in the Gaussian approximation both share, return
      // the conditioned errors: report the largest deviation of each. Since
      // 2026-09-15b the companion IS the quoted stat error, so this is the
      // check that the conditioned (refit-free) answer agrees with it.
      std::cout << "[extract-simfit] " << fitName << "_statonlyfit vs conditioned covariance: max |r_fit - r| = "
                << Form("%.4f", maxDr)
                << (maxDr > 0.01 ? "  WARN (> 0.01: nuisances not frozen at the post-fit values?)" : "")
                << ", max |rErr_fit / rErr_conditioned - 1| = " << Form("%.4f", maxDe)
                << (maxDe > 0.02 ? "  WARN (> 2%: non-Gaussian likelihood or a mis-classified parameter)" : "") << "\n";
    }
    scsv << fitName << ",r_Z," << rz.v << "," << rz.e << ",,,\n";
    scsv << fitName << ",fit_status," << fr->status() << "," << fr->covQual() << ",,,\n";

    // ---- 24x24 fitted-yield covariance from the fit correlation matrix ------
    // cov(Y_a, Y_b) = S_a S_b rho_ab sigma_a sigma_b, fixed order
    // [Wp_y0..11, Wm_y0..11] (axis labels = "Wp_y0" ... for readability).
    fy->cd();
    const char *covName = (ib == 0) ? "h_cov_yield" : "h_cov_yield_FB";
    TH2D *hcov = new TH2D(covName, Form("%s (fitted-yield covariance);index;index", covName),
                          NPOI, 0, NPOI, NPOI, 0, NPOI);
    for (int a = 0; a < NPOI; ++a) {
      TString lab = pois[a]; lab.ReplaceAll("r_", "");
      hcov->GetXaxis()->SetBinLabel(a + 1, lab.Data());
      hcov->GetYaxis()->SetBinLabel(a + 1, lab.Data());
    }
    for (int a = 0; a < NPOI; ++a)
      for (int b = 0; b < NPOI; ++b) {
        if (!ok[a] || !ok[b]) continue;
        const double rho = (a == b) ? 1.0 : fr->correlation(pois[a].Data(), pois[b].Data());
        hcov->SetBinContent(a + 1, b + 1, S[a] * S[b] * rho * re[a] * re[b]);
      }
    hcov->SetDirectory(fy);
    fy->WriteTObject(hcov, covName, "Overwrite");

    // the STATISTICAL part of the same matrix (2026-09-15): h_cov_yield_stat /
    // h_cov_yield_FB_stat = S_a S_b cov_stat(r_a, r_b) from the conditioned
    // covariance -- the inner error bars of xsec_fiducial_comb (and any
    // downstream stat/syst split)
    if (statAvail()) {
      const TString covNameS = TString(covName) + "_stat";
      TH2D *hcovs = new TH2D(covNameS, Form("%s (fitted-yield covariance, stat component: constrained nuisances conditioned out);index;index", covNameS.Data()),
                             NPOI, 0, NPOI, NPOI, 0, NPOI);
      for (int a = 0; a < NPOI; ++a) {
        TString lab = pois[a]; lab.ReplaceAll("r_", "");
        hcovs->GetXaxis()->SetBinLabel(a + 1, lab.Data());
        hcovs->GetYaxis()->SetBinLabel(a + 1, lab.Data());
      }
      for (int a = 0; a < NPOI; ++a)
        for (int b = 0; b < NPOI; ++b) {
          if (!ok[a] || !ok[b]) continue;
          hcovs->SetBinContent(a + 1, b + 1, S[a] * S[b] * statCov(pois[a], pois[b]));
        }
      hcovs->SetDirectory(fy);
      fy->WriteTObject(hcovs, covNameS, "Overwrite");
    }

    // ---- 25x25 POI covariance INCLUDING r_Z (2026-09-15) --------------------
    // h_cov_poi[_FB][_stat]: cov(p_a, p_b) in PARAMETER space (not yield
    // space), order [r_Wp_y0..11, r_Wm_y0..11, r_Z], axis labels = the POI
    // names. h_cov_yield deliberately covers the 24 W POIs only, so the
    // sigma_W <-> sigma_Z cross term had no home; this matrix is what the
    // (sigma_W, sigma_Z) covariance ellipse needs -- sigma_W = Sum_i r_i
    // sigma_gen,i and sigma_Z = r_Z sigma_gen,Z are both LINEAR in these
    // parameters, so their 2x2 covariance is J V J^T with J the gen sigmas.
    // Written unconditionally (it needs no prefit template integrals, unlike
    // the yield matrix), so it is also the cleaner input for any future r
    // propagation.
    {
      std::vector<TString> pn(pois.begin(), pois.begin() + NPOI);
      pn.push_back("r_Z");
      const int NP = NPOI + 1;
      std::vector<double> pe(NP, 0);
      std::vector<bool> pok(NP, false);
      for (int a = 0; a < NPOI; ++a) { pe[a] = re[a]; pok[a] = ok[a]; }
      pe[NPOI] = rz.e; pok[NPOI] = rz.ok;
      const char *pcName = (ib == 0) ? "h_cov_poi" : "h_cov_poi_FB";
      for (int pass = 0; pass < 2; ++pass) {         // 0 = total, 1 = stat
        if (pass == 1 && !statAvail()) continue;
        const TString nm = TString(pcName) + (pass ? "_stat" : "");
        TH2D *hp = new TH2D(nm, Form("%s (POI covariance%s);;", nm.Data(),
                                     pass ? ", stat component: constrained nuisances conditioned out" : ""),
                            NP, 0, NP, NP, 0, NP);
        for (int a = 0; a < NP; ++a) {
          hp->GetXaxis()->SetBinLabel(a + 1, pn[a].Data());
          hp->GetYaxis()->SetBinLabel(a + 1, pn[a].Data());
        }
        for (int a = 0; a < NP; ++a)
          for (int b = 0; b < NP; ++b) {
            if (!pok[a] || !pok[b]) continue;
            const double v = pass ? statCov(pn[a], pn[b])
                                  : ((a == b) ? 1.0 : fr->correlation(pn[a].Data(), pn[b].Data())) * pe[a] * pe[b];
            hp->SetBinContent(a + 1, b + 1, v);
          }
        hp->SetDirectory(fy);
        fy->WriteTObject(hp, nm, "Overwrite");
      }
      if (rz.ok && ok[0])
        std::cout << "[extract-simfit] " << fitName << ": wrote " << pcName
                  << " (" << NP << "x" << NP << ", r_Z included); corr(r_Wp_y0, r_Z) = "
                  << Form("%+.3f", fr->correlation("r_Wp_y0", "r_Z")) << "\n";
    }

    // ---- covariance-propagated inclusive sums (diagnostics / AN numbers) ----
    // (total from the nominal fit's correlation matrix; stat from the
    // conditioned covariance; the companion fit's own numbers as the cross-check)
    auto sumWithCov = [&](const std::function<double(int, int)> &covr, const std::vector<double> &rr,
                          const std::vector<bool> &kk, int lo, int hi, double &val, double &err) { // [lo,hi)
      val = 0; double var = 0;
      for (int a = lo; a < hi; ++a) {
        if (!kk[a]) continue;
        val += rr[a] * S[a];
        for (int b = lo; b < hi; ++b) {
          if (!kk[b]) continue;
          var += S[a] * S[b] * covr(a, b);
        }
      }
      err = (var > 0) ? std::sqrt(var) : 0.0;
    };
    auto covTot = [&](int a, int b) {
      return ((a == b) ? 1.0 : fr->correlation(pois[a].Data(), pois[b].Data())) * re[a] * re[b];
    };
    double vWp, eWp, vWm, eWm, vW, eW;
    sumWithCov(covTot, rv, ok, 0, NB, vWp, eWp);
    sumWithCov(covTot, rv, ok, NB, NPOI, vWm, eWm);
    sumWithCov(covTot, rv, ok, 0, NPOI, vW, eW);
    scsv << fitName << ",Wp_sum_yield," << vWp << "," << eWp << ",,,\n";
    scsv << fitName << ",Wm_sum_yield," << vWm << "," << eWm << ",,,\n";
    scsv << fitName << ",W_sum_yield,"  << vW  << "," << eW  << ",,,\n";
    if (statAvail()) {
      auto covStat = [&](int a, int b) { return statCov(pois[a], pois[b]); };
      double vWps, eWps, vWms, eWms, vWs, eWs;
      sumWithCov(covStat, rv, ok, 0, NB, vWps, eWps);
      sumWithCov(covStat, rv, ok, NB, NPOI, vWms, eWms);
      sumWithCov(covStat, rv, ok, 0, NPOI, vWs, eWs);
      scsv << fitName << "_stat,Wp_sum_yield," << vWps << "," << eWps << ",,,\n";
      scsv << fitName << "_stat,Wm_sum_yield," << vWms << "," << eWms << ",,,\n";
      scsv << fitName << "_stat,W_sum_yield,"  << vWs  << "," << eWs  << ",,,\n";
      std::cout << "[extract-simfit] " << fitName << ": W sum yield "
                << Form("%.1f +/- %.1f total, +/- %.1f stat -> syst %.1f", vW, eW, eWs,
                        std::sqrt(std::max(0.0, eW * eW - eWs * eWs))) << "\n";
    }
    if (frs) {
      auto covFit = [&](int a, int b) {
        return ((a == b) ? 1.0 : frs->correlation(pois[a].Data(), pois[b].Data())) * ref[a] * ref[b];
      };
      double vWpf, eWpf, vWmf, eWmf, vWf, eWf;
      sumWithCov(covFit, rvf, okf, 0, NB, vWpf, eWpf);
      sumWithCov(covFit, rvf, okf, NB, NPOI, vWmf, eWmf);
      sumWithCov(covFit, rvf, okf, 0, NPOI, vWf, eWf);
      scsv << fitName << "_statonlyfit,Wp_sum_yield," << vWpf << "," << eWpf << ",,,\n";
      scsv << fitName << "_statonlyfit,Wm_sum_yield," << vWmf << "," << eWmf << ",,,\n";
      scsv << fitName << "_statonlyfit,W_sum_yield,"  << vWf  << "," << eWf  << ",,,\n";
    }

    if (fds) { fds->Close(); delete fds; fds = nullptr; }

    // ---- Asimov closure (present only when run with --asimov) ---------------
    TFile *fda = nullptr;
    RooFitResult *fra = OpenFitS(TString::Format("%s/%s/fitDiagnostics_%s_asimov.root",
                                                 fitsDir, fitName.Data(), fitName.Data()), fda);
    if (fra) {
      double worst = -1.0; TString worstName = "(none)";
      for (int a = 0; a < NPOI; ++a) {
        Par p = GetPar(fra, pois[a]);
        if (!p.ok) continue;
        scsv << fitName << "_asimov," << pois[a] << "," << p.v << "," << p.e << ",,,\n";
        if (std::fabs(p.v - 1.0) > worst) { worst = std::fabs(p.v - 1.0); worstName = pois[a]; }
      }
      Par pz = GetPar(fra, "r_Z");
      if (pz.ok) {
        scsv << fitName << "_asimov,r_Z," << pz.v << "," << pz.e << ",,,\n";
        if (std::fabs(pz.v - 1.0) > worst) { worst = std::fabs(pz.v - 1.0); worstName = "r_Z"; }
      }
      // constrained nuisances must come back at theta = 0 on the prefit Asimov
      const char *nuisN[5] = {"qcd_rate_mu_Wp", "qcd_rate_mu_Wm",
                              "qcd_rate_ele_Wp", "qcd_rate_ele_Wm", "lumi"};
      for (int in2 = 0; in2 < 5; ++in2) {
        Par t = GetPar(fra, nuisN[in2]);
        if (!t.ok) continue; // absent in free/legacy cards
        scsv << fitName << "_asimov," << nuisN[in2] << "_theta," << t.v << "," << t.e << ",,,\n";
        if (std::fabs(t.v) > worst) { worst = std::fabs(t.v); worstName = nuisN[in2]; }
      }
      // LHE shape nuisances (2026-09-07): the prefit Asimov is generated at
      // theta = 0, so every shape nuisance must come back at 0 too
      for (size_t il = 0; il < lheNames.size(); ++il) {
        Par t = GetPar(fra, lheNames[il]);
        if (!t.ok) continue;
        scsv << fitName << "_asimov," << lheNames[il] << "_theta," << t.v << "," << t.e << ",,,\n";
        if (std::fabs(t.v) > worst) { worst = std::fabs(t.v); worstName = lheNames[il]; }
      }
      // abcd-mode CR scales must come back at 1 (the CR templates hold the
      // prefit counts, so the Asimov is generated at scale 1; absent in
      // lnN/free cards and silently skipped there)
      const char *sclN[3] = {"B", "C", "D"};
      for (int ifl = 0; ifl < 2; ++ifl)
        for (int ic2 = 0; ic2 < 2; ++ic2)
          for (int is2 = 0; is2 < 3; ++is2) {
            const TString nn = TString::Format("qcd_s%s_%s_%s", sclN[is2],
                                               (ifl == 0) ? "mu" : "ele", charges[ic2]);
            Par s = GetPar(fra, nn);
            if (!s.ok) continue;
            scsv << fitName << "_asimov," << nn << "," << s.v << "," << s.e << ",,,\n";
            if (std::fabs(s.v - 1.0) > worst) { worst = std::fabs(s.v - 1.0); worstName = nn; }
          }
      const bool pass = (worst >= 0.0 && worst < 0.01);
      std::cout << "[asimov] " << fitName << " closure: max |POI-1| = " << Form("%.4f", worst)
                << " (" << worstName << ") => " << (pass ? "PASS" : "FAIL (tolerance 0.01)") << "\n";
      if (!pass)
        std::cout << "[asimov] WARN closure failed -- model/extraction inconsistency, inspect fit_asimov.log\n";
      fda->Close(); delete fda;
    }

    fd->Close(); delete fd;
  }

  csv.close(); scsv.close();
  fy->Close(); delete fy;
  if (wmu) { wmu->Close(); delete wmu; }
  if (wel) { wel->Close(); delete wel; }
  std::cout << "[extract-simfit] wrote " << yieldsRoot << "\n";
  std::cout << "[extract-simfit] wrote " << outDir << "/" << outTag << "_W_yields.csv, "
            << outTag << "_summary.csv\n";
}
