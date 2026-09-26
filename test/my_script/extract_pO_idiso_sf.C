// =============================================================================
// extract_pO_idiso_sf.C -- extraction of the ELECTRON ID+ISO SF FIT
// (2026-09-24; driver run_pO_idiso_sf.sh, cards make_pO_idiso_sf_cards.sh).
// SEPARATE from the nominal fit stream: reads only this fit's files.
//
// The efficiency of eleMVAIdWP90 && eleMVAIsoWP90 on W -> e nu events, from the
// POST-FIT W counts of the pass and fail channels (disjoint samples, one
// likelihood; total = pass + fail):
//   N_pass = r_pass x S_pass ,  N_fail = r_fail x S_fail   (S = the W signal
//            template integral over the fitted bins, from the input file)
//   eps_data = N_pass / (N_pass + N_fail) ,  eps_MC = S_pass / (S_pass + S_fail)
//   SF = eps_data / eps_MC
// per charge and per coarse bin k, charge-combined per bin, and combined over
// all bins of the scheme ("all"). Errors: the fit covariance of the r's
// (fr->correlation x errors), propagated through the exact gradient of eps --
// the pattern of extract_pO_simfit.C::AbcdScale. The MC-statistical error of
// S is neglected (MC samples ~100x the data).
// The INCLUSIVE Z TAG-AND-PROBE (2026-09-25) is its own fit (input scheme
// "ztnp"), extracted into a single row bin = "ztnp": a Z_PP event holds two
// passing probes, a Z_PF event one failing probe, so
//   eps = 2 r_ZPP S_PP / (2 r_ZPP S_PP + r_ZPF S_PF),  eps_MC = 2 S_PP / (2 S_PP + S_PF)
// (its CSV N_pass / S_pass count PROBES, i.e. 2x the Z_PP events).
// QCD (the nominal in-fit ABCD): per SR channel the multiplier kappa^theta sB sC
// / sD with its error is printed; bkg_norm_* too.
// Asimov closure (when the _asimov fitDiagnostics exists): every r and every
// ABCD scale / rateParam must return 1 and every lnN theta 0, hence SF = 1 in
// every bin -- printed as PASS/FAIL.
//
// Outputs in <outDir>:
//   idiso_sf_<tag>.csv          one row per (bin, charge) + the charge-combined
//                               rows (charge = all) + the scheme-combined row
//                               (bin = all): columns below
//   idiso_sf_<tag>_summary.csv  every floating parameter of fit_s (+ Asimov)
// Run by run_pO_idiso_sf.sh, which tees the console to
// summary/extract_idiso_sf_<tag>.log (the record: fit quality, closure, WARNs).
//
//   root -b -q 'extract_pO_idiso_sf.C("fitDiag.root","fitDiag_asimov.root|none",
//               "combine_input_idiso_<tag>.root","..._meta.txt","summary","<tag>"[,kappa])'
//   kappa = the residual lnN on the SR qcd the card was built with (default 1.15)
// =============================================================================
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
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

// nullptr when the file does not exist / has no fit_s; caller owns *f.
RooFitResult *OpenFitS(const TString &path, TFile *&f) {
  f = nullptr;
  if (path == "none" || gSystem->AccessPathName(path)) return nullptr;
  f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) { if (f) { f->Close(); delete f; f = nullptr; } return nullptr; }
  return (RooFitResult *)f->Get("fit_s");
}

// var(sum_a g_a x_a) over parameters with errors, correlations from the fit
double Var(const RooFitResult *fr, const std::vector<TString> &names,
           const std::vector<double> &grad, const std::vector<double> &err) {
  double v = 0.0;
  for (size_t a = 0; a < names.size(); ++a) {
    if (err[a] <= 0.0) continue;
    for (size_t b = 0; b < names.size(); ++b) {
      if (err[b] <= 0.0) continue;
      const double rho = (a == b) ? 1.0 : fr->correlation(names[a].Data(), names[b].Data());
      v += grad[a] * grad[b] * rho * err[a] * err[b];
    }
  }
  return v > 0.0 ? v : 0.0;
}

// one (r, S) term of a pass or a fail sum
struct Term { TString name; double r, rErr, S; bool pass; };

// eps = sum_pass r S / sum_all r S, with its propagated error
void Eps(const RooFitResult *fr, const std::vector<Term> &t, double &eps, double &err,
         double &nPass, double &nFail, double &epsMC) {
  nPass = nFail = 0.0; double sP = 0.0, sF = 0.0;
  for (const auto &x : t) { (x.pass ? nPass : nFail) += x.r * x.S; (x.pass ? sP : sF) += x.S; }
  const double tot = nPass + nFail;
  eps = tot > 0 ? nPass / tot : 0.0;
  epsMC = (sP + sF) > 0 ? sP / (sP + sF) : 0.0;
  std::vector<TString> n; std::vector<double> g, e;
  for (const auto &x : t) {
    n.push_back(x.name); e.push_back(x.rErr);
    // d eps / d r = S nFail / tot^2 (pass term) or -S nPass / tot^2 (fail term)
    g.push_back(tot > 0 ? (x.pass ? x.S * nFail : -x.S * nPass) / (tot * tot) : 0.0);
  }
  err = std::sqrt(Var(fr, n, g, e));
}

std::string MetaValue(const std::string &meta, const std::string &key) {
  std::ifstream in(meta);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string k; ss >> k;
    if (k == key) { std::string rest; std::getline(ss, rest); return rest.empty() ? "" : rest.substr(1); }
  }
  return "";
}

} // namespace

void extract_pO_idiso_sf(const char *fitDiag, const char *fitDiagAsimov, const char *inputFile,
                         const char *metaFile, const char *outDir, const char *tag, double kappa = 1.15)
{
  const std::string scheme = MetaValue(metaFile, "scheme"), disc = MetaValue(metaFile, "disc"),
                    sbwin = MetaValue(metaFile, "sbwin");
  const int nb = std::atoi(MetaValue(metaFile, "nbins").c_str());
  if (nb <= 0) { std::cerr << "[idiso-sf] FAIL no nbins in " << metaFile << "\n"; return; }
  std::vector<double> lo(nb, 0.0), hi(nb, 0.0);
  {
    std::ifstream in(metaFile);
    std::string line;
    while (std::getline(in, line)) {
      std::istringstream ss(line);
      std::string k; int ib; double a, b;
      if ((ss >> k) && k == "bin" && (ss >> ib >> a >> b) && ib >= 0 && ib < nb) { lo[ib] = a; hi[ib] = b; }
    }
  }

  TFile *fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) { std::cerr << "[idiso-sf] FAIL cannot open " << inputFile << "\n"; return; }
  TFile *ff = nullptr, *fa = nullptr;
  RooFitResult *fr = OpenFitS(fitDiag, ff);
  if (!fr) { std::cerr << "[idiso-sf] FAIL no fit_s in " << fitDiag << "\n"; return; }
  RooFitResult *fra = OpenFitS(fitDiagAsimov, fa);

  const char *chg[2] = {"Wp", "Wm"}, *cat[2] = {"pass", "fail"};
  std::cout << "[idiso-sf] " << tag << ": scheme " << scheme << ", disc " << disc << ", QCD sideband " << sbwin
            << ", " << nb << " bin(s); fit status " << fr->status() << ", covQual " << fr->covQual() << "\n";
  if (fr->status() != 0 || fr->covQual() < 3)
    std::cout << "[idiso-sf] WARN fit not clean (status " << fr->status() << ", covQual " << fr->covQual() << ")\n";

  // S and r per (charge, category, bin)
  auto S = [&](int q, int c, int k) {
    TH1D *h = (TH1D *)fin->Get(TString::Format("%s_%s_k%d/signal", chg[q], cat[c], k));
    return h ? h->Integral(1, h->GetNbinsX()) : 0.0; // the fitted bins only (Combine ignores under/overflow)
  };
  auto T = [&](int q, int c, int k) {
    Term t; t.name = TString::Format("r_%s_%s_k%d", chg[q], cat[c], k);
    const Par p = GetPar(fr, t.name);
    if (!p.ok) std::cout << "[idiso-sf] WARN " << t.name << " not in fit_s -> taken as 1 +- 0\n";
    t.r = p.ok ? p.v : 1.0; t.rErr = p.ok ? p.e : 0.0; t.S = S(q, c, k); t.pass = (c == 0);
    return t;
  };

  gSystem->mkdir(outDir, kTRUE);
  const TString csv = TString::Format("%s/idiso_sf_%s.csv", outDir, tag);
  std::ofstream o(csv.Data());
  o << "tag,scheme,disc,sbwin,bin,lo,hi,charge,S_pass,S_fail,eps_mc,r_pass,r_pass_err,r_fail,r_fail_err,"
       "rho_pass_fail,N_pass,N_fail,eps_data,eps_data_err,sf,sf_err,status,covQual\n";
  auto row = [&](const std::string &bin, double l, double h, const std::string &q, const std::vector<Term> &t) {
    double eps, err, nP, nF, eMC;
    Eps(fr, t, eps, err, nP, nF, eMC);
    double sP = 0, sF = 0;
    for (const auto &x : t) (x.pass ? sP : sF) += x.S;
    // single (pass, fail) pair: report its r's and their correlation
    double rp = -1, rpe = 0, rf = -1, rfe = 0, rho = 0;
    if (t.size() == 2) {
      rp = t[0].r; rpe = t[0].rErr; rf = t[1].r; rfe = t[1].rErr;
      if (rpe > 0 && rfe > 0) rho = fr->correlation(t[0].name.Data(), t[1].name.Data());
    }
    const double sf = eMC > 0 ? eps / eMC : 0.0, sfe = eMC > 0 ? err / eMC : 0.0;
    o << tag << "," << scheme << "," << disc << "," << sbwin << "," << bin << "," << l << "," << h << "," << q << ","
      << sP << "," << sF << "," << eMC << "," << rp << "," << rpe << "," << rf << "," << rfe << "," << rho << ","
      << nP << "," << nF << "," << eps << "," << err << "," << sf << "," << sfe << ","
      << fr->status() << "," << fr->covQual() << "\n";
    std::cout << TString::Format("[idiso-sf] %-6s bin %-3s %-3s  eps_data %.4f +- %.4f  eps_MC %.4f  SF %.4f +- %.4f"
                                 "   (N_pass %.0f, N_fail %.0f)\n",
                                 tag, bin.c_str(), q.c_str(), eps, err, eMC, sf, sfe, nP, nF).Data();
  };

  if (scheme == "ztnp") {
    // the Z tag-and-probe fit (its own input and card): eps = 2 N_PP / (2 N_PP + N_PF)
    auto SZ = [&](const char *dir) {
      TH1D *h = (TH1D *)fin->Get(TString::Format("%s/signal", dir));
      return h ? h->Integral(1, h->GetNbinsX()) : 0.0;
    };
    const double sPP = SZ("Z_PP"), sPF = SZ("Z_PF");
    if (!(sPP > 0 && sPF > 0)) { std::cerr << "[idiso-sf] FAIL no Z_PP/Z_PF signal in " << inputFile << "\n"; return; }
    const Par pp = GetPar(fr, "r_ZPP"), pf = GetPar(fr, "r_ZPF");
    if (!pp.ok || !pf.ok) std::cout << "[idiso-sf] WARN r_ZPP or r_ZPF not in fit_s -> taken as 1 +- 0\n";
    Term tp, tf;
    tp.name = "r_ZPP"; tp.r = pp.ok ? pp.v : 1.0; tp.rErr = pp.ok ? pp.e : 0.0; tp.S = 2.0 * sPP; tp.pass = true;
    tf.name = "r_ZPF"; tf.r = pf.ok ? pf.v : 1.0; tf.rErr = pf.ok ? pf.e : 0.0; tf.S = sPF;       tf.pass = false;
    row("ztnp", 0, 0, "all", {tp, tf});
  } else {
    std::vector<Term> all;
    for (int k = 0; k < nb; ++k) {
      std::vector<Term> both;
      for (int q = 0; q < 2; ++q) {
        std::vector<Term> one = {T(q, 0, k), T(q, 1, k)};
        row(std::to_string(k), lo[k], hi[k], chg[q], one);
        both.insert(both.end(), one.begin(), one.end());
      }
      row(std::to_string(k), lo[k], hi[k], "all", both);
      all.insert(all.end(), both.begin(), both.end());
    }
    if (nb > 1) row("all", lo[0], hi[nb - 1], "all", all);
  }
  o.close();

  // QCD: the in-fit ABCD multiplier of every SR channel, M = kappa^theta sB sC / sD
  // (the formula rateParam qcd_abcd_ele_R is a RooFormulaVar, absent from
  // floatParsFinal -- evaluated here from its floating constituents with the
  // full correlation, the pattern of extract_pO_simfit.C::AbcdScale)
  if (scheme != "ztnp") {
    for (int k = 0; k < nb; ++k)
      for (int q = 0; q < 2; ++q)
        for (int c = 0; c < 2; ++c) {
          const TString R = TString::Format("%s_%s_k%d", chg[q], cat[c], k);
          const std::vector<TString> n = {TString::Format("qcd_rate_ele_%s_%s", chg[q], cat[c]), "qcd_sB_ele_" + R,
                                          "qcd_sC_ele_" + R, "qcd_sD_ele_" + R};
          std::vector<Par> p;
          for (const auto &x : n) p.push_back(GetPar(fr, x));
          if (!p[1].ok || !p[2].ok || !p[3].ok || p[3].v <= 0) {
            std::cout << "[idiso-sf] WARN no in-fit ABCD scales for " << R << " in fit_s\n";
            continue;
          }
          const double th = p[0].ok ? p[0].v : 0.0;
          const double M = std::pow(kappa, th) * p[1].v * p[2].v / p[3].v;
          const std::vector<double> g = {M * std::log(kappa), M / p[1].v, M / p[2].v, -M / p[3].v};
          const std::vector<double> e = {p[0].ok ? p[0].e : 0.0, p[1].e, p[2].e, p[3].e};
          std::cout << TString::Format("[idiso-sf] %s QCD in-fit ABCD %-12s sB %.3f +- %.3f  sC %.3f +- %.3f  sD %.3f +- %.3f"
                                       "  theta %+.2f +- %.2f  -> SR multiplier %.3f +- %.3f\n",
                                       tag, R.Data(), p[1].v, p[1].e, p[2].v, p[2].e, p[3].v, p[3].e, th,
                                       p[0].ok ? p[0].e : 0.0, M, std::sqrt(Var(fr, n, g, e))).Data();
        }
  }
  for (auto *a : fr->floatParsFinal()) {
    RooRealVar *v = (RooRealVar *)a;
    const TString n = v->GetName();
    if (n.BeginsWith("bkg_norm_"))
      std::cout << TString::Format("[idiso-sf] %s %-28s = %.3f +- %.3f\n", tag, n.Data(), v->getVal(), v->getError()).Data();
  }

  // every floating parameter (+ the Asimov closure)
  const TString sum = TString::Format("%s/idiso_sf_%s_summary.csv", outDir, tag);
  std::ofstream os(sum.Data());
  os << "fit,param,value,error\n";
  os << "nominal,status," << fr->status() << ",0\n" << "nominal,covQual," << fr->covQual() << ",0\n";
  for (auto *a : fr->floatParsFinal()) {
    RooRealVar *v = (RooRealVar *)a;
    os << "nominal," << v->GetName() << "," << v->getVal() << "," << v->getError() << "\n";
  }
  const Par rz = GetPar(fr, "r_Z");
  if (rz.ok) std::cout << TString::Format("[idiso-sf] %s r_Z = %.4f +- %.4f (Z_PP, both legs pass)\n", tag, rz.v, rz.e).Data();
  if (fra) {
    // every POI / rateParam returns 1, every lnN theta 0
    double maxDev = 0.0; TString worst;
    for (auto *a : fra->floatParsFinal()) {
      RooRealVar *v = (RooRealVar *)a;
      os << "asimov," << v->GetName() << "," << v->getVal() << "," << v->getError() << "\n";
      const TString n = v->GetName();
      double dev = -1.0;
      if (n.BeginsWith("r_") || n.BeginsWith("qcd_s") || n.BeginsWith("bkg_norm_")) dev = std::fabs(v->getVal() - 1.0);
      else if (n.BeginsWith("qcd_rate_")) dev = std::fabs(v->getVal());
      if (dev > maxDev) { maxDev = dev; worst = n; }
    }
    std::cout << TString::Format("[asimov] %s: max |POI - 1|, |rateParam - 1|, |theta| = %.2e (%s) -> %s   (status %d, covQual %d)\n",
                                 tag, maxDev, worst.Data(), maxDev < 1e-3 ? "PASS" : "FAIL", fra->status(), fra->covQual()).Data();
    // the EXPECTED precision: eps and SF from the Asimov fit (every r = 1, so SF = 1 +- its expected error).
    // On data a poorly constrained r_fail can sit on its 0 boundary, where the fitted error means
    // nothing -- this is the number that says how well the channel can measure the efficiency.
    auto expected = [&](const std::string &bin, std::vector<Term> t) {
      for (auto &x : t) {
        const Par p = GetPar(fra, x.name);
        x.r = p.ok ? p.v : 1.0; x.rErr = p.ok ? p.e : 0.0;
      }
      double eps, err, nP, nF, eMC;
      Eps(fra, t, eps, err, nP, nF, eMC);
      std::cout << TString::Format("[expected] %-6s bin %-4s all  eps_MC %.4f  expected eps error %.4f -> SF 1 +- %.4f\n",
                                   tag, bin.c_str(), eMC, err, eMC > 0 ? err / eMC : 0.0).Data();
    };
    if (scheme == "ztnp") {
      auto SZ = [&](const char *dir) {
        TH1D *h = (TH1D *)fin->Get(TString::Format("%s/signal", dir));
        return h ? h->Integral(1, h->GetNbinsX()) : 0.0;
      };
      Term tp, tf;
      tp.name = "r_ZPP"; tp.S = 2.0 * SZ("Z_PP"); tp.pass = true;
      tf.name = "r_ZPF"; tf.S = SZ("Z_PF");       tf.pass = false;
      expected("ztnp", {tp, tf});
    } else {
      std::vector<Term> all;
      for (int k = 0; k < nb; ++k) {
        std::vector<Term> both = {T(0, 0, k), T(0, 1, k), T(1, 0, k), T(1, 1, k)};
        expected(std::to_string(k), both);
        all.insert(all.end(), both.begin(), both.end());
      }
      if (nb > 1) expected("all", all);
    }
  } else {
    std::cout << "[asimov] " << tag << ": no Asimov fit (run with --asimov for the closure)\n";
  }
  os.close();
  std::cout << "[idiso-sf] wrote " << csv << " and " << sum << "\n";
  if (fa) fa->Close();
  ff->Close();
  fin->Close();
}
