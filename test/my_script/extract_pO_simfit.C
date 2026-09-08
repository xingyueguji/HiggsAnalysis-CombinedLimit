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
// Nuisance thetas (pulls) are also dumped to comb_summary.csv (<name>_theta
// rows) -- THE check that the ABCD prediction and its assigned uncertainty are
// consistent with the data (|pull| ~> 1 means kappa too small or template biased).
// A final sweep prints any floating parameter of fit_s not reported above, so
// a new nuisance can never be silently invisible.
//
// It writes into <outDir> (both variants into the same files):
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
#include <cmath>
#include <fstream>
#include <iostream>
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
                       const char *lheSysts = "") // comma list of the LHE shape
                                                  //  nuisances in the cards ("" = none)
{
  gSystem->mkdir(outDir, kTRUE);

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

  TFile *wmu = TFile::Open(muWFile, "READ");
  TFile *wel = TFile::Open(eleWFile, "READ");
  if (!wmu || wmu->IsZombie() || !wel || wel->IsZombie()) {
    std::cerr << "[extract-simfit] WARN cannot open W input file(s) -- aborting extraction\n";
    return;
  }
  auto sigPrefit = [&](TFile *fin, const TString &region) -> double {
    TH1 *h = (TH1 *)fin->Get(region + "/signal");
    return h ? h->Integral() : -1.0; // bins 1..N, matches Combine rate convention
  };

  TString yieldsRoot = TString::Format("%s/comb_fitted_yields.root", outDir);
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

  std::ofstream csv(TString::Format("%s/comb_W_yields.csv", outDir).Data());
  csv << "region,charge,binning,ybin,r,rErr,signal_prefit,fitted_yield,fitted_yield_err,"
         "qcd_norm_mu,qcd_norm_muErr,qcd_norm_ele,qcd_norm_eleErr,r_Z,r_ZErr,lumi,lumiErr,"
         "qcd_model\n";

  std::ofstream scsv(TString::Format("%s/comb_summary.csv", outDir).Data());
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
    std::vector<double> rv(NPOI, 0), re(NPOI, 0), S(NPOI, 0);
    std::vector<bool> ok(NPOI, false);
    double rmin = 1e30, rmax = -1e30;
    for (int ic = 0; ic < 2; ++ic)
      for (int iy = 0; iy < NB; ++iy) {
        const int k = ic * NB + iy;
        pois.push_back(TString::Format("r_%s_y%d", charges[ic], iy));
        regs.push_back(TString::Format("%s_%s_y%d", charges[ic], B, iy));
        Par p = GetPar(fr, pois[k]);
        const double smu = sigPrefit(wmu, regs[k]), sel = sigPrefit(wel, regs[k]);
        if (smu < 0 || sel < 0)
          std::cerr << "[extract-simfit] WARN missing prefit signal for " << regs[k] << "\n";
        S[k]  = (smu > 0 ? smu : 0) + (sel > 0 ? sel : 0);
        ok[k] = p.ok && S[k] > 0;
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
        csv << regs[k] << "," << charges[ic] << "," << B << "," << iy << ","
            << rv[k] << "," << re[k] << "," << S[k] << "," << y << "," << e << ","
            << qmu.v << "," << qmu.e << "," << qel.v << "," << qel.e << ","
            << rz.v << "," << rz.e << "," << lum.v << "," << lum.e << ","
            << qModel << "\n";
        scsv << fitName << "," << pois[k] << "," << rv[k] << "," << re[k] << ","
             << S[k] << "," << y << "," << e << "\n";
      }
    std::cout << "[extract-simfit] " << fitName << ": 24 x r in ["
              << Form("%.3f, %.3f", rmin, rmax) << "]\n";
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

    // ---- covariance-propagated inclusive sums (diagnostics / AN numbers) ----
    auto sumWithCov = [&](int lo, int hi, double &val, double &err) { // [lo,hi)
      val = 0; double var = 0;
      for (int a = lo; a < hi; ++a) {
        if (!ok[a]) continue;
        val += rv[a] * S[a];
        for (int b = lo; b < hi; ++b) {
          if (!ok[b]) continue;
          const double rho = (a == b) ? 1.0 : fr->correlation(pois[a].Data(), pois[b].Data());
          var += S[a] * S[b] * rho * re[a] * re[b];
        }
      }
      err = (var > 0) ? std::sqrt(var) : 0.0;
    };
    double vWp, eWp, vWm, eWm, vW, eW;
    sumWithCov(0, NB, vWp, eWp);
    sumWithCov(NB, NPOI, vWm, eWm);
    sumWithCov(0, NPOI, vW, eW);
    scsv << fitName << ",Wp_sum_yield," << vWp << "," << eWp << ",,,\n";
    scsv << fitName << ",Wm_sum_yield," << vWm << "," << eWm << ",,,\n";
    scsv << fitName << ",W_sum_yield,"  << vW  << "," << eW  << ",,,\n";

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
  wmu->Close(); delete wmu;
  wel->Close(); delete wel;
  std::cout << "[extract-simfit] wrote " << yieldsRoot << "\n";
  std::cout << "[extract-simfit] wrote " << outDir << "/comb_W_yields.csv, comb_summary.csv\n";
}
