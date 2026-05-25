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

static double GetFitPar(RooFitResult* fr, const char* name, double& err) {
  err = 0.0;
  if (!fr)
    return 0.0;

  RooArgList pars = fr->floatParsFinal();
  RooRealVar* v = dynamic_cast<RooRealVar*>(pars.find(name));
  if (!v) {
    std::cerr << "[WARN] Cannot find fit parameter: " << name << "\n";
    return 0.0;
  }

  err = v->getError();
  return v->getVal();
}

static TH1* RemapToReferenceBinning(TH1* hIn, TH1* hRef, const char* newName) {
  if (!hIn || !hRef)
    return nullptr;

  TH1* hOut = dynamic_cast<TH1*>(hRef->Clone(newName));
  if (!hOut)
    return nullptr;

  hOut->Reset();
  hOut->SetDirectory(nullptr);

  const int nIn = hIn->GetNbinsX();
  const int nRef = hRef->GetNbinsX();

  if (nIn != nRef) {
    std::cerr << "[WARN] Bin count mismatch in RemapToReferenceBinning: " << nIn << " vs " << nRef << "\n";
  }

  const int nCopy = std::min(nIn, nRef);
  for (int ib = 1; ib <= nCopy; ++ib) {
    hOut->SetBinContent(ib, hIn->GetBinContent(ib));
    hOut->SetBinError(ib, hIn->GetBinError(ib));
  }

  return hOut;
}

static TH1D *MakePullHist(TH1 *hData, TH1 *hPost, const char *name = "h_pull")
{
    if (!hData || !hPost)
        return nullptr;

    TH1D *hPull = dynamic_cast<TH1D *>(hData->Clone(name));
    if (!hPull)
        return nullptr;

    hPull->Reset();
    hPull->SetDirectory(nullptr);

    const int nb = std::min(hData->GetNbinsX(), hPost->GetNbinsX());
    for (int ib = 1; ib <= nb; ++ib)
    {
        const double d  = hData->GetBinContent(ib);
        const double m  = hPost->GetBinContent(ib);
        double err      = hData->GetBinError(ib);

        if (err <= 0.0)
        {
            // fallback for empty/zero-error bins
            err = (d > 0.0) ? std::sqrt(d) : 1.0;
        }

        const double pull = (d - m) / err;
        hPull->SetBinContent(ib, pull);
        hPull->SetBinError(ib, 0.0);
    }

    return hPull;
}

static void SaveNicePlot1D_WithBkg_Postfit_Pull(
    TH1 *hData,
    const std::vector<TH1 *> &bkgs,
    const std::vector<std::string> &bkgNames,
    TH1 *hTotal,
    const std::string &outPathNoExt,
    const std::string &xTitle,
    const std::string &yTitle,
    const std::string &mainTitle,
    const std::string &subTitle1,
    const std::string &subTitle2,
    const std::vector<std::string> &boxLines,
    const PlotStyle &ps = PlotStyle(),
    PlotTuner tuner = nullptr)
{
    if (!hData || !hTotal)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas *c = new TCanvas(Form("c_%s_postfit_pull", hData->GetName()), "", ps.w, ps.h);

    // ---------------- top / bottom pads ----------------
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.25);

    pad1->SetLeftMargin(ps.lm);
    pad1->SetRightMargin(ps.rm);
    pad1->SetTopMargin(ps.tm);
    pad1->SetBottomMargin(0.02);
    pad1->SetTicks(1, 1);
    pad1->SetLogy(ps.logy);

    pad2->SetLeftMargin(ps.lm);
    pad2->SetRightMargin(ps.rm);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);
    pad2->SetTicks(1, 1);

    c->cd();
    pad1->Draw();
    pad2->Draw();

    // ---------------- top pad ----------------
    pad1->cd();

    THStack *hs = new THStack("hs_postfit_pull", "");

    std::vector<int> colors = {
        kAzure - 9,
        kOrange - 3,
        kGreen + 2,
        kMagenta - 3,
        kCyan + 1};

    for (int i = static_cast<int>(bkgs.size()) - 1; i >= 0; --i)
    {
        TH1 *b = bkgs[i];
        if (!b)
            continue;

        b->SetFillColor(colors[i % colors.size()]);
        b->SetLineColor(kBlack);
        b->SetLineWidth(1);
        hs->Add(b);
    }

    ApplyHistStyle(hData, ps, xTitle, yTitle);
    hData->GetXaxis()->SetLabelSize(0.0);
    hData->GetXaxis()->SetTitleSize(0.0);

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.2);
    hData->SetLineColor(kBlack);
    hData->SetMarkerColor(kBlack);

    hData->Draw("E");
    hs->Draw("HIST SAME");

    hTotal->SetLineColor(kRed + 1);
    hTotal->SetLineWidth(3);
    hTotal->SetFillStyle(0);
    hTotal->Draw("HIST SAME");

    hData->Draw("E SAME");

    double ymax = hData->GetMaximum();
    ymax = std::max(ymax, hTotal->GetMaximum());
    for (auto *b : bkgs)
    {
        if (b)
            ymax = std::max(ymax, b->GetMaximum());
    }
    hData->SetMaximum(ps.logy ? 10.25 * ymax : 1.35 * ymax);
    if (ps.logy)
        hData->SetMinimum(1.0);

    TLegend *leg = new TLegend(ps.boxX1 + 0.4, ps.boxY1 - 0.22, ps.boxX2 + 0.4, ps.boxY2 - 0.12);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.032);

    leg->AddEntry(hData, "Data", "lep");
    for (size_t i = 0; i < bkgs.size(); ++i)
    {
        if (bkgs[i])
            leg->AddEntry(bkgs[i], bkgNames[i].c_str(), "f");
    }
    leg->AddEntry(hTotal, "Postfit total", "l");
    leg->Draw();

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    if (tuner)
        tuner((TCanvas *)pad1, hData);

    CMS_lumi(pad1, 13, 10);

    // ---------------- bottom pad ----------------
    pad2->cd();

    TH1D *hPull = MakePullHist(hData, hTotal, "h_pull");
    if (!hPull)
        return;

    hPull->SetTitle("");
    hPull->SetMarkerStyle(20);
    hPull->SetMarkerSize(0.9);
    hPull->SetLineColor(kBlack);
    hPull->SetMarkerColor(kBlack);

    hPull->GetYaxis()->SetTitle("Pull");
    hPull->GetXaxis()->SetTitle(xTitle.c_str());

    hPull->GetXaxis()->SetTitleFont(42);
    hPull->GetYaxis()->SetTitleFont(42);
    hPull->GetXaxis()->SetLabelFont(42);
    hPull->GetYaxis()->SetLabelFont(42);

    hPull->GetXaxis()->SetTitleSize(0.12);
    hPull->GetYaxis()->SetTitleSize(0.10);
    hPull->GetXaxis()->SetLabelSize(0.10);
    hPull->GetYaxis()->SetLabelSize(0.09);

    hPull->GetXaxis()->SetTitleOffset(1.15);
    hPull->GetYaxis()->SetTitleOffset(0.55);

    hPull->GetYaxis()->CenterTitle(true);
    hPull->GetYaxis()->SetNdivisions(505);
    hPull->GetYaxis()->SetRangeUser(-5.0, 5.0);

    hPull->Draw("EP");

    TLine *l0 = new TLine(hPull->GetXaxis()->GetXmin(), 0.0,
                          hPull->GetXaxis()->GetXmax(), 0.0);
    l0->SetLineStyle(2);
    l0->SetLineWidth(2);
    l0->Draw("SAME");

    TLine *lp2 = new TLine(hPull->GetXaxis()->GetXmin(), 2.0,
                           hPull->GetXaxis()->GetXmax(), 2.0);
    TLine *lm2 = new TLine(hPull->GetXaxis()->GetXmin(), -2.0,
                           hPull->GetXaxis()->GetXmax(), -2.0);
    lp2->SetLineStyle(3);
    lm2->SetLineStyle(3);
    lp2->SetLineColor(kRed + 1);
    lm2->SetLineColor(kRed + 1);
    lp2->Draw("SAME");
    lm2->Draw("SAME");

    c->cd();
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

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