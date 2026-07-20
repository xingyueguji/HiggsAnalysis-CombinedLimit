#include "CMS_lumi.C"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TStyle.h"
#include "TString.h"
#include <algorithm>
#include <cmath>

// Global plot style: thicken the frame (the box / "bezel" around every plot).
// This runs once when the header is loaded, so it applies to ALL pads and
// canvases — including the ratio sub-pads in SaveDataMCRatio and any plot that
// doesn't go through ApplyCanvasStyle. Change the number here to taste;
// PlotStyle::FrameLineWidth below mirrors it for the per-canvas path.
namespace
{
struct ApplyGlobalPlotStyle
{
    ApplyGlobalPlotStyle()
    {
        if (gStyle)
            gStyle->SetFrameLineWidth(3);
    }
};
const ApplyGlobalPlotStyle kApplyGlobalPlotStyle;
} // namespace

// -----------------------------
// 1) Generic style/options blob
// -----------------------------
struct PlotStyle
{
    // canvas
    int w = 800, h = 800;
    double lm = 0.14, rm = 0.05, bm = 0.12, tm = 0.07;
    bool ticks = true;

    // fonts
    int font = 42;
    double titleSize = 0.050;
    double subSize = 0.040;
    double boxTextSize = 0.034;

    // header positions (NDC)
    double headerX = 0.6;
    double headerY = 0.93;  // main title top
    double headerDy = 0.06; // spacing between lines

    // box (NDC)
    double boxX1 = 0.55, boxY1 = 0.42, boxX2 = 0.93, boxY2 = 0.88;

    // legend rectangle (NDC) for SaveNicePlot1D_WithBkg; <0 -> auto = derived
    // from the info box (boxX1, boxY1-0.2, boxX2, boxY2-0.2). Set explicitly to
    // decouple the legend from the info box (e.g. legend lower-right on postfit).
    double legX1 = -1, legY1 = -1, legX2 = -1, legY2 = -1;

    // axis sizes
    double xTitleSize = 0.045, yTitleSize = 0.045;
    double xLabelSize = 0.040, yLabelSize = 0.040;
    double xTitleOffset = 1.10, yTitleOffset = 1.35;
    double FrameLineWidth = 3; // mirror of the global gStyle->SetFrameLineWidth above

    // draw options
    std::string drawOpt = "E"; // e.g. "E", "hist", "E1", etc.
    std::string drawOptGraph = "AP";
    bool logy = false;

    // stat box
    bool showStats = false;

    // background/MC normalization in SaveNicePlot1D_WithBkg:
    //   true  -> scale the MC stack to the data integral (shape comparison)
    //   false -> draw the MC at its absolute (already-scaled) yield
    bool normBkgToData = true;

    // pull sub-pad in SaveNicePlot1D_WithBkg: per-bin (data - MC_total)/sigma,
    // sigma^2 = MC_total + sigma_MCstat^2 (Poisson variance of the data around
    // the prediction + MC statistical error), drawn as zero-anchored bars
    // under the main pad. The canvas stays NEAR-SQUARE: its height grows only
    // by pullCanvasScale (not by the full pull-pad fraction), so the main pad
    // compresses a little vertically. All fonts/positions are pad-relative, so
    // the main-plot layout scales consistently; only the aspect changes.
    bool pullPad = false;      // enable the pull sub-pad
    double pullSplit = 0.30;   // fraction of the canvas given to the pull pad
    double pullCanvasScale = 1.125; // canvas-height multiplier when pullPad is on
    double pullYRange = -1.0;  // fixed symmetric |pull| axis; <0 -> auto (max(3, 1.15*max|pull|))
};

// -----------------------------
// 2) Basic helpers
// -----------------------------
static void ApplyCanvasStyle(TCanvas *c, const PlotStyle &ps)
{
    c->SetLeftMargin(ps.lm);
    c->SetRightMargin(ps.rm);
    c->SetBottomMargin(ps.bm);
    c->SetTopMargin(ps.tm);
    c->SetFrameLineWidth(ps.FrameLineWidth);
    if (ps.ticks)
        c->SetTicks(1, 1);
    c->SetLogy(ps.logy);
}

static void ApplyHistStyle(TH1 *h, const PlotStyle &ps,
                           const std::string &xtitle,
                           const std::string &ytitle)
{
    h->SetTitle("");

    h->GetXaxis()->SetTitle(xtitle.c_str());
    h->GetYaxis()->SetTitle(ytitle.c_str());

    h->GetXaxis()->SetTitleFont(ps.font);
    h->GetYaxis()->SetTitleFont(ps.font);
    h->GetXaxis()->SetLabelFont(ps.font);
    h->GetYaxis()->SetLabelFont(ps.font);

    h->GetXaxis()->SetTitleSize(ps.xTitleSize);
    h->GetYaxis()->SetTitleSize(ps.yTitleSize);
    h->GetXaxis()->SetLabelSize(ps.xLabelSize);
    h->GetYaxis()->SetLabelSize(ps.yLabelSize);

    h->GetXaxis()->SetTitleOffset(ps.xTitleOffset);
    h->GetYaxis()->SetTitleOffset(ps.yTitleOffset);
}

static void DrawHeader(const PlotStyle &ps,
                       const std::string &mainTitle,
                       const std::string &sub1,
                       const std::string &sub2)
{
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(ps.font);
    lat.SetTextAlign(13); // left-top

    lat.SetTextSize(ps.titleSize);
    lat.DrawLatex(ps.headerX, ps.headerY, mainTitle.c_str());

    lat.SetTextSize(ps.subSize);
    if (!sub1.empty())
        lat.DrawLatex(ps.headerX, ps.headerY - ps.headerDy, sub1.c_str());
    if (!sub2.empty())
        lat.DrawLatex(ps.headerX, ps.headerY - 2.0 * ps.headerDy, sub2.c_str());
}

static void DrawInfoBox(const PlotStyle &ps,
                        const std::vector<std::string> &lines)
{
    if (lines.empty())
        return;

    // Size the box to its content (fixed per-line spacing), centered in the
    // [boxY1,boxY2] band, so 2+ lines pack tightly instead of stretching evenly
    // across a tall fixed box (a TPaveText spaces lines as boxHeight/nLines).
    const double lineH = 0.05; // NDC line spacing (boxTextSize is 0.034)
    const double yC = 0.5 * (ps.boxY1 + ps.boxY2);
    const double yHalf = 0.5 * lineH * (double)lines.size();
    TPaveText *p = new TPaveText(ps.boxX1, yC - yHalf, ps.boxX2, yC + yHalf, "NDC");
    p->SetFillStyle(0);
    p->SetBorderSize(0);
    p->SetTextFont(ps.font);
    p->SetTextSize(ps.boxTextSize);
    p->SetTextAlign(12);
    for (const auto &s : lines)
        p->AddText(s.c_str());
    p->Draw();
}

// ---------------------------------------------------------
// 3) The ONE generic plot function (works for any TH1*)
//    - h can be something you read from file
//    - optional "tuner" lets you do plot-specific tweaks
// ---------------------------------------------------------
using PlotTuner = std::function<void(TCanvas *, TH1 *)>;
using GraphTuner = std::function<void(TCanvas *, TGraphErrors *)>;

/* Example of how to use a lambda version of type defined Plot Tuner
PlotTuner mtTuner = [](TCanvas* c, TH1* h) {
    c->SetLogy();
    h->GetXaxis()->SetRangeUser(40, 200);
};
*/

// The following main function takes genearal input of
// Saving Path
// Axis Title
// ... etc check it yourself

static void SaveNicePlot1D(TH1 *h,
                           const std::string &outPathNoExt, // e.g. "./mT/foo"
                           const std::string &xTitle,
                           const std::string &yTitle,
                           const std::string &mainTitle,
                           const std::string &subTitle1,
                           const std::string &subTitle2,
                           const std::vector<std::string> &boxLines,
                           const PlotStyle &ps = PlotStyle(),
                           PlotTuner tuner = nullptr)
{
    if (!h)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);

    c->cd();

    ApplyHistStyle(h, ps, xTitle, yTitle);

    h->Draw(ps.drawOpt.c_str());

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    if (tuner)
        tuner(c, h); // e.g. change ranges, add extra lines, etc.

    // ensure parent directory exists
    // outPathNoExt could include dirs, so you can mkdir manually outside too.

    CMS_lumi(c, 13, 10);
    c->RedrawAxis(); // redraw frame + ticks on top of the histograms/fills
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

static std::string RangeLabel(const char *var, int iy, const double *edges, int nY)
{
    const double a = edges[iy];
    const double b = edges[iy + 1];

    // Use [a, b) except last bin [a, b]
    if (iy == nY - 1)
        return Form("%s #in [%.1f, %.1f]", var, a, b);
    return Form("%s #in [%.1f, %.1f)", var, a, b);
}

static void ApplyGraphStyle(TGraphErrors *g, const PlotStyle &ps,
                            const std::string &xtitle,
                            const std::string &ytitle)
{
    g->SetTitle("");
    g->GetXaxis()->SetTitle(xtitle.c_str());
    g->GetYaxis()->SetTitle(ytitle.c_str());

    g->GetXaxis()->SetTitleFont(ps.font);
    g->GetYaxis()->SetTitleFont(ps.font);
    g->GetXaxis()->SetLabelFont(ps.font);
    g->GetYaxis()->SetLabelFont(ps.font);

    g->GetXaxis()->SetTitleSize(ps.xTitleSize);
    g->GetYaxis()->SetTitleSize(ps.yTitleSize);
    g->GetXaxis()->SetLabelSize(ps.xLabelSize);
    g->GetYaxis()->SetLabelSize(ps.yLabelSize);

    g->GetXaxis()->SetTitleOffset(ps.xTitleOffset);
    g->GetYaxis()->SetTitleOffset(ps.yTitleOffset);
}

static void SaveNiceGraph(TGraphErrors *g,
                          const std::string &outPathNoExt,
                          const std::string &xTitle,
                          const std::string &yTitle,
                          const std::string &mainTitle,
                          const std::string &subTitle1,
                          const std::string &subTitle2,
                          const std::vector<std::string> &boxLines,
                          const PlotStyle &ps = PlotStyle(),
                          GraphTuner tuner = nullptr,
                          TGraphErrors *g1 = nullptr,
                          TGraphErrors *g2 = nullptr,
                          TGraphErrors *g3 = nullptr,
                          TGraphErrors *g4 = nullptr)
{
    if (!g)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas *c = new TCanvas(Form("c_%s", g->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    ApplyGraphStyle(g, ps, xTitle, yTitle);

    g->Draw(ps.drawOptGraph.c_str());

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    if (tuner)
        tuner(c, g);

    CMS_lumi(c, 13, 10);

    if (g1)
    {
        if (tuner)
            tuner(c, g1);

        g1->SetMarkerColor(kRed);
        g1->Draw("P SAME");
    }
    if (g2)
    {
        if (tuner)
            tuner(c, g2);

        g2->SetMarkerColor(kBlue);
        g2->Draw("P SAME");
    }
    if (g3)
    {
        if (tuner)
            tuner(c, g3);

        g3->SetMarkerColor(kGreen);
        g3->Draw("P SAME");
    }

    if (g4)
    {
        if (tuner)
            tuner(c, g4);

        g4->SetMarkerColor(kPink);
        g4->Draw("P SAME");
    }
    c->RedrawAxis(); // redraw frame + ticks on top
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

static void SaveNiceGraph_ErrorBand(TGraphErrors *g,
                                    const std::string &outPathNoExt,
                                    const std::string &xTitle,
                                    const std::string &yTitle,
                                    const std::string &mainTitle,
                                    const std::string &subTitle1,
                                    const std::string &subTitle2,
                                    const std::vector<std::string> &boxLines,
                                    const PlotStyle &ps = PlotStyle(),
                                    GraphTuner tuner = nullptr,
                                    TGraphErrors *g1 = nullptr,
                                    TGraphErrors *g2 = nullptr,
                                    TGraphErrors *g3 = nullptr,
                                    TGraphErrors *g4 = nullptr)
{
    if (!g)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas *c = new TCanvas(Form("c_%s", g->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    ApplyGraphStyle(g, ps, xTitle, yTitle);

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->Draw(ps.drawOptGraph.c_str());

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    if (tuner)
        tuner(c, g);

    CMS_lumi(c, 13, 10);

    // --- all four nPDF model predictions as translucent error bands ---
    // (order matches observables.C: g1=EPPS21, g2=nCTEQ15HQ, g3=nNNPDF3.0, g4=TUJU21nlo)
    TGraphErrors *models[4] = {g1, g2, g3, g4};
    const int     mcol[4]   = {kAzure - 2, kGreen + 2, kOrange + 1, kViolet + 1};
    const char   *mname[4]  = {"EPPS21", "nCTEQ15HQ", "nNNPDF3.0", "TUJU21nlo"};
    for (int i = 0; i < 4; ++i)
    {
        TGraphErrors *gm = models[i];
        if (!gm) continue;
        gm->SetFillColorAlpha(mcol[i], 0.35); // translucent band
        gm->SetFillStyle(1001);
        gm->SetLineWidth(0); // band ONLY -- no central line, no error bars
        gm->SetMarkerSize(0);
        gm->SetMarkerStyle(0);
        gm->Draw("3");       // filled band only
    }

    g->Draw("P SAME"); // data points back on top of the bands

    // >>>>>>>>>>>>>>>> TEMPORARY: "Projection with Electrons" band >>>>>>>>>>>>>>>>
    // Overlays the data with shrunk error bars (errors/1.4 ~ errors/sqrt(2), i.e.
    // the stat reach with ~2x the data). To REMOVE later: set kShowProjection to
    // false (one line), or delete this whole block AND its legend entry below.
    const bool kShowProjection = false;
    TGraphErrors *g_stat2x = nullptr;
    if (kShowProjection)
    {
        g_stat2x = (TGraphErrors *)g->Clone(Form("%s_stat2x", g->GetName()));
        const double scale = 1.0 / 1.4; // ~1/sqrt(2)
        for (int i = 0; i < g_stat2x->GetN(); ++i)
            g_stat2x->SetPointError(i, 0, g_stat2x->GetErrorY(i) * scale);
        g_stat2x->SetLineColor(kRed + 1);
        g_stat2x->SetMarkerSize(0);
        g_stat2x->SetLineWidth(3);
        g_stat2x->Draw("E SAME");
    }
    // <<<<<<<<<<<<<<<<<<<< end temporary projection block <<<<<<<<<<<<<<<<<<<<<<<<<

    TLegend *leg = new TLegend(0.20, 0.15, 0.45, 0.40);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.033);
    leg->AddEntry(g, "Data", "p");
    if (kShowProjection && g_stat2x) // TEMPORARY -- remove with the projection block above
        leg->AddEntry(g_stat2x, "Projection with Electrons", "l");
    for (int i = 0; i < 4; ++i)
        if (models[i]) leg->AddEntry(models[i], mname[i], "f");
    leg->Draw();

    c->RedrawAxis(); // redraw frame + ticks on top of the filled error bands
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

// Overlay TWO data graphs (e.g. muon vs electron) on the same axes, with the
// shared nPDF theory error bands. g1..g4 are the (flavour-universal) model
// graphs in the same order as SaveNiceGraph_ErrorBand; pass them null for a
// data-only overlay (e.g. the charge asymmetry, which has no theory bands here).
static void SaveNiceGraph_ErrorBand_TwoData(TGraphErrors *gA, const std::string &labelA,
                                            TGraphErrors *gB, const std::string &labelB,
                                            const std::string &outPathNoExt,
                                            const std::string &xTitle,
                                            const std::string &yTitle,
                                            const std::string &mainTitle,
                                            const std::string &subTitle1,
                                            const std::string &subTitle2,
                                            const std::vector<std::string> &boxLines,
                                            const PlotStyle &ps = PlotStyle(),
                                            GraphTuner tuner = nullptr,
                                            TGraphErrors *g1 = nullptr,
                                            TGraphErrors *g2 = nullptr,
                                            TGraphErrors *g3 = nullptr,
                                            TGraphErrors *g4 = nullptr)
{
    if (!gA && !gB)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);
    TGraphErrors *base = gA ? gA : gB; // defines the frame
    TCanvas *c = new TCanvas(Form("c_ovl_%s", base->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    ApplyGraphStyle(base, ps, xTitle, yTitle);
    // data A (e.g. muon): black circles ; data B (e.g. electron): red squares
    if (gA) { gA->SetMarkerStyle(20); gA->SetMarkerSize(1.2); gA->SetMarkerColor(kBlack); gA->SetLineColor(kBlack); }
    if (gB) { gB->SetMarkerStyle(21); gB->SetMarkerSize(1.2); gB->SetMarkerColor(kRed + 1); gB->SetLineColor(kRed + 1); }

    base->Draw(ps.drawOptGraph.c_str()); // axes + first data graph
    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);
    if (tuner)
        tuner(c, base);
    CMS_lumi(c, 13, 10);

    // shared theory bands (same order/colours/names as SaveNiceGraph_ErrorBand)
    TGraphErrors *models[4] = {g1, g2, g3, g4};
    const int     mcol[4]   = {kAzure - 2, kGreen + 2, kOrange + 1, kViolet + 1};
    const char   *mname[4]  = {"EPPS21", "nCTEQ15HQ", "nNNPDF3.0", "TUJU21nlo"};
    for (int i = 0; i < 4; ++i)
    {
        TGraphErrors *gm = models[i];
        if (!gm) continue;
        gm->SetFillColorAlpha(mcol[i], 0.35);
        gm->SetFillStyle(1001);
        gm->SetLineWidth(0); // band ONLY -- no central line, no error bars
        gm->SetMarkerSize(0);
        gm->SetMarkerStyle(0);
        gm->Draw("3");
    }

    // data points on top of the bands
    if (gA) gA->Draw("P SAME");
    if (gB) gB->Draw("P SAME");

    TLegend *leg = new TLegend(0.20, 0.15, 0.47, 0.42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.032);
    if (gA) leg->AddEntry(gA, labelA.c_str(), "p");
    if (gB) leg->AddEntry(gB, labelB.c_str(), "p");
    for (int i = 0; i < 4; ++i)
        if (models[i]) leg->AddEntry(models[i], mname[i], "f");
    leg->Draw();

    c->RedrawAxis();
    c->Modified();
    c->Update();
    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());
    delete c;
}

static void SaveNicePlot1D_twoplots(TH1 *h, TH1 *h2,
                                    const std::string &outPathNoExt, // e.g. "./mT/foo"
                                    const std::string &xTitle,
                                    const std::string &yTitle,
                                    const std::string &mainTitle,
                                    const std::string &subTitle1,
                                    const std::string &subTitle2,
                                    const std::vector<std::string> &boxLines,
                                    const PlotStyle &ps = PlotStyle(), const PlotStyle &ps1 = PlotStyle(),
                                    PlotTuner tuner = nullptr)
{
    if (!h)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);

    c->cd();

    ApplyHistStyle(h, ps, xTitle, yTitle);
    ApplyHistStyle(h2, ps1, xTitle, yTitle);

    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.2);
    h->SetLineColor(kBlue);
    h->SetMarkerColor(kBlue);

    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(1.2);
    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);

    h->Draw(ps.drawOpt.c_str());
    h2->Draw((ps1.drawOpt + "SAME").c_str());

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    TLegend *leg1 = new TLegend(0.65, 0.45, 0.88, 0.65);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(h, "Muon", "l");
    leg1->AddEntry(h2, "Electron", "l");
    leg1->Draw();

    if (tuner)
        tuner(c, h); // e.g. change ranges, add extra lines, etc.

    if (tuner)
        tuner(c, h2); // e.g. change ranges, add extra lines, etc.

    // ensure parent directory exists
    // outPathNoExt could include dirs, so you can mkdir manually outside too.

    CMS_lumi(c, 13, 10);
    c->RedrawAxis(); // redraw frame + ticks on top of the histograms/fills
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

static void SaveNicePlot1D_WithBkg(
    TH1 *h,                                   // data / signal
    const std::vector<TH1 *> &bkgs,           // background histograms
    const std::vector<std::string> &bkgNames, // legend names (same size as bkgs)

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
    if (!h)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    // With the pull pad the canvas grows only slightly (pullCanvasScale), so
    // the overall shape stays near-square; the main pad absorbs the rest by
    // compressing a little (fonts/positions are pad-relative and follow).
    const int canvasH = ps.pullPad
                            ? (int)((double)ps.h * ps.pullCanvasScale + 0.5)
                            : ps.h;
    TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", ps.w, canvasH);
    ApplyCanvasStyle(c, ps);

    TPad *pTop = nullptr, *pBot = nullptr;
    if (ps.pullPad)
    {
        pTop = new TPad(Form("pTop_%s", h->GetName()), "", 0.0, ps.pullSplit, 1.0, 1.0);
        pBot = new TPad(Form("pBot_%s", h->GetName()), "", 0.0, 0.0, 1.0, ps.pullSplit);
        pTop->SetLeftMargin(ps.lm);  pTop->SetRightMargin(ps.rm);
        pTop->SetTopMargin(ps.tm);   pTop->SetBottomMargin(0.025);
        pBot->SetLeftMargin(ps.lm);  pBot->SetRightMargin(ps.rm);
        pBot->SetTopMargin(0.045);   pBot->SetBottomMargin(0.34);
        pTop->SetFrameLineWidth(ps.FrameLineWidth);
        pBot->SetFrameLineWidth(ps.FrameLineWidth);
        if (ps.ticks) { pTop->SetTicks(1, 1); pBot->SetTicks(1, 1); }
        pTop->SetLogy(ps.logy);
        c->cd();
        pTop->Draw();
        pBot->Draw();
        pTop->cd();
    }
    else
        c->cd();

    //--------------------------------------------------
    // 1️⃣ Compute normalization
    //--------------------------------------------------
    double dataInt = h->Integral("width");

    double bkgSumInt = 0.0;
    for (auto b : bkgs)
    {
        if (!b)
            continue;
        bkgSumInt += b->Integral("width");
    }

    double scale = 1.0;
    if (ps.normBkgToData && bkgSumInt > 0)
        scale = dataInt / bkgSumInt;

    //--------------------------------------------------
    // 2️⃣ Prepare stack
    //--------------------------------------------------
    THStack *hs = new THStack("hs", "");

    std::vector<int> colors = {
        kAzure - 9,
        kOrange - 3,
        kGreen + 2,
        kMagenta - 3,
        kCyan + 1,
        kRed - 7,
        kGray + 1}; // 6th/7th added so the 6 MC samples + ABCD QCD each get a distinct color

    for (int i = static_cast<int>(bkgs.size()) - 1; i >= 0; --i)
    {
        TH1 *b = bkgs[i];
        if (!b)
            continue;

        b->Scale(scale);

        const int col = colors[i % colors.size()];
        b->SetFillColorAlpha(col, 0.65); // translucent fill (matches SaveDataMCRatio look)
        b->SetLineColor(col);            // color-matched outline
        b->SetLineWidth(3);              // thick edge

        hs->Add(b);
    }

    //--------------------------------------------------
    // 4️⃣ Draw data on top
    //--------------------------------------------------
    ApplyHistStyle(h, ps, xTitle, yTitle);
    if (ps.pullPad)
    {
        // x axis moves to the pull pad
        h->GetXaxis()->SetLabelSize(0.0);
        h->GetXaxis()->SetTitleSize(0.0);
    }

    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.2);
    h->SetLineColor(kBlack);

    h->Draw();

    //--------------------------------------------------
    // 3️⃣ Draw stack first
    //--------------------------------------------------
    hs->Draw("HIST SAME");

    h->Draw("E SAME");

    //--------------------------------------------------
    // 5️⃣ Legend
    //--------------------------------------------------
    const bool legAuto = (ps.legX1 < 0);
    TLegend *leg = new TLegend(legAuto ? ps.boxX1 : ps.legX1,
                               legAuto ? ps.boxY1 - 0.2 : ps.legY1,
                               legAuto ? ps.boxX2 : ps.legX2,
                               legAuto ? ps.boxY2 - 0.2 : ps.legY2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h, "Data", "lep");

    for (size_t i = 0; i < bkgs.size(); i++)
    {
        if (bkgs[i])
            leg->AddEntry(bkgs[i], bkgNames[i].c_str(), "lf"); // line + block (like SaveDataMCRatio)
    }

    leg->Draw();

    //--------------------------------------------------
    // 6️⃣ Header + info
    //--------------------------------------------------
    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    //--------------------------------------------------
    // 7️⃣ Optional tuner
    //--------------------------------------------------
    if (tuner)
        tuner(c, h);

    TPad *mainPad = ps.pullPad ? pTop : (TPad *)c;
    CMS_lumi(mainPad, 13, 10);
    mainPad->RedrawAxis(); // redraw frame + ticks on top of the stacked histograms

    //--------------------------------------------------
    // 8️⃣ Pull sub-pad: (data - MC_total)/sigma, zero-anchored bars
    //--------------------------------------------------
    if (ps.pullPad)
    {
        // MC total with full MC-stat errors (bkgs are already scaled above;
        // Add() propagates the Sumw2 arrays).
        TH1 *mcTot = nullptr;
        for (auto b : bkgs)
        {
            if (!b)
                continue;
            if (!mcTot)
            {
                mcTot = (TH1 *)b->Clone(Form("%s_pullMCtot", h->GetName()));
                mcTot->SetDirectory(nullptr);
            }
            else
                mcTot->Add(b);
        }

        TH1 *hPull = (TH1 *)h->Clone(Form("%s_pull", h->GetName()));
        hPull->SetDirectory(nullptr);
        hPull->Reset();

        double maxAbs = 0.0;
        if (mcTot)
        {
            for (int i = 1; i <= h->GetNbinsX(); ++i)
            {
                const double d = h->GetBinContent(i);
                const double m = mcTot->GetBinContent(i), em = mcTot->GetBinError(i);
                if (d == 0.0 && m == 0.0)
                    continue; // nothing in this bin at all
                // Poisson pull under H0 "data fluctuates around the prediction":
                // sigma^2 = m (expected Poisson variance of the data) + MC-stat^2.
                // Using the OBSERVED data error instead (sqrt(d), 0 for empty
                // bins) blows up empty-data bins against a small-error MC tail.
                const double s2 = std::max(m, 0.0) + em * em;
                if (s2 <= 0.0)
                    continue; // no prediction at all -> pull undefined
                const double pull = (d - m) / std::sqrt(s2);
                hPull->SetBinContent(i, pull);
                hPull->SetBinError(i, 0.0);
                maxAbs = std::max(maxAbs, std::fabs(pull));
            }
        }

        pBot->cd();

        const double yr = (ps.pullYRange > 0.0) ? ps.pullYRange
                                                : std::max(3.0, 1.15 * maxAbs);
        hPull->SetTitle("");
        hPull->SetMinimum(-yr);
        hPull->SetMaximum(yr);
        hPull->SetFillColorAlpha(kAzure + 1, 0.7);
        hPull->SetLineColor(kAzure + 2);
        hPull->SetLineWidth(1);
        hPull->SetStats(false);

        // Pad-relative font sizes: scale by the pad-height ratio so the pull
        // pad text matches the main pad visually.
        const double sf = (1.0 - ps.pullSplit) / ps.pullSplit;
        hPull->GetXaxis()->SetTitle(xTitle.c_str());
        hPull->GetYaxis()->SetTitle("(Data#minusMC)/#sigma");
        hPull->GetXaxis()->SetTitleFont(ps.font);
        hPull->GetYaxis()->SetTitleFont(ps.font);
        hPull->GetXaxis()->SetLabelFont(ps.font);
        hPull->GetYaxis()->SetLabelFont(ps.font);
        hPull->GetXaxis()->SetTitleSize(ps.xTitleSize * sf);
        hPull->GetYaxis()->SetTitleSize(ps.yTitleSize * sf * 0.85);
        hPull->GetXaxis()->SetLabelSize(ps.xLabelSize * sf);
        hPull->GetYaxis()->SetLabelSize(ps.yLabelSize * sf);
        hPull->GetXaxis()->SetTitleOffset(1.05);
        hPull->GetYaxis()->SetTitleOffset(ps.yTitleOffset / sf);
        hPull->GetYaxis()->SetNdivisions(505);
        hPull->GetYaxis()->CenterTitle(true);

        hPull->Draw("B"); // bars anchored at zero (handles negative pulls)

        const double x1 = hPull->GetXaxis()->GetXmin();
        const double x2 = hPull->GetXaxis()->GetXmax();
        TLine *l0 = new TLine(x1, 0.0, x2, 0.0);
        l0->SetLineStyle(2);
        l0->SetLineColor(kRed + 1);
        l0->Draw("SAME");
        for (const double g : {-2.0, 2.0}) // +-2 sigma guides
        {
            if (g <= -yr || g >= yr)
                continue;
            TLine *lg = new TLine(x1, g, x2, g);
            lg->SetLineStyle(3);
            lg->SetLineColor(kGray + 2);
            lg->Draw("SAME");
        }
        pBot->RedrawAxis();
        c->cd();
    }

    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

// ---------------------------------------------------------
// 7b) Data vs a SINGLE signal-MC template, overlaid, with the signal
//     normalized to the DATA PEAK (max-bin content) rather than the integral.
//       - hSignal drawn as a filled/line histogram, hData as points on top.
//       - Use for the inclusive m_T shape check where we only want to compare
//         the signal shape against data (no backgrounds, no stack).
//     Inputs are cloned, so the caller's histograms are untouched.
// ---------------------------------------------------------
static void SaveNicePlot1D_DataSignalPeak(
    TH1 *hData, TH1 *hSignal,
    const std::string &outPathNoExt,
    const std::string &xTitle,
    const std::string &yTitle,
    const std::string &mainTitle,
    const std::string &subTitle1,
    const std::string &subTitle2,
    const std::vector<std::string> &boxLines,
    const PlotStyle &ps = PlotStyle(),
    PlotTuner tuner = nullptr,
    const std::string &dataLabel = "Data",
    const std::string &signalLabel = "W+/W- simulation")
{
    if (!hData || !hSignal)
        return;

    gStyle->SetOptStat(ps.showStats ? 1110 : 0);

    TCanvas *c = new TCanvas(Form("c_%s_dsp", hData->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);
    c->cd();

    // Clone so we never mutate the file-owned / caller-owned histograms.
    TH1 *hd = (TH1 *)hData->Clone(Form("%s_dspD", hData->GetName()));
    TH1 *hs = (TH1 *)hSignal->Clone(Form("%s_dspS", hSignal->GetName()));
    hd->SetDirectory(nullptr);
    hs->SetDirectory(nullptr);

    // Normalize signal MC so its maximum bin matches the data maximum bin
    // (peak normalization, NOT integral normalization). Computed before any
    // SetMaximum from the tuner so GetMaximum() returns true bin contents.
    const double dPeak = hd->GetMaximum();
    const double sPeak = hs->GetMaximum();
    if (sPeak > 0.0 && dPeak > 0.0)
        hs->Scale(dPeak / sPeak);

    // signal: filled histogram (defines the frame); data: points on top.
    ApplyHistStyle(hs, ps, xTitle, yTitle);
    hs->SetLineColor(kRed + 1);
    hs->SetLineWidth(2);
    hs->SetFillColorAlpha(kRed - 9, 0.5);
    hs->SetMarkerSize(0);

    hd->SetMarkerStyle(20);
    hd->SetMarkerSize(1.2);
    hd->SetLineColor(kBlack);
    hd->SetMarkerColor(kBlack);

    hs->Draw("HIST");
    hd->Draw("E SAME");

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    TLegend *leg = new TLegend(ps.boxX1, ps.boxY1 - 0.2, ps.boxX2, ps.boxY2 - 0.2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(ps.font);
    leg->AddEntry(hd, dataLabel.c_str(), "lep");
    leg->AddEntry(hs, signalLabel.c_str(), "lf");
    leg->Draw();

    // Tuner acts on the frame-defining histogram (the signal).
    if (tuner)
        tuner(c, hs);

    CMS_lumi(c, 13, 10);
    c->RedrawAxis(); // redraw frame + ticks on top of the histograms/fills
    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}

// ---------------------------------------------------------
// 8) Data vs MC overlay with a ratio pad (data / MC).
//    Added for the correction/ Data-vs-MC kinematic (boson-pT) check.
//      - hData drawn as points, hMC as a filled/line histogram.
//      - normToData == true (default): MC is scaled to the data integral so the
//        comparison is shape-only (right for a pT shape check; backgrounds are
//        negligible in the Z peak). Pass false to compare absolute yields.
//      - Bottom pad shows hData/hMC via TH1::Divide (binomial-free error prop).
//    Inputs are cloned, so the caller's histograms are untouched.
// ---------------------------------------------------------
static void SaveDataMCRatio(TH1 *hData, TH1 *hMC,
                            const std::string &outPathNoExt,
                            const std::string &xTitle,
                            const std::string &yTitle,
                            const std::string &mainTitle,
                            const std::string &subTitle1,
                            const std::string &subTitle2,
                            const PlotStyle &ps = PlotStyle(),
                            bool normToData = true,
                            const std::string &dataLabel = "Data",
                            const std::string &mcLabel = "Signal MC")
{
    if (!hData || !hMC)
        return;

    gStyle->SetOptStat(0);

    TH1 *hd = (TH1 *)hData->Clone(Form("%s_dClone", hData->GetName()));
    TH1 *hm = (TH1 *)hMC->Clone(Form("%s_mClone", hMC->GetName()));
    hd->SetDirectory(nullptr);
    hm->SetDirectory(nullptr);

    if (normToData)
    {
        const double iD = hd->Integral();
        const double iM = hm->Integral();
        if (iM > 0.0 && iD > 0.0)
            hm->Scale(iD / iM);
    }

    TCanvas *c = new TCanvas(Form("c_ratio_%s", hData->GetName()), "", ps.w, ps.h);

    // Two pads: top (main) ~70%, bottom (ratio) ~30%.
    TPad *pTop = new TPad("pTop", "", 0.0, 0.30, 1.0, 1.0);
    TPad *pBot = new TPad("pBot", "", 0.0, 0.00, 1.0, 0.30);
    pTop->SetTopMargin(ps.tm);
    pTop->SetBottomMargin(0.02);
    pTop->SetLeftMargin(ps.lm);
    pTop->SetRightMargin(ps.rm);
    pBot->SetTopMargin(0.04);
    pBot->SetBottomMargin(0.35);
    pBot->SetLeftMargin(ps.lm);
    pBot->SetRightMargin(ps.rm);
    if (ps.ticks) { pTop->SetTicks(1, 1); pBot->SetTicks(1, 1); }
    pTop->SetLogy(ps.logy);
    pTop->Draw();
    pBot->Draw();

    // ---------------- top pad: overlay ----------------
    pTop->cd();
    ApplyHistStyle(hm, ps, "", yTitle); // MC defines the frame
    hm->GetXaxis()->SetLabelSize(0.0);  // hide x labels on the top pad
    hm->GetXaxis()->SetTitleSize(0.0);

    hm->SetLineColor(kRed + 1);
    hm->SetLineWidth(2);
    hm->SetFillColorAlpha(kRed - 9, 0.5);
    hm->SetMarkerSize(0);

    hd->SetMarkerStyle(20);
    hd->SetMarkerSize(0.9);
    hd->SetLineColor(kBlack);
    hd->SetMarkerColor(kBlack);

    const double ymax = std::max(hd->GetMaximum(), hm->GetMaximum());
    hm->SetMaximum(1.4 * ymax);
    if (!ps.logy)
        hm->SetMinimum(0.0);

    hm->Draw("HIST");
    hd->Draw("E1 SAME");

    TLegend *leg = new TLegend(0.62, 0.70, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(ps.font);
    leg->SetTextSize(ps.boxTextSize);
    leg->AddEntry(hd, dataLabel.c_str(), "lep");
    leg->AddEntry(hm, mcLabel.c_str(), "lf");
    leg->Draw();

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    CMS_lumi(pTop, 13, 10);
    pTop->RedrawAxis(); // redraw frame + ticks on top of the filled MC histogram

    // ---------------- bottom pad: ratio ----------------
    pBot->cd();
    TH1 *hr = (TH1 *)hd->Clone(Form("%s_ratio", hData->GetName()));
    hr->SetDirectory(nullptr);
    hr->Divide(hm); // data / MC with error propagation
    hr->SetTitle("");
    hr->SetMarkerStyle(20);
    hr->SetMarkerSize(0.9);
    hr->SetLineColor(kBlack);
    hr->SetMarkerColor(kBlack);

    hr->GetYaxis()->SetTitle("Data / MC");
    hr->GetXaxis()->SetTitle(xTitle.c_str());

    // The bottom pad is ~0.3 of the height, so scale fonts up to match the top.
    const double sf = 0.7 / 0.3;
    hr->GetXaxis()->SetTitleFont(ps.font); hr->GetYaxis()->SetTitleFont(ps.font);
    hr->GetXaxis()->SetLabelFont(ps.font); hr->GetYaxis()->SetLabelFont(ps.font);
    hr->GetXaxis()->SetTitleSize(ps.xTitleSize * sf);
    hr->GetYaxis()->SetTitleSize(ps.yTitleSize * sf);
    hr->GetXaxis()->SetLabelSize(ps.xLabelSize * sf);
    hr->GetYaxis()->SetLabelSize(ps.yLabelSize * sf);
    hr->GetXaxis()->SetTitleOffset(1.0);
    hr->GetYaxis()->SetTitleOffset(ps.yTitleOffset / sf);
    hr->GetYaxis()->SetNdivisions(505);

    hr->SetMinimum(0.5);
    hr->SetMaximum(1.5);
    hr->Draw("E1");

    TLine *l1 = new TLine(hr->GetXaxis()->GetXmin(), 1.0,
                          hr->GetXaxis()->GetXmax(), 1.0);
    l1->SetLineStyle(2);
    l1->SetLineColor(kRed + 1);
    l1->Draw("SAME");
    pBot->RedrawAxis(); // redraw frame + ticks on top of the ratio points

    c->cd();
    c->Modified();
    c->Update();
    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}
