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

    g->Draw(ps.drawOptGraph.c_str());

    // --- Clone graph to show 2x statistics (errors / sqrt(2)) ---
    TGraphErrors *g_stat2x = (TGraphErrors *)g->Clone(Form("%s_stat2x", g->GetName()));

    const double scale = 1.0 / 1.4; // ~1/sqrt(2)

    for (int i = 0; i < g_stat2x->GetN(); ++i)
    {
        double ex = g_stat2x->GetErrorX(i);
        double ey = g_stat2x->GetErrorY(i);

        g_stat2x->SetPointError(i, 0, ey * scale);
    }

    // style: only draw error bars, different color
    g_stat2x->SetLineColor(kRed + 1);
    g_stat2x->SetMarkerSize(0); // hide markers
    g_stat2x->SetLineWidth(3);

    // draw only the error bars
    g_stat2x->Draw("E SAME");

    DrawHeader(ps, mainTitle, subTitle1, subTitle2);
    DrawInfoBox(ps, boxLines);

    if (tuner)
        tuner(c, g);

    CMS_lumi(c, 13, 10);

    if (g1)
    {

        g1->SetFillColorAlpha(kBlue, 0.35); // translucent band
        g1->SetFillStyle(1001);             // solid fill
        g1->SetMarkerSize(0);               // hide points
        g1->SetMarkerStyle(0);              // hide markers
        g1->SetLineWidth(0);                // hide error bar stems

        g1->Draw("3");      // A = draw axes, 3 = filled error band
        g1->Draw("L SAME"); // draw central line on top
    }
    if (g2)
    {

        g2->SetFillColorAlpha(kGreen, 0.35); // translucent band
        g2->SetFillStyle(1001);              // solid fill
        g2->SetMarkerSize(0);                // hide points
        g2->SetMarkerStyle(0);               // hide markers
        g2->SetLineWidth(0);                 // hide error bar stems

        // nCTEQ15HQ: kept available but not drawn (see legend below)
        // g2->Draw("3");      // A = draw axes, 3 = filled error band
        // g2->Draw("L SAME"); // draw central line on top
    }
    if (g3)
    {

        g3->SetFillColorAlpha(kOrange, 0.35); // translucent band
        g3->SetFillStyle(1001);               // solid fill
        g3->SetMarkerSize(0);                 // hide points
        g3->SetMarkerStyle(0);                // hide markers
        g3->SetLineWidth(0);                  // hide error bar stems

        g3->Draw("3");      // A = draw axes, 3 = filled error band
        g3->Draw("L SAME"); // draw central line on top
    }

    if (g4)
    {

        g4->SetFillColorAlpha(kRed, 0.35); // translucent band
        g4->SetFillStyle(1001);            // solid fill
        g4->SetMarkerSize(0);              // hide points
        g4->SetMarkerStyle(0);             // hide markers
        g4->SetLineWidth(0);               // hide error bar stems

        // TUJU21nlo: kept available but not drawn (see legend below)
        // g4->Draw("3");      // A = draw axes, 3 = filled error band
        // g4->Draw("L SAME"); // draw central line on top
    }

    TLegend *leg = new TLegend(0.20, 0.15, 0.45, 0.38);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    // --- Data: points ---
    leg->AddEntry(g, "Data", "p");
    leg->AddEntry(g_stat2x, "Projection with Electrons", "l");

    // --- Models: error bands ---
    leg->AddEntry(g1, "EPPS21", "f");
    // leg->AddEntry(g2, "nCTEQ15HQ", "f");
    leg->AddEntry(g3, "nNNPDF3.0", "f");
    // leg->AddEntry(g4, "TUJU21nlo", "f");

    leg->Draw();

    c->RedrawAxis(); // redraw frame + ticks on top of the filled error bands
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

    TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", ps.w, ps.h);
    ApplyCanvasStyle(c, ps);

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
    TLegend *leg = new TLegend(ps.boxX1, ps.boxY1 - 0.2, ps.boxX2, ps.boxY2 - 0.2);
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

    CMS_lumi(c, 13, 10);

    c->RedrawAxis(); // redraw frame + ticks on top of the stacked histograms
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
