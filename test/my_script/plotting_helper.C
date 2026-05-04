#include "CMS_lumi.C"

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
    double boxTextSize = 0.025;

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
    double FrameLineWidth = 3;

    // draw options
    std::string drawOpt = "E"; // e.g. "E", "hist", "E1", etc.
    std::string drawOptGraph = "AP";
    bool logy = false;

    // stat box
    bool showStats = false;
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

    TPaveText *p = new TPaveText(ps.boxX1, ps.boxY1, ps.boxX2, ps.boxY2, "NDC");
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
    if (bkgSumInt > 0)
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
        kCyan + 1};

    for (int i = static_cast<int>(bkgs.size()) - 1; i >= 0; --i)
    {
        TH1 *b = bkgs[i];
        if (!b)
            continue;

        b->Scale(scale);

        b->SetFillColor(colors[i % colors.size()]);
        b->SetLineColor(kBlack);
        b->SetLineWidth(1);

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
            leg->AddEntry(bkgs[i], bkgNames[i].c_str(), "f");
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

    c->Modified();
    c->Update();

    c->SaveAs((outPathNoExt + ".png").c_str());
    c->SaveAs((outPathNoExt + ".pdf").c_str());

    delete c;
}