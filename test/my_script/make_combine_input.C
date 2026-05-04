#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TClass.h"
#include "TKey.h"
#include "TF1.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

// ------------------------------------------------------------
// Your Rayleigh-like function
// p[0] = overall normalization (not important for shape template)
// p[1], p[2], p[3] = shape parameters
// ------------------------------------------------------------
double QCDRayleighLike(double *x, double *p)
{
    const double xx = x[0];
    if (xx <= 0)
        return 0.0;

    const double xt = xx / 50.0 - 1.0;
    const double sigma = p[1] + p[2] * xt + p[3] * (2.0 * xt * xt - 1.0);

    if (sigma <= 0)
        return 0.0;

    return p[0] * xx * std::exp(-(xx * xx) / (2.0 * sigma * sigma));
}

// ------------------------------------------------------------
// Helper: copy one histogram from input file into output file
// ------------------------------------------------------------
bool CopyHist(TFile *fin, TFile *fout, const char *srcName, const char *dstName)
{
    if (!fin || fin->IsZombie())
    {
        std::cerr << "[ERROR] Bad input file while copying " << srcName << "\n";
        return false;
    }

    TH1 *h = dynamic_cast<TH1 *>(fin->Get(srcName));
    if (!h)
    {
        std::cerr << "[ERROR] Histogram not found: " << srcName
                  << " in file " << fin->GetName() << "\n";
        return false;
    }

    fout->cd();
    TH1 *hc = dynamic_cast<TH1 *>(h->Clone(dstName));
    if (!hc)
    {
        std::cerr << "[ERROR] Failed to clone histogram " << srcName << "\n";
        return false;
    }

    hc->SetDirectory(fout);
    hc->Write(dstName, TObject::kOverwrite);
    return true;
}

// ------------------------------------------------------------
// Main macro
// ------------------------------------------------------------
void make_combine_input()
{
    const char *outFileName = "combine_input.root";

    // Better to avoid "~" here
    const char *file_all = "/afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots/combine_input_inclusive.root";

    const char *hist_data = "data_obs";
    const char *hist_sig  = "signal";
    const char *hist_bkg1 = "z";
    const char *hist_bkg2 = "ztau";
    const char *hist_bkg3 = "wtau";

    TFile *fin = TFile::Open(file_all, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open " << file_all << "\n";
        return;
    }

    TFile *fout = TFile::Open(outFileName, "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "[ERROR] Cannot create " << outFileName << "\n";
        return;
    }

    if (!CopyHist(fin, fout, hist_data, "data_obs")) return;
    if (!CopyHist(fin, fout, hist_sig,  "signal"))   return;
    if (!CopyHist(fin, fout, hist_bkg1, "z"))        return;
    if (!CopyHist(fin, fout, hist_bkg2, "ztau"))     return;
    if (!CopyHist(fin, fout, hist_bkg3, "wtau"))     return;

    TH1 *h_data_out = dynamic_cast<TH1 *>(fout->Get("data_obs"));
    TH1 *h_sig_out  = dynamic_cast<TH1 *>(fout->Get("signal"));
    TH1 *h_z_out    = dynamic_cast<TH1 *>(fout->Get("z"));
    TH1 *h_ztau_out = dynamic_cast<TH1 *>(fout->Get("ztau"));
    TH1 *h_wtau_out = dynamic_cast<TH1 *>(fout->Get("wtau"));

    if (!h_data_out || !h_sig_out || !h_z_out || !h_ztau_out || !h_wtau_out)
    {
        std::cerr << "[ERROR] Failed to retrieve one or more written histograms from output file.\n";
        return;
    }

    const double data_int = h_data_out->Integral("width");
    const double mc_sum =
        h_sig_out->Integral("width") +
        h_z_out->Integral("width") +
        h_ztau_out->Integral("width") +
        h_wtau_out->Integral("width");

    if (data_int <= 0.0)
    {
        std::cerr << "[ERROR] Data integral (with width) is <= 0, cannot normalize.\n";
        return;
    }
    if (mc_sum <= 0.0)
    {
        std::cerr << "[ERROR] Total MC integral (with width) is <= 0, cannot normalize.\n";
        return;
    }

    const double scale_to_data = data_int / mc_sum;

    h_sig_out->Scale(scale_to_data);
    h_z_out->Scale(scale_to_data);
    h_ztau_out->Scale(scale_to_data);
    h_wtau_out->Scale(scale_to_data);

    fout->cd();
    h_sig_out->Write("signal", TObject::kOverwrite);
    h_z_out->Write("z", TObject::kOverwrite);
    h_ztau_out->Write("ztau", TObject::kOverwrite);
    h_wtau_out->Write("wtau", TObject::kOverwrite);

    std::cout << "[INFO] Data integral (width)     = " << data_int << "\n";
    std::cout << "[INFO] Total MC integral (width) = " << mc_sum << "\n";
    std::cout << "[INFO] Applied common MC scale   = " << scale_to_data << "\n";

    // ------------------------------------------------------------
    // Create histogram template from the old PDF shape
    // ------------------------------------------------------------
    TH1 *h_ref = dynamic_cast<TH1 *>(fout->Get("data_obs"));
    if (!h_ref)
    {
        std::cerr << "[ERROR] Could not retrieve data_obs from output file.\n";
        return;
    }

    TH1D *h_pdfbkg = dynamic_cast<TH1D *>(h_ref->Clone("pdfbkg"));
    if (!h_pdfbkg)
    {
        std::cerr << "[ERROR] Failed to clone reference histogram for pdfbkg.\n";
        return;
    }
    h_pdfbkg->Reset();
    h_pdfbkg->SetDirectory(nullptr);

    TF1 fQCD("fQCD", QCDRayleighLike,
             h_ref->GetXaxis()->GetXmin(),
             h_ref->GetXaxis()->GetXmax(), 4);

    // p0 is only an overall scale for building the template
    fQCD.SetParameter(0, 1.0);
    fQCD.SetParameter(1, 17.2223);
    fQCD.SetParameter(2, 10.57);
    fQCD.SetParameter(3, 0.6);

    for (int ib = 1; ib <= h_pdfbkg->GetNbinsX(); ++ib)
    {
        const double xlo = h_pdfbkg->GetXaxis()->GetBinLowEdge(ib);
        const double xhi = h_pdfbkg->GetXaxis()->GetBinUpEdge(ib);

        // Fill each bin by integrating the function over the bin
        const double val = fQCD.Integral(xlo, xhi);
        h_pdfbkg->SetBinContent(ib, val);
        h_pdfbkg->SetBinError(ib, 0.0);
    }

    // Optional: normalize this template to area 1 so its amplitude is purely from rateParam
    const double pdf_int = h_pdfbkg->Integral("width");
    if (pdf_int > 0.0)
    {
        h_pdfbkg->Scale(1.0 / pdf_int);
    }
    else
    {
        std::cerr << "[WARNING] pdfbkg histogram has non-positive integral.\n";
    }

    fout->cd();
    h_pdfbkg->Write("pdfbkg", TObject::kOverwrite);

    fout->Close();
    fin->Close();

    std::cout << "[INFO] Wrote output file: " << outFileName << "\n";
    std::cout << "[INFO] It contains TH1 templates:\n";
    std::cout << "       - data_obs\n";
    std::cout << "       - signal\n";
    std::cout << "       - z\n";
    std::cout << "       - ztau\n";
    std::cout << "       - wtau\n";
    std::cout << "       - pdfbkg\n";
}