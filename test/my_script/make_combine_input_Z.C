#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

#include <iostream>

// ------------------------------------------------------------
// Helper: copy one histogram from input file into output file
// ------------------------------------------------------------
bool CopyHist(TFile *fin, TFile *fout, const char *srcName, const char *dstName)
{
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Bad input file while copying " << srcName << "\n";
        return false;
    }
    TH1 *h = dynamic_cast<TH1 *>(fin->Get(srcName));
    if (!h) {
        std::cerr << "[ERROR] Histogram not found: " << srcName
                  << " in " << fin->GetName() << "\n";
        return false;
    }
    fout->cd();
    TH1 *hc = dynamic_cast<TH1 *>(h->Clone(dstName));
    if (!hc) return false;
    hc->SetDirectory(fout);
    hc->Write(dstName, TObject::kOverwrite);
    return true;
}

// ------------------------------------------------------------
// Main macro
//   isElec = 0 -> Z -> mu mu
//   isElec = 1 -> Z -> e e
// ------------------------------------------------------------
void make_combine_input_Z(bool isElec = 0)
{
    const char *channel = isElec ? "ee" : "mumu";

    const char *file_all = isElec
        ? "/afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots/Elec/combine_input_dilepton.root"
        : "/afs/cern.ch/user/z/zheng/pO_analysis/plotting/plots/combine_input_dilepton.root";

    const TString outFileName = TString::Format("combine_input_Z%s.root", channel);

    TFile *fin = TFile::Open(file_all, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open " << file_all << "\n";
        return;
    }

    TFile *fout = TFile::Open(outFileName, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create " << outFileName << "\n";
        return;
    }

    // Copy the templates (names already match what dileptonpeak.C writes)
    if (!CopyHist(fin, fout, "data_obs", "data_obs")) return;
    if (!CopyHist(fin, fout, "signal",   "signal"))   return;   // DY -> ll
    if (!CopyHist(fin, fout, "w",        "w"))        return;   // W+/W-
    if (!CopyHist(fin, fout, "wtau",     "wtau"))     return;   // W+/W- tau
    if (!CopyHist(fin, fout, "ztau",     "ztau"))     return;   // DY tau

    TH1 *h_data = dynamic_cast<TH1 *>(fout->Get("data_obs"));
    TH1 *h_sig  = dynamic_cast<TH1 *>(fout->Get("signal"));
    TH1 *h_w    = dynamic_cast<TH1 *>(fout->Get("w"));
    TH1 *h_wtau = dynamic_cast<TH1 *>(fout->Get("wtau"));
    TH1 *h_ztau = dynamic_cast<TH1 *>(fout->Get("ztau"));
    if (!h_data || !h_sig || !h_w || !h_wtau || !h_ztau) {
        std::cerr << "[ERROR] Failed to retrieve one or more written histograms.\n";
        return;
    }

    // Common MC -> data normalization (mirrors the W macro)
    const double data_int = h_data->Integral("width");
    const double mc_sum =
        h_sig ->Integral("width") +
        h_w   ->Integral("width") +
        h_wtau->Integral("width") +
        h_ztau->Integral("width");

    if (data_int <= 0 || mc_sum <= 0) {
        std::cerr << "[ERROR] Data or MC integral non-positive (data="
                  << data_int << ", mc=" << mc_sum << ")\n";
        return;
    }

    const double scale = data_int / mc_sum;
    h_sig ->Scale(scale);
    h_w   ->Scale(scale);
    h_wtau->Scale(scale);
    h_ztau->Scale(scale);

    fout->cd();
    h_sig ->Write("signal", TObject::kOverwrite);
    h_w   ->Write("w",      TObject::kOverwrite);
    h_wtau->Write("wtau",   TObject::kOverwrite);
    h_ztau->Write("ztau",   TObject::kOverwrite);

    std::cout << "[INFO] Channel                   : Z -> " << (isElec ? "e+e-" : "mu+mu-") << "\n";
    std::cout << "[INFO] Data integral (width)     = " << data_int << "\n";
    std::cout << "[INFO] Total MC integral (width) = " << mc_sum   << "\n";
    std::cout << "[INFO] Applied common MC scale   = " << scale    << "\n";

    fout->Close();
    fin ->Close();

    std::cout << "[INFO] Wrote " << outFileName << "\n"
              << "       Templates: data_obs, signal, w, wtau, ztau\n";
}