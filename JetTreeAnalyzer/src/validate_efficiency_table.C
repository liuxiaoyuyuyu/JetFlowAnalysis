#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

namespace {

const double kCrossSectionsPb[] = {
    552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236
};
const int kQHatBins = 7;

bool closeEnough(double actual, double expected) {
    return std::fabs(actual - expected) <=
           1.0e-8 * std::max(1.0, std::fabs(expected));
}

}  // namespace

int main(int argc, const char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: validate_efficiency_table <global-qhat-counts.root> "
                     "<merged-correction.root>\n";
        return 2;
    }

    TFile countsFile(argv[1], "read");
    TFile correctionFile(argv[2], "read");
    if (countsFile.IsZombie() || correctionFile.IsZombie()) {
        std::cerr << "Cannot open one or both input ROOT files.\n";
        return 1;
    }

    TH1D* counts = dynamic_cast<TH1D*>(countsFile.Get("qhatCounts"));
    TH1D* weighted = dynamic_cast<TH1D*>(correctionFile.Get("qhatWeighted"));
    TH2D* reco = dynamic_cast<TH2D*>(
        correctionFile.Get("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All"));
    TH2D* gen = dynamic_cast<TH2D*>(
        correctionFile.Get("h2_Dau_Gen_Pt_Eta_Lab_All"));
    if (!counts || !weighted || !reco || !gen) {
        std::cerr << "A required validation histogram is missing.\n";
        return 1;
    }
    if (counts->GetNbinsX() != kQHatBins || weighted->GetNbinsX() != kQHatBins) {
        std::cerr << "Unexpected q-hat histogram binning.\n";
        return 1;
    }

    bool valid = true;
    for (int bin = 1; bin <= kQHatBins; ++bin) {
        if (counts->GetBinContent(bin) <= 0.0) {
            std::cerr << "q-hat bin " << bin << " has no generated events.\n";
            valid = false;
        }
        if (!closeEnough(weighted->GetBinContent(bin), kCrossSectionsPb[bin - 1])) {
            std::cerr << "q-hat bin " << bin << " integrates to "
                      << weighted->GetBinContent(bin) << " pb, expected "
                      << kCrossSectionsPb[bin - 1] << " pb.\n";
            valid = false;
        }
    }
    if (reco->GetNbinsX() != gen->GetNbinsX() ||
        reco->GetNbinsY() != gen->GetNbinsY()) {
        std::cerr << "Reco and gen pT-eta histograms have different binning.\n";
        return 1;
    }

    std::unique_ptr<TH2D> efficiency(
        dynamic_cast<TH2D*>(reco->Clone("validation_efficiency")));
    efficiency->SetDirectory(0);
    efficiency->Divide(gen);

    double minimum = 1.0e100;
    double maximum = -1.0e100;
    int populatedBins = 0;
    int emptyDenominatorBins = 0;
    int suspiciousBins = 0;
    for (int xbin = 1; xbin <= gen->GetNbinsX(); ++xbin) {
        const double pt = gen->GetXaxis()->GetBinCenter(xbin);
        if (pt < 0.3 || pt > 100.0) continue;
        for (int ybin = 1; ybin <= gen->GetNbinsY(); ++ybin) {
            const double eta = gen->GetYaxis()->GetBinCenter(ybin);
            if (std::fabs(eta) > 2.4) continue;
            const double denominator = gen->GetBinContent(xbin, ybin);
            if (denominator <= 0.0) {
                ++emptyDenominatorBins;
                continue;
            }
            const double value = efficiency->GetBinContent(xbin, ybin);
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            ++populatedBins;
            if (!std::isfinite(value) || value <= 0.0 || value > 1.2) {
                ++suspiciousBins;
            }
        }
    }

    std::cout << "Efficiency validation region: " << populatedBins
              << " populated bins, " << emptyDenominatorBins
              << " empty denominator bins, range [" << minimum << ", "
              << maximum << "], suspicious bins=" << suspiciousBins << '\n';
    if (populatedBins == 0 || suspiciousBins != 0) valid = false;
    if (!valid) {
        std::cerr << "Validation failed; do not use this correction table.\n";
        return 1;
    }
    std::cout << "Validation passed.\n";
    return 0;
}
