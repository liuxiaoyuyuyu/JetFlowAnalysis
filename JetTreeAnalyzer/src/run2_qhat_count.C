#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TBranch.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

namespace {

const double kQHatBoundaries[] = {
    470.0, 600.0, 800.0, 1000.0, 1400.0,
    1800.0, 2400.0, 3200.0
};
const int kQHatBins = 7;
const char* kTreePath = "analyzer/trackTree";

bool supportedEra(const std::string& era) {
    return era == "2016APV" || era == "2016" || era == "2017" || era == "2018";
}

}  // namespace

int main(int argc, const char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: run2_qhat_count <file-list> <production-era> <output.root>\n";
        return 2;
    }

    const std::string fileListPath = argv[1];
    const std::string productionEra = argv[2];
    const std::string outputPath = argv[3];
    if (!supportedEra(productionEra)) {
        std::cerr << "Unsupported production era '" << productionEra
                  << "'. Choose 2016APV, 2016, 2017, or 2018.\n";
        return 1;
    }

    std::ifstream fileList(fileListPath.c_str());
    if (!fileList) {
        std::cerr << "Cannot open input file list: " << fileListPath << '\n';
        return 1;
    }

    std::vector<std::string> inputPaths;
    std::string inputPath;
    while (fileList >> inputPath) {
        inputPaths.push_back(inputPath);
    }
    if (inputPaths.empty()) {
        std::cerr << "Input file list is empty: " << fileListPath << '\n';
        return 1;
    }

    TH1D counts("qhatCounts", ";#hat{p}_{T} (GeV);generated events",
                kQHatBins, kQHatBoundaries);
    Long64_t totalEntries = 0;
    for (std::size_t fileIndex = 0; fileIndex < inputPaths.size(); ++fileIndex) {
        TFile* inputFile = TFile::Open(inputPaths[fileIndex].c_str(), "read");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Cannot open input file: " << inputPaths[fileIndex] << '\n';
            delete inputFile;
            return 1;
        }
        TTree* tree = dynamic_cast<TTree*>(inputFile->Get(kTreePath));
        if (!tree) {
            std::cerr << "Missing tree '" << kTreePath << "' in "
                      << inputPaths[fileIndex] << '\n';
            inputFile->Close();
            delete inputFile;
            return 1;
        }
        TBranch* qhatBranch = tree->GetBranch("genQScale");
        if (!qhatBranch) {
            std::cerr << "Missing genQScale branch in "
                      << inputPaths[fileIndex] << '\n';
            inputFile->Close();
            delete inputFile;
            return 1;
        }

        Float_t genQScale = 0.0f;
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("genQScale", 1);
        if (tree->SetBranchAddress("genQScale", &genQScale) < 0) {
            std::cerr << "Cannot bind genQScale branch in "
                      << inputPaths[fileIndex] << '\n';
            inputFile->Close();
            delete inputFile;
            return 1;
        }
        const Long64_t entries = tree->GetEntries();
        if (entries <= 0) {
            std::cerr << "No entries found at tree path '" << kTreePath
                      << "' in " << inputPaths[fileIndex] << '\n';
            inputFile->Close();
            delete inputFile;
            return 1;
        }
        for (Long64_t entry = 0; entry < entries; ++entry) {
            if (tree->GetEntry(entry) <= 0) {
                std::cerr << "Failed to read entry " << entry << " in "
                          << inputPaths[fileIndex] << '\n';
                inputFile->Close();
                delete inputFile;
                return 1;
            }
            counts.Fill(genQScale);
        }
        totalEntries += entries;
        inputFile->Close();
        delete inputFile;
    }

    if (counts.GetBinContent(0) != 0.0 ||
        counts.GetBinContent(kQHatBins + 1) != 0.0) {
        std::cerr << "Found q-hat values outside the configured range: underflow="
                  << counts.GetBinContent(0) << ", overflow="
                  << counts.GetBinContent(kQHatBins + 1) << ".\n";
        return 1;
    }

    TFile output(outputPath.c_str(), "recreate");
    if (output.IsZombie()) {
        std::cerr << "Cannot create output file: " << outputPath << '\n';
        return 1;
    }
    counts.Write();
    output.Close();

    std::cout << "Counted " << totalEntries << " events from " << inputPaths.size()
              << " files into " << outputPath << '\n';
    return 0;
}
