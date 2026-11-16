#define MyClass_cxx
#include "MyTrimPythia.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TTree.h>

namespace {

const double kQHatBoundaries[] = {
    470.0, 600.0, 800.0, 1000.0, 1400.0,
    1800.0, 2400.0, 3200.0
};
const double kCrossSectionsPb[] = {
    552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236
};
const int kQHatBins = 7;
const double kJetEtaMax = 1.6;
const double kJetPtMin = 550.0;
const double kMatchDeltaRMax = 0.4;  // AK8 radius / 2, as in the revised workflow.
const double kTrackEtaMax = 2.4;
const char* kTreePath = "analyzer/trackTree";

std::string gQHatCountsPath;
std::string gOutputPath;
bool gSucceeded = false;

bool isSupportedEra(const std::string& era) {
    return era == "2016APV" || era == "2016" ||
           era == "2017" || era == "2018";
}

struct MatchCandidate {
    std::size_t reco;
    std::size_t gen;
    double deltaR;
};

std::vector<int> matchJets(const std::vector<float>& recoEta,
                           const std::vector<float>& recoPhi,
                           const std::vector<float>& genEta,
                           const std::vector<float>& genPhi) {
    std::vector<MatchCandidate> candidates;
    for (std::size_t reco = 0; reco < recoEta.size(); ++reco) {
        for (std::size_t gen = 0; gen < genEta.size(); ++gen) {
            const double deltaEta = recoEta[reco] - genEta[gen];
            const double deltaPhi = std::remainder(
                recoPhi[reco] - genPhi[gen], 2.0 * TMath::Pi());
            const double deltaR = std::hypot(deltaEta, deltaPhi);
            if (deltaR < kMatchDeltaRMax) {
                MatchCandidate candidate = {reco, gen, deltaR};
                candidates.push_back(candidate);
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const MatchCandidate& lhs, const MatchCandidate& rhs) {
                  return lhs.deltaR < rhs.deltaR;
              });

    std::vector<int> recoToGen(recoEta.size(), -1);
    std::vector<bool> genUsed(genEta.size(), false);
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        const MatchCandidate& candidate = candidates[index];
        if (recoToGen[candidate.reco] >= 0 || genUsed[candidate.gen]) continue;
        recoToGen[candidate.reco] = static_cast<int>(candidate.gen);
        genUsed[candidate.gen] = true;
    }
    return recoToGen;
}

// Common Run 2/Run 3 selection from the analysis note.  The relative-pT-error
// and DCA requirements apply to every selected reconstructed charged particle.
bool passesRecoTrackSelection(const MyClass& reader,
                              std::size_t jet, std::size_t track) {
    if ((*reader.dau_chg)[jet][track] == 0) return false;
    const double pt = (*reader.dau_pt)[jet][track];
    if (pt <= 0.0) return false;
    if (std::fabs((*reader.dau_eta)[jet][track]) > kTrackEtaMax) return false;
    if ((*reader.dau_ptError)[jet][track] / pt > 0.1) return false;
    if (std::fabs((*reader.dau_XYDCAsig)[jet][track]) > 3.0) return false;
    if (std::fabs((*reader.dau_ZDCAsig)[jet][track]) > 3.0) return false;
    return true;
}

bool passesGenTrackSelection(const MyClass& reader,
                             std::size_t jet, std::size_t track) {
    return (*reader.genDau_chg)[jet][track] != 0 &&
           std::fabs((*reader.genDau_eta)[jet][track]) <= kTrackEtaMax;
}

bool hasRequiredBranches(TTree* tree) {
    const char* required[] = {
        "genQScale", "jetEta", "jetPt", "jetPhi", "dau_chg", "dau_pt",
        "dau_ptError", "dau_XYDCAsig", "dau_ZDCAsig", "dau_eta", "dau_phi",
        "genJetEta", "genJetPt", "genJetPhi", "genDau_chg", "genDau_pt",
        "genDau_eta", "genDau_phi"
    };
    const std::size_t count = sizeof(required) / sizeof(required[0]);
    for (std::size_t index = 0; index < count; ++index) {
        if (!tree->GetBranch(required[index])) {
            std::cerr << "Missing required branch '" << required[index] << "'.\n";
            return false;
        }
    }
    return true;
}

bool collectionsAreConsistent(const MyClass& reader) {
    const std::size_t recoJets = reader.jetPt->size();
    const std::size_t genJets = reader.genJetPt->size();
    if (reader.jetEta->size() != recoJets || reader.jetPhi->size() != recoJets ||
        reader.dau_chg->size() != recoJets || reader.dau_pt->size() != recoJets ||
        reader.dau_ptError->size() != recoJets ||
        reader.dau_XYDCAsig->size() != recoJets ||
        reader.dau_ZDCAsig->size() != recoJets ||
        reader.dau_eta->size() != recoJets || reader.dau_phi->size() != recoJets ||
        reader.genJetEta->size() != genJets || reader.genJetPhi->size() != genJets ||
        reader.genDau_chg->size() != genJets || reader.genDau_pt->size() != genJets ||
        reader.genDau_eta->size() != genJets || reader.genDau_phi->size() != genJets) {
        return false;
    }
    for (std::size_t jet = 0; jet < recoJets; ++jet) {
        const std::size_t tracks = (*reader.dau_pt)[jet].size();
        if ((*reader.dau_chg)[jet].size() != tracks ||
            (*reader.dau_ptError)[jet].size() != tracks ||
            (*reader.dau_XYDCAsig)[jet].size() != tracks ||
            (*reader.dau_ZDCAsig)[jet].size() != tracks ||
            (*reader.dau_eta)[jet].size() != tracks ||
            (*reader.dau_phi)[jet].size() != tracks) return false;
    }
    for (std::size_t jet = 0; jet < genJets; ++jet) {
        const std::size_t tracks = (*reader.genDau_pt)[jet].size();
        if ((*reader.genDau_chg)[jet].size() != tracks ||
            (*reader.genDau_eta)[jet].size() != tracks ||
            (*reader.genDau_phi)[jet].size() != tracks) return false;
    }
    return true;
}

}  // namespace

void MyClass::Loop(int job, std::string fileListPath) {
    (void)job;
    (void)fileListPath;
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const double ptBins[] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
        1.8, 1.9, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0,
        15.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0
    };
    const int ptBinCount = sizeof(ptBins) / sizeof(ptBins[0]) - 1;
    // Keep the expanded table axes used in the revised workflow.  The selected
    // tracks occupy |eta| <= 2.4, but the wider axis makes under/overflow easy
    // to diagnose and preserves compatibility with future selections.
    const int etaBins = 80;
    const double etaMax = 4.0;
    const int phiBins = 80;

    TH1D recoPt("h_MAIN_DCA_Dau_Reco_Pt_Lab_All",
                ";reco track p_{T} (GeV);weighted tracks", ptBinCount, ptBins);
    TH1D recoEta("h_MAIN_DCA_Dau_Reco_Eta_Lab_All",
                 ";reco track #eta;weighted tracks", etaBins, -etaMax, etaMax);
    TH1D recoPhi("h_MAIN_DCA_Dau_Reco_Phi_Lab_All",
                 ";reco track #phi;weighted tracks", phiBins,
                 -TMath::Pi(), TMath::Pi());
    TH2D recoPtEta("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All",
                   ";reco track p_{T} (GeV);reco track #eta", ptBinCount, ptBins,
                   etaBins, -etaMax, etaMax);
    TH2D recoPhiEta("h2_MAIN_DCA_Dau_Reco_Phi_Eta_Lab_All",
                    ";reco track #phi;reco track #eta", phiBins,
                    -TMath::Pi(), TMath::Pi(), etaBins, -etaMax, etaMax);
    TH2D recoPhiPt("h2_MAIN_DCA_Dau_Reco_Phi_Pt_Lab_All",
                   ";reco track #phi;reco track p_{T} (GeV)", phiBins,
                   -TMath::Pi(), TMath::Pi(), ptBinCount, ptBins);

    TH1D genPt("h_Dau_Gen_Pt_Lab_All",
               ";gen charged-particle p_{T} (GeV);weighted particles",
               ptBinCount, ptBins);
    TH1D genEta("h_Dau_Gen_Eta_Lab_All",
                ";gen charged-particle #eta;weighted particles",
                etaBins, -etaMax, etaMax);
    TH1D genPhi("h_Dau_Gen_Phi_Lab_All",
                ";gen charged-particle #phi;weighted particles", phiBins,
                -TMath::Pi(), TMath::Pi());
    TH2D genPtEta("h2_Dau_Gen_Pt_Eta_Lab_All",
                  ";gen charged-particle p_{T} (GeV);gen charged-particle #eta",
                  ptBinCount, ptBins, etaBins, -etaMax, etaMax);
    TH2D genPhiEta("h2_Dau_Gen_Phi_Eta_Lab_All",
                   ";gen charged-particle #phi;gen charged-particle #eta",
                   phiBins, -TMath::Pi(), TMath::Pi(), etaBins, -etaMax, etaMax);
    TH2D genPhiPt("h2_Dau_Gen_Phi_Pt_Lab_All",
                  ";gen charged-particle #phi;gen charged-particle p_{T} (GeV)",
                  phiBins, -TMath::Pi(), TMath::Pi(), ptBinCount, ptBins);

    TH1D qhatWeighted("qhatWeighted", ";#hat{p}_{T} (GeV);cross section (pb)",
                      kQHatBins, kQHatBoundaries);
    TH1D matchedDeltaR("matchedJetDeltaR", ";#DeltaR(reco,gen);matched jets",
                       80, 0.0, kMatchDeltaRMax);
    TH1D matchStatus("jetMatchStatus", ";status;jets", 3, 0.5, 3.5);
    matchStatus.GetXaxis()->SetBinLabel(1, "selected reco");
    matchStatus.GetXaxis()->SetBinLabel(2, "matched");
    matchStatus.GetXaxis()->SetBinLabel(3, "unmatched");

    TFile* qhatFile = TFile::Open(gQHatCountsPath.c_str(), "read");
    TH1D* qhatCounts = (!qhatFile || qhatFile->IsZombie())
        ? 0 : dynamic_cast<TH1D*>(qhatFile->Get("qhatCounts"));
    if (!qhatCounts || qhatCounts->GetNbinsX() != kQHatBins) {
        std::cerr << "Invalid global q-hat count file: " << gQHatCountsPath << '\n';
        return;
    }
    for (int bin = 1; bin <= kQHatBins; ++bin) {
        if (qhatCounts->GetXaxis()->GetBinLowEdge(bin) != kQHatBoundaries[bin - 1] ||
            qhatCounts->GetXaxis()->GetBinUpEdge(bin) != kQHatBoundaries[bin] ||
            qhatCounts->GetBinContent(bin) <= 0.0) {
            std::cerr << "Invalid q-hat count bin " << bin << ".\n";
            return;
        }
    }

    for (std::size_t fileIndex = 0; fileIndex < fileList.size(); ++fileIndex) {
        fFile = TFile::Open(fileList[fileIndex].c_str(), "read");
        if (!fFile || fFile->IsZombie()) {
            std::cerr << "Cannot open input file: " << fileList[fileIndex] << '\n';
            return;
        }
        TTree* tree = dynamic_cast<TTree*>(fFile->Get(kTreePath));
        if (!tree) {
            std::cerr << "Missing tree '" << kTreePath << "' in "
                      << fileList[fileIndex] << '\n';
            fFile->Close();
            return;
        }
        if (!hasRequiredBranches(tree)) {
            std::cerr << "Required branch check failed for " << fileList[fileIndex]
                      << '\n';
            fFile->Close();
            return;
        }
        Init(tree);

        const Long64_t entries = tree->GetEntries();
        for (Long64_t entry = 0; entry < entries; ++entry) {
            if (LoadTree(entry) < 0 || fChain->GetEntry(entry) <= 0) {
                std::cerr << "Failed to read entry " << entry << " in "
                          << fileList[fileIndex] << '\n';
                return;
            }
            if (!collectionsAreConsistent(*this)) {
                std::cerr << "Inconsistent collection sizes at entry " << entry
                          << " in " << fileList[fileIndex] << '\n';
                return;
            }
            const int qhatBin = qhatCounts->FindFixBin(genQScale);
            if (qhatBin < 1 || qhatBin > kQHatBins) {
                std::cerr << "q-hat value outside configured range: "
                          << genQScale << '\n';
                return;
            }
            const double qhatWeight =
                kCrossSectionsPb[qhatBin - 1] / qhatCounts->GetBinContent(qhatBin);
            const double eventWeight = qhatWeight;
            qhatWeighted.Fill(genQScale, qhatWeight);

            // Mirror the revised Run 3 event logic: require at least one
            // reconstructed jet above threshold.  Run 2 has no JER/jet-ID
            // inputs in these trees, so the stored corrected jet pT is used.
            if (jetPt->empty() ||
                *std::max_element(jetPt->begin(), jetPt->end()) < kJetPtMin) continue;
            const std::vector<int> recoToGen =
                matchJets(*jetEta, *jetPhi, *genJetEta, *genJetPhi);

            for (std::size_t recoJet = 0; recoJet < jetPt->size(); ++recoJet) {
                if (std::fabs((*jetEta)[recoJet]) > kJetEtaMax ||
                    (*jetPt)[recoJet] < kJetPtMin) continue;
                matchStatus.Fill(1.0, eventWeight);
                const int genJet = recoToGen[recoJet];
                if (genJet < 0) {
                    matchStatus.Fill(3.0, eventWeight);
                    continue;
                }
                matchStatus.Fill(2.0, eventWeight);
                const double deltaEta = (*jetEta)[recoJet] - (*genJetEta)[genJet];
                const double deltaPhi = std::remainder(
                    (*jetPhi)[recoJet] - (*genJetPhi)[genJet], 2.0 * TMath::Pi());
                matchedDeltaR.Fill(std::hypot(deltaEta, deltaPhi), eventWeight);

                for (std::size_t track = 0;
                     track < (*dau_pt)[recoJet].size(); ++track) {
                    if (!passesRecoTrackSelection(*this, recoJet, track)) continue;
                    const double pt = (*dau_pt)[recoJet][track];
                    const double eta = (*dau_eta)[recoJet][track];
                    const double phi = (*dau_phi)[recoJet][track];
                    recoPt.Fill(pt, eventWeight);
                    recoEta.Fill(eta, eventWeight);
                    recoPhi.Fill(phi, eventWeight);
                    recoPtEta.Fill(pt, eta, eventWeight);
                    recoPhiEta.Fill(phi, eta, eventWeight);
                    recoPhiPt.Fill(phi, pt, eventWeight);
                }
                for (std::size_t track = 0;
                     track < (*genDau_pt)[genJet].size(); ++track) {
                    if (!passesGenTrackSelection(*this, genJet, track)) continue;
                    const double pt = (*genDau_pt)[genJet][track];
                    const double eta = (*genDau_eta)[genJet][track];
                    const double phi = (*genDau_phi)[genJet][track];
                    genPt.Fill(pt, eventWeight);
                    genEta.Fill(eta, eventWeight);
                    genPhi.Fill(phi, eventWeight);
                    genPtEta.Fill(pt, eta, eventWeight);
                    genPhiEta.Fill(phi, eta, eventWeight);
                    genPhiPt.Fill(phi, pt, eventWeight);
                }
            }
        }
        fChain = 0;
        fFile->Close();
        delete fFile;
        fFile = 0;
    }

    TFile output(gOutputPath.c_str(), "recreate");
    if (output.IsZombie()) {
        std::cerr << "Cannot create output file: " << gOutputPath << '\n';
        return;
    }
    recoPt.Write(); recoEta.Write(); recoPhi.Write();
    recoPtEta.Write(); recoPhiEta.Write(); recoPhiPt.Write();
    genPt.Write(); genEta.Write(); genPhi.Write();
    genPtEta.Write(); genPhiEta.Write(); genPhiPt.Write();
    qhatWeighted.Write(); matchedDeltaR.Write(); matchStatus.Write();
    output.Close();
    gSucceeded = true;
}

int main(int argc, const char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: run2_efficiency_correction <file-list> <production-era> "
                     "<global-qhat-counts.root> <output.root>\n";
        return 2;
    }
    const std::string fileListPath = argv[1];
    const std::string productionEra = argv[2];
    if (!isSupportedEra(productionEra)) {
        std::cerr << "Unsupported production era '" << productionEra
                  << "'. Choose 2016APV, 2016, 2017, or 2018.\n";
        return 1;
    }
    gQHatCountsPath = argv[3];
    gOutputPath = argv[4];

    std::ifstream input(fileListPath.c_str());
    if (!input) {
        std::cerr << "Cannot open input file list: " << fileListPath << '\n';
        return 1;
    }
    std::vector<std::string> files;
    std::string path;
    while (input >> path) files.push_back(path);
    if (files.empty()) {
        std::cerr << "Input file list is empty: " << fileListPath << '\n';
        return 1;
    }

    MyClass reader(files);
    reader.Loop(0, fileListPath);
    return gSucceeded ? 0 : 1;
}
