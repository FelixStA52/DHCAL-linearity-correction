#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

using namespace std;

pair<double, double> calculateStats(const vector<double>& values) {
    if (values.empty()) return {0, 0};
    double sum = 0, sumSq = 0;
    for (double v : values) {
        sum += v;
        sumSq += v * v;
    }
    double mean = sum / values.size();
    double stdev = sqrt((sumSq - sum * sum / values.size()) / (values.size() - 1));
    double sem = stdev / sqrt(values.size());
    return {mean, sem};
}

pair<double, double> calculateMeanAndStdev(const vector<double>& values) {
    if (values.empty()) return {0, 0};
    double sum = 0.0, sumSq = 0.0;
    for (double v : values) {
        sum += v;
        sumSq += v * v;
    }
    double mean = sum / values.size();
    if (values.size() == 1) return {mean, 0.0};
    double variance = (sumSq - (sum * sum) / values.size()) / (values.size() - 1);
    return {mean, sqrt(variance)};
}

vector<double> filterOutliers(const vector<double>& data, double sigmaThreshold = 1000.0) {
    if (data.empty()) return data;
    auto [mean, stdev] = calculateMeanAndStdev(data);
    vector<double> filtered;
    for (double v : data) {
        if (abs(v - mean) <= sigmaThreshold * stdev) filtered.push_back(v);
    }
    return filtered;
}

void ProcessParticle(const string& particle,
        const map<double, vector<double>>& energyHits,
        const map<double, vector<double>>& energySigma,
        const map<double, vector<double>>& energyCorrected,
        const map<double, vector<double>>& energySigmaCorrected) {
    vector<double> energies, meanHits, meanCorrected;
    vector<double> semHits, semCorrected;
    vector<double> resolutionHits, resolutionCorrected;
    vector<double> resHitsErrors, resCorrectedErrors;

    // Process each energy point
    for (const auto& [e, hits] : energyHits) {
        auto itSigma = energySigma.find(e);
        auto itCorr = energyCorrected.find(e);
        auto itSigmaCorr = energySigmaCorrected.find(e);

        if (itSigma == energySigma.end() || 
            itCorr == energyCorrected.end() ||
            itSigmaCorr == energySigmaCorrected.end()) continue;

        vector<double> filteredHits = filterOutliers(hits);
        vector<double> filteredSigma = filterOutliers(itSigma->second);
        vector<double> filteredCorr = filterOutliers(itCorr->second);
        vector<double> filteredSigmaCorr = filterOutliers(itSigmaCorr->second);

        if (filteredHits.empty() || filteredSigma.empty() || 
            filteredCorr.empty() || filteredSigmaCorr.empty()) continue;

        auto [meanH, semH] = calculateStats(filteredHits);
        auto [sigmaH, sigmaH_err] = calculateStats(filteredSigma);
        auto [meanC, semC] = calculateStats(filteredCorr);
        auto [sigmaC, sigmaC_err] = calculateStats(filteredSigmaCorr);

        if (meanH > 0 && meanC > 0 && sigmaH > 0 && sigmaC > 0) {
            energies.push_back(e);
            meanHits.push_back(meanH);
            semHits.push_back(semH);
            meanCorrected.push_back(meanC);
            semCorrected.push_back(semC);

            double resH = sigmaH / meanH;
            double resC = sigmaC / meanC;
            double errResH = resH * sqrt(pow(sigmaH_err/sigmaH, 2) + pow(semH/meanH, 2));
            double errResC = resC * sqrt(pow(sigmaC_err/sigmaC, 2) + pow(semC/meanC, 2));

            resolutionHits.push_back(resH);
            resolutionCorrected.push_back(resC);
            resHitsErrors.push_back(errResH);
            resCorrectedErrors.push_back(errResC);
        }
    }

    // Create graphs
    TGraphErrors* grResHits = new TGraphErrors(energies.size(), &energies[0], 
                                &resolutionHits[0], nullptr, &resHitsErrors[0]);
    TGraphErrors* grResCorrected = new TGraphErrors(energies.size(), &energies[0],
                                    &resolutionCorrected[0], nullptr, &resCorrectedErrors[0]);

    // Setup fit function
    TF1* fRes = new TF1("fRes", "sqrt([0]^2/x + [1]^2/(x*x) + [2]^2)", 10, 150);
    if (particle == "Pions") {
        fRes->SetParameters(0.65, 0.0, 0.1);
        fRes->FixParameter(1, 0.0);
    } else {
        fRes->SetParameters(0.65, 0.1, 0.1);
    }

    // Create canvas and set styles
    TCanvas* cRes = new TCanvas(Form("cRes_%s", particle.c_str()), 
                Form("%s Energy Resolution", particle.c_str()), 1200, 800);
    cRes->SetGrid();

    // Calculate axis limits
    double maxX = *max_element(energies.begin(), energies.end()) * 1.1;
    double ymax = max(*max_element(resolutionHits.begin(), resolutionHits.end()),
                    *max_element(resolutionCorrected.begin(), resolutionCorrected.end())) * 1.2;

    // Create frame and set titles
    TH1F* frame = cRes->DrawFrame(0, 0, maxX, ymax);
    frame->SetTitle(Form("%s Resolution;Energy [GeV];#sigma/#mu", particle.c_str()));
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.04);

    // Set graph styles
    grResHits->SetTitle("");
    grResHits->SetMarkerStyle(20);
    grResHits->SetMarkerColor(kBlue);
    grResHits->SetLineColor(kBlue);
    grResHits->SetMarkerSize(1.5);

    grResCorrected->SetTitle("");
    grResCorrected->SetMarkerStyle(22);
    grResCorrected->SetMarkerColor(kRed);
    grResCorrected->SetLineColor(kRed);
    grResCorrected->SetMarkerSize(1.5);

    // Set fit styles
    fRes->SetLineColor(kBlue);
    fRes->SetLineStyle(2);
    fRes->SetLineWidth(2);

    TF1* fResCorr = (TF1*)fRes->Clone("fResCorr");
    fResCorr->SetLineColor(kRed);
    fResCorr->SetLineStyle(2);
    fResCorr->SetLineWidth(2);

    // Draw everything in correct order
    grResHits->Draw("P");
    grResCorrected->Draw("P SAME");
    grResHits->Fit(fRes, "R");
    fRes->Draw("same");
    grResCorrected->Fit(fResCorr, "R+");
    fResCorr->Draw("same");

    // Add legend
    TLegend* leg = new TLegend(0.15, 0.12, 0.68, 0.40);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetMargin(0.08);
    leg->SetEntrySeparation(0.008);
    leg->AddEntry(grResHits, "Original Hits", "pe");
    leg->AddEntry(fRes, Form("Original Fit: a=%.3f, b=%.3f, c=%.3f", 
            fRes->GetParameter(0), fRes->GetParameter(1), fRes->GetParameter(2)), "l");
    leg->AddEntry(grResCorrected, "Corrected Hits", "pe");
    leg->AddEntry(fResCorr, Form("Corrected Fit: a=%.3f, b=%.3f, c=%.3f", 
                fResCorr->GetParameter(0), fResCorr->GetParameter(1), fResCorr->GetParameter(2)), "l");
    leg->Draw();

    cRes->SaveAs(Form("%s_Resolution.png", particle.c_str()));

    // Cleanup
    delete grResHits;
    delete grResCorrected;
    delete fRes;
    delete fResCorr;
    delete cRes;
    }

void PlotRes() {
    vector<string> filenames = {
        "Fit_Params_Electron_DHCAL_Jun2011_20250316.txt",
        "Fit_Params_Muons_DHCAL_Jun2011_20250316.txt",
        "Fit_Params_Pions_DHCAL_Jun2011_20250316.txt"
    };

    // Maps to store μ (p1), σ (p2), and their errors
    map<string, map<double, vector<double>>> energyHitsMap; // μ (p1) for hits
    map<string, map<double, vector<double>>> energySigmaMap; // σ (p2) for hits
    map<string, map<double, vector<double>>> energyCorrectedMap; // μ (p1) for corrected
    map<string, map<double, vector<double>>> energySigmaCorrectedMap; // σ (p2) for corrected

    for (const auto& filename : filenames) {
        string particle;
        if (filename.find("Electron") != string::npos) particle = "Electron";
        else if (filename.find("Muons") != string::npos) particle = "Muons";
        else if (filename.find("Pions") != string::npos) particle = "Pions";
        else continue;

        ifstream file(filename);
        string line;
        while (getline(file, line)) {
            if (line.empty()) continue;
            stringstream ss(line);
            string part, period, run;
            double e, p0, p0e, p1, p1e, p2, p2e;
            char comma;

            getline(ss, part, ',');
            part.erase(0, part.find_first_not_of(" "));

            getline(ss, period, ',');
            ss >> e >> comma >> run >> comma >> p0 >> comma >> p0e >> comma >> p1 >> comma >> p1e >> comma >> p2 >> comma >> p2e;

            if (p1 <= 0 || p1e <= 0 || p2 <= 0 || p2e <= 0) continue;

            if (part == particle) {
                energyHitsMap[particle][e].push_back(p1); // Store μ (p1)
                energySigmaMap[particle][e].push_back(p2); // Store σ (p2)
            } else if (part == particle + "_Corrected") {
                energyCorrectedMap[particle][e].push_back(p1);
                energySigmaCorrectedMap[particle][e].push_back(p2);
            }
        }
        file.close();
    }

    // Process each particle using stored μ and σ
    for (const auto& [particle, hits] : energyHitsMap) {
        ProcessParticle(
            particle,
            hits,                  // μ for hits
            energySigmaMap[particle], // σ for hits
            energyCorrectedMap[particle], // μ for corrected
            energySigmaCorrectedMap[particle] // σ for corrected
        );
    }
}