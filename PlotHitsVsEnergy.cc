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
    if (values.size() == 1) return {mean, 0.0}; // Avoid division by zero
    double variance = (sumSq - (sum * sum) / values.size()) / (values.size() - 1);
    return {mean, sqrt(variance)};
}

vector<double> filterOutliers(const vector<double>& data, double sigmaThreshold = 2.0) {
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
                     const map<double, vector<double>>& energyCorrected) {
    vector<double> energies, meanHits, meanCorrected, semHits, semCorrected;
    for (const auto& [e, hits] : energyHits) {
        auto it = energyCorrected.find(e);
        if (it == energyCorrected.end()) continue;

        // Filter outliers for this energy group (2σ threshold)
        vector<double> filteredHits = filterOutliers(hits, 2.0);
        vector<double> filteredCorrected = filterOutliers(it->second, 2.0);

        if (filteredHits.empty() || filteredCorrected.empty()) continue; // Skip if no data

        auto [mh, sh] = calculateStats(filteredHits);
        auto [mc, sc] = calculateStats(filteredCorrected);

        if (mh > 0 && mc > 0) {
            energies.push_back(e);
            meanHits.push_back(mh);
            semHits.push_back(sh);
            meanCorrected.push_back(mc);
            semCorrected.push_back(sc);
        }
    }

    TGraphErrors* grHits = new TGraphErrors(energies.size(), &energies[0], &meanHits[0], nullptr, &semHits[0]);
    TGraphErrors* grCorrected = new TGraphErrors(energies.size(), &energies[0], &meanCorrected[0], nullptr, &semCorrected[0]);

    // Define fit functions with distinct colors
    TF1* fHits = new TF1("fHits", "[0]*pow(x,[1]) + [2]", 10, 150);
    TF1* fCorrectedFree = new TF1("fCorrFree", "[0]*pow(x,[1]) + [2]", 10, 150);
    TF1* fCorrectedFixed = new TF1("fCorrFixed", "[0]*x + [1]", 10, 150);
    
    // Set line colors (changed here)
    fHits->SetLineColor(kBlue);          // Original fit: blue
    fCorrectedFree->SetLineColor(kRed); // Corrected free fit: magenta
    fCorrectedFixed->SetLineColor(kGreen+2); // Corrected fixed fit: green

    // Style data points (unchanged)
    grHits->SetMarkerStyle(20);
    grHits->SetMarkerColor(kBlue);
    grHits->SetLineColor(kBlue);
    grCorrected->SetMarkerStyle(22);
    grCorrected->SetMarkerColor(kRed);
    grCorrected->SetLineColor(kRed);

    TCanvas* c = new TCanvas(Form("c_%s", particle.c_str()), 
           Form("%s Hits vs Energy", particle.c_str()), 1200, 800);
    c->SetGrid();
    
    grHits->Fit(fHits, "R");
    grCorrected->Fit(fCorrectedFree, "R+");
    grCorrected->Fit(fCorrectedFixed, "R+");

    double maxY = 0;
    for (size_t i = 0; i < meanHits.size(); i++) {
        if (meanHits[i] + semHits[i] > maxY) maxY = meanHits[i] + semHits[i];
        if (meanCorrected[i] + semCorrected[i] > maxY) maxY = meanCorrected[i] + semCorrected[i];
    }
    grHits->SetMinimum(0); // Start y-axis at 0
    grHits->SetMaximum(maxY * 1.2); // Add 20% margin
    grHits->GetXaxis()->SetLimits(0, *max_element(energies.begin(), energies.end()) * 1.1);

    // --- Move legend to bottom-right ---
    TLegend* leg = new TLegend(0.5, 0.15, 0.85, 0.4); // (x1, y1, x2, y2)
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);

    grHits->SetTitle(Form("Number of hits per event as a function of the incident beam's energy for %s;Beam Energy [GeV];Number of Hits", particle.c_str()));
    grHits->GetXaxis()->SetLimits(0, 130);
    grHits->Draw("AP");
    grCorrected->Draw("P SAME");

    leg->AddEntry(grHits, "Original Hits", "pe");
    leg->AddEntry(fHits, Form("Fit: (%.1f#pm%.1f)E^{%.2f#pm%.2f}+(%.1f#pm%.1f)", 
                            fHits->GetParameter(0), fHits->GetParError(0),
                            fHits->GetParameter(1), fHits->GetParError(1),
                            fHits->GetParameter(2), fHits->GetParError(2)), "l");
    leg->AddEntry(grCorrected, "Corrected Hits", "pe");
    leg->AddEntry(fCorrectedFree, Form("Free: (%.1f#pm%.1f)E^{%.2f#pm%.2f}+(%.1f#pm%.1f)", 
                                     fCorrectedFree->GetParameter(0), fCorrectedFree->GetParError(0),
                                     fCorrectedFree->GetParameter(1), fCorrectedFree->GetParError(1),
                                     fCorrectedFree->GetParameter(2), fCorrectedFree->GetParError(2)), "l");
    leg->AddEntry(fCorrectedFixed, Form("#gamma=1: (%.1f#pm%.1f)E+(%.1f#pm%.1f)", 
                                      fCorrectedFixed->GetParameter(0), fCorrectedFixed->GetParError(0),
                                      fCorrectedFixed->GetParameter(1), fCorrectedFixed->GetParError(1)), "l");
    leg->Draw();

    /*TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.028);
    tex.DrawLatex(0.6, 0.25, Form("Original #chi^{2}/NDF = %.1f/%d", fHits->GetChisquare(), fHits->GetNDF()));
    tex.DrawLatex(0.6, 0.20, Form("Corrected (Free) #chi^{2}/NDF = %.1f/%d", fCorrectedFree->GetChisquare(), fCorrectedFree->GetNDF()));
    tex.DrawLatex(0.6, 0.15, Form("Corrected (#gamma=1) #chi^{2}/NDF = %.1f/%d", fCorrectedFixed->GetChisquare(), fCorrectedFixed->GetNDF()));
    */
    c->SaveAs(Form("%s_HitsVsEnergy.png", particle.c_str()));

    vector<double> ratioMeans, ratioErrors;
    vector<double> ratioEnergies;

    for (const auto& [e, hits] : energyHits) {
        auto it = energyCorrected.find(e);
        if (it == energyCorrected.end()) continue;

        // Filter outliers for both raw and corrected hits
        vector<double> filteredHits = filterOutliers(hits);
        vector<double> filteredCorrected = filterOutliers(it->second);
        
        if (filteredHits.empty() || filteredCorrected.empty()) continue;

        // Calculate means and SEMs
        auto [meanRaw, semRaw] = calculateStats(filteredHits);
        auto [meanCorr, semCorr] = calculateStats(filteredCorrected);
        
        if (meanRaw == 0) continue; // Avoid division by zero

        // Calculate ratio of means with error propagation
        double ratio = meanCorr / meanRaw;
        double ratioError = ratio * sqrt(pow(semRaw/meanRaw, 2) + pow(semCorr/meanCorr, 2));

        ratioEnergies.push_back(e);
        ratioMeans.push_back(ratio);
        ratioErrors.push_back(ratioError);
    }

    // Create and plot ratio of means
    TCanvas* c2 = new TCanvas(Form("c2_%s", particle.c_str()), 
             Form("%s Mean Corrected/Raw Ratio", particle.c_str()), 1200, 800);
    TGraphErrors* grRatio = new TGraphErrors(ratioEnergies.size(), &ratioEnergies[0],
                                           &ratioMeans[0], nullptr, &ratioErrors[0]);
    
//     // Formatting
//     grRatio->SetTitle(Form("%s: Ratio of mean number of corrected hits to mean number of detected hits;Energy [GeV];Ratio of Means", 
//                         particle.c_str()));
//     grRatio->SetMarkerStyle(24);
//     grRatio->SetMarkerColor(kBlue);
//     grRatio->SetLineColor(kBlue);
    
//     // Axis setup
//     double xMax = *max_element(ratioEnergies.begin(), ratioEnergies.end()) * 1.1;
//     grRatio->GetXaxis()->SetLimits(0, xMax);
//     grRatio->SetMinimum(1); // Force y-axis start at 0
//     grRatio->SetMaximum(*max_element(ratioMeans.begin(), ratioMeans.end()) * 1.1);

//     grRatio->Draw("AP");
//     c2->SaveAs(Form("%s_MeanRatio.png", particle.c_str()));
// }


    // Get parameters from original fits
    double a_corr = fCorrectedFree->GetParameter(0);
    double gamma_corr = fCorrectedFree->GetParameter(1);
    double c_corr = fCorrectedFree->GetParameter(2);
    double a_uncorr = fHits->GetParameter(0);
    double gamma_uncorr = fHits->GetParameter(1);
    double c_uncorr = fHits->GetParameter(2);

    // Define the ratio function using parameters from the original fits
    TF1 *fRatioFit = new TF1("fRatioFit", 
        Form("( (%.5f)*pow(x,%.5f) + (%.5f) ) / ( (%.5f)*pow(x,%.5f) + (%.5f) )", 
            a_corr, gamma_corr, c_corr, 
            a_uncorr, gamma_uncorr, c_uncorr), 
        1, 150);
    fRatioFit->SetLineColor(kBlue);
    fRatioFit->SetLineWidth(2);

    // Formatting
    grRatio->SetTitle(Form("%s: Ratio of mean number of corrected hits to mean number of detected hits;Energy [GeV];Ratio of Means", 
                        particle.c_str()));
    grRatio->SetMarkerStyle(24);
    grRatio->SetMarkerColor(kBlue);
    grRatio->SetLineColor(kBlue);
    
    // Axis setup
    double xMax = *max_element(ratioEnergies.begin(), ratioEnergies.end()) * 1.1;
    grRatio->GetXaxis()->SetLimits(0, xMax);
    grRatio->SetMinimum(1); // Force y-axis start at 0
    grRatio->SetMaximum(*max_element(ratioMeans.begin(), ratioMeans.end()) * 1.1);

    grRatio->Draw("AP");
    fRatioFit->Draw("same");

    // Add legend for the ratio plot
    TLegend *legRatio = new TLegend(0.15, 0.75, 0.4, 0.85);
    legRatio->SetBorderSize(0);
    legRatio->AddEntry(grRatio, "Data", "pe");
    legRatio->AddEntry(fRatioFit, "Fit Ratio", "l");
    legRatio->Draw();

    c2->SaveAs(Form("%s_MeanRatio.png", particle.c_str()));
}

void PlotHitsVsEnergy() {
    vector<string> filenames = {
        "Fit_Params_Electron_DHCAL_Jun2011_20250316.txt",
        "Fit_Params_Muons_DHCAL_Jun2011_20250316.txt",
        "Fit_Params_Pions_DHCAL_Jun2011_20250316.txt"
    };

    map<string, map<double, vector<double>>> energyHitsMap, energyCorrectedMap;

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

            if (p1 <= 0 || p1e <= 0) continue;

            if (part == particle) 
                energyHitsMap[particle][e].push_back(p1);
            else if (part == particle + "_Corrected") 
                energyCorrectedMap[particle][e].push_back(p1);
        }
        file.close();
    }

    for (const auto& [particle, hits] : energyHitsMap) 
        ProcessParticle(particle, hits, energyCorrectedMap[particle]);
}