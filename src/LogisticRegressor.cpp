#include "LogisticRegressor.h"
#include "AMRDatabase.h"
#include <cmath>
#include <algorithm>
#include <unordered_set>

LogisticRegressor::LogisticRegressor() {}

double LogisticRegressor::computeDotProduct(const std::vector<double>& weights, const std::vector<int>& featureVector, double bias) const {
    double z = bias;
    size_t size = std::min(weights.size(), featureVector.size());
    for (size_t i = 0; i < size; ++i) {
        if (featureVector[i] > 0) {
            z += weights[i] * featureVector[i];
        }
    }
    return z;
}

double LogisticRegressor::sigmoid(double z) const {
    // Clamped sigmoid to prevent numerical overflow
    if (z > 20.0) return 1.0;
    if (z < -20.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-z));
}

double LogisticRegressor::predictProbability(const std::string& drugName, const std::vector<int>& featureVector) const {
    auto it = drugModels.find(drugName);
    if (it != drugModels.end()) {
        double z = computeDotProduct(it->second.weights, featureVector, it->second.bias);
        return sigmoid(z);
    }

    // Default prior for unknown drug classes: susceptible baseline
    double z = -3.5;
    return sigmoid(z);
}

void LogisticRegressor::calibrateFromDatabase(const AMRDatabase& db) {
    drugModels.clear();

    const auto& records = db.getRecords();
    const size_t numFeatures = records.size();
    if (numFeatures == 0) return;

    // Collect all unique drug classes from authentic database records
    std::unordered_set<std::string> drugClasses;
    for (const auto& rec : records) {
        if (!rec.targetDrug.empty()) {
            drugClasses.insert(rec.targetDrug);
        }
    }

    // Helper to evaluate cross-class beta-lactam homology
    auto isRelatedDrug = [](const std::string& targetDrug, const std::string& modelDrug) -> bool {
        if (targetDrug == modelDrug) return true;
        if (targetDrug == "Beta-Lactams") {
            return (modelDrug == "Penicillins" || modelDrug == "Cephalosporins" || modelDrug == "Carbapenems");
        }
        if (modelDrug == "Beta-Lactams") {
            return (targetDrug == "Penicillins" || targetDrug == "Cephalosporins" || targetDrug == "Carbapenems");
        }
        return false;
    };

    for (const auto& drug : drugClasses) {
        DrugModelWeights model;
        model.drugName = drug;
        // Prior baseline bias = -3.5 (P_0 = sigmoid(-3.5) approx 2.9% background wild-type resistance)
        model.bias = -3.5;
        model.weights.assign(numFeatures, 0.0);

        for (size_t i = 0; i < numFeatures; ++i) {
            const auto& rec = records[i];
            if (rec.recordType == "TargetGene") {
                // Presence of intact wild-type molecular target confers susceptibility
                if (rec.targetDrug == drug) {
                    model.weights[i] = -2.5;
                }
            } else if (rec.targetDrug == drug) {
                // Direct validated resistance determinant for target drug class
                model.weights[i] = 5.5; // Single hit pushes z from -3.5 to +2.0 (P = 88.1% resistant)
            } else if (isRelatedDrug(rec.targetDrug, drug)) {
                // Secondary homologous cross-resistance determinant
                model.weights[i] = 3.2;
            }
        }

        drugModels[drug] = std::move(model);
    }
}

void LogisticRegressor::initializePretrainedWeights(size_t featureVectorSize) {
    // If models are already calibrated from real database records, retain them
    if (!drugModels.empty()) return;

    static const char* defaultDrugs[] = {
        "Beta-Lactams", "Carbapenems", "Cephalosporins", "Penicillins",
        "Fluoroquinolones", "Aminoglycosides", "Tetracyclines", "Macrolides",
        "Glycopeptides", "Rifamycins", "Phenicols", "Sulfonamides", "Lincosamides"
    };

    for (const char* drug : defaultDrugs) {
        DrugModelWeights model;
        model.drugName = drug;
        model.bias = -3.5;
        model.weights.assign(featureVectorSize, 0.0);
        drugModels[drug] = std::move(model);
    }
}
