#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

struct MarkerEntry {
    std::string id;
    std::string sequence;
    std::string targetDrug;
    std::string recordType;
};

static inline void trim(std::string& s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }
}

int main(int argc, char* argv[]) {
    std::string inputPath = (argc > 1) ? argv[1] : "data/markers_card_top.csv";
    std::string outputPath = (argc > 2) ? argv[2] : "data/markers_verified.csv";

    std::cout << "========================================================\n";
    std::cout << "   GENOME FIREWALL - NATIVE C++ CARD DATASET EXTRACTOR   \n";
    std::cout << "========================================================\n";
    std::cout << "[*] Input  : " << inputPath << "\n";
    std::cout << "[*] Output : " << outputPath << "\n\n";

    std::ifstream inFile(inputPath);
    if (!inFile.is_open()) {
        std::cerr << "[!] Error opening input file: " << inputPath << "\n";
        return 1;
    }

    std::vector<MarkerEntry> markers;
    std::string line;
    std::getline(inFile, line); // Header

    while (std::getline(inFile, line)) {
        trim(line);
        if (line.empty()) continue;

        std::stringstream ss(line);
        MarkerEntry m;
        std::getline(ss, m.id, ',');
        std::getline(ss, m.sequence, ',');
        std::getline(ss, m.targetDrug, ',');
        std::getline(ss, m.recordType, ',');

        trim(m.id);
        trim(m.sequence);
        trim(m.targetDrug);
        trim(m.recordType);

        if (!m.id.empty() && !m.sequence.empty()) {
            if (m.recordType.empty()) m.recordType = "CARD_Homolog";
            markers.push_back(m);
        }
    }
    inFile.close();

    std::cout << "[+] Parsed " << markers.size() << " authentic resistance markers.\n";

    // Write sanitized output CSV
    std::ofstream outFile(outputPath);
    if (!outFile.is_open()) {
        std::cerr << "[!] Error opening output file: " << outputPath << "\n";
        return 1;
    }

    outFile << "markerId,sequence,targetDrug,recordType\n";
    std::unordered_map<std::string, int> drugCounts;

    for (const auto& m : markers) {
        outFile << m.id << "," << m.sequence << "," << m.targetDrug << "," << m.recordType << "\n";
        drugCounts[m.targetDrug]++;
    }
    outFile.close();

    std::cout << "[+] Successfully exported to " << outputPath << "\n";
    std::cout << "[+] Classes indexed: " << drugCounts.size() << "\n";
    for (const auto& pair : drugCounts) {
        std::cout << "    - " << pair.first << ": " << pair.second << " determinants\n";
    }

    return 0;
}
