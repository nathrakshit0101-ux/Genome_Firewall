#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <cmath>
#include "AMRBloomFilter.h"
#include "TrieArena.h"
#include "AMRTrie.h"
#include "AMRDatabase.h"
#include "FastaMmapParser.h"
#include "LogisticRegressor.h"
#include "ScoringEngine.h"
#include "ProteomeTranslator.h"
#include "nucleotide_bitpack.h"

void testBloomFilter() {
    std::cout << "[*] Testing AMRBloomFilter..." << std::endl;
    AMRBloomFilter bf(1048576, 3);
    bf.add("ATGAGTATTCAACATTTCCGT");
    bf.add("CGTAACCGCCAGTGTCAGTCG");

    assert(bf.contains("ATGAGTATTCAACATTTCCGT") == true);
    assert(bf.contains("CGTAACCGCCAGTGTCAGTCG") == true);
    // Unrelated k-mer should be rejected by fast Bloom path
    assert(bf.contains("CCCCCCCCCCCCCCCCCCCCC") == false);
    std::cout << "    -> AMRBloomFilter tests passed." << std::endl;
}

void testAMRTrie() {
    std::cout << "[*] Testing AMRTrie..." << std::endl;
    AMRTrie trie;
    std::string kmer = "ATGAGTATTCAACATTTCCGT";
    trie.insertMarker(kmer, 42);

    int outIdx = -1;
    // Exact match
    assert(trie.searchExact(kmer, &outIdx) == true);
    assert(outIdx == 42);

    // Single SNP mismatch (A -> G at position 5)
    std::string mutatedKmer = "ATGAGTGTTCAACATTTCCGT";
    outIdx = -1;
    assert(trie.searchExact(mutatedKmer, &outIdx) == false);
    assert(trie.searchWithMismatch(mutatedKmer, 1, &outIdx) == true);
    assert(outIdx == 42);

    // Double SNP mismatch should fail with maxMismatches = 1
    std::string doubleMutated = "ATGAGTGTTCAACATTTCCAT";
    outIdx = -1;
    assert(trie.searchWithMismatch(doubleMutated, 1, &outIdx) == false);
    std::cout << "    -> AMRTrie exact & mutation tolerance tests passed." << std::endl;
}

void testLogisticRegressor() {
    std::cout << "[*] Testing LogisticRegressor..." << std::endl;
    LogisticRegressor reg;
    std::vector<int> X = {0, 0, 0};
    double pZero = reg.predictProbability("Carbapenems", X);
    // Baseline prior with no features should be approx 0.029 (susceptible)
    assert(pZero < 0.15);

    double zVal = reg.computeDotProduct({5.5, 0.0, 0.0}, {1, 0, 0}, -3.5);
    assert(std::abs(zVal - 2.0) < 1e-4);
    double pRes = reg.sigmoid(zVal);
    assert(pRes > 0.85); // High resistance probability
    std::cout << "    -> LogisticRegressor mathematical tests passed." << std::endl;
}

void testProteomeTranslator() {
    std::cout << "[*] Testing ProteomeTranslator (64-Codon NCBI Standard)..." << std::endl;
    // ATG -> M, GCT -> A, TAA -> * (Stop)
    std::string protein = ProteomeTranslator::translateToProtein("ATGGCTTAA", 0);
    assert(protein == "MA*");

    // Reverse complement of ATGGCTTAA -> TTAAGCCAT
    std::string rc = ProteomeTranslator::reverseComplement("ATGGCTTAA");
    assert(rc == "TTAAGCCAT");

    // 6-frame translation
    auto frames = ProteomeTranslator::translateAllFrames("ATGGCTTAA");
    assert(frames.size() == 6);
    std::cout << "    -> ProteomeTranslator translation & 6-frame tests passed." << std::endl;
}

void testNucleotideBitpack() {
    std::cout << "[*] Testing C Nucleotide Bitpacker & O(1) Hamming Distance..." << std::endl;
    const char* seq1 = "ACGTACGTACGT";
    uint64_t packed1 = pack_kmer_2bit(seq1, 12);

    char decoded[32];
    unpack_kmer_2bit(packed1, 12, decoded);
    assert(std::string(decoded) == "ACGTACGTACGT");

    // Single SNP mismatch: A -> G at position 0
    const char* seq2 = "GCGTACGTACGT";
    uint64_t packed2 = pack_kmer_2bit(seq2, 12);
    int dist = hamming_distance_2bit(packed1, packed2, 12);
    assert(dist == 1);

    // Identical k-mers have distance 0
    assert(hamming_distance_2bit(packed1, packed1, 12) == 0);
    std::cout << "    -> C Nucleotide Bitpack & Hamming distance tests passed." << std::endl;
}

int main() {
    std::cout << "========================================================\n";
    std::cout << "        GENOME FIREWALL v2.0 - UNIT TEST SUITE          \n";
    std::cout << "========================================================\n";

    testBloomFilter();
    testAMRTrie();
    testLogisticRegressor();
    testProteomeTranslator();
    testNucleotideBitpack();

    std::cout << "\n[+] ALL UNIT TESTS COMPLETED SUCCESSFULLY (100% PASS RATE)!\n\n";
    return 0;
}
