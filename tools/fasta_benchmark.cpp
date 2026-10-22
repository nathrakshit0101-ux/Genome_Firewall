#include <iostream>
#include <chrono>
#include <iomanip>
#include "FastaMmapParser.h"
#include "AMRTrie.h"
#include "AMRDatabase.h"
#include "nucleotide_bitpack.h"

int main(int argc, char* argv[]) {
    std::string fastaPath = (argc > 1) ? argv[1] : "data/sample_real.fasta";
    std::string dbPath = (argc > 2) ? argv[2] : "data/markers_card_top.csv";

    std::cout << "========================================================\n";
    std::cout << "     GENOME FIREWALL - NATIVE C/C++ BENCHMARK HARNESS    \n";
    std::cout << "========================================================\n";
    std::cout << "[*] Target FASTA : " << fastaPath << "\n";
    std::cout << "[*] Reference DB : " << dbPath << "\n\n";

    // 1. Benchmark 2-bit packing speed
    std::cout << "[1] Benchmarking 2-bit Nucleotide Bitpacking & O(1) Hamming Distance...\n";
    const char* testKmer1 = "ATGAGTATTCAACATTTCCGT";
    const char* testKmer2 = "ATGAGTGTTCAACATTTCCGT";

    auto t0 = std::chrono::high_resolution_clock::now();
    uint64_t p1 = pack_kmer_2bit(testKmer1, 21);
    uint64_t p2 = pack_kmer_2bit(testKmer2, 21);
    int dist = hamming_distance_2bit(p1, p2, 21);

    char unpacked[32];
    unpack_kmer_2bit(p1, 21, unpacked);
    auto t1 = std::chrono::high_resolution_clock::now();
    double bitpackNs = std::chrono::duration<double, std::nano>(t1 - t0).count();

    std::cout << "    - K-mer 1: " << testKmer1 << " (Packed: 0x" << std::hex << p1 << std::dec << ")\n";
    std::cout << "    - K-mer 2: " << testKmer2 << " (Packed: 0x" << std::hex << p2 << std::dec << ")\n";
    std::cout << "    - Decoded: " << unpacked << "\n";
    std::cout << "    - Bitwise Hamming Distance: " << dist << " (SNP match)\n";
    std::cout << "    - Bitpack + Hamming Latency: " << bitpackNs << " ns\n\n";

    // 2. Benchmark Trie & Bloom Filter
    std::cout << "[2] Benchmarking 4-ary Prefix Trie & Bloom Ingestion...\n";
    AMRTrie trie;
    AMRDatabase db;
    auto dbStart = std::chrono::high_resolution_clock::now();
    db.loadFromFile(dbPath, trie);
    auto dbEnd = std::chrono::high_resolution_clock::now();
    double dbMs = std::chrono::duration<double, std::milli>(dbEnd - dbStart).count();

    std::cout << "    - Database Ingestion: " << db.getMarkerCount() << " markers in " << dbMs << " ms\n";
    std::cout << "    - Arena Memory: " << (trie.getMemoryUsage() / 1024) << " KB\n\n";

    // 3. Benchmark Memory-Mapped Zero-Copy Streaming
    std::cout << "[3] Benchmarking Zero-Copy Streaming Throughput...\n";
    try {
        FastaMmapParser parser(fastaPath);
        std::string chunk;
        size_t totalBytes = 0;
        size_t chunks = 0;

        auto streamStart = std::chrono::high_resolution_clock::now();
        while (parser.getNextChunk(chunk, 4096, 20)) {
            totalBytes += chunk.size();
            chunks++;
        }
        auto streamEnd = std::chrono::high_resolution_clock::now();
        double streamMs = std::chrono::duration<double, std::milli>(streamEnd - streamStart).count();
        double mbScanned = (double)totalBytes / (1024.0 * 1024.0);
        double mbPerSec = (streamMs > 0) ? (mbScanned / (streamMs / 1000.0)) : 0.0;

        std::cout << "    - Scanned: " << totalBytes << " bytes across " << chunks << " chunks\n";
        std::cout << "    - Streaming Latency: " << streamMs << " ms (" << mbPerSec << " MB/s)\n";
    } catch (const std::exception& e) {
        std::cerr << "[!] Error: " << e.what() << "\n";
    }

    std::cout << "\n[+] Benchmark completed successfully.\n";
    return 0;
}
