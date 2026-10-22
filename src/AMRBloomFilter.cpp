#include "AMRBloomFilter.h"
#include <algorithm>

AMRBloomFilter::AMRBloomFilter(size_t bitSize, size_t hashCount)
    : bits(bitSize, false), hashCount(hashCount) {}

void AMRBloomFilter::add(const std::string& kmer) {
    add(kmer.data(), kmer.size());
}

void AMRBloomFilter::add(const char* kmer, size_t len) {
    if (!kmer || len == 0 || bits.empty()) return;

    size_t sz = bits.size();
    uint64_t h1 = hash1(kmer, len) % sz;
    uint64_t h2 = hash2(kmer, len) % sz;
    bits[h1] = true;
    bits[h2] = true;

    if (hashCount >= 3) {
        uint64_t h3 = hash3(kmer, len) % sz;
        bits[h3] = true;
    }
}

bool AMRBloomFilter::contains(const std::string& kmer) const {
    return contains(kmer.data(), kmer.size());
}

bool AMRBloomFilter::contains(const char* kmer, size_t len) const {
    if (!kmer || len == 0 || bits.empty()) return false;

    size_t sz = bits.size();
    uint64_t h1 = hash1(kmer, len) % sz;
    if (!bits[h1]) return false;

    uint64_t h2 = hash2(kmer, len) % sz;
    if (!bits[h2]) return false;

    if (hashCount >= 3) {
        uint64_t h3 = hash3(kmer, len) % sz;
        if (!bits[h3]) return false;
    }

    return true; // Might be a false positive, never a false negative
}

void AMRBloomFilter::clear() {
    std::fill(bits.begin(), bits.end(), false);
}
