#include "nucleotide_bitpack.h"

uint64_t pack_kmer_2bit(const char* seq, size_t length) {
    if (!seq || length == 0 || length > 32) return 0;

    uint64_t packed = 0ULL;
    for (size_t i = 0; i < length; ++i) {
        uint64_t code = (uint64_t)base_to_2bit(seq[i]);
        packed = (packed << 2) | code;
    }
    return packed;
}

void unpack_kmer_2bit(uint64_t packed, size_t length, char* out_seq) {
    if (!out_seq || length == 0 || length > 32) return;

    for (size_t i = length; i > 0; --i) {
        uint8_t code = (uint8_t)(packed & 0x03ULL);
        out_seq[i - 1] = two_bit_to_base(code);
        packed >>= 2;
    }
    out_seq[length] = '\0';
}

int hamming_distance_2bit(uint64_t kmer_a, uint64_t kmer_b, size_t length) {
    if (length == 0 || length > 32) return 0;

    // XOR reveals bit differences: any mismatch in 2-bit symbol produces non-zero in that 2-bit pair
    uint64_t xor_diff = kmer_a ^ kmer_b;
    if (length < 32) {
        uint64_t mask = (1ULL << (length * 2)) - 1ULL;
        xor_diff &= mask;
    }

    // Mask to 2-bit symbols
    // If either bit in a 2-bit pair is 1, that nucleotide was mismatched
    uint64_t odd_bits = (xor_diff & 0xAAAAAAAAAAAAAAAAULL) >> 1;
    uint64_t even_bits = (xor_diff & 0x5555555555555555ULL);
    uint64_t symbol_diff = odd_bits | even_bits;

    // Each non-zero 2-bit position in symbol_diff has a 1 at its even position (0, 2, 4...)
    // The population count directly equals the Hamming distance!
    #if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(symbol_diff);
    #else
    uint64_t v = symbol_diff;
    int count = 0;
    while (v) {
        v &= (v - 1);
        count++;
    }
    return count;
    #endif
}
