#ifndef NUCLEOTIDE_BITPACK_H
#define NUCLEOTIDE_BITPACK_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * 2-bit nucleotide encoding:
 * A = 0b00 (0)
 * C = 0b01 (1)
 * G = 0b10 (2)
 * T = 0b11 (3)
 *
 * Encodes up to 32 nucleotides into a single 64-bit unsigned integer (uint64_t).
 */

// Encode a DNA string of length up to 32 into a 64-bit packed word
uint64_t pack_kmer_2bit(const char* seq, size_t length);

// Decode a 64-bit packed word back into a DNA string buffer (must be at least length + 1)
void unpack_kmer_2bit(uint64_t packed, size_t length, char* out_seq);

// Fast bitwise Hamming distance between two 2-bit encoded k-mers of identical length
// Evaluates 2-bit symbol differences using parallel bit operations in O(1) CPU cycles
int hamming_distance_2bit(uint64_t kmer_a, uint64_t kmer_b, size_t length);

// Check if a character is a valid canonical nucleotide (A, C, G, T, or lowercase)
static inline bool is_valid_canonical_base(char c) {
    switch (c) {
        case 'A': case 'a':
        case 'C': case 'c':
        case 'G': case 'g':
        case 'T': case 't':
            return true;
        default:
            return false;
    }
}

// Map nucleotide character to 2-bit integer code
static inline uint8_t base_to_2bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0x00;
        case 'C': case 'c': return 0x01;
        case 'G': case 'g': return 0x02;
        case 'T': case 't': return 0x03;
        default:            return 0x00;
    }
}

// Map 2-bit integer code back to uppercase nucleotide character
static inline char two_bit_to_base(uint8_t code) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    return BASES[code & 0x03];
}

#ifdef __cplusplus
}
#endif

#endif // NUCLEOTIDE_BITPACK_H
