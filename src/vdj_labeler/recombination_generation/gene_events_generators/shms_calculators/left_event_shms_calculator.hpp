#include "shm_calculator.hpp"

class LeftEventSHMsCalculator : public SHMsCalculator {
    int ComputeNumberCleavedSHMs(IgGeneAlignmentPtr gene_alignment, size_t left_cleavage_length);

    int ComputeNumberPalindromeSHMs(IgGeneAlignmentPtr gene_alignment, size_t left_palindrome_length);

public:
    int ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                          int left_cleavage_length,
                          int right_cleavage_length);
};