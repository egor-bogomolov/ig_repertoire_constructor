#pragma once

#include "shm_calculator.hpp"

class RightEventSHMsCalculator : public SHMsCalculator {
    int ComputeNumberCleavedSHMs(IgGeneAlignmentPtr gene_alignment, size_t cleavage_length);

    int ComputeNumberPalindromeSHMs(IgGeneAlignmentPtr gene_alignment, size_t palindrome_length);

public:
    int ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                          int left_cleavage_length,
                          int right_cleavage_length);

    int ComputeNumberSHMsForLeftEvent(IgGeneAlignmentPtr, int) { return 0; }

    int ComputeNumberSHMsForRightEvent(IgGeneAlignmentPtr gene_alignment,
                                       int right_cleavage_length);
};