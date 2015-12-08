#include "right_event_shms_calculator.hpp"

int RightEventSHMsCalculator::ComputeNumberCleavedSHMs(IgGeneAlignmentPtr gene_alignment, size_t cleavage_length) {
    assert(false);
    return -1;
}

int RightEventSHMsCalculator::ComputeNumberPalindromeSHMs(IgGeneAlignmentPtr gene_alignment, size_t palindrome_length) {
    assert(false);
    return 1;
}

int RightEventSHMsCalculator::ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                                                int left_cleavage_length,
                                                int right_cleavage_length) {
    if(right_cleavage_length == 0)
        return 0;
    // if gene was cleaved, number of SHMs would 0 or negative
    // since some unaligned nucleotides were cleaved
    if(right_cleavage_length > 0)
        return ComputeNumberCleavedSHMs(gene_alignment, size_t(right_cleavage_length));
    // if gene contains palindrome, number of SHMs would be 0 or positive
    // since some nucleotides in palindrome are mutated
    return ComputeNumberPalindromeSHMs(gene_alignment, size_t(right_cleavage_length * -1));
}