#include "versatile_shms_calculator.hpp"

int VersatileGeneSHMsCalculator::ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                                                   int left_cleavage_length,
                                                   int right_cleavage_length) {
    return int(gene_alignment->SHMsNumber()) +
            left_shms_calculator_.ComputeNumberSHMs(gene_alignment, left_cleavage_length, right_cleavage_length) +
            right_shms_calculator_.ComputeNumberSHMs(gene_alignment, left_cleavage_length, right_cleavage_length);
}

int VersatileGeneSHMsCalculator::ComputeNumberSHMsForLeftEvent(IgGeneAlignmentPtr gene_alignment,
                                                               int left_cleavage_length) {
    return left_shms_calculator_.ComputeNumberSHMsForLeftEvent(gene_alignment, left_cleavage_length);
}

int VersatileGeneSHMsCalculator::ComputeNumberSHMsForRightEvent(IgGeneAlignmentPtr gene_alignment,
                                                               int right_cleavage_length) {
    return left_shms_calculator_.ComputeNumberSHMsForRightEvent(gene_alignment, right_cleavage_length);
}