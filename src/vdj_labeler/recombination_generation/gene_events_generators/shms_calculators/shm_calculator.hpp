#pragma once

#include "../../../vdj_alignments/alignment_structs.hpp"

class SHMsCalculator {
public:
    virtual int ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                                  int left_cleavage_length,
                                  int right_cleavage_length) = 0;

    virtual int ComputeNumberSHMsForLeftEvent(IgGeneAlignmentPtr gene_alignment,
                                              int left_cleavage_length) = 0;

    virtual int ComputeNumberSHMsForRightEvent(IgGeneAlignmentPtr gene_alignment,
                                              int right_cleavage_length) = 0;

    virtual ~SHMsCalculator() { }
};