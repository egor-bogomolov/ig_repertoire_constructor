#pragma once

#include "gene_segment_aligner.hpp"

class RightVTailAligner : public GeneSegmentAligner {
    size_t right_shift_;

    void RefineAlignmentPositions(IgGeneAlignmentPtr alignment_ptr);

public:
    RightVTailAligner(size_t right_shift = 0) :
        right_shift_(right_shift) { }

    IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions);
};
