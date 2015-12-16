#pragma once

#include "gene_segment_aligner.hpp"

class RightVTailAligner : public GeneSegmentAligner {
    size_t right_shift_;

    void RefineAlignmentPositionsForEmptyAlign(IgGeneAlignmentPtr alignment_ptr);

    IgGeneAlignmentPositions RefineAlignmentPositions(IgGeneAlignmentPositions alignment_positions,
                                                      seqan::Align<Dna5String> &alignment);

public:
    RightVTailAligner(size_t right_shift = 0) :
        right_shift_(right_shift) { }

    IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions);
};
