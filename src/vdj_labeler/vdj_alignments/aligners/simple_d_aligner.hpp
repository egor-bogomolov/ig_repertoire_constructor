#pragma once

#include "gene_segment_aligner.hpp"

class SimpleDAligner : public GeneSegmentAligner {

    void RefineAlignmentPositions(IgGeneAlignmentPtr d_gene_alignment);

public:
    IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions);
};