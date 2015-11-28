#pragma once

#include "gene_segment_aligner.hpp"

class SimpleDAligner : public GeneSegmentAligner {
public:
    IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions);
};