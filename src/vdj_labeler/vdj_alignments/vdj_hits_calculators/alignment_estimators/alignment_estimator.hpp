#pragma once

#include "../../alignment_structs.hpp"

class AlignmentEstimator {
public:
    virtual bool AlignmentIsGood(IgGeneAlignmentPtr ig_gene_alignment) = 0;
    virtual ~AlignmentEstimator() { }
};