#pragma once

#include "ig_gene_hits_calculator.hpp"
#include "../vj_alignment_info.hpp"

class InfoBasedVHitsCalculator : public IgGeneHitsCalculator {
    const VJAlignmentInfo& vj_alignment_info_;

public:
    InfoBasedVHitsCalculator(const VJAlignmentInfo& vj_alignment_info) :
        vj_alignment_info_(vj_alignment_info) { }

    IgGeneSegmentHitsPtr ComputeHits(ReadPtr read_ptr);
};