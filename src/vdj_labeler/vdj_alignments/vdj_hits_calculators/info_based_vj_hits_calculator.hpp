#pragma once

#include "ig_gene_hits_calculator.hpp"
#include "../vj_alignment_info.hpp"
#include "../aligners/gene_segment_aligner.hpp"

class InfoBasedVJHitsCalculator : public IgGeneHitsCalculator {
    IgGeneType gene_type_;
    const FastqReadArchive& read_archive_;
    const VJAlignmentInfo& vj_alignment_info_;
    GeneSegmentAligner& gene_aligner_;

public:
    InfoBasedVJHitsCalculator(IgGeneType gene_type,
                             const FastqReadArchive& read_archive,
                             const VJAlignmentInfo& vj_alignment_info,
                             GeneSegmentAligner& gene_aligner) :
        gene_type_(gene_type),
        read_archive_(read_archive),
        vj_alignment_info_(vj_alignment_info),
        gene_aligner_(gene_aligner) {
        assert(gene_type_ == IgGeneType::variable_gene or gene_type_ == IgGeneType::join_gene);
    }

    IgGeneSegmentHitsPtr ComputeHits(ReadPtr read_ptr);
};