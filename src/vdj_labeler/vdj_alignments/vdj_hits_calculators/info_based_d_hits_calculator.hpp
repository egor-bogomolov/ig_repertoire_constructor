#pragma once

#include "ig_gene_hits_calculator.hpp"
#include "../vj_alignment_info.hpp"
#include "../aligners/gene_segment_aligner.hpp"
#include "alignment_estimators/alignment_estimator.hpp"

class InfoBasedDHitsCalculator : public IgGeneHitsCalculator {
    const FastqReadArchive& read_archive_;
    const VJAlignmentInfo& vj_alignment_info_;
    const IgGeneDatabase& d_gene_database_;
    GeneSegmentAligner& d_gene_aligner_;
    AlignmentEstimator& estimator_;

    IgGeneAlignmentPositions ComputeDAlignmentPositions(IgGeneAlignmentPositions v_positions,
                                                        IgGeneAlignmentPositions j_positions,
                                                        IgGenePtr gene_ptr,
                                                        ReadPtr read_ptr);

public:
    InfoBasedDHitsCalculator(const FastqReadArchive& read_archive,
                             const VJAlignmentInfo& vj_alignment_info,
                             const IgGeneDatabase& d_gene_database,
                             GeneSegmentAligner& d_gene_aligner,
                             AlignmentEstimator& estimator) :
            read_archive_(read_archive),
            vj_alignment_info_(vj_alignment_info),
            d_gene_database_(d_gene_database),
            d_gene_aligner_(d_gene_aligner),
            estimator_(estimator) { }

    IgGeneSegmentHitsPtr ComputeHits(ReadPtr read_ptr);
};