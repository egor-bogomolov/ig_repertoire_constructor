#include "info_based_d_hits_calculator.hpp"

IgGeneAlignmentPositions InfoBasedDHitsCalculator::ComputeDAlignmentPositions(IgGeneAlignmentPositions v_positions,
                                                    IgGeneAlignmentPositions j_positions) {
    assert(false);
    return IgGeneAlignmentPositions();
}

bool InfoBasedDHitsCalculator::DAlignmentIsGood(IgGeneAlignmentPtr d_alignment) {
    assert(false);
    return false;
}

IgGeneSegmentHitsPtr InfoBasedDHitsCalculator::ComputeHits(ReadPtr read_ptr) {
    size_t read_index = read_archive_.GetIndexByRead(read_ptr);
    IgGeneAlignmentPositions v_alignment_positions = vj_alignment_info_.GetVAlignmentByIndex(read_index);
    IgGeneAlignmentPositions j_alignment_positions = vj_alignment_info_.GetJAlignmentByIndex(read_index);
    IgGeneAlignmentPositions d_aligment_positions = ComputeDAlignmentPositions(v_alignment_positions,
                                                                               j_alignment_positions);
    IgGeneSegmentHitsPtr d_hits_ptr(new IgGeneSegmentHits(IgGeneType::diversity_gene, read_ptr));
    for(auto d_gene = d_gene_database_.cbegin(); d_gene != d_gene_database_.cend(); d_gene++) {
        d_aligment_positions.ig_gene = *d_gene;
        auto d_alignment = d_gene_aligner_.ComputeAlignment(d_aligment_positions);
        if(DAlignmentIsGood(d_alignment))
            d_hits_ptr->AddHit(d_alignment);
    }
    return d_hits_ptr;
}