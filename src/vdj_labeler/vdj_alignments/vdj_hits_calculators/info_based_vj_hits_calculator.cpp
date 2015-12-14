
#include "info_based_vj_hits_calculator.hpp"

IgGeneSegmentHitsPtr InfoBasedVJHitsCalculator::ComputeHits(ReadPtr read_ptr) {
    size_t read_index = read_archive_.GetIndexByRead(read_ptr);
    assert(read_index == read_ptr->id);
    IgGeneAlignmentPositions aligment_positions;
    if(gene_type_ == IgGeneType::variable_gene)
        aligment_positions = vj_alignment_info_.GetVAlignmentByReadIndex(read_index);
    else
        aligment_positions = vj_alignment_info_.GetJAlignmentByReadIndex(read_index);
    IgGeneAlignmentPtr alignment = gene_aligner_.ComputeAlignment(aligment_positions);
    IgGeneSegmentHitsPtr hits_ptr(new IgGeneSegmentHits(gene_type_, read_ptr));
    hits_ptr->AddHit(alignment);
    return hits_ptr;
}