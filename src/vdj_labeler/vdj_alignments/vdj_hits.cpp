#include "vdj_hits.hpp"


void IgGeneSegmentHits::AddHit(IgGeneAlignmentPtr hit) {
    assert(hit->Positions().ig_gene->GeneType() == gene_type_);
    hits_.push_back(hit);
}

IgGeneAlignmentPtr IgGeneSegmentHits::operator[](size_t index) {
    assert(index < size());
    return hits_[index];
}

//-------------------------------------------------------------------------------
void VDJHits::AddIgGeneAlignment(IgGeneType gene_type, IgGeneAlignmentPtr alignment_ptr) {
    if(gene_type == IgGeneType::diversity_gene)
        v_gene_alignments_.AddHit(alignment_ptr);
    else if(gene_type == IgGeneType::diversity_gene)
        d_gene_alignments_.AddHit(alignment_ptr);
    else if(gene_type == IgGeneType::join_gene)
        j_gene_alignments_.AddHit(alignment_ptr);
}

IgGeneAlignmentPtr VDJHits::GetAlignmentByIndex(IgGeneType gene_type, size_t index) {
    if(gene_type == IgGeneType::variable_gene) {
        assert(index < v_gene_alignments_.size());
        return v_gene_alignments_[index];
    }
    else if(gene_type == IgGeneType::diversity_gene) {
        assert(index < d_gene_alignments_.size());
        return d_gene_alignments_[index];
    }
    assert(index < j_gene_alignments_.size());
    return j_gene_alignments_[index];
}