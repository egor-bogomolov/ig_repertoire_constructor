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
void VDJHits::AddIgGeneAlignment(IgGeneAlignmentPtr alignment_ptr) {
    IgGeneType gene_type = alignment_ptr->GeneType();
    AddIgGeneAlignment(gene_type, alignment_ptr);
}

void VDJHits::AddIgGeneAlignment(IgGeneType gene_type, IgGeneAlignmentPtr alignment_ptr) {
    if(gene_type == IgGeneType::variable_gene)
        v_hits_.AddHit(alignment_ptr);
    else if(gene_type == IgGeneType::diversity_gene)
        d_hits_.AddHit(alignment_ptr);
    else if(gene_type == IgGeneType::join_gene)
        j_hits_.AddHit(alignment_ptr);
}

IgGeneAlignmentPtr VDJHits::GetAlignmentByIndex(IgGeneType gene_type, size_t index) {
    if(gene_type == IgGeneType::variable_gene) {
        assert(index < v_hits_.size());
        return v_hits_[index];
    }
    else if(gene_type == IgGeneType::diversity_gene) {
        assert(index < d_hits_.size());
        return d_hits_[index];
    }
    assert(index < j_hits_.size());
    return j_hits_[index];
}