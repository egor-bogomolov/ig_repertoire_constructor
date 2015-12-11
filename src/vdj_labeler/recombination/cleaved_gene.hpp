#pragma once

#include "../vdj_alignments/alignment_structs.hpp"

class CleavedIgGeneAlignment {
    IgGeneAlignmentPtr gene_alignment_ptr_;
    int left_cleavage_length_;
    int right_cleavage_length_;
    size_t num_shms_;

public:
    CleavedIgGeneAlignment(IgGeneAlignmentPtr gene_alignment_ptr,
                           int left_cleavage_length,
                           int right_cleavage_length,
                           size_t num_shms = 0) :
            gene_alignment_ptr_(gene_alignment_ptr),
            left_cleavage_length_(left_cleavage_length),
            right_cleavage_length_(right_cleavage_length),
            num_shms_(num_shms) { }

    CleavedIgGeneAlignment(const CleavedIgGeneAlignment& cleaved_gene) {
        gene_alignment_ptr_ = cleaved_gene.gene_alignment_ptr_;
        left_cleavage_length_ = cleaved_gene.left_cleavage_length_;
        right_cleavage_length_ = cleaved_gene.right_cleavage_length_;
        num_shms_ = cleaved_gene.num_shms_;
    }

    IgGeneAlignmentPtr GeneAlignment() const { return gene_alignment_ptr_; }

    // negative length of cleavage means existence of the left palindrome of this length
    // positive length of cleavage shows length of gene cleavage
    int LeftCleavageLength() const { return left_cleavage_length_; }

    int RightCleavageLength() const { return right_cleavage_length_; }

    size_t SHMsNumber() const { return num_shms_ + gene_alignment_ptr_->SHMsNumber(); }

    size_t GeneId() const { return gene_alignment_ptr_->GeneId(); }

    size_t StartReadPosition() const { return size_t(int(gene_alignment_ptr_->Positions().ReadStartPos()) +
                                                             left_cleavage_length_); }

    size_t EndReadPosition() const { return size_t(int(gene_alignment_ptr_->Positions().ReadEndPos()) +
                                                           right_cleavage_length_ * -1); }
};

std::ostream& operator<<(std::ostream &out, const CleavedIgGeneAlignment& cleaved_gene);