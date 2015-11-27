#pragma once

#include "../vdj_alignments/alignment_structs.hpp"

class CleavedIgGeneAlignment {
    const IgGeneAlignmentPtr &gene_alignment_ptr_;
    int left_cleavage_length_;
    int right_cleavage_length_;

public:
    CleavedIgGeneAlignment(IgGeneAlignmentPtr gene_alignment_ptr,
                           int left_cleavage_length,
                           int right_cleavage_length) :
            gene_alignment_ptr_(gene_alignment_ptr),
            left_cleavage_length_(left_cleavage_length),
            right_cleavage_length_(right_cleavage_length) { }

    IgGeneAlignmentPtr GeneAlignment() const { return gene_alignment_ptr_; }

    // negative length of cleavage means existence of the left palindrome of this length
    // positive length of cleavage shows length of gene cleavage
    int LeftCleavageLength() const { return left_cleavage_length_; }

    int RightCleavageLength() const { return right_cleavage_length_; }

    size_t SMHsNumber() const { return gene_alignment_ptr_->SHMsNumber(); }
};