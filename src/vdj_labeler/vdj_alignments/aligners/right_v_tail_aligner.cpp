//
// Created by yana on 17.11.15.
//

#include "logger/logger.hpp"
#include "right_v_tail_aligner.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;
using namespace seqan;

void RightVTailAligner::RefineAlignmentPositionsForEmptyAlign(IgGeneAlignmentPtr alignment_ptr) {
    if(alignment_ptr->IsEmpty())
        return;
    size_t old_read_second = alignment_ptr->Positions().alignment.query_pos.second + 1;
    size_t old_gene_second = length(alignment_ptr->Positions().ig_gene->seq());
    alignment_ptr->RefineAlignmentPositions(AlignmentPositions(
            make_pair(old_read_second, old_read_second),
            make_pair(old_gene_second, old_gene_second)));
}

IgGeneAlignmentPositions RightVTailAligner::RefineAlignmentPositions(IgGeneAlignmentPositions alignment_positions,
                                                                     seqan::Align<Dna5String> &alignment) {
    TRACE("Refinement of V alignment positions");
    auto row1 = row(alignment, 0);
    auto row2 = row(alignment, 1);
    // computation of alignment start positions
    size_t start_gene_pos = alignment_positions.ReadEndPos() + 1 + toSourcePosition(row1, 0);
    size_t start_read_pos = alignment_positions.GeneEndPos() + 1 + toSourcePosition(row2, 0);
    // computation of alignment end positions
    // we will find the first position from the right that is not gap on gene
    size_t end_read_pos = 0;
    size_t end_gene_pos = 0;
    for (int i = int(seqan::length(row2) - 1); i >= 0; i--)
        if (row2[i] != '-') {
            end_read_pos = alignment_positions.ReadEndPos() + 1 + toSourcePosition(row1, i);
            end_gene_pos = alignment_positions.GeneEndPos() + 1 + toSourcePosition(row2, i);
            break;
        }
    TRACE("Refined read positions: " << start_read_pos << " - " << end_read_pos);
    TRACE("Refined gene positions: " << start_gene_pos << " - " << end_gene_pos);
    return IgGeneAlignmentPositions(AlignmentPositions(make_pair(start_read_pos, end_read_pos),
                                                       make_pair(start_gene_pos, end_gene_pos)),
                                    alignment_positions.Gene(), alignment_positions.Read());
}

IgGeneAlignmentPtr RightVTailAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    TRACE("Computation of V alignment for positions: " << alignment_positions);
    Align<Dna5String> align;
    resize(rows(align), 2);
    if(length(alignment_positions.ig_gene->seq()) == alignment_positions.alignment.subject_pos.second)
        return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align, -1));
    size_t tail_length = length(alignment_positions.ig_gene->seq()) -
            alignment_positions.alignment.subject_pos.second - 1;
    if(tail_length + right_shift_ == 0) {
        IgGeneAlignmentPtr v_alignment(new IgGeneAlignment(alignment_positions, align, -1));
        RefineAlignmentPositionsForEmptyAlign(v_alignment);
        return v_alignment;
    }
    auto read_segment = infixWithLength(alignment_positions.read->seq,
                                        alignment_positions.ReadEndPos() + 1,
                                        tail_length + right_shift_);
    TRACE("Read segment (" << length(read_segment) << "): " << read_segment);
    auto gene_segment = suffix(alignment_positions.ig_gene->seq(),
                               alignment_positions.alignment.subject_pos.second + 1);
    TRACE("Gene segment (" << length(gene_segment) << "): " << gene_segment);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), gene_segment);
    int score = globalAlignment(align, Score<int, Simple>(2, -1, -10, -10));
    TRACE(align);
    auto refined_positions = RefineAlignmentPositions(alignment_positions, align);
    IgGeneAlignmentPtr v_alignment(new IgGeneAlignment(refined_positions, align, score));
    return v_alignment;
}
