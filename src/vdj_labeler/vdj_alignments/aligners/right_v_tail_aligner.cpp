//
// Created by yana on 17.11.15.
//

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

void RightVTailAligner::RefineAlignmentPositions(IgGeneAlignmentPtr alignment_ptr) {
    // todo: debug me!
    auto row1 = row(alignment_ptr->Alignment(), 0);
    auto row2 = row(alignment_ptr->Alignment(), 1);
    size_t alignment_length = length(row1);
    // alignment can contain gaps =>
    // alignment length and number of read or gene nucleotides participating in alignment can be different
    // refinement of read positions
    size_t read_alignment_length = toSourcePosition(row1, alignment_length - 1) + 1;
    size_t start_read_pos = alignment_ptr->Positions().alignment.query_pos.second + 1;
    size_t end_read_pos = start_read_pos + toSourcePosition(row1, read_alignment_length - 1);
    // refinement of gene positions
    size_t gene_alignment_length = toSourcePosition(row2, alignment_length - 1) + 1;
    size_t start_gene_pos = alignment_ptr->Positions().alignment.subject_pos.second + 1;
    size_t end_gene_pos = start_gene_pos + toSourcePosition(row2, gene_alignment_length - 1);
    alignment_ptr->RefineAlignmentPositions(AlignmentPositions(
            make_pair(start_read_pos, end_read_pos),
            make_pair(start_gene_pos, end_gene_pos)));
}

IgGeneAlignmentPtr RightVTailAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    cout << alignment_positions << endl;
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
    auto read_segment = prefix(
            suffix(alignment_positions.read->seq, alignment_positions.alignment.query_pos.second + 1),
            tail_length + right_shift_);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), suffix(alignment_positions.ig_gene->seq(),
                                       alignment_positions.alignment.subject_pos.second + 1));
    int score = globalAlignment(align, Score<int, Simple>(2, -1, -3, -2));
    IgGeneAlignmentPtr v_alignment(new IgGeneAlignment(alignment_positions, align, score));
    RefineAlignmentPositions(v_alignment);
    return v_alignment;
}
