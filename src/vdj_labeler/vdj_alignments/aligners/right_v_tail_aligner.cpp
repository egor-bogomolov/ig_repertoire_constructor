//
// Created by yana on 17.11.15.
//

#include "right_v_tail_aligner.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;
using namespace seqan;

void RightVTailAligner::RefineAlignmentPositions(IgGeneAlignmentPtr alignment_ptr) {
    // todo: debug me!
    auto row1 = row(alignment_ptr->Alignment(), 0);
    auto row2 = row(alignment_ptr->Alignment(), 1);
    size_t alignment_length = length(row1);
    size_t start_read_pos = alignment_ptr->Positions().alignment.query_pos.second + 1;
    size_t end_read_pos = start_read_pos + toSourcePosition(row1, alignment_length - 1);
    size_t start_gene_pos = alignment_ptr->Positions().alignment.subject_pos.second + 1;
    size_t end_gene_pos = start_gene_pos + toSourcePosition(row2, alignment_length - 1);
    alignment_ptr->RefineAlignmentPositions(AlignmentPositions(
            make_pair(start_read_pos, end_read_pos),
            make_pair(start_gene_pos, end_gene_pos)));
}

IgGeneAlignmentPtr RightVTailAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    Align<Dna5String> align;
    resize(rows(align), 2);
    //cout << length(alignment_positions.ig_gene->seq()) << " " << alignment_positions.alignment.query_pos.second << endl;
    if(length(alignment_positions.ig_gene->seq()) == alignment_positions.alignment.subject_pos.second)
        return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align, -1));
    size_t tail_length = length(alignment_positions.ig_gene->seq()) - alignment_positions.alignment.subject_pos.second;
    //cout << "Tail length: " << tail_length << endl;
    auto read_segment = prefix(
            suffix(alignment_positions.read->seq, alignment_positions.alignment.query_pos.second),
            tail_length + right_shift_);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), suffix(alignment_positions.ig_gene->seq(),
                                       alignment_positions.alignment.subject_pos.second));
    int score = globalAlignment(align, Score<int, Simple>(2, -1, -3, -2));
    IgGeneAlignmentPtr v_alignment(new IgGeneAlignment(alignment_positions, align, score));
    RefineAlignmentPositions(v_alignment);
    return v_alignment;
}
