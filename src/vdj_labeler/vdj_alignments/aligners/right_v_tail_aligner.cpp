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
    // todo: compute alignment positions
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
    return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align, score));
}
