#include "simple_d_aligner.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;
using namespace seqan;

void SimpleDAligner::RefineAlignmentPositions(IgGeneAlignmentPtr d_gene_alignment) {
    auto row1 = row(d_gene_alignment->Alignment(), 0);
    auto row2 = row(d_gene_alignment->Alignment(), 1);
    size_t alignment_length = length(row1);
    //cout << "Read. Start pos: " << toSourcePosition(row1, 0) << ", end pos: " <<
    //        toSourcePosition(row1, alignment_length - 1) << endl;
    //cout << "Gene. Start pos: " << toSourcePosition(row2, 0) << ", end pos: " <<
    //        toSourcePosition(row2, alignment_length - 1) << endl;
    size_t start_read_pos = toSourcePosition(row1, 0) + d_gene_alignment->Positions().alignment.query_pos.first;
    size_t end_read_pos = toSourcePosition(row1, alignment_length - 1) +
            d_gene_alignment->Positions().alignment.query_pos.first;
    d_gene_alignment->RefineAlignmentPositions(AlignmentPositions(
            make_pair(start_read_pos, end_read_pos),
            make_pair(toSourcePosition(row2, 0), toSourcePosition(row2, alignment_length - 1))));
}

IgGeneAlignmentPtr SimpleDAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    Align<Dna5String> align;
    resize(rows(align), 2);
    //cout << length(alignment_positions.ig_gene->seq()) << " " << alignment_positions.alignment.query_pos.second << endl;
    //if(length(alignment_positions.ig_gene->seq()) == alignment_positions.alignment.subject_pos.second)
    //    return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align));
    auto read_segment = suffix(
            prefix(alignment_positions.read->seq, alignment_positions.alignment.query_pos.second - 1),
            alignment_positions.alignment.query_pos.first);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), alignment_positions.ig_gene->seq());
    //cout << "Ig gene: " << alignment_positions << endl;
    //cout << "Read segment (" << length(read_segment) << "): " << read_segment << endl;
    int score = localAlignment(align, Score<int, Simple>(2, -1, -3, -2));
    //cout << "Score: " << score << endl;
    //cout << align << endl;
    IgGeneAlignmentPtr d_gene_alignment(new IgGeneAlignment(alignment_positions, align, score));
    RefineAlignmentPositions(d_gene_alignment);
    return d_gene_alignment;
}