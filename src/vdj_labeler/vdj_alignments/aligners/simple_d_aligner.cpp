#include "simple_d_aligner.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;
using namespace seqan;

IgGeneAlignmentPtr SimpleDAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    Align<Dna5String> align;
    resize(rows(align), 2);
    //cout << length(alignment_positions.ig_gene->seq()) << " " << alignment_positions.alignment.query_pos.second << endl;
    //if(length(alignment_positions.ig_gene->seq()) == alignment_positions.alignment.subject_pos.second)
    //    return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align));
    auto read_segment = prefix(
            prefix(alignment_positions.read->seq, alignment_positions.alignment.query_pos.second),
            alignment_positions.alignment.query_pos.second - alignment_positions.alignment.query_pos.first);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), alignment_positions.ig_gene->seq());
    cout << "Ig gene: " << alignment_positions << endl;
    cout << "Read segment: " << read_segment << endl;
    Score<int, Simple> scoringScheme(2, -1, -2, -1);
    localAlignment(align, scoringScheme);
    cout << align << endl;
    return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align));
}