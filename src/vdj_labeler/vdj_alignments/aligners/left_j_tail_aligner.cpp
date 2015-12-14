#include "logger/logger.hpp"
#include "left_j_tail_aligner.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;
using namespace seqan;

void LeftJTailAligner::RefineAlignmentPositionsForEmptyAlign(IgGeneAlignmentPtr alignment_ptr) {
    if(alignment_ptr->IsEmpty())
        return;
    size_t old_read_first = alignment_ptr->Positions().alignment.query_pos.first - 1;
    alignment_ptr->RefineAlignmentPositions(AlignmentPositions(
            make_pair(old_read_first, old_read_first),
            make_pair(size_t(-1), size_t(-1))));
}

void LeftJTailAligner::RefineAlignmentPositions(IgGeneAlignmentPtr alignment_ptr) {
    // todo: debug me!
    auto row1 = row(alignment_ptr->Alignment(), 0);
    auto row2 = row(alignment_ptr->Alignment(), 1);
    size_t alignment_length = length(row1);
    // alignment can contain gaps =>
    // alignment length and number of read or gene nucleotides participating in alignment can be different
    // refinement of read position
    size_t read_alignment_length = toSourcePosition(row1, alignment_length - 1) + 1;
    size_t end_read_pos = alignment_ptr->Positions().alignment.query_pos.first - 1;
    size_t start_read_pos = end_read_pos - read_alignment_length + 1;
    // refinement of gene position
    size_t gene_alignment_length = toSourcePosition(row2, alignment_length - 1) + 1;
    size_t end_gene_pos = alignment_ptr->Positions().alignment.subject_pos.first - 1;
    size_t start_gene_pos = end_gene_pos - gene_alignment_length + 1;
    alignment_ptr->RefineAlignmentPositions(AlignmentPositions(
            make_pair(start_read_pos, end_read_pos),
            make_pair(start_gene_pos, end_gene_pos)));
}

IgGeneAlignmentPtr LeftJTailAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    INFO("Computation of J alignment for positions: " << alignment_positions);
    Align<Dna5String> align;
    resize(rows(align), 2);
    if(alignment_positions.alignment.query_pos.first == 0)
        return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align, -1));
    size_t tail_length = alignment_positions.alignment.subject_pos.first;
    auto read_segment = prefix(
            suffix(alignment_positions.read->seq, alignment_positions.alignment.query_pos.first -
            tail_length - left_shift_), tail_length + left_shift_);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), prefix(alignment_positions.ig_gene->seq(), tail_length));
    int score = globalAlignment(align, Score<int, Simple>(2, -1, -3, -2));
    IgGeneAlignmentPtr j_alignment(new IgGeneAlignment(alignment_positions, align, score));
    RefineAlignmentPositions(j_alignment);
    return j_alignment;
}