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
    auto row1 = row(alignment_ptr->Alignment(), 0);
    auto row2 = row(alignment_ptr->Alignment(), 1);
    size_t alignment_length = length(row1);
    // computation of alignment end positions
    size_t end_read_pos = alignment_ptr->Positions().ReadStartPos() - 1;
    size_t end_gene_pos = alignment_ptr->Positions().GeneStartPos() - 1;
    // computation of alignment start positions
    // we need to skip all gapped positions
    size_t start_gene_pos = size_t(-1);
    size_t start_read_pos = size_t(-1);
    for(size_t i = 0; i < length(row1); i++)
        if(row2[i] != '-') {
            start_read_pos = end_read_pos - toSourcePosition(row1, alignment_length - 1) + toSourcePosition(row1, i);
            start_gene_pos = end_gene_pos - toSourcePosition(row2, alignment_length - 1) + toSourcePosition(row2, i);
            break;
        }
    alignment_ptr->RefineAlignmentPositions(AlignmentPositions(
            make_pair(start_read_pos, end_read_pos),
            make_pair(start_gene_pos, end_gene_pos)));
    TRACE("Refined read positions: " << start_read_pos << " - " << end_read_pos);
    TRACE("Refined gene positions: " << start_gene_pos << " - " << end_gene_pos);
}

IgGeneAlignmentPtr LeftJTailAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    TRACE("Computation of J alignment for positions: " << alignment_positions);
    Align<Dna5String> align;
    resize(rows(align), 2);
    if(alignment_positions.alignment.subject_pos.first == 0)
        return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align, -1));
    size_t tail_length = alignment_positions.GeneStartPos();
    auto read_segment = infixWithLength(alignment_positions.read->seq,
                                        alignment_positions.ReadStartPos() - tail_length - left_shift_,
                                        tail_length + left_shift_);
    TRACE("Read segment (" << length(read_segment) << "): " << read_segment);
    auto gene_segment = prefix(alignment_positions.ig_gene->seq(), tail_length);
    TRACE("Gene segment (" << length(gene_segment) << "): " << gene_segment);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), gene_segment);
    int score = globalAlignment(align, Score<int, Simple>(2, -1, -10, -10));
    TRACE(align);
    IgGeneAlignmentPtr j_alignment(new IgGeneAlignment(alignment_positions, align, score));
    RefineAlignmentPositions(j_alignment);
    return j_alignment;
}