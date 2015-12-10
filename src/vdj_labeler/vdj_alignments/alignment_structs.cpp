#include "alignment_structs.hpp"

#include "seqan/sequence.h"
#include <seqan/stream.h>

using namespace std;

ostream& operator<<(ostream &out, const AlignmentPositions& obj) {
    out << "Query pos: " << obj.query_pos.first << " - " << obj.query_pos.second << ". Subject pos: " <<
            obj.subject_pos.first << " - " << obj.subject_pos.second;
    return out;
}

//-----------------------------------------------------------------------------

ostream& operator<<(ostream& out, const IgGeneAlignmentPositions& obj) {
    out << "Read. " << *obj.read << endl;
    out << "Gene. " << *obj.ig_gene << endl;
    out << obj.alignment << endl;
    return out;
}

//-----------------------------------------------------------------------------

void IgGeneAlignment::ComputeSHMsNumber() {
    num_shms_ = 0;
    typedef seqan::Row<DnaAlignment>::Type DnaAlignmentRow;
    DnaAlignmentRow &row1 = seqan::row(alignment_, 0);
    DnaAlignmentRow &row2 = seqan::row(alignment_, 1);
    assert(length(row1) == length(row2));
    if(length(row1) == 0)
        return;
    size_t row_length = length(row1);
    for(size_t i = 0; i < row_length; i++)
        if(row1[i] != row2[i])
            num_shms_++;
}

void IgGeneAlignment::ComputeNormalizedScore() {
    auto alignment_row = seqan::row(alignment_, 0);
    normalized_score_ = double(score_) / double(seqan::length(alignment_row));
}

void IgGeneAlignment::ComputeGapsNumber() {
    num_gaps_ = 0;
    typedef seqan::Row<DnaAlignment>::Type DnaAlignmentRow;
    DnaAlignmentRow &row1 = seqan::row(alignment_, 0);
    DnaAlignmentRow &row2 = seqan::row(alignment_, 1);
    assert(length(row1) == length(row2));
    if(length(row1) == 0)
        return;
    size_t row_length = length(row1);
    for(size_t i = 0; i < row_length; i++)
        if(row1[i] == '-' or row2[i] == '-')
            num_gaps_++;
}

void IgGeneAlignment::ComputeAlignmentLength() {
    auto row1 = seqan::row(alignment_, 0);
    auto row2 = seqan::row(alignment_, 1);
    alignment_length_ = length(row1);
    read_alignment_length_ = seqan::toSourcePosition(row1, alignment_length_ - 1) -
            seqan::toSourcePosition(row1, 0) + 1;
    gene_alignment_length_ = seqan::toSourcePosition(row2, alignment_length_ - 1) -
            seqan::toSourcePosition(row2, 0) + 1;
}

void IgGeneAlignment::RefineAlignmentPositions(AlignmentPositions alignment_positions) {
    positions_.alignment.query_pos = alignment_positions.query_pos;
    positions_.alignment.subject_pos = alignment_positions.subject_pos;
}

ostream& operator<<(ostream &out, const IgGeneAlignment& ig_gene_alignment) {
    out << ig_gene_alignment.Positions() << endl;
    out << ig_gene_alignment.ConstAlignment();
    out << "# SHMs: " << ig_gene_alignment.SHMsNumber() << ", # gaps: " << ig_gene_alignment.GapsNumber() << endl;
    out << "Score: " << ig_gene_alignment.Score() << ", normalized score: " <<
            ig_gene_alignment.NormalizedScore() << endl;
    return out;
}

//-----------------------------------------------------------------

seqan::Dna getRevCompl(seqan::Dna const & nucleotide) {
    if (nucleotide == 'A')
        return 'T';
    if (nucleotide == 'T')
        return 'A';
    if (nucleotide == 'C')
        return 'G';
    return 'C';
}