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
    std::cout << row1 << endl;
    std::cout << row2 << endl;
    size_t row_length = length(row1);
    for(size_t i = 0; i < row_length; i++)
        if(row1[i] != row2[i])
            num_shms_++;
}

size_t IgGeneAlignment::SHMsNumber() {
    if(num_shms_ == size_t(-1))
        ComputeSHMsNumber();
    return num_shms_;
}

ostream& operator<<(ostream &out, const IgGeneAlignment& ig_gene_alignment) {
    out << ig_gene_alignment.Positions() << endl;
    out << ig_gene_alignment.Alignment();
    return out;
}