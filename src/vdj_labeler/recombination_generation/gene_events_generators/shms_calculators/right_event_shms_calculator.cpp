#include "right_event_shms_calculator.hpp"
#include "seqan/modifier.h"

using std::cout;
using std::endl;

int RightEventSHMsCalculator::ComputeNumberCleavedSHMs(IgGeneAlignmentPtr gene_alignment, size_t cleavage_length) {
    size_t alignment_cleavage = gene_alignment->GeneLength() - 1 - gene_alignment->Positions().GeneEndPos();
    assert(cleavage_length >= alignment_cleavage);
    cout << "Computation of # SHMs in right cleavage of length " << cleavage_length << endl;
    cout << *gene_alignment << endl;
    size_t rel_cleavage_length = cleavage_length - alignment_cleavage;
    size_t abs_start_cleavage_pos = gene_alignment->Positions().ReadEndPos() - rel_cleavage_length + 1;
    size_t rel_start_cleavage_pos = abs_start_cleavage_pos - gene_alignment->Positions().ReadStartPos();
    cout << "Cleavage length: " << cleavage_length << ", rel cleavage length: " << rel_cleavage_length <<
        ", abs start cleavage position: " << abs_start_cleavage_pos <<
        ", rel start cleavage position: " << rel_start_cleavage_pos << endl;
    auto alignment = gene_alignment->Alignment();
    typedef seqan::Row<IgGeneAlignment::DnaAlignment>::Type DnaAlignmentRow;
    DnaAlignmentRow &row1 = seqan::row(alignment, 0);
    DnaAlignmentRow &row2 = seqan::row(alignment, 1);
    size_t start_alignment_cleavage_pos = seqan::toViewPosition(row2, rel_start_cleavage_pos);
    int num_shms = 0;
    for(size_t i = start_alignment_cleavage_pos; i < seqan::length(row2); i++)
        if(row1[i] != row2[i])
            num_shms++;
    cout << "End alignment position: " << start_alignment_cleavage_pos << ", # shms: " << num_shms << endl;
    cout << "------------------------------" << endl;
    return -1 * num_shms;
}

int RightEventSHMsCalculator::ComputeNumberPalindromeSHMs(IgGeneAlignmentPtr gene_alignment,
                                                          size_t palindrome_length) {
    // if gene has alignment to read with gaps at the end, we can not compute
    assert(gene_alignment->Positions().GeneEndPos() == gene_alignment->GeneLength() - 1);
    //cout << "Computation of #SHMs in right palindrome of length " << palindrome_length << endl;
    //cout << *gene_alignment << endl;
    int num_shms = 0;
    for(size_t i = 0; i < palindrome_length; i++) {
        size_t gene_pos = gene_alignment->GeneLength() - 1 - i;
        size_t read_pos = gene_alignment->Positions().ReadEndPos() + 1 + i;
        if(getRevCompl(gene_alignment->GeneSeq()[gene_pos]) != gene_alignment->ReadSeq()[read_pos])
            num_shms++;
    }
    //cout << "#SHMs: " << num_shms << endl;
    //cout << "------------------------------" << endl;
    return num_shms;
}

int RightEventSHMsCalculator::ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                                                int,
                                                int right_cleavage_length) {
    if(right_cleavage_length == 0)
        return 0;
    // if gene was cleaved, number of SHMs would 0 or negative
    // since some unaligned nucleotides were cleaved
    if(right_cleavage_length > 0)
        return ComputeNumberCleavedSHMs(gene_alignment, size_t(right_cleavage_length));
    // if gene contains palindrome, number of SHMs would be 0 or positive
    // since some nucleotides in palindrome are mutated
    return ComputeNumberPalindromeSHMs(gene_alignment, size_t(right_cleavage_length * -1));
}