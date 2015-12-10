#include "left_event_shms_calculator.hpp"

using std::cout;
using std::endl;

int LeftEventSHMsCalculator::ComputeNumberCleavedSHMs(IgGeneAlignmentPtr gene_alignment, size_t cleavage_length) {
    cout << "Computation of # SHMs in left cleavage" << endl;
    cout << *gene_alignment << endl;
    size_t abs_end_cleavage_position = gene_alignment->Positions().ReadStartPos() + cleavage_length - 1;
    size_t rel_end_cleavage_position = cleavage_length - 1;
    cout << "Cleavage length: " << cleavage_length <<
            ", abs end cleavage position: " << abs_end_cleavage_position <<
            ", rel end cleavage position: " << rel_end_cleavage_position << endl;
    auto alignment = gene_alignment->Alignment();
    typedef seqan::Row<IgGeneAlignment::DnaAlignment>::Type DnaAlignmentRow;
    DnaAlignmentRow &row1 = seqan::row(alignment, 0);
    DnaAlignmentRow &row2 = seqan::row(alignment, 1);
    size_t end_alignment_cleavage_position = seqan::toViewPosition(row2, rel_end_cleavage_position);
    int num_shms = 0;
    for(size_t i = 0; i <= end_alignment_cleavage_position; i++)
        if(row1[i] != row2[i])
            num_shms++;
    cout << "End alignment position: " << end_alignment_cleavage_position << ", # shms: " << num_shms << endl;
    cout << "------------------------------" << endl;
    return -1 * num_shms;
}

int LeftEventSHMsCalculator::ComputeNumberPalindromeSHMs(IgGeneAlignmentPtr gene_alignment, size_t palindrome_length) {
    // if gene has alignment to read with gaps at the end, we can not compute
    assert(gene_alignment->Positions().GeneStartPos() == 0);
    cout << "Computation of #SHMs in left palindrome of length " << palindrome_length << endl;
    cout << *gene_alignment << endl;
    int num_shms = 0;
    for(size_t i = 0; i < palindrome_length; i++) {
        size_t gene_pos = i;
        size_t read_pos = gene_alignment->Positions().ReadStartPos() - 1 - i;
        if(getRevCompl(gene_alignment->GeneSeq()[gene_pos]) == gene_alignment->ReadSeq()[read_pos])
            num_shms++;
    }
    cout << "#SHMs: " << num_shms << endl;
    cout << "------------------------------" << endl;
    return num_shms;
}

int LeftEventSHMsCalculator::ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                                               int left_cleavage_length,
                                               int) {
    if(left_cleavage_length == 0)
        return 0;
    // if gene was cleaved, number of SHMs would 0 or negative
    // since some unaligned nucleotides were cleaved
    if(left_cleavage_length > 0)
        return ComputeNumberCleavedSHMs(gene_alignment, size_t(left_cleavage_length));
    // if gene contains palindrome, number of SHMs would be 0 or positive
    // since some nucleotides in palindrome are mutated
    return ComputeNumberPalindromeSHMs(gene_alignment, size_t(left_cleavage_length * -1));
}