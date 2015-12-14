#include "logger/logger.hpp"
#include "left_event_shms_calculator.hpp"

using std::cout;
using std::endl;

int LeftEventSHMsCalculator::ComputeNumberCleavedSHMs(IgGeneAlignmentPtr gene_alignment, size_t cleavage_length) {
    assert(cleavage_length >= gene_alignment->Positions().GeneStartPos());
    INFO("Computation of # SHMs in left cleavage of length " << cleavage_length);
    size_t rel_cleavage_len = cleavage_length - gene_alignment->Positions().GeneStartPos();
    size_t rel_end_cleavage_position = rel_cleavage_len;
    INFO("Cleavage length: " << cleavage_length << ", rel cleavage len: " << rel_cleavage_len <<
                 ", rel end cleavage position: " << rel_end_cleavage_position);
    auto alignment = gene_alignment->Alignment();
    typedef seqan::Row<IgGeneAlignment::DnaAlignment>::Type DnaAlignmentRow;
    DnaAlignmentRow &row1 = seqan::row(alignment, 0);
    DnaAlignmentRow &row2 = seqan::row(alignment, 1);
    int num_shms = 0;
    for(size_t i = 0; i < rel_end_cleavage_position; i++)
        if(row1[i] != row2[i])
            num_shms++;
    INFO("#SHMs: -" << num_shms);
    return -1 * num_shms;
}

int LeftEventSHMsCalculator::ComputeNumberPalindromeSHMs(IgGeneAlignmentPtr gene_alignment, size_t palindrome_length) {
    // if gene has alignment to read with gaps at the end, we can not compute
    assert(gene_alignment->Positions().GeneStartPos() == 0);
    INFO("Computation of #SHMs in left palindrome of length " << palindrome_length);
    int num_shms = 0;
    for(size_t i = 0; i < palindrome_length; i++) {
        size_t gene_pos = i;
        size_t read_pos = gene_alignment->Positions().ReadStartPos() - 1 - i;
        if(getRevCompl(gene_alignment->GeneSeq()[gene_pos]) != gene_alignment->ReadSeq()[read_pos])
            num_shms++;
    }
    INFO("#SHMs: +" << num_shms);
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