#include "logger/logger.hpp"
#include "right_event_shms_calculator.hpp"
#include "seqan/modifier.h"

using std::cout;
using std::endl;

int RightEventSHMsCalculator::ComputeNumberCleavedSHMs(IgGeneAlignmentPtr gene_alignment, size_t cleavage_length) {
    size_t alignment_cleavage = gene_alignment->GeneLength() - 1 - gene_alignment->Positions().GeneEndPos();
    assert(cleavage_length >= alignment_cleavage);
    INFO("Computation of # SHMs in right cleavage of length " << cleavage_length);
    size_t rel_cleavage_length = cleavage_length - alignment_cleavage;
    if(rel_cleavage_length == 0)
        return 0;
    auto alignment = gene_alignment->Alignment();
    typedef seqan::Row<IgGeneAlignment::DnaAlignment>::Type DnaAlignmentRow;
    DnaAlignmentRow &row1 = seqan::row(alignment, 0);
    DnaAlignmentRow &row2 = seqan::row(alignment, 1);
    int num_shms = 0;
    int cur_cleavage = 0;
    for(int i = int(seqan::length(row1)) - 1; i >= 0; i--) {
        if(row1[i] != '-')
            cur_cleavage++;
        if(row1[i] != row2[i])
            num_shms++;
        if(cur_cleavage == rel_cleavage_length)
            break;
    }
    INFO("Cleavage length: " << cleavage_length << ", rel cleavage length: " << rel_cleavage_length);
    INFO("#SHMs: -" << num_shms);
    return -1 * num_shms;
}

int RightEventSHMsCalculator::ComputeNumberPalindromeSHMs(IgGeneAlignmentPtr gene_alignment,
                                                          size_t palindrome_length) {
    // if gene has alignment to read with gaps at the end, we can not compute
    assert(gene_alignment->Positions().GeneEndPos() == gene_alignment->GeneLength() - 1);
    INFO("Computation of #SHMs in right palindrome of length " << palindrome_length);
    int num_shms = 0;
    for(size_t i = 0; i < palindrome_length; i++) {
        size_t gene_pos = gene_alignment->GeneLength() - 1 - i;
        size_t read_pos = gene_alignment->Positions().ReadEndPos() + 1 + i;
        if(getRevCompl(gene_alignment->GeneSeq()[gene_pos]) != gene_alignment->ReadSeq()[read_pos])
            num_shms++;
    }
    INFO("#SHMs: +" << num_shms);
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