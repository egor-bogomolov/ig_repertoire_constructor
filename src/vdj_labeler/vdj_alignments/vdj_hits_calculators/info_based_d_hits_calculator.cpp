#include "info_based_d_hits_calculator.hpp"

using namespace std;

IgGeneAlignmentPositions InfoBasedDHitsCalculator::ComputeDAlignmentPositions(IgGeneAlignmentPositions v_positions,
                                                                              IgGeneAlignmentPositions j_positions,
                                                                              IgGenePtr gene_ptr,
                                                                              ReadPtr read_ptr) {
    //cout << "V alignment positions: " << v_positions << endl;
    //cout << "J alignment positions: " << j_positions << endl;
    AlignmentPositions alignment_pos(make_pair(v_positions.alignment.query_pos.second,
                                               j_positions.alignment.query_pos.first),
                                    make_pair(0, length(gene_ptr->seq())));
    return IgGeneAlignmentPositions(alignment_pos, gene_ptr, read_ptr);
}

bool InfoBasedDHitsCalculator::DAlignmentIsGood(IgGeneAlignmentPtr d_alignment) {
    //assert(false);
    return d_alignment != NULL;
}

IgGeneSegmentHitsPtr InfoBasedDHitsCalculator::ComputeHits(ReadPtr read_ptr) {
    size_t read_index = read_archive_.GetIndexByRead(read_ptr);
    IgGeneAlignmentPositions v_alignment_positions = vj_alignment_info_.GetVAlignmentByIndex(read_index);
    IgGeneAlignmentPositions j_alignment_positions = vj_alignment_info_.GetJAlignmentByIndex(read_index);
    IgGeneSegmentHitsPtr d_hits_ptr(new IgGeneSegmentHits(IgGeneType::diversity_gene, read_ptr));
    for(auto d_gene = d_gene_database_.cbegin(); d_gene != d_gene_database_.cend(); d_gene++) {
        IgGeneAlignmentPositions d_alignment_pos = ComputeDAlignmentPositions(v_alignment_positions,
                                                                              j_alignment_positions,
                                                                              *d_gene,
                                                                              read_ptr);
        auto d_alignment = d_gene_aligner_.ComputeAlignment(d_alignment_pos);
        if(DAlignmentIsGood(d_alignment))
            d_hits_ptr->AddHit(d_alignment);
    }
    assert(false);
    return d_hits_ptr;
}