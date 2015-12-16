#include "logger/logger.hpp"
#include "info_based_d_hits_calculator.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;

AlignmentPositions InfoBasedDHitsCalculator::ComputeDPositions(IgGeneAlignmentPositions v_positions,
                                                               IgGeneAlignmentPositions j_positions) {
    return AlignmentPositions(make_pair(v_positions.alignment.query_pos.second + 1,
                                        j_positions.alignment.query_pos.first - 1),
                              make_pair(0, size_t(-1)));
}

IgGeneAlignmentPositions InfoBasedDHitsCalculator::CreateDAlignmentPositions(AlignmentPositions d_alignment_positions,
                                                   IgGenePtr gene_ptr,
                                                   ReadPtr read_ptr) {
    d_alignment_positions.subject_pos.second = seqan::length(gene_ptr->seq()) - 1;
    return IgGeneAlignmentPositions(d_alignment_positions, gene_ptr, read_ptr);
}

// todo: refactor magic number!
bool InfoBasedDHitsCalculator::DAlignmentPositionsAreGood(AlignmentPositions d_alignment_positions) {
    return d_alignment_positions.QueryAlignmentLength() >= 5;
}

IgGeneSegmentHitsPtr InfoBasedDHitsCalculator::ComputeHits(ReadPtr read_ptr) {
    size_t read_index = read_archive_.GetIndexByRead(read_ptr);
    assert(read_index == read_ptr->id);
    IgGeneAlignmentPositions v_alignment_positions = vj_alignment_info_.GetVAlignmentByReadIndex(read_index);
    IgGeneAlignmentPositions j_alignment_positions = vj_alignment_info_.GetJAlignmentByReadIndex(read_index);
    AlignmentPositions d_positions = ComputeDPositions(v_alignment_positions, j_alignment_positions);
    IgGeneSegmentHitsPtr d_hits_ptr(new IgGeneSegmentHits(IgGeneType::diversity_gene, read_ptr));
    if(!DAlignmentPositionsAreGood(d_positions)) {
        INFO("D positions are too short to generate alignment");
        INFO(d_positions);
        // add single empty alignment and return hits storage
        seqan::Align<Dna5String> align;
        seqan::resize(seqan::rows(align), 2);
        d_hits_ptr->AddHit(IgGeneAlignmentPtr(new IgGeneAlignment(
                CreateDAlignmentPositions(d_positions,
                                          *d_gene_database_.cbegin(),
                                          read_ptr),
                align, -1)));
        return d_hits_ptr;
    }
    vector<IgGeneAlignmentPtr> d_alignments;
    int max_score = 0;
    for(auto d_gene = d_gene_database_.cbegin(); d_gene != d_gene_database_.cend(); d_gene++) {
        IgGeneAlignmentPositions d_alignment_pos = CreateDAlignmentPositions(d_positions,
                                                                              *d_gene,
                                                                              read_ptr);
        auto d_alignment = d_gene_aligner_.ComputeAlignment(d_alignment_pos);
        d_alignments.push_back(d_alignment);
        max_score = max<int>(max_score, d_alignment->Score());
        //if(estimator_.AlignmentIsGood(d_alignment))
        //    d_hits_ptr->AddHit(d_alignment);
    }
    // just a stub: we select the only alignment with the best score
    for(auto it = d_alignments.begin(); it != d_alignments.end(); it++)
        if((*it)->Score() == max_score) {
            d_hits_ptr->AddHit(*it);
            break;
        }
    return d_hits_ptr;
}