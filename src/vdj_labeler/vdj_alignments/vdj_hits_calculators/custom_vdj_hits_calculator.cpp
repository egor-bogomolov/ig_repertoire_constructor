#include "custom_vdj_hits_calculator.hpp"

void CustomVDJHitsCalculator::AddHits(VDJHitsPtr vdj_hits, IgGeneSegmentHitsPtr ig_gene_hits) {
    for(auto it = ig_gene_hits->cbegin(); it != ig_gene_hits->cend(); it++)
        vdj_hits->AddIgGeneAlignment(*it);
}

VDJHitsStoragePtr CustomVDJHitsCalculator::ComputeHits() {
    VDJHitsStoragePtr vdj_hits_storage(new VDJHitsStorage());
    for(auto it = read_archive_.cbegin(); it != read_archive_.cend(); it++) {
        IgGeneSegmentHitsPtr v_hits = v_hits_calculator_.ComputeHits(*it);
        IgGeneSegmentHitsPtr d_hits = d_hits_calculator_.ComputeHits(*it);
        IgGeneSegmentHitsPtr j_hits = j_hits_calculator_.ComputeHits(*it);
        VDJHitsPtr vdj_hits_ptr(new VDJHits(*it));
        AddHits(vdj_hits_ptr, v_hits);
        AddHits(vdj_hits_ptr, d_hits);
        AddHits(vdj_hits_ptr, j_hits);
        vdj_hits_storage->AddVDJHits(vdj_hits_ptr);
    }
    return vdj_hits_storage;
}