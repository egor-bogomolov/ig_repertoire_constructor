#include "custom_vdj_hits_calculator.hpp"

VDJHitsStoragePtr CustomVDJHitsCalculator::ComputeHits() {
    VDJHitsStoragePtr vdj_hits_storage(new VDJHitsStorage());
    for(auto it = read_archive_.cbegin(); it != read_archive_.cend(); it++) {
        IgGeneSegmentHitsPtr v_hits = v_hits_calculator_.ComputeHits(*it);
        IgGeneSegmentHitsPtr d_hits = d_hits_calculator_.ComputeHits(*it);
        IgGeneSegmentHitsPtr j_hits = j_hits_calculator_.ComputeHits(*it);
        VDJHitsPtr vdj_hits_ptr(new VDJHits(*it));
        vdj_hits_storage->AddVDJHits(vdj_hits_ptr);
    }
    return vdj_hits_storage;
}