#include "custom_hc_recombination_generator.hpp"

using namespace std;

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::CreateRecombinations(
        HcRecombinationStoragePtr recombination_storage,
        IgGeneRecombinationEventStoragePtr v_events,
        IgGeneRecombinationEventStoragePtr d_events,
        IgGeneRecombinationEventStoragePtr j_events) {
    return recombination_storage;
}

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::ComputeRecombinations(VDJHitsPtr vdj_hits) {
    size_t num_v_hits = vdj_hits->VHitsNumber();
    size_t num_d_hits = vdj_hits->DHitsNumber();
    size_t num_j_hits = vdj_hits->JHitsNumber();
    HcRecombinationStoragePtr recombination_storage(new HcRecombinationStorage(vdj_hits->Read()));
    for(size_t vi = 0; vi < num_v_hits; vi++) {
        auto v_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::variable_gene, vi);
        auto v_events = v_events_generator_.ComputeEvents(v_alignment);
        for(size_t di = 0; di < num_d_hits; di++) {
            auto d_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::diversity_gene, di);
            auto d_events = d_events_generator_.ComputeEvents(d_alignment);
            for(size_t ji = 0; ji < num_j_hits; ji++) {
                auto j_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::join_gene, ji);
                auto j_events = j_events_generator_.ComputeEvents(j_alignment);
                cout << "V. " << *v_alignment << endl;
                cout << "D. " << *d_alignment << endl;
                cout << "J. " << *j_alignment << endl;
                cout << "-----------------------------" << endl;
            }
        }
    }
    return recombination_storage;
}