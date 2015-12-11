#include "custom_hc_recombination_generator.hpp"

using namespace std;

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::CreateRecombinations(
        HcRecombinationStoragePtr recombination_storage,
        CleavedIgGeneAlignment v_gene,
        CleavedIgGeneAlignment d_gene,
        CleavedIgGeneAlignment j_gene,
        InsertionEventStoragePtr vd_insertions,
        InsertionEventStoragePtr dj_insertions) {
    for(auto vd_it = vd_insertions->cbegin(); vd_it != vd_insertions->cend(); vd_it++)
        for(auto dj_it = dj_insertions->cbegin(); dj_it != dj_insertions->cend(); dj_it++) {
            HCRecombination recombination(recombination_storage->Read(), v_gene, d_gene, j_gene, *vd_it, *dj_it);
            cout << recombination << endl;
            cout << "---------------------------" << endl;
            if(recombination.Valid())
                recombination_storage->AddRecombination(recombination);
        }
    return recombination_storage;
}

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::CreateRecombinations(
        HcRecombinationStoragePtr recombination_storage,
        IgGeneRecombinationEventStoragePtr v_events,
        IgGeneRecombinationEventStoragePtr d_events,
        IgGeneRecombinationEventStoragePtr j_events) {
    cout << v_events->size() << " V events were computed" << endl;
    cout << d_events->size() << " D events were computed" << endl;
    cout << j_events->size() << " J events were computed" << endl;
    for(auto vit = v_events->cbegin(); vit != v_events->cend(); vit++)
        for(auto dit = d_events->cbegin(); dit != d_events->cend(); dit++) {
            auto vd_insertions = vd_insertion_generator_.ComputeInsertionEvents(*vit, *dit);
            for(auto jit = j_events->cbegin(); jit != j_events->cend(); jit++) {
                auto dj_insertions = dj_insertion_generator_.ComputeInsertionEvents(*dit, *jit);
                recombination_storage = CreateRecombinations(recombination_storage, *vit, *dit, *jit,
                                                             vd_insertions, dj_insertions);
            }
        }
    return recombination_storage;
}

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::ComputeRecombinations(VDJHitsPtr vdj_hits) {
    HcRecombinationStoragePtr recombination_storage(new HcRecombinationStorage(vdj_hits->Read()));
    size_t num_v_hits = vdj_hits->VHitsNumber();
    size_t num_d_hits = vdj_hits->DHitsNumber();
    size_t num_j_hits = vdj_hits->JHitsNumber();
    cout << "Generation of recombinations for read " << *(vdj_hits->Read()) << endl;
    for(size_t vi = 0; vi < num_v_hits; vi++) {
        auto v_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::variable_gene, vi);
        auto v_events = v_events_generator_.ComputeEvents(v_alignment);
        for(size_t di = 0; di < num_d_hits; di++) {
            auto d_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::diversity_gene, di);
            auto d_events = d_events_generator_.ComputeEvents(d_alignment);
            for(size_t ji = 0; ji < num_j_hits; ji++) {
                auto j_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::join_gene, ji);
                auto j_events = j_events_generator_.ComputeEvents(j_alignment);
                //cout << "V. " << *v_alignment << endl;
                //cout << "D. " << *d_alignment << endl;
                //cout << "J. " << *j_alignment << endl;
                cout << "-----------------------------" << endl;
                recombination_storage = CreateRecombinations(recombination_storage, v_events, d_events, j_events);
            }
        }
    }
    return recombination_storage;
}