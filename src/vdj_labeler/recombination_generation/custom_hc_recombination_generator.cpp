#include "logger/logger.hpp"
#include "custom_hc_recombination_generator.hpp"

using namespace std;

void CustomHeavyChainRecombinationGenerator::Clear() {
    v_storages_.clear();
    d_storages_.clear();
    j_storages_.clear();
}

void CustomHeavyChainRecombinationGenerator::ComputeVEventStorages(VDJHitsPtr vdj_hits) {
    INFO("Computation of event vector for V hits starts...");
    size_t v_events_num = 0;
    for(size_t vi = 0; vi < vdj_hits->VHitsNumber(); vi++) {
        auto v_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::variable_gene, vi);
        auto v_events = v_events_generator_.ComputeEvents(v_alignment);
        v_events_num += v_events->size();
        v_storages_.push_back(v_events);
    }
    INFO(v_events_num << " events were computed for " << vdj_hits->VHitsNumber() << " V hits");
}

void CustomHeavyChainRecombinationGenerator::ComputeDEventStorages(VDJHitsPtr vdj_hits) {
    INFO("Computation of events vector for D hits starts...");
    size_t d_events_num = 0;
    for(size_t di = 0; di < vdj_hits->DHitsNumber(); di++) {
        auto d_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::diversity_gene, di);
        auto d_events = d_events_generator_.ComputeEvents(d_alignment);
        d_events_num += d_events->size();
        d_storages_.push_back(d_events);
    }
    INFO(d_events_num << " events were computed for " << vdj_hits->DHitsNumber() << " D hits");
    if(d_events_num == 0) {
        assert(vdj_hits->DHitsNumber() == 1);
        d_storages_[0]->AddEvent(CleavedIgGeneAlignment(
                vdj_hits->GetAlignmentByIndex(IgGeneType::diversity_gene, 0), 0, 0));
        INFO("Himeric event was added to D event storage");
    }
}

void CustomHeavyChainRecombinationGenerator::ComputeJEventStorages(VDJHitsPtr vdj_hits) {
    size_t j_events_num = 0;
    INFO("Computation of events vector for J segments starts...");
    for(size_t ji = 0; ji < vdj_hits->JHitsNumber(); ji++) {
        auto j_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::join_gene, ji);
        auto j_events = j_events_generator_.ComputeEvents(j_alignment);
        j_storages_.push_back(j_events);
        j_events_num += j_events->size();
    }
    INFO(j_events_num << " events were computed for " << vdj_hits->JHitsNumber() << " J hits");
}

std::pair<NongenomicInsertion, NongenomicInsertion> CustomHeavyChainRecombinationGenerator::RefineNongenomicInsertions(
        NongenomicInsertion vd_insertion,
        NongenomicInsertion dj_insertion) {
    NongenomicInsertion vd_2(vd_insertion.StartPosition(), dj_insertion.EndPosition());
    NongenomicInsertion dj_2(dj_insertion.EndPosition(), dj_insertion.EndPosition() - 1);
    return make_pair(vd_2, dj_2);
}

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::CreateRecombinations(
        HcRecombinationStoragePtr recombination_storage,
        CleavedIgGeneAlignment v_gene,
        CleavedIgGeneAlignment d_gene,
        CleavedIgGeneAlignment j_gene,
        InsertionEventStoragePtr vd_insertions,
        InsertionEventStoragePtr dj_insertions) {
    for(auto vd_it = vd_insertions->cbegin(); vd_it != vd_insertions->cend(); vd_it++)
        for(auto dj_it = dj_insertions->cbegin(); dj_it != dj_insertions->cend(); dj_it++) {
            NongenomicInsertion vd_insersion = *vd_it;
            NongenomicInsertion dj_insersion = *dj_it;
            if(d_gene.GeneAlignment()->IsEmpty()) {
                auto refined_insersions = RefineNongenomicInsertions(*vd_it, *dj_it);
                vd_insersion = refined_insersions.first;
                dj_insersion = refined_insersions.second;
            }
            HCRecombination recombination(recombination_storage->Read(),
                                          v_gene, d_gene, j_gene,
                                          vd_insersion, dj_insersion);
            INFO(recombination);
            INFO("-------")
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
    INFO(v_events->size() << " V events were computed");
    INFO(d_events->size() << " D events were computed");
    INFO(j_events->size() << " J events were computed");
    for(auto vit = v_events->cbegin(); vit != v_events->cend(); vit++) {
        size_t v_end_position = (*vit).EndReadPosition();
        auto d_range = d_events->GetIndexRangeFromLeftPosition(v_end_position + 1);
        INFO("D range: " << d_range.first << " - " << d_range.second);
        for(size_t di = d_range.first; di <= d_range.second; di++) {
            auto vd_insertions = vd_insertion_generator_.ComputeInsertionEvents(*vit, (*d_events)[di]);
            size_t d_end_position = (*d_events)[di].EndReadPosition();
            auto j_range = j_events->GetIndexRangeFromLeftPosition(d_end_position + 1);
            INFO("J range: " << j_range.first << " - " << j_range.second);
            for(size_t ji = j_range.first; ji <= j_range.second; ji++) {
                auto dj_insertions = dj_insertion_generator_.ComputeInsertionEvents((*d_events)[di], (*j_events)[ji]);
                recombination_storage = CreateRecombinations(recombination_storage,
                                                             *vit, (*d_events)[di], (*j_events)[ji],
                                                             vd_insertions, dj_insertions);
            }
        }
    }
    return recombination_storage;
}

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::ComputeRecombinations(VDJHitsPtr vdj_hits) {
    Clear();
    HcRecombinationStoragePtr recombination_storage(new HcRecombinationStorage(vdj_hits->Read()));
    INFO("Generation of recombinations for read " << vdj_hits->Read()->id);
    ComputeVEventStorages(vdj_hits);
    ComputeDEventStorages(vdj_hits);
    ComputeJEventStorages(vdj_hits);
    for(size_t vi = 0; vi < vdj_hits->VHitsNumber(); vi++) {
        auto v_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::variable_gene, vi);
        IgGeneRecombinationEventStoragePtr v_events = v_storages_[vi];
        for(size_t di = 0; di < vdj_hits->DHitsNumber(); di++) {
            auto d_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::diversity_gene, di);
            IgGeneRecombinationEventStoragePtr d_events = d_storages_[di];
            for(size_t ji = 0; ji < vdj_hits->JHitsNumber(); ji++) {
                auto j_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::join_gene, ji);
                IgGeneRecombinationEventStoragePtr j_events = j_storages_[ji];
                //cout << "V. " << *v_alignment << endl;
                //cout << "D. " << *d_alignment << endl;
                //cout << "J. " << *j_alignment << endl;
                recombination_storage = CreateRecombinations(recombination_storage, v_events, d_events, j_events);
            }
        }
    }
    INFO(recombination_storage->size() << " valid recombinations were found");
    return recombination_storage;
}