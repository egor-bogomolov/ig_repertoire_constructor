#include "recombination_event_storage.hpp"

void IgGeneRecombinationEventStorage::AddEvent(CleavedIgGeneAlignment new_gene_event) {
    if(CheckConsistency(new_gene_event))
        recombination_events_.push_back(new_gene_event);
}

CleavedIgGeneAlignment IgGeneRecombinationEventStorage::GetByIndex(size_t index) {
    assert(index < size());
    return recombination_events_[index];
}