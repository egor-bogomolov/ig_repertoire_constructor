#pragma once

#include "../recombination/cleaved_gene.hpp"

class IgGeneRecombinationEventStorage {
    IgGeneType gene_type_;
    std::vector<CleavedIgGeneAlignment> recombination_events_;

    bool CheckConsistency(CleavedIgGeneAlignment new_gene_event) {
        return new_gene_event.GeneAlignment()->GeneType() == gene_type_;
    }

public:
    IgGeneRecombinationEventStorage(IgGeneType gene_type) :
            gene_type_(gene_type) { }

    void AddNewEvent(CleavedIgGeneAlignment new_gene_event);

    typedef std::vector<CleavedIgGeneAlignment>::const_iterator gene_segment_event_iterator;

    gene_segment_event_iterator cbegin() const { return recombination_events_.cbegin(); }

    gene_segment_event_iterator cend() const { return recombination_events_.cend(); }

    size_t size() const { return recombination_events_.size(); }

    CleavedIgGeneAlignment GetByIndex(size_t index);
};

typedef std::shared_ptr<IgGeneRecombinationEventStorage > IgGeneRecombinationEventStoragePtr;