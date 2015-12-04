#pragma once

#include "../insertion_event_storage.hpp"
#include "../../recombination/cleaved_gene.hpp"

class InsertionEventGenerator {
public:
    virtual InsertionEventStoragePtr ComputeInsertionEvents(CleavedIgGeneAlignment left_gene_alignment,
                                                            CleavedIgGeneAlignment right_gene_alignment) = 0;

    virtual ~InsertionEventGenerator() { }
};