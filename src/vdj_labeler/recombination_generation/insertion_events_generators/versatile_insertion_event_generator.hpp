#pragma once

#include "insertion_event_generator.hpp"

class VersatileInsertionGenerator : public InsertionEventGenerator {
public:
    InsertionEventStoragePtr ComputeInsertionEvents(CleavedIgGeneAlignment left_gene_alignment,
                                                    CleavedIgGeneAlignment right_gene_alignment);
};