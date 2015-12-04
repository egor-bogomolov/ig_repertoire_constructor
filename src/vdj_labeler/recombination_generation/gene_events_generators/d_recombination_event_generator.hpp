#pragma once

#include "ig_gene_recombination_event_generator.hpp"

class DRecombinationEventGenerator : public IgGeneRecombinationEventsGenerator {

public:
    IgGeneRecombinationEventStoragePtr ComputeEvents(IgGeneAlignmentPtr gene_segment_alignment);
};