#include "ig_gene_recombination_event_generator.hpp"

class JRecombinationEventGenerator : public IgGeneRecombinationEventsGenerator {

public:
    IgGeneRecombinationEventStoragePtr ComputeEvents(IgGeneAlignmentPtr gene_segment_alignment);
};