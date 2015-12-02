#include "j_recombination_event_generator.hpp"

IgGeneRecombinationEventStoragePtr JRecombinationEventGenerator::ComputeEvents(
        IgGeneAlignmentPtr gene_segment_alignment) {
    assert(false);
    return IgGeneRecombinationEventStoragePtr(new IgGeneRecombinationEventStorage(IgGeneType::join_gene));
}