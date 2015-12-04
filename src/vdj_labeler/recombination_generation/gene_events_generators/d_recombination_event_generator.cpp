#include "d_recombination_event_generator.hpp"

IgGeneRecombinationEventStoragePtr DRecombinationEventGenerator::ComputeEvents(
        IgGeneAlignmentPtr gene_segment_alignment) {
    assert(false);
    return IgGeneRecombinationEventStoragePtr(new IgGeneRecombinationEventStorage(IgGeneType::diversity_gene));
}