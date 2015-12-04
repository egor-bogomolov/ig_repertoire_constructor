#include "v_recombination_event_generator.hpp"

IgGeneRecombinationEventStoragePtr VRecombinationEventGenerator::ComputeEvents(IgGeneAlignmentPtr v_alignment) {
    assert(false);
    return IgGeneRecombinationEventStoragePtr(new IgGeneRecombinationEventStorage(IgGeneType::variable_gene));
}