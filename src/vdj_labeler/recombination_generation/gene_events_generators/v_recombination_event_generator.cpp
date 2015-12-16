#include "logger/logger.hpp"
#include "v_recombination_event_generator.hpp"

using namespace std;

// universal calculator of number of somatic hypermutations for gene alignment with fixed cleavage length
// positive cleavage_length means cleavage
// negative cleavage_length means palindromic insertion
size_t VRecombinationEventGenerator::ComputeSHMsNumber(IgGeneAlignmentPtr v_alignment, int cleavage_length) {
    return shms_calculator_.ComputeNumberSHMs(v_alignment, 0, cleavage_length);
}

CleavedIgGeneAlignment VRecombinationEventGenerator::GenerateCleavageEvent(IgGeneAlignmentPtr v_alignment,
                                                                           size_t cleavage_length) {
    return CleavedIgGeneAlignment(v_alignment, 0, int(cleavage_length),
                                  ComputeSHMsNumber(v_alignment, int(cleavage_length)));
}

void VRecombinationEventGenerator::GenerateCleavageEvents(IgGeneAlignmentPtr v_alignment,
                                                          IgGeneRecombinationEventStoragePtr v_events) {
    // if alignment is empty cleavage can not be observed
    if(v_alignment->IsEmpty())
        return;
    size_t min_cleavage = v_alignment->GeneLength() - v_alignment->Positions().GeneEndPos() - 1;
    size_t max_cleavage = min<size_t>(max_cleavage_, v_alignment->ReadAlignmentLength());
    INFO("Min cleavage: " <<  min_cleavage << ", max cleavage length in V: " << max_cleavage);
    for(size_t clen = min_cleavage; clen <= max_cleavage; clen++)
        v_events->AddEvent(GenerateCleavageEvent(v_alignment, clen));
}

CleavedIgGeneAlignment VRecombinationEventGenerator::GeneratePalindromicEvent(IgGeneAlignmentPtr v_alignment,
                                                                              size_t palindrome_length) {
    int event_length = int(palindrome_length) * - 1;
    return CleavedIgGeneAlignment(v_alignment, 0, event_length,
                                  ComputeSHMsNumber(v_alignment, event_length));
}

void VRecombinationEventGenerator::GeneratePalindromicEvents(IgGeneAlignmentPtr v_alignment,
                                                             IgGeneRecombinationEventStoragePtr v_events) {
    if(v_alignment->Positions().GeneEndPos() != v_alignment->GeneLength() - 1)
        return;
    for(size_t plen = 1; plen <= max_palindrome_; plen++)
        v_events->AddEvent(GeneratePalindromicEvent(v_alignment, plen));
}

IgGeneRecombinationEventStoragePtr VRecombinationEventGenerator::ComputeEvents(IgGeneAlignmentPtr v_alignment) {
    IgGeneRecombinationEventStoragePtr v_events(new IgGeneRecombinationEventStorage(IgGeneType::variable_gene));
    INFO(*v_alignment);
    // generation of cleavage events
    INFO("Generation of cleavage events");
    GenerateCleavageEvents(v_alignment, v_events);
    // generation of zero event
    INFO("Generation of zero event");
    v_events->AddEvent(CleavedIgGeneAlignment(v_alignment, 0, 0, ComputeSHMsNumber(v_alignment, 0)));
    // generation of palindromic events
    INFO("Generartion of palindromic events");
    GeneratePalindromicEvents(v_alignment, v_events);
    return v_events;
}