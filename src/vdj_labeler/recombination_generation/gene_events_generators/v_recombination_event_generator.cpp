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
    size_t cleavage_length = min<size_t>(max_cleavage_, v_alignment->ReadAlignmentLength());
    cout << "Max cleavage length in V: " << cleavage_length << endl;
    for(size_t clen = 1; clen <= cleavage_length; clen++)
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
    for(size_t plen = 1; plen < max_palindrome_; plen++)
        v_events->AddEvent(GeneratePalindromicEvent(v_alignment, plen));
}

IgGeneRecombinationEventStoragePtr VRecombinationEventGenerator::ComputeEvents(IgGeneAlignmentPtr v_alignment) {
    IgGeneRecombinationEventStoragePtr v_events(new IgGeneRecombinationEventStorage(IgGeneType::variable_gene));
    // generation of cleavage events
    GenerateCleavageEvents(v_alignment, v_events);
    // generation of zero event
    v_events->AddEvent(CleavedIgGeneAlignment(v_alignment, 0, 0, ComputeSHMsNumber(v_alignment, 0)));
    // generation of palindromic events
    GeneratePalindromicEvents(v_alignment, v_events);
    return v_events;
}