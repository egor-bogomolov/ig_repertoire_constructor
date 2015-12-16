#include "logger/logger.hpp"
#include "j_recombination_event_generator.hpp"

using namespace std;

size_t JRecombinationEventGenerator::ComputeSHMsNumber(IgGeneAlignmentPtr j_alignment, int cleavage_length) {
    return shms_calculator_.ComputeNumberSHMs(j_alignment, cleavage_length, 0);
}

CleavedIgGeneAlignment JRecombinationEventGenerator::GenerateCleavageEvent(IgGeneAlignmentPtr j_alignment,
                                                                           size_t cleavage_length) {
    return CleavedIgGeneAlignment(j_alignment, int(cleavage_length), 0,
                                  ComputeSHMsNumber(j_alignment, int(cleavage_length)));
}

void JRecombinationEventGenerator::GenerateCleavageEvents(IgGeneAlignmentPtr j_alignment,
                                                          IgGeneRecombinationEventStoragePtr j_events) {
    // if alignment is empty, cleavage can not be observed
    if(j_alignment->IsEmpty())
        return;
    size_t min_cleavage_length = j_alignment->Positions().GeneStartPos();
    size_t max_cleavage_length = min<size_t>(max_cleavage_, j_alignment->ReadAlignmentLength());
    INFO("Min J cleavage:" << min_cleavage_length << ", max J cleavage: " << max_cleavage_length);
    for(size_t clen = min_cleavage_length; clen <= max_cleavage_length; clen++)
        j_events->AddEvent(GenerateCleavageEvent(j_alignment, clen));
}

CleavedIgGeneAlignment JRecombinationEventGenerator::GeneratePalindromicEvent(IgGeneAlignmentPtr j_alignment,
                                                                              size_t palindrome_length) {
    int event_length = int(palindrome_length) * -1;
    return CleavedIgGeneAlignment(j_alignment, event_length, 0,
                                  ComputeSHMsNumber(j_alignment, event_length));
}

void JRecombinationEventGenerator::GeneratePalindromicEvents(IgGeneAlignmentPtr j_alignment,
                                                             IgGeneRecombinationEventStoragePtr j_events) {
    if(j_alignment->Positions().GeneStartPos() != 0)
        return;
    for(int plen = int(max_palindrome_); plen > 0; plen--)
        j_events->AddEvent(GeneratePalindromicEvent(j_alignment, size_t(plen)));
}

IgGeneRecombinationEventStoragePtr JRecombinationEventGenerator::ComputeEvents(IgGeneAlignmentPtr j_alignment) {
    IgGeneRecombinationEventStoragePtr j_events(new IgGeneRecombinationEventStorage(IgGeneType::join_gene));
    INFO(*j_alignment);
    // generation of palindromic events
    INFO("Generarion of palindromic events");
    GeneratePalindromicEvents(j_alignment, j_events);
    // generation of zero event
    INFO("Generation of zero event");
    j_events->AddEvent(CleavedIgGeneAlignment(j_alignment, 0, 0, ComputeSHMsNumber(j_alignment, 0)));
    // generation of cleavage events
    INFO("Generation of cleavage events");
    GenerateCleavageEvents(j_alignment, j_events);
    return j_events;
}