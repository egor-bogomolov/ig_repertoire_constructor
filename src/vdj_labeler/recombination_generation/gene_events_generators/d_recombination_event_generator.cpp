#include "d_recombination_event_generator.hpp"

using namespace std;

size_t DRecombinationEventGenerator::ComputeMaxRightConsistentCleavage(IgGeneAlignmentPtr d_alignment,
                                                                    int left_event_size) {
    int left_positive_len = max<int>(0, left_event_size);
    return size_t(min<int>(max_cleavage_, d_alignment->GeneLength() - left_positive_len));
}

void DRecombinationEventGenerator::GenerateRightConsistentEvents(IgGeneAlignmentPtr d_alignment, int left_event_size,
                                                                 IgGeneRecombinationEventStoragePtr d_events) {
    int max_right_palindrome = int(max_palindrome_) * -1;
    int max_right_cleavage = int(ComputeMaxRightConsistentCleavage(d_alignment, left_event_size));
    for(int relen = max_right_palindrome; relen <= max_right_cleavage; relen++)
        d_events->AddEvent(CleavedIgGeneAlignment(d_alignment, left_event_size, relen,
                                                  shm_calculator_.ComputeNumberSHMs(d_alignment,
                                                                                    left_event_size, relen)));
}

IgGeneRecombinationEventStoragePtr DRecombinationEventGenerator::ComputeEvents(IgGeneAlignmentPtr d_alignment) {
    IgGeneRecombinationEventStoragePtr d_events(new IgGeneRecombinationEventStorage(IgGeneType::diversity_gene));
    int max_left_palindrome = int(max_palindrome_) * -1;
    int max_left_cleavage = int(min<size_t>(d_alignment->GeneLength(), max_cleavage_));
    // we iterate from max allowed palindrome to max allowed cleavage and
    // consider that this event occurred at the start of D segment
    for(int elen = max_left_palindrome; elen <= max_left_cleavage; elen++) {
        GenerateRightConsistentEvents(d_alignment, elen, d_events);
    }
    return d_events;
}