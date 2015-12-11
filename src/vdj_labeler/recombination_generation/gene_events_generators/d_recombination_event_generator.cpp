#include "d_recombination_event_generator.hpp"

using namespace std;

int DRecombinationEventGenerator::ComputeLeftEvent(IgGeneAlignmentPtr d_alignment) {
    if(d_alignment->Positions().GeneStartPos() != 0)
        return 1;
    return int(max_palindrome_) * -1;
}

size_t DRecombinationEventGenerator::ComputeMaxRightConsistentCleavage(IgGeneAlignmentPtr d_alignment,
                                                                    int left_event_size) {
    int left_positive_len = max<int>(0, left_event_size);
    return size_t(min<int>(int(max_cleavage_), int(d_alignment->ReadAlignmentLength()) - left_positive_len));
}

void DRecombinationEventGenerator::GenerateRightConsistentEvents(IgGeneAlignmentPtr d_alignment, int left_event_size,
                                                                 IgGeneRecombinationEventStoragePtr d_events) {
    if(d_alignment->Positions().GeneEndPos() != d_alignment->GeneLength() - 1)
        return;
    int max_right_palindrome = int(max_palindrome_) * -1;
    int max_right_cleavage = int(ComputeMaxRightConsistentCleavage(d_alignment, left_event_size));
    //cout << "Max right cleavage of D segment: " << max_right_cleavage << endl;
    for(int relen = max_right_palindrome; relen <= max_right_cleavage; relen++) {
        //cout << "Left event size: " << left_event_size << ", right event size: " << relen << endl;
        d_events->AddEvent(CleavedIgGeneAlignment(d_alignment, left_event_size, relen,
                                                  shm_calculator_.ComputeNumberSHMs(d_alignment,
                                                                                    left_event_size, relen)));
    }
}

IgGeneRecombinationEventStoragePtr DRecombinationEventGenerator::ComputeEvents(IgGeneAlignmentPtr d_alignment) {
    IgGeneRecombinationEventStoragePtr d_events(new IgGeneRecombinationEventStorage(IgGeneType::diversity_gene));
    if(d_alignment->IsEmpty())
        return d_events;
    int left_bound = ComputeLeftEvent(d_alignment);
    int max_left_cleavage = int(min<size_t>(d_alignment->ReadAlignmentLength(), max_cleavage_));
    //cout << " Max left cleavage of D segment: " << max_left_cleavage << endl;
    // we iterate from max allowed palindrome to max allowed cleavage and
    // consider that this event occurred at the start of D segment
    for(int elen = left_bound; elen <= max_left_cleavage; elen++) {
        //cout << "============================" << endl;
        GenerateRightConsistentEvents(d_alignment, elen, d_events);
        //cout << "============================" << endl;
    }
    return d_events;
}