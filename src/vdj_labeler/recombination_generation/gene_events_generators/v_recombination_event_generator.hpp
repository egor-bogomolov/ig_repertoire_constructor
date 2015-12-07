#pragma once

#include "ig_gene_recombination_event_generator.hpp"

class VRecombinationEventGenerator : public IgGeneRecombinationEventsGenerator {
    size_t max_cleavage_;
    size_t max_palindrome_;

    size_t ComputeSHMsNumber(IgGeneAlignmentPtr v_alignment, int cleavage_length);

    CleavedIgGeneAlignment GenerateCleavageEvent(IgGeneAlignmentPtr v_alignment, size_t cleavage_length);

    void GenerateCleavageEvents(IgGeneAlignmentPtr v_alignment,
                                IgGeneRecombinationEventStoragePtr v_events);

    CleavedIgGeneAlignment GeneratePalindromicEvent(IgGeneAlignmentPtr v_alignment, size_t palindrome_length);

    void GeneratePalindromicEvents(IgGeneAlignmentPtr v_alignment,
                                   IgGeneRecombinationEventStoragePtr v_events);

public:
    VRecombinationEventGenerator(size_t max_cleavage, size_t max_palindrome) :
            max_cleavage_(max_cleavage),
            max_palindrome_(max_palindrome) { }

    IgGeneRecombinationEventStoragePtr ComputeEvents(IgGeneAlignmentPtr gene_segment_alignment);
};