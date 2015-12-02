#pragma once

#include "../vdj_alignments/vdj_hits.hpp"

#include "recombination_generator.hpp"
#include "../recombination/recombination.hpp"
#include "recombination_storage.hpp"
#include "gene_events_generators/ig_gene_recombination_event_generator.hpp"

class CustomHeavyChainRecombinationGenerator :
        public AbstractRecombinationGenerator<HCRecombination, HcRecombinationStoragePtr> {
    IgGeneRecombinationEventsGenerator &v_events_generator_;
    IgGeneRecombinationEventsGenerator &d_events_generator_;
    IgGeneRecombinationEventsGenerator &j_events_generator_;

public:
    CustomHeavyChainRecombinationGenerator(IgGeneRecombinationEventsGenerator &v_events_generator,
                                           IgGeneRecombinationEventsGenerator &d_events_generator,
                                           IgGeneRecombinationEventsGenerator &j_events_generator) :
            v_events_generator_(v_events_generator),
            d_events_generator_(d_events_generator),
            j_events_generator_(j_events_generator) { }

    HcRecombinationStoragePtr ComputeRecombinations(VDJHitsPtr vdj_hits);
};