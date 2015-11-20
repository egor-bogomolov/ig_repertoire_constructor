#pragma once

#include "../vdj_alignments/vdj_hits.hpp"

#include "recombination_generator.hpp"
#include "../recombination/recombination.hpp"
#include "recombination_storage.hpp"
#include "gene_events_generators/ig_gene_recombination_event_generator.hpp"
#include "insertion_events_generators/insertion_event_generator.hpp"

class CustomHeavyChainRecombinationGenerator :
        public AbstractRecombinationGenerator<HCRecombination, HcRecombinationStoragePtr> {
    IgGeneRecombinationEventsGenerator &v_events_generator_;
    IgGeneRecombinationEventsGenerator &d_events_generator_;
    IgGeneRecombinationEventsGenerator &j_events_generator_;
    InsertionEventGenerator& vd_insertion_generator_;
    InsertionEventGenerator& dj_insertion_generator_;

    // inner structures
    std::vector<IgGeneRecombinationEventStoragePtr> v_storages_;
    std::vector<IgGeneRecombinationEventStoragePtr> d_storages_;
    std::vector<IgGeneRecombinationEventStoragePtr> j_storages_;

    void Clear();

    void ComputeVEventStorages(VDJHitsPtr vdj_hits);

    void ComputeDEventStorages(VDJHitsPtr vdj_hits);

    void ComputeJEventStorages(VDJHitsPtr vdj_hits);

    // if d alignment is empty, we will assign length from V end to J start to the first insertion
    // and zero length to the second insertion
    // todo: extract to separate class
    std::pair<NongenomicInsertion, NongenomicInsertion>
            RefineNongenomicInsertions(NongenomicInsertion vd_insertion, NongenomicInsertion dj_insertion);

    HcRecombinationStoragePtr CreateRecombinations(HcRecombinationStoragePtr recombination_storage,
                                                   CleavedIgGeneAlignment v_gene,
                                                   CleavedIgGeneAlignment d_gene,
                                                   CleavedIgGeneAlignment j_gene,
                                                   InsertionEventStoragePtr vd_insertions,
                                                   InsertionEventStoragePtr dj_insertions);

    HcRecombinationStoragePtr CreateRecombinations(HcRecombinationStoragePtr recombination_storage,
                                                   IgGeneRecombinationEventStoragePtr v_events,
                                                   IgGeneRecombinationEventStoragePtr d_events,
                                                   IgGeneRecombinationEventStoragePtr j_events);

public:
    CustomHeavyChainRecombinationGenerator(IgGeneRecombinationEventsGenerator &v_events_generator,
                                           IgGeneRecombinationEventsGenerator &d_events_generator,
                                           IgGeneRecombinationEventsGenerator &j_events_generator,
                                           InsertionEventGenerator& vd_insertion_generator,
                                           InsertionEventGenerator& dj_insertion_generator) :
            v_events_generator_(v_events_generator),
            d_events_generator_(d_events_generator),
            j_events_generator_(j_events_generator),
            vd_insertion_generator_(vd_insertion_generator),
            dj_insertion_generator_(dj_insertion_generator) { }

    HcRecombinationStoragePtr ComputeRecombinations(VDJHitsPtr vdj_hits);
};