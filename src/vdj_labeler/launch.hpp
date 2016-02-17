#pragma once

#include "vdj_config.hpp"

#include "fastq_read_archive.hpp"
#include "vdj_alignments/gene_database.hpp"
#include "vdj_alignments/vj_alignment_info.hpp"

#include "vdj_alignments/aligners/right_v_tail_aligner.hpp"
#include "vdj_alignments/aligners/left_j_tail_aligner.hpp"
#include "vdj_alignments/aligners/simple_d_aligner.hpp"

#include "vdj_alignments/vdj_hits_calculators/custom_vdj_hits_calculator.hpp"
#include "vdj_alignments/vdj_hits_calculators/info_based_vj_hits_calculator.hpp"
#include "vdj_alignments/vdj_hits_calculators/info_based_d_hits_calculator.hpp"
#include "vdj_alignments/vdj_hits_calculators/alignment_estimators/threshold_alignment_estimator.hpp"

#include "recombination_generation/gene_events_generators/shms_calculators/left_event_shms_calculator.hpp"
#include "recombination_generation/gene_events_generators/shms_calculators/right_event_shms_calculator.hpp"
#include "recombination_generation/gene_events_generators/shms_calculators/versatile_shms_calculator.hpp"
#include "recombination_generation/custom_hc_recombination_generator.hpp"
#include "recombination_generation/gene_events_generators/v_recombination_event_generator.hpp"
#include "recombination_generation/gene_events_generators/d_recombination_event_generator.hpp"
#include "recombination_generation/gene_events_generators/j_recombination_event_generator.hpp"
#include "recombination_generation/insertion_events_generators/versatile_insertion_event_generator.hpp"

#include "recombination_estimators/hc_recombination_estimator.hpp"
#include "recombination_calculator/hc_model_based_recombination_calculator.hpp"

class VDJLabeler {
    void TestRecombinationCalculator(const FastqReadArchive& reads_archive,
                                     VDJHitsStoragePtr &hits_storage);

public:
    void Run(const vdj_config::io_params &io,
             const vdj_config::run_params &rp);
};