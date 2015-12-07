#include "standard.hpp"
#include "logger/log_writers.hpp"

#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"
#include "segfault_handler.hpp"

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

#include "recombination_generation/custom_hc_recombination_generator.hpp"
#include "recombination_generation/gene_events_generators/v_recombination_event_generator.hpp"
#include "recombination_generation/gene_events_generators/d_recombination_event_generator.hpp"
#include "recombination_generation/gene_events_generators/j_recombination_event_generator.hpp"
#include "recombination_generation/insertion_events_generators/versatile_insertion_event_generator.hpp"

#include "recombination_calculator/hc_model_based_recombination_calculator.hpp"

void create_console_logger() {
    using namespace logging;
    string log_props_file = "";
    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

void TestRecombinationCalculator(const FastqReadArchive& reads_archive, VDJHitsStoragePtr hits_storage) {
    size_t read_index = 3;
    ReadPtr read_3 = reads_archive[read_index];
    VDJHitsPtr hits_3 = (*hits_storage)[read_index];
    INFO("Read 3. #V: " << hits_3->VHitsNumber() <<
         ", #D: " << hits_3->DHitsNumber() <<
         ", #J: " << hits_3->JHitsNumber());

    auto v_alignment = hits_3->GetAlignmentByIndex(IgGeneType::variable_gene, 0);
    CleavedIgGeneAlignment v_event_0(v_alignment, 0, 0, 0);
    CleavedIgGeneAlignment v_event_1(v_alignment, 0, -1, 0);
    CleavedIgGeneAlignment v_event_2(v_alignment, 0, -2, 0);
    CleavedIgGeneAlignment v_event_3(v_alignment, 0, -3, 1);

    auto d_alignment = hits_3->GetAlignmentByIndex(IgGeneType::diversity_gene, 0);
    CleavedIgGeneAlignment d_event_0(d_alignment, 1, 8, 0);

    auto j_alignment = hits_3->GetAlignmentByIndex(IgGeneType::join_gene, 0);
    CleavedIgGeneAlignment j_event_0(j_alignment, 0, 0, 1);
    CleavedIgGeneAlignment j_event_1(j_alignment, 1, 0, 0);

    NongenomicInsertion vd_insertion_0(425, 441);
    NongenomicInsertion vd_insertion_1(426, 441);
    NongenomicInsertion vd_insertion_2(427, 441);
    NongenomicInsertion vd_insertion_3(428, 441);

    NongenomicInsertion dj_insertion_0(453, 452);
    NongenomicInsertion dj_insertion_1(453, 453);

    RecombinationStorage<HCRecombination> recombination_storage(read_3);
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_0, d_event_0, j_event_0,
                                                           vd_insertion_0, dj_insertion_0));
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_1, d_event_0, j_event_0,
                                                           vd_insertion_1, dj_insertion_0));
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_2, d_event_0, j_event_0,
                                                           vd_insertion_2, dj_insertion_0));
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_3, d_event_0, j_event_0,
                                                           vd_insertion_3, dj_insertion_0));
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_0, d_event_0, j_event_1,
                                                           vd_insertion_0, dj_insertion_1));
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_1, d_event_0, j_event_1,
                                                           vd_insertion_1, dj_insertion_1));
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_2, d_event_0, j_event_1,
                                                           vd_insertion_2, dj_insertion_1));
    recombination_storage.AddRecombination(HCRecombination(read_3, v_event_3, d_event_0, j_event_1,
                                                           vd_insertion_3, dj_insertion_1));
    INFO(recombination_storage.size() << " recombinaions were generated");
}

int main(int, char**) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();

    std::string fastq_reads_fname = "src/vdj_labeler/test/vdj_labeling.fastq";
    std::string vj_alignment_fname = "src/vdj_labeler/test/vdj_labeling.csv";
    std::string v_germline_genes_fname = "src/fast_ig_tools/germline/human/IGHV.fa";
    std::string d_germline_genes_fname = "src/fast_ig_tools/germline/human/IGHD.fa";
    std::string j_germline_genes_fname = "src/fast_ig_tools/germline/human/IGHJ.fa";

    INFO("VDJ labeler starts");

    FastqReadArchive reads_archive(fastq_reads_fname);
    INFO(reads_archive.size() << " reads were extracted from " << fastq_reads_fname);

    HC_GenesDatabase hc_db;
    hc_db.AddGenesFromFile(IgGeneType::variable_gene, v_germline_genes_fname);
    hc_db.AddGenesFromFile(IgGeneType::diversity_gene, d_germline_genes_fname);
    hc_db.AddGenesFromFile(IgGeneType::join_gene, j_germline_genes_fname);

    INFO(hc_db.VariableGenes().size() << " variable genes were extracted from " << v_germline_genes_fname);
    INFO(hc_db.DiversityGenes().size() << " diversity genes were extracted from " << d_germline_genes_fname);
    INFO(hc_db.JoinGenes().size() << " join genes were extracted from " << j_germline_genes_fname);

    VJAlignmentInfo vj_alignment_info(hc_db.VariableGenes(), hc_db.JoinGenes(), reads_archive);
    vj_alignment_info.ExtractAlignment(vj_alignment_fname);
    INFO(vj_alignment_info.size() << " alignment lines were extracted from " << vj_alignment_fname);

    INFO("Best VDJ hits alignment calculation starts");
    RightVTailAligner v_aligner;
    InfoBasedVJHitsCalculator v_hits_calc(IgGeneType::variable_gene, reads_archive, vj_alignment_info, v_aligner);
    SimpleDAligner d_aligner;
    ThresholdAlignmentEstimator d_estimator(1.0);
    InfoBasedDHitsCalculator d_hits_calc(reads_archive,
                                         vj_alignment_info,
                                         hc_db.DiversityGenes(),
                                         d_aligner, d_estimator);
    LeftJTailAligner j_aligner;
    InfoBasedVJHitsCalculator j_hits_calc(IgGeneType::join_gene, reads_archive, vj_alignment_info, j_aligner);
    CustomVDJHitsCalculator vdj_hits_calc(reads_archive, v_hits_calc, d_hits_calc, j_hits_calc);
    auto hits_storage = vdj_hits_calc.ComputeHits();
    INFO("Best VDJ hits alignment calculation ends");

    size_t max_cleavage = 6;
    size_t max_palindrome = 7;
    VRecombinationEventGenerator v_generator(max_cleavage, max_palindrome);
    DRecombinationEventGenerator d_generator;
    JRecombinationEventGenerator j_generator;
    VersatileInsertionGenerator insertion_generator;
    CustomHeavyChainRecombinationGenerator recombination_generator(v_generator,
                                                                   d_generator,
                                                                   j_generator,
                                                                   insertion_generator,
                                                                   insertion_generator);
    INFO("Generator of VDJ recombinations starts");
    for(auto it = hits_storage->cbegin(); it != hits_storage->cend(); it++) {
        cout << "Read " << (*it)->Read()->name << endl;
        auto recombination_storage = recombination_generator.ComputeRecombinations(*it);
        cout << endl << endl;
    }
    INFO("Generator of VDJ recombinations ends");

    // test recombinations for Andrey
    INFO("Test recombination for Andrey");
    TestRecombinationCalculator(reads_archive, hits_storage);

    INFO("VDJ labeler ends");
    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    INFO("Running time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    return 0;
}