#include "logger/logger.hpp"
#include "launch.hpp"

void VDJLabeler::TestRecombinationCalculator(const FastqReadArchive &reads_archive,
                                             VDJHitsStoragePtr &hits_storage) {
    size_t read_index = 3;
    ReadPtr read_3 = reads_archive[read_index];
    VDJHitsPtr hits_3 = (*hits_storage)[read_index];
    INFO("Read 3. #V: " << hits_3->VHitsNumber() <<
         ", #D: " << hits_3->DHitsNumber() <<
         ", #J: " << hits_3->JHitsNumber());

    auto v_alignment = hits_3->GetAlignmentByIndex(IgGeneType::variable_gene, 0);
    CleavedIgGeneAlignment v_event_0(v_alignment, 0, 0, 0, 0);
    CleavedIgGeneAlignment v_event_1(v_alignment, 0, -1, 0, 0);
    CleavedIgGeneAlignment v_event_2(v_alignment, 0, -2, 0, 0);
    CleavedIgGeneAlignment v_event_3(v_alignment, 0, -3, 0, 1);

    auto d_alignment = hits_3->GetAlignmentByIndex(IgGeneType::diversity_gene, 0);
    CleavedIgGeneAlignment d_event_0(d_alignment, 1, 8, 0, 0);

    auto j_alignment = hits_3->GetAlignmentByIndex(IgGeneType::join_gene, 0);
    CleavedIgGeneAlignment j_event_0(j_alignment, 0, 0, 1, 0);
    CleavedIgGeneAlignment j_event_1(j_alignment, 1, 0, 0, 0);

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

void VDJLabeler::Run(const vdj_config::io_params &io, const vdj_config::run_params &rp) {
    INFO("VDJ labeler starts");

    std::string fastq_reads_fname = io.input.input_sequences;
    std::string vj_alignment_fname = io.input.vj_alignment_file;
    std::string v_germline_genes_fname = io.input.germlines.igh.variable_genes;
    std::string d_germline_genes_fname = io.input.germlines.igh.diversity_genes;
    std::string j_germline_genes_fname = io.input.germlines.igh.join_genes;

    FastqReadArchive reads_archive(fastq_reads_fname);
    INFO(reads_archive.size() << " reads were extracted from " << fastq_reads_fname);

    HC_GenesDatabase hc_db;
    hc_db.AddGenesFromFile(IgGeneType::variable_gene, v_germline_genes_fname);
    hc_db.AddGenesFromFile(IgGeneType::diversity_gene, d_germline_genes_fname);
    hc_db.AddGenesFromFile(IgGeneType::join_gene, j_germline_genes_fname);

    INFO(hc_db.VariableGenes().size() << " variable genes were extracted from " << v_germline_genes_fname);
    INFO(hc_db.DiversityGenes().size() << " diversity genes were extracted from " << d_germline_genes_fname);
    INFO(hc_db.JoinGenes().size() << " join genes were extracted from " << j_germline_genes_fname);

    // VJAlignmentInfo vj_alignment_info(hc_db.VariableGenes(), hc_db.JoinGenes(), reads_archive);
    // vj_alignment_info.ExtractAlignment(vj_alignment_fname);
    // INFO(vj_alignment_info.size() << " alignment lines were extracted from " << vj_alignment_fname);

    // INFO("Best VDJ hits alignment calculation starts");
    // RightVTailAligner v_aligner;
    // InfoBasedVJHitsCalculator v_hits_calc(IgGeneType::variable_gene, reads_archive, vj_alignment_info, v_aligner);
    // SimpleDAligner d_aligner;
    // ThresholdAlignmentEstimator d_estimator(1.0);
    // InfoBasedDHitsCalculator d_hits_calc(reads_archive,
    //                                      vj_alignment_info,
    //                                      hc_db.DiversityGenes(),
    //                                      d_aligner, d_estimator);
    // LeftJTailAligner j_aligner;
    // InfoBasedVJHitsCalculator j_hits_calc(IgGeneType::join_gene, reads_archive, vj_alignment_info, j_aligner);
    // CustomVDJHitsCalculator vdj_hits_calc(reads_archive, vj_alignment_info, v_hits_calc, d_hits_calc, j_hits_calc);
    // auto hits_storage = vdj_hits_calc.ComputeHits();
    // INFO("Best VDJ hits alignment calculation ends");

    // size_t max_cleavage = 20;
    // size_t max_palindrome = 7;
    // LeftEventSHMsCalculator left_shms_calculator;
    // RightEventSHMsCalculator right_shms_calculator;
    // VersatileGeneSHMsCalculator shms_calculator(left_shms_calculator, right_shms_calculator);
    // VRecombinationEventGenerator v_generator(shms_calculator, max_cleavage, max_palindrome);
    // DRecombinationEventGenerator d_generator(shms_calculator, max_cleavage, max_palindrome);
    // JRecombinationEventGenerator j_generator(shms_calculator, max_cleavage, max_palindrome);
    // VersatileInsertionGenerator insertion_generator;
    // CustomHeavyChainRecombinationGenerator recombination_generator(v_generator,
    //                                                                d_generator,
    //                                                                j_generator,
    //                                                                insertion_generator,
    //                                                                insertion_generator);
    // HcRecombinationEstimator recombination_estimator;
    // INFO("Generator of VDJ recombinations starts");
    // for(auto it = hits_storage->cbegin(); it != hits_storage->cend(); it++) {
    //     auto recombination_storage = recombination_generator.ComputeRecombinations(*it);
    //     recombination_estimator.Update(recombination_storage);
    // }
    // INFO("Generator of VDJ recombinations ends");
    // recombination_estimator.OutputRecombinationNumber();
    // recombination_estimator.OutputSHMsDistribution();
    // recombination_estimator.OutputRecombinationEvents();

    {
        std::ifstream in("src/vdj_labeler/test/blank_model.csv");
        // HCProbabilityRecombinationModel model(in, hc_db);
        // cout << model;
        IgGeneProbabilityModel model_V(in, hc_db.VariableGenes());
        //cout << model_V;
        IgGeneProbabilityModel model_D(in, hc_db.DiversityGenes());
        // cout << model_D;
        IgGeneProbabilityModel model_J(in, hc_db.JoinGenes());
        // cout << model_J;
        // NongenomicInsertionModel modelVD(in);
        // NongenomicInsertionModel modelDJ(in);
        // cout << modelVD;
        // cout << modelDJ;
        // PalindromeDeletionModel modelDelV(in, hc_db.VariableGenes());
        // cout << modelDelV;
        // PalindromeDeletionModel modelDelJ(in, hc_db.JoinGenes());
        // PalindromeDeletionModel modelDelDL(in, hc_db.DiversityGenes());
        // PalindromeDeletionModel modelDelDR(in, hc_db.DiversityGenes());

        // HCModelBasedRecombinationCalculator recombination_calculator(model);
    }
    INFO("VDJ labeler ends");
}
