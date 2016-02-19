#pragma once

#include "alignment_estimator.hpp"

class MatchThresholdAlignmentEstimator : public AlignmentEstimator {
    /*
     * The default value is 5 following IgBlast
     * link: https://nar.oxfordjournals.org/content/early/2013/05/11/nar.gkt382.full.pdf :
     * Identifying the V, D and J gene hits:
     * "The default word size is 5, which requires a minimum of five consecutive nucleotide
     * matches for a D gene to be found."
     */
    int normalized_score_threshold_ = 5;

public:
    MatchThresholdAlignmentEstimator(int normalized_score_threshold) :
            normalized_score_threshold_(normalized_score_threshold) { }

    bool AlignmentIsGood(IgGeneAlignmentPtr ig_gene_alignment);
};
