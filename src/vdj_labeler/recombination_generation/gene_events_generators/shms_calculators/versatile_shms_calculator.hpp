#include "shm_calculator.hpp"

class VersatileGeneSHMsCalculator : public SHMsCalculator {
    SHMsCalculator& left_shms_calculator_;
    SHMsCalculator& right_shms_calculator_;

public:
    VersatileGeneSHMsCalculator(SHMsCalculator &left_shms_calculator, SHMsCalculator &right_shms_calculator) :
            left_shms_calculator_(left_shms_calculator),
            right_shms_calculator_(right_shms_calculator) { }

    int ComputeNumberSHMs(IgGeneAlignmentPtr gene_alignment,
                          int left_cleavage_length,
                          int right_cleavage_length);
};