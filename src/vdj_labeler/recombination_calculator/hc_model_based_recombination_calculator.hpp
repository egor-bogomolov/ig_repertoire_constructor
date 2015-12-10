#pragma once

#include "abstract_recombination_calculator.hpp"

class HCModelBasedRecombinationCalculator : public AbstractRecombinationCalculator<HCRecombination> {
    HCProbabilityRecombinationModel model_;

public:
    HCModelBasedRecombinationCalculator(const HCProbabilityRecombinationModel& model) :
            model_(model) { }

    double ComputeAssemblyProbability(const HCRecombination& recombination) const;
};
