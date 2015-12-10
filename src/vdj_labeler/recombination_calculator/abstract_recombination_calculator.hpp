#pragma once

#include "../model/recombination_model.hpp"
#include "../recombination/recombination.hpp"

template<class Recombination>
class AbstractRecombinationCalculator {
public:
    virtual double ComputeAssemblyProbability(const Recombination& recombination) const = 0;
    virtual ~AbstractRecombinationCalculator() { }
};
