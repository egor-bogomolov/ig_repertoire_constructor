#pragma once

#include "../vdj_alignments/vdj_hits.hpp"

#include "recombination_generator.hpp"
#include "../recombination/recombination.hpp"

typedef std::vector<HCRecombination>::iterator HCRecombinationIterator;

class BaseHCRecombinationGenerator : public AbstractRecombinationGenerator<HCRecombination, HCRecombinationIterator> {
    VDJHitsPtr vdj_hits_;
    std::vector<HCRecombination> recombinations_;

public:
    BaseHCRecombinationGenerator(VDJHitsPtr vdj_hits) : vdj_hits_(vdj_hits) { }

    void ComputeRecombinations();

    HCRecombinationIterator begin() { return recombinations_.begin(); }

    HCRecombinationIterator end() { return recombinations_.end(); }

    // std::size_t size() const { return recombinations_.size(); }
};