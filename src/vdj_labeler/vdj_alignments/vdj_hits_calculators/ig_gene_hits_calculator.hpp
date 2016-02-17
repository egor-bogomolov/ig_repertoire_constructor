#pragma once

#include "../vdj_hits.hpp"

class IgGeneHitsCalculator {
public:
    virtual IgGeneSegmentHitsPtr ComputeHits(ReadPtr read_ptr) = 0;

    virtual ~IgGeneHitsCalculator() { }
};