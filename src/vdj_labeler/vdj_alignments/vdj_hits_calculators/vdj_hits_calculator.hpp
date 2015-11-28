#pragma once

#include "../vdj_hits_storage.hpp"
# include "../../fastq_read_archive.hpp"

class VDJHitsCalculator {
protected:
    const FastqReadArchive& read_archive_;

public:
    VDJHitsCalculator(const FastqReadArchive& read_archive) :
            read_archive_(read_archive) { }

    virtual VDJHitsStoragePtr ComputeHits() = 0;

    virtual ~VDJHitsCalculator() { }
};