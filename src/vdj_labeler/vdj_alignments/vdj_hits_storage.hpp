#pragma once

#include "vdj_hits.hpp"

class VDJHitsStorage {
    std::vector<VDJHitsPtr> vdj_hits_;

public:
    void AddVDJHits(VDJHitsPtr vdj_hits_ptr) {
        vdj_hits_.push_back(vdj_hits_ptr);
    }

    size_t size() const { return vdj_hits_.size(); }

    typedef std::vector<VDJHitsPtr>::const_iterator vdj_hits_citerator;

    vdj_hits_citerator cbegin() const { return vdj_hits_.cbegin(); }

    vdj_hits_citerator cend() const { return vdj_hits_.cend(); }
};

typedef std::shared_ptr<VDJHitsStorage> VDJHitsStoragePtr;