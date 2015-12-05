#pragma once

#include "../recombination/recombination.hpp"

template<class Recombination>
class RecombinationStorage {
    ReadPtr read_ptr_;
    std::vector<Recombination> recombinations_;

    bool CheckConsistency(Recombination recombination) {
        assert(false);
        return false;
    }

public:
    RecombinationStorage(ReadPtr read_ptr) :
            read_ptr_(read_ptr) { }

    void AddRecombination(Recombination recombination) {
        if(CheckConsistency(recombination))
            recombinations_.push_back(recombination);
    }

    typedef typename std::vector<Recombination>::const_iterator recombination_iterator;

    recombination_iterator cbegin() const { return recombinations_.cbegin(); }

    recombination_iterator cend() const { return recombinations_.cend(); }

    size_t size() const { return recombinations_.size(); }

    ReadPtr Read() const { return read_ptr_; }
};

typedef RecombinationStorage<HCRecombination> HcRecombinationStorage;

typedef std::shared_ptr<HcRecombinationStorage> HcRecombinationStoragePtr;