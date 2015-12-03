#pragma once

#include "../recombination/nongenomic_insertion.hpp"

class InsertionEventStorage {
    std::vector<NongenomicInsertion> insertions_;

public:
    void AddInsertion(NongenomicInsertion insertion) {
        insertions_.push_back(insertion);
    }

    typedef std::vector<NongenomicInsertion>::const_iterator insertion_event_iterator;

    insertion_event_iterator cbegin() const { return insertions_.cbegin(); }

    insertion_event_iterator cend() const { return insertions_.cend(); }

    size_t size() const { return insertions_.size(); }
};

typedef std::shared_ptr<InsertionEventStorage> InsertionEventStoragePtr;