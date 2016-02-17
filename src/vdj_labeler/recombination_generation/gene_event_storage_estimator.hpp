#pragma once

#include "recombination_event_storage.hpp"

class GeneEventRecombinationStorageEsimator {
    std::vector<size_t> max_zero_v_palindromes_;
    std::vector<size_t> max_zero_dl_palindromes_;
    std::vector<size_t> max_zero_dr_palindromes_;
    std::vector<size_t> max_zero_j_palindromes_;

    void UpdateVariableGeneStatistics(IgGeneRecombinationEventStoragePtr storage_ptr);

    void UpdateDiversityGeneStatistics(IgGeneRecombinationEventStoragePtr storage_ptr);

    void UpdateJoinGeneStatistics(IgGeneRecombinationEventStoragePtr storage_ptr);

public:
    void Update(IgGeneRecombinationEventStoragePtr storage_ptr, IgGeneType gene_type);
};