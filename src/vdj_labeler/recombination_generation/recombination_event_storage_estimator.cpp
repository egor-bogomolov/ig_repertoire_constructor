#include "gene_event_storage_estimator.hpp"

void GeneEventRecombinationStorageEsimator::Update(IgGeneRecombinationEventStoragePtr storage_ptr,
                                                   IgGeneType gene_type) {
    size_t max_left_palindrome = 0;
    size_t max_right_pslindrome = 0;
    for(auto it = storage_ptr->cbegin(); it != storage_ptr->cend(); it++) {
        //max_left_palindrome = (*it)
    }
}