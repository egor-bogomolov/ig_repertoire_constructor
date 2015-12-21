#include "gene_event_storage_estimator.hpp"

void GeneEventRecombinationStorageEsimator::UpdateVariableGeneStatistics(
        IgGeneRecombinationEventStoragePtr storage_ptr) {
    for(auto it = storage_ptr->cbegin(); it != storage_ptr->cend(); it++)
        if(it->RightEventIsPalindrome() and it->NumberRightSHMs() == 0)
            max_zero_v_palindromes_.push_back(it->RightCleavageLength() * -1);
}

void GeneEventRecombinationStorageEsimator::UpdateDiversityGeneStatistics(
        IgGeneRecombinationEventStoragePtr storage_ptr) {
    for(auto it = storage_ptr->cbegin(); it != storage_ptr->cend(); it++) {
        if(it->LeftEventIsPalindrome() and it->NumberLeftSHMs() == 0)
            max_zero_dl_palindromes_.push_back(it->LeftCleavageLength() * -1);
        if(it->RightEventIsPalindrome() and it->NumberRightSHMs() == 0)
            max_zero_dr_palindromes_.push_back(it->RightCleavageLength() * -1);
    }
}

void GeneEventRecombinationStorageEsimator::UpdateJoinGeneStatistics(
        IgGeneRecombinationEventStoragePtr storage_ptr) {
    for(auto it = storage_ptr->cbegin(); it != storage_ptr->cend(); it++)
        if(it->LeftEventIsPalindrome() and it->NumberLeftSHMs() == 0)
            max_zero_j_palindromes_.push_back(it->LeftCleavageLength() * -1);
}

void GeneEventRecombinationStorageEsimator::Update(IgGeneRecombinationEventStoragePtr storage_ptr,
                                                   IgGeneType gene_type) {
    if(gene_type == IgGeneType::variable_gene)
        UpdateVariableGeneStatistics(storage_ptr);
    if(gene_type == IgGeneType::diversity_gene)
        UpdateDiversityGeneStatistics(storage_ptr);
    if(gene_type == IgGeneType::join_gene)
        UpdateJoinGeneStatistics(storage_ptr);
}