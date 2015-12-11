#include "versatile_insertion_event_generator.hpp"

InsertionEventStoragePtr VersatileInsertionGenerator::ComputeInsertionEvents(
        CleavedIgGeneAlignment left_gene_alignment, CleavedIgGeneAlignment right_gene_alignment) {
    InsertionEventStoragePtr insertions(new InsertionEventStorage());
    NongenomicInsertion insertion(left_gene_alignment.EndReadPosition() + 1,
            right_gene_alignment.StartReadPosition() - 1);
    insertions->AddInsertion(insertion);
    return insertions;
}