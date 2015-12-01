#include "base_hc_recombination_generator.hpp"

using namespace std;

void BaseHCRecombinationGenerator::ComputeRecombinations() {
    size_t num_v_hits = vdj_hits_->VHitsNumber();
    size_t num_d_hits = vdj_hits_->DHitsNumber();
    size_t num_j_hits = vdj_hits_->JHitsNumber();
    for(size_t vi = 0; vi < num_v_hits; vi++)
        for(size_t di = 0; di < num_d_hits; di++)
            for(size_t ji = 0; ji < num_j_hits; ji++) {
                auto v_alignment = vdj_hits_->GetAlignmentByIndex(IgGeneType::variable_gene, vi);
                auto d_alignment = vdj_hits_->GetAlignmentByIndex(IgGeneType::diversity_gene, di);
                auto j_alignment = vdj_hits_->GetAlignmentByIndex(IgGeneType::join_gene, ji);
                cout << "V. " << *v_alignment << endl;
                cout << "D. " << *d_alignment << endl;
                cout << "J. " << *j_alignment << endl;
                cout << "-----------------------------" << endl;
            }
}