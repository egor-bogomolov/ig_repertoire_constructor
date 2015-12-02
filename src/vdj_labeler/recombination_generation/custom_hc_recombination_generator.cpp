#include "custom_hc_recombination_generator.hpp"

using namespace std;

HcRecombinationStoragePtr CustomHeavyChainRecombinationGenerator::ComputeRecombinations(VDJHitsPtr vdj_hits) {
    size_t num_v_hits = vdj_hits->VHitsNumber();
    size_t num_d_hits = vdj_hits->DHitsNumber();
    size_t num_j_hits = vdj_hits->JHitsNumber();
    for(size_t vi = 0; vi < num_v_hits; vi++) {
        auto v_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::variable_gene, vi);
        for(size_t di = 0; di < num_d_hits; di++) {
            auto d_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::diversity_gene, di);
            for(size_t ji = 0; ji < num_j_hits; ji++) {
                auto j_alignment = vdj_hits->GetAlignmentByIndex(IgGeneType::join_gene, ji);
                cout << "V. " << *v_alignment << endl;
                cout << "D. " << *d_alignment << endl;
                cout << "J. " << *j_alignment << endl;
                cout << "-----------------------------" << endl;
            }
        }
    }
}