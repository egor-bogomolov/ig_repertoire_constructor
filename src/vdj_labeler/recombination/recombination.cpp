#include "recombination.hpp"

bool HCRecombination::Valid() const {
    return vd_insertion_.Valid() and dj_insertion_.Valid();
}