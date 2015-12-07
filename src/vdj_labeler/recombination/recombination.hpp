#pragma once

#include "cleaved_gene.hpp"
#include "nongenomic_insertion.hpp"
#include "../fastq_read_archive.hpp"

class HCRecombination {
    ReadPtr read_ptr_;
    const CleavedIgGeneAlignment& v_gene_;
    const CleavedIgGeneAlignment& d_gene_;
    const CleavedIgGeneAlignment& j_gene_;
    NongenomicInsertion vd_insertion_;
    NongenomicInsertion dj_insertion_;

public:
    HCRecombination(ReadPtr read_ptr,
                    const CleavedIgGeneAlignment& v_gene,
                    const CleavedIgGeneAlignment& d_gene,
                    const CleavedIgGeneAlignment& j_gene,
                    NongenomicInsertion vd_insertion,
                    NongenomicInsertion dj_insertion) :
            read_ptr_(read_ptr),
            v_gene_(v_gene),
            d_gene_(d_gene),
            j_gene_(j_gene),
            vd_insertion_(vd_insertion),
            dj_insertion_(dj_insertion) { }

    ReadPtr Read() const { return read_ptr_; }

    const CleavedIgGeneAlignment& V() const { return v_gene_; }

    const CleavedIgGeneAlignment& D() const { return d_gene_; }

    const CleavedIgGeneAlignment& J() const { return j_gene_; }

    NongenomicInsertion VDInsertion() const { return vd_insertion_; }

    NongenomicInsertion DJInsertion() const { return dj_insertion_; }

    size_t SHMsNumber() const { return v_gene_.SHMsNumber() + d_gene_.SHMsNumber() + j_gene_.SHMsNumber(); }

    bool Valid() const;
};