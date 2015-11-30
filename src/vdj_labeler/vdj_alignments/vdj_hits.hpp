#pragma once

#include "alignment_structs.hpp"
#include "../fastq_read_archive.hpp"

typedef std::vector<IgGeneAlignmentPtr>::const_iterator hits_citerator;

class IgGeneSegmentHits {
    IgGeneType gene_type_;
    ReadPtr read_ptr_;
    std::vector<IgGeneAlignmentPtr> hits_;

public:
    IgGeneSegmentHits(IgGeneType gene_type, ReadPtr read_ptr) :
            gene_type_(gene_type),
            read_ptr_(read_ptr) { }

    void AddHit(IgGeneAlignmentPtr hit);

    size_t size() const { return hits_.size(); }

    hits_citerator cbegin() const  { return hits_.cbegin(); }

    hits_citerator cend() const { return hits_.cend(); }

    IgGeneAlignmentPtr operator[](size_t index);

    IgGeneType GeneType() const { return gene_type_; }
};

typedef std::shared_ptr<IgGeneSegmentHits> IgGeneSegmentHitsPtr;

//------------------------------------------------------------

class VDJHits {
    ReadPtr read_ptr_;
    IgGeneSegmentHits v_hits_;
    IgGeneSegmentHits d_hits_;
    IgGeneSegmentHits j_hits_;

public:
    VDJHits(ReadPtr read_ptr) :
            read_ptr_(read_ptr),
            v_hits_(IgGeneType::variable_gene, read_ptr),
            d_hits_(IgGeneType::diversity_gene, read_ptr),
            j_hits_(IgGeneType::join_gene, read_ptr) { }

    //VDJHits(ReadPtr read_ptr, IgGeneSegmentHits v_hits, IgGeneSegmentHits j_hits, IgGeneSegmentHits d_hits) :
    //    read_ptr_(read_ptr),
    //    v_hits_(v_hits),
    //    d_hits_(d_hits),
    //    j_hits_(j_hits) { }

    void AddIgGeneAlignment(IgGeneType gene_type, IgGeneAlignmentPtr alignment_ptr);

    void AddIgGeneAlignment(IgGeneAlignmentPtr alignment_ptr);

    size_t VHitsNumber() const { return v_hits_.size(); }

    size_t DHitsNumber() const { return d_hits_.size(); }

    size_t JHitsNumber() const { return j_hits_.size(); }

    IgGeneAlignmentPtr GetAlignmentByIndex(IgGeneType gene_type, size_t index);
};

typedef std::shared_ptr<VDJHits> VDJHitsPtr;