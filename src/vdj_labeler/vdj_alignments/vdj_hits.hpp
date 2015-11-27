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

    hits_citerator cend() const { return hits_.end(); }

    IgGeneAlignmentPtr operator[](size_t index);
};

typedef std::shared_ptr<IgGeneSegmentHits> IgGeneSegmentHitsPtr;

//------------------------------------------------------------

class VDJHits {
    ReadPtr read_ptr_;
    IgGeneSegmentHits v_gene_alignments_;
    IgGeneSegmentHits d_gene_alignments_;
    IgGeneSegmentHits j_gene_alignments_;

public:
    VDJHits(ReadPtr read_ptr) :
            read_ptr_(read_ptr),
            v_gene_alignments_(IgGeneType::variable_gene, read_ptr),
            d_gene_alignments_(IgGeneType::diversity_gene, read_ptr),
            j_gene_alignments_(IgGeneType::join_gene, read_ptr) { }

    void AddIgGeneAlignment(IgGeneType gene_type, IgGeneAlignmentPtr alignment_ptr);

    size_t VHitsNumber() const { return v_gene_alignments_.size(); }

    size_t DHitsNumber() const { return d_gene_alignments_.size(); }

    size_t JHitsNumber() const { return j_gene_alignments_.size(); }

    IgGeneAlignmentPtr GetAlignmentByIndex(IgGeneType gene_type, size_t index);
};

typedef std::shared_ptr<VDJHits> VDJHitsPtr;