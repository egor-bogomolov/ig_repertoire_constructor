#pragma once

#include <seqan/align.h>

#include "gene_database.hpp"
#include "../fastq_read_archive.hpp"

// in our context, query is a reads, subject is a gene segment
struct AlignmentPositions {
    std::pair<size_t, size_t> query_pos;
    std::pair<size_t, size_t> subject_pos;

    AlignmentPositions() :
        query_pos(std::make_pair(size_t(-1), size_t(-1))),
        subject_pos(std::make_pair(size_t(-1), size_t(-1))) { }

    AlignmentPositions(std::pair<size_t, size_t> new_query_pos,
              std::pair<size_t, size_t> new_subject_pos) :
            query_pos(new_query_pos),
            subject_pos(new_subject_pos) { }
};

std::ostream& operator<<(std::ostream &out, const AlignmentPositions& obj);

//-----------------------------------------------------------

struct IgGeneAlignmentPositions {
    AlignmentPositions alignment;
    IgGenePtr ig_gene;
    ReadPtr read;

    IgGeneAlignmentPositions() :
        alignment(),
        ig_gene(NULL),
        read(NULL) { }

    IgGeneAlignmentPositions(AlignmentPositions new_alignment,
                             IgGenePtr new_ig_gene,
                             ReadPtr new_read) :
            alignment(new_alignment),
            ig_gene(new_ig_gene),
            read(new_read) { }

    size_t GeneStartPos() const { return alignment.subject_pos.first; }

    size_t GeneEndPos() const { return alignment.subject_pos.second; }

    size_t ReadStartPos() const { return alignment.query_pos.first; }

    size_t ReadEndPos() const { return alignment.query_pos.second; }

    bool IsEmpty() const { return ReadStartPos() > ReadEndPos(); }

};

std::ostream& operator<<(std::ostream& out, const IgGeneAlignmentPositions& obj);

//-----------------------------------------------------------

class IgGeneAlignment {
    IgGeneAlignmentPositions positions_;
    seqan::Align<Dna5String> alignment_;
    int score_;
    double normalized_score_;
    size_t num_shms_;
    size_t num_gaps_;
    size_t alignment_length_;
    size_t read_alignment_length_;
    size_t gene_alignment_length_;

    void ComputeNormalizedScore();

    void ComputeSHMsNumber();

    void ComputeGapsNumber();

    void ComputeAlignmentLength();

public:
    IgGeneAlignment(IgGeneAlignmentPositions new_positions,
                    seqan::Align<Dna5String> new_alignment,
                    int score) :
            positions_(new_positions),
            alignment_(new_alignment),
            score_(score),
            normalized_score_(-1),
            num_shms_(size_t(-1)),
            num_gaps_(size_t(-1)) {
        ComputeNormalizedScore();
        ComputeSHMsNumber();
        ComputeGapsNumber();
        ComputeAlignmentLength();
    }

    IgGeneAlignmentPositions Positions() const { return positions_; }

    typedef seqan::Align<Dna5String> DnaAlignment;

    DnaAlignment& Alignment() { return alignment_; }

    const DnaAlignment& ConstAlignment() const { return alignment_; }

    IgGeneType GeneType() const { return positions_.ig_gene->GeneType(); }

    size_t GeneId() const { return positions_.ig_gene->id(); }

    size_t GeneLength() const { return positions_.ig_gene->length(); }

    void RefineAlignmentPositions(AlignmentPositions alignment_positions);

    size_t SHMsNumber() const { return num_shms_; }

    int Score() const { return score_; }

    double NormalizedScore() const { return normalized_score_; }

    size_t GapsNumber() const { return num_gaps_; }

    size_t AlignmentLength() const { return alignment_length_; }

    size_t ReadAlignmentLength() const { return read_alignment_length_; }

    size_t GeneAlignmentLength() const { return gene_alignment_length_; }

    const Dna5String& GeneSeq() const { return positions_.ig_gene->seq(); }

    const Dna5String& ReadSeq() const { return positions_.read->seq; }

    bool IsEmpty() const { return AlignmentLength() == 0; }
};

std::ostream& operator<<(std::ostream &out, const IgGeneAlignment& ig_gene_alignment);

typedef std::shared_ptr<IgGeneAlignment> IgGeneAlignmentPtr;

//----------------------------------------------------

seqan::Dna getRevCompl(seqan::Dna const & nucleotide);