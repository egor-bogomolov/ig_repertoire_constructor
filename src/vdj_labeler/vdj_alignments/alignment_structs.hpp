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
};

std::ostream& operator<<(std::ostream& out, const IgGeneAlignmentPositions& obj);

//-----------------------------------------------------------

class IgGeneAlignment {
    IgGeneAlignmentPositions positions_;
    seqan::Align<Dna5String> alignment_;
    int score_;
    size_t num_shms_;
    double normalized_score_;

    void ComputeSHMsNumber();

    void ComputeNormalizedScore();

public:
    IgGeneAlignment(IgGeneAlignmentPositions new_positions,
                    seqan::Align<Dna5String> new_alignment,
                    int score) :
            positions_(new_positions),
            alignment_(new_alignment),
            score_(score),
            num_shms_(size_t(-1)) {
        ComputeNormalizedScore();
        ComputeSHMsNumber();
    }

    IgGeneAlignmentPositions Positions() const { return positions_; }

    typedef seqan::Align<Dna5String> DnaAlignment;

    DnaAlignment& Alignment() { return alignment_; }

    const DnaAlignment& ConstAlignment() const { return alignment_; }

    IgGeneType GeneType() const { return positions_.ig_gene->GeneType(); }

    void RefineAlignmentPositions(AlignmentPositions alignment_positions);

    size_t SHMsNumber() const { return num_shms_; }

    int Score() const { return score_; }

    double NormalizedScore() const { return normalized_score_; }
};

std::ostream& operator<<(std::ostream &out, const IgGeneAlignment& ig_gene_alignment);

typedef std::shared_ptr<IgGeneAlignment> IgGeneAlignmentPtr;