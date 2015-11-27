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
    size_t num_shms_;

    void ComputeSHMsNumber();

public:
    IgGeneAlignment(IgGeneAlignmentPositions new_positions,
                    seqan::Align<Dna5String> new_alignment) :
            positions_(new_positions),
            alignment_(new_alignment),
            num_shms_(size_t(-1)) { }

    size_t SHMsNumber();

    IgGeneAlignmentPositions Positions() const { return positions_; }

    typedef seqan::Align<Dna5String> DnaAlignment;

    const DnaAlignment& Alignment() const { return alignment_; }
};

std::ostream& operator<<(std::ostream &out, const IgGeneAlignment& ig_gene_alignment);

typedef std::shared_ptr<IgGeneAlignment> IgGeneAlignmentPtr;