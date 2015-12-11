#pragma once

#include "../fastq_read_archive.hpp"

class NongenomicInsertion {
    size_t start_position_;
    size_t end_position_;

public:
    NongenomicInsertion() :
        start_position_(size_t(-1)),
        end_position_(size_t(-1)) { }

    NongenomicInsertion(size_t start_position,
                        size_t end_position) :
            start_position_(start_position),
            end_position_(end_position) { }

    size_t StartPosition() const { return start_position_; }

    size_t EndPosition() const { return end_position_; }

    size_t length() const {
        if(end_position_ < start_position_)
            return 0;
        return end_position_ - start_position_ + 1;
    }

    bool Valid() const { return start_position_ <= end_position_ or
                start_position_ - end_position_ <= 1; }
};