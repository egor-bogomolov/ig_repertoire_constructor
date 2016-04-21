#pragma once
#include <utils/fastq_reader.hpp>
#include <utils/include_me.hpp>
#include <utils/sequence_tools.hpp>
#include <utils/string_tools.hpp>

#include "logger/log_writers.hpp"

#include <fstream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>

namespace reads_merger {

using namespace seqan;
using namespace std;

const size_t INFTY = size_t(-1);
const char worst_quality = char(33);

class Read {
public:
	CharString id, qual; 
	Dna5String seq;

	Read() = default;
	Read(const CharString &new_id, const CharString &new_qual, const Dna5String &new_seq) {
		id = new_id;
		qual = new_qual;
		seq = new_seq;
	}
	void read(SeqFileIn &in) {	
		readRecord(id, seq, qual, in);
	} 
	inline void reverseRead() {
		reverseComplement(seq);
		reverse(qual);
	}	
	void write(SeqFileOut &out) {
		writeRecord(out, id, seq, qual);
	}
};

struct merger_setting {
    size_t min_overlap;
    double max_mismatch_rate;
    bool simulated_mode;

    merger_setting() :
        min_overlap(50),
        max_mismatch_rate(.1),
        simulated_mode(false) { }
    void print() {
        INFO("Min overlap size: " << min_overlap);
        INFO("Max mismatch rate: " << max_mismatch_rate);
        if(simulated_mode) {
            INFO("Simulated mode is ON");
		}
    }
};

struct PairerInfo {
	size_t position;
	PairerInfo(): position(INFTY) {}
	PairerInfo(size_t pos): position(pos) {}
};

// вместо простого можно k-меры
class SimplePairer {
private:
	merger_setting settings;
	inline size_t hamming_distance(const Dna5String &s1, const Dna5String &s2, size_t pos) const {
		size_t result = 0;		
		for (size_t i = pos; i < length(s1) && i - pos < length(s2); ++i) {
			if (s1[i] != s2[i - pos]) {
				++result;
			}
		} 
		return result;
	} 
public:
	SimplePairer(const merger_setting &settings) : settings(settings) {}
	PairerInfo find_best_overlap(const Dna5String &left, const Dna5String &right) const {
		size_t best = INFTY, pos = INFTY; 
		size_t llen = length(left), rlen = length(right);
		for (size_t i = 0; i < llen; ++i) {
			size_t overlap = min(llen - i, rlen);
			if (overlap < settings.min_overlap) {
				break;
			}
			size_t dist = hamming_distance(left, right, i);
			if (double(dist) / double(overlap) <= settings.max_mismatch_rate && best > dist) {
				best = dist, pos = i;
			}
		}
		return PairerInfo(pos);
	}
};

class SequenceMerger {
private:
	inline CharString merge_ids(const CharString &id, size_t index) {
		string result;
		result.resize(length(id));
		for (size_t i = 0; i < length(id); ++i) {
			if (id[i] == ' ') {
				result[i] = '_';
			} else {
				result[i] = id[i];
			}
		}
		return CharString(to_string(index) + "_merged_read_" + result);
	}
	inline Dna5String merge_seqs(const Dna5String &lseq, const Dna5String &rseq, 
								 const CharString &lq, const CharString &rq, size_t pos) {
		Dna5String result = lseq;
		size_t llen = length(lseq), rlen = length(rseq);
		for (size_t i = pos; i < min(llen, pos + rlen); ++i) {
			if (result[i] == 'N') {
				result[i] = rseq[i - pos];
				continue;
			}
			if (lq[i] < rq[i - pos] && rq[i - pos] != 'N') {
				result[i] = rseq[i - pos];
			}
		}
		append(result, suffix(rseq, min(llen - pos, rlen)));
		return result;
	}
	inline CharString merge_quals(const Dna5String &lseq, const Dna5String &rseq, 
								  const CharString &lqual, const CharString &rqual, size_t pos) {
		CharString result = lqual;
		size_t llen = length(lseq), rlen = length(rseq);
		for (size_t i = pos; i < min(llen, pos + rlen); ++i) {
			if (lseq[i] == 'N' && rseq[i - pos] == 'N') {
				result[i] = worst_quality;
			} else if (lseq[i] == 'N') {
				result[i] = rqual[i - pos];
			} else if (rseq[i - pos] == 'N') {
				result[i] = lqual[i];
			} else if (lseq[i] == rseq[i - pos]) {
				result[i] = char(lqual[i] - worst_quality + rqual[i - pos]);
			} else if (lqual[i] < rqual[i - pos]) { 
				result[i] = rqual[i - pos]; 
			}
		}
		append(result, suffix(rqual, min(llen - pos, rlen)));
		return result;
	}
public:
	template <class Pairer>
	pair <Read, bool> Merge(const Read &left, const Read &right, const Pairer &pairer, size_t index) {
		PairerInfo info = pairer.find_best_overlap(left.seq, right.seq);
		if (info.position == INFTY) {
			return (make_pair(Read(), false));
		}
		CharString new_id = merge_ids(left.id, index);
		CharString new_qual = merge_quals(left.seq, right.seq, left.qual, right.qual, info.position);
		Dna5String new_seq = merge_seqs(left.seq, right.seq, left.qual, right.qual, info.position);
		assert(length(new_qual) == length(new_seq));
		return make_pair(Read(new_id, new_qual, new_seq), true);
	}
};

}
