#pragma once
#include <utils/sequence_tools.hpp>

#include "logger/log_writers.hpp"

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>

namespace reads_merger {

	using namespace seqan;
	using namespace std;

	const size_t INFTY = size_t(-1);

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
			void reverseRead() {
				reverseComplement(seq);
				reverse(qual);
			}	
			void write(SeqFileOut &out) {
				writeRecord(out, id, seq, qual);
			}
	};

	struct merger_setting {
		size_t min_overlap;
		size_t base_quality;
		double max_mismatch_rate;

		merger_setting() = default;
		
		void print() {
			INFO("Min overlap size: " << min_overlap);
			INFO("Max mismatch rate: " << max_mismatch_rate);
			INFO("Base quality: " << base_quality);
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
			size_t hamming_distance(const Dna5String &s1, const Dna5String &s2, size_t pos) const {
				size_t result = 0;
				size_t length1 = length(s1), length2 = length(s2);		
				for (size_t i = pos; i < length1 && i - pos < length2; ++i) {
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
				size_t left_len = length(left), right_len = length(right);
				for (size_t i = 0; i < left_len; ++i) {
					size_t overlap = min(left_len - i, right_len);
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
			size_t quality() const {
				return settings.base_quality;
			}
	};
	
	class SuffArrayPairer {
		private:
			merger_setting settings;
			vector<int> col, num, p, p2, lcp; 
			vector<vector<int>> sparse_table;
			
			void build_suff_array(const Dna5String &s, int len) {
				int ma = max(len, 256);
				col.resize(len);
				p.resize(len);
				p2.resize(len);
				num.resize(len);
				for (int i = 0; i < len; ++i) {
					col[i] = s[i];
					p[i] = i;
				}	
				for (int k2 = 1; k2 / 2 < len; k2 *= 2) {
					int k = k2 / 2;
					num.assign(len, 0);
					for (int i = 0; i < len; ++i) {
						++num[col[i] + 1];
					}
					for (int i = 0; i < ma; ++i) {
						num[i + 1] += num[i];
					}
					for (int i = 0; i < len; ++i) {
						p2[num[col[(p[i] - k + len) % len]]++] = (p[i] - k + len) % len;
					}
					
					int cc = 0;
					for (int i = 0; i < len; ++i) {
						if (i && (col[p2[i]] != col[p2[i - 1]] || col[(p2[i] + k) % len] != col[(p2[i - 1] + k) % len])) {
							cc++;
						}
						num[p2[i]] = cc;
					}
					for (int i = 0; i < len; ++i) {
						p[i] = p2[i];
						col[i] = num[i];
					}
				}
				num.assign(len, 0);
				for (int i = 0; i < len; ++i) {
					++num[col[i] + 1];
				}
				for (int i = 0; i < ma; ++i) {
					num[i + 1] += num[i];
				}
				for (int i = 0; i < len; ++i) {
					p2[num[col[i]]] = i;
					++num[col[i]];
				}
				for (int i = 0; i < len; ++i) {
					p[i] = p2[i];
				}
				for (int i = 0; i < len; ++i) {
					p2[p[i]] = i;
				}
			}
			
			void build_lcp(const Dna5String &s, int len) {
				lcp.resize(len);
				int now = 0;
				for (int i = 0; i < len; ++i) {
					int j = p2[i];
					now = max(now - 1, 0);
					if (j != len - 1) {
						while (now < len && s[(p[j] + now) % len] == s[(p[j + 1] + now) % len]) {
							++now;
						}
					} 
					lcp[j] = now;
					if (j != len - 1 && p[j + 1] == len - 1) {
						now = 0;
					}
				}
			}
			/*
			void build_sparse(int len) {
				int log = 8;
				while
				sparse.assign(
			}*/
			
		public:
			SuffArrayPairer(const merger_setting &settings) : settings(settings) {}
			PairerInfo find_best_overlap(const Dna5String &left, const Dna5String &right) {
				size_t best = INFTY, pos = INFTY;
				size_t left_len = length(left), right_len = length(right);
				build_suff_array(left, int(left_len + 1 + right_len));
				build_lcp(left, int(left_len + 1 + right_len));
				return PairerInfo(pos);
			}
	};

	class SequenceMerger {
		private:
			CharString merge_ids(const CharString &id, size_t index) {
				string result;
				size_t length_id = length(id);
				result.resize(length_id);
				for (size_t i = 0; i < length_id; ++i) {
					if (id[i] == ' ') {
						result[i] = '_';
					} else {
						result[i] = id[i];
					}
				}
				return CharString(to_string(index) + "_merged_read_" + result);
			}
			Dna5String merge_seqs(const Dna5String &lseq, const Dna5String &rseq, 
					const CharString &lq, const CharString &rq, size_t pos) {
				Dna5String result = lseq;
				size_t left_len = length(lseq), right_len = length(rseq);
				for (size_t i = pos; i < min(left_len, pos + right_len); ++i) {
					if (result[i] == 'N') {
						result[i] = rseq[i - pos];
						continue;
					}
					if (lq[i] < rq[i - pos] && rq[i - pos] != 'N') {
						result[i] = rseq[i - pos];
					}
				}
				append(result, suffix(rseq, min(left_len - pos, right_len)));
				return result;
			}
			CharString merge_quals(const Dna5String &lseq, const Dna5String &rseq, 
					const CharString &lqual, const CharString &rqual, size_t pos, size_t base_quality) {
				CharString result = lqual;
				size_t left_len = length(lseq), right_len = length(rseq);
				for (size_t i = pos; i < min(left_len, pos + right_len); ++i) {
					if (lseq[i] == 'N' && rseq[i - pos] == 'N') {
						result[i] = char(base_quality);
					} else if (lseq[i] == 'N') {
						result[i] = rqual[i - pos];
					} else if (rseq[i - pos] == 'N') {
						result[i] = lqual[i];
					} else if (lseq[i] == rseq[i - pos]) {
						result[i] = char(lqual[i] - base_quality + rqual[i - pos]);
					} else if (lqual[i] < rqual[i - pos]) { 
						result[i] = rqual[i - pos]; 
					}
				}
				append(result, suffix(rqual, min(left_len - pos, right_len)));
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
					CharString new_qual = merge_quals(left.seq, right.seq, left.qual, right.qual, info.position, pairer.quality());
					Dna5String new_seq = merge_seqs(left.seq, right.seq, left.qual, right.qual, info.position);
					assert(length(new_qual) == length(new_seq));
					return make_pair(Read(new_id, new_qual, new_seq), true);
				}
	};
}
