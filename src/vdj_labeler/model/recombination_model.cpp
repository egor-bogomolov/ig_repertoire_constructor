#include "recombination_model.hpp"

#include <cassert>

#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include <utility>

#include <boost/tokenizer.hpp>
#include <seqan/sequence.h>

using boost::tokenizer;
using boost::escaped_list_separator;
using Tokenizer = tokenizer<escaped_list_separator<char>>;

IgGeneProbabilityModel::IgGeneProbabilityModel(const IgGeneProbabilityVector& ig_gene_probabilities,
                                               const IgGeneDatabasePtrConst& ig_gene_database) :
        ig_gene_probabilities_(ig_gene_probabilities),
        ig_gene_database_(ig_gene_database) { }

size_t IgGeneProbabilityModel::size() const { return ig_gene_probabilities_.size(); }

IgGeneProbabilityModel::IgGeneProbabilityModel(std::ifstream& in,
                                               const IgGeneDatabasePtrConst& ig_gene_database) :
        ig_gene_database_(ig_gene_database) {
    assert(in.is_open());
    string line;
    vector<string> parsed_vector;
    ig_gene_probabilities_.resize(ig_gene_database_ -> size());
    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == 2);
        size_t index_of_current_gene = ig_gene_database -> GetIndexByName(parsed_vector.front());
        ig_gene_probabilities_[index_of_current_gene] = std::stod(parsed_vector.back());
    }
}

IgGeneProbabilityModel::IgGeneProbabilityModel(std::ifstream& in, const IgGeneDatabase& database) :
    IgGeneProbabilityModel(in, make_shared<const IgGeneDatabase>(database)) { }

std::ostream& operator<<(std::ostream& out, const IgGeneProbabilityModel& model) {
    for (size_t i = 0; i < model.size(); ++i)
        out << "Gene_id: " << model.GetIgGeneDatabase() -> GetByIndex(i) -> name()
            << ", " << "Gene probability: " << model.GetIgGeneProbabilities().at(i) << "\n";
    return out;
}


/************************************************************************************************/

NongenomicInsertionModel::NongenomicInsertionModel(
            vector<double> insertion_probabilities,
            NongenomicInsertionMatrix transition_matrix) :
        insertion_probabilities_(insertion_probabilities),
        transition_matrix_(transition_matrix) { }

NongenomicInsertionModel::NongenomicInsertionModel(std::ifstream& in) {
    assert(in.is_open());
    string line;
    vector<string> parsed_vector;
    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == 2);
        insertion_probabilities_.push_back(std::stod(parsed_vector[1]));
    }

    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        transition_matrix_.emplace_back(std::distance(tokenizer.begin(), tokenizer.end()));
        std::transform(tokenizer.begin(), tokenizer.end(),
                       transition_matrix_.back().begin(),
                       [] (string str) { return std::stod(str); });
    }
}

double NongenomicInsertionModel::GetInsertionProbabilityByLength(
        const unsigned int ins_length) const {
    assert(ins_length < insertion_probabilities_.size());
    return insertion_probabilities_[ins_length];
}

double NongenomicInsertionModel::GetTransitionProbability(char in, char out) const {
    using Alphabet = seqan::Dna5;
    size_t in_pos = static_cast<size_t>(seqan::ordValue(Alphabet(in)));
    size_t out_pos = static_cast<size_t>(seqan::ordValue(Alphabet(out)));
    assert(in_pos < 4);
    assert(out_pos < 4);
    return transition_matrix_[in_pos][out_pos];
}


std::ostream& operator<<(std::ostream& out, const NongenomicInsertionModel& model) {
    out << "Length, Probability" << "\n";
    for (auto it = model.GetInsertionProbabilities().cbegin();
            it != model.GetInsertionProbabilities().cend(); ++it)
        out << it - model.GetInsertionProbabilities().cbegin() << ", " << *it << "\n";
    out << "\n";

    out << "Markov chain transition matrix" << "\n";
    for (auto it1 = model.GetTransitionMatrix().cbegin();
            it1 != model.GetTransitionMatrix().cend(); ++it1) {
        for (auto it2 = it1->cbegin(); it2 != it1->cend(); ++it2)
            out << *it2 << " ";
        out << "\n";
    }
    return out;
}

double NongenomicInsertionModel::GetTransitionProbability(
        const std::pair<char, char>& transition) const {
    return GetTransitionProbability(transition.first, transition.second);
}

/**************************************************************************************************/

PalindromeDeletionModel::PalindromeDeletionModel(const DeletionTableVector& deletion_table,
                                                 const vector<int>& deletion_length,
                                                 const IgGeneDatabasePtrConst ig_gene_database) :
        deletion_table_(deletion_table),
        deletion_length_(deletion_length),
        ig_gene_database_(ig_gene_database) { }


PalindromeDeletionModel::PalindromeDeletionModel(std::ifstream& in,
                                                 const IgGeneDatabasePtrConst& ig_gene_database) :
        ig_gene_database_(ig_gene_database) {
    INFO("OK");
    assert(in.is_open());
    deletion_table_.resize(ig_gene_database -> size());

    string line;
    vector<string> parsed_vector;
    getline(in, line);
    Tokenizer tokenizer(line);
    size_t palidrome_len_diversity = std::distance(tokenizer.begin(), tokenizer.end()) - 1;
    vector<int> palindrome_lengths(palidrome_len_diversity);

    auto tokenizer_it = tokenizer.begin();
    tokenizer_it++;

    for (auto it = palindrome_lengths.begin();
             it != palindrome_lengths.end();
             ++it, ++tokenizer_it) {
        *it = std::stoi(*tokenizer_it);
    }

    for (auto it = deletion_table_.begin(); it != deletion_table_.end(); ++it)
        it -> resize(palidrome_len_diversity);

    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == palidrome_len_diversity + 1);
        auto palidrome_len_it = palindrome_lengths.begin();
        size_t index_of_current_gene = ig_gene_database -> GetIndexByName(parsed_vector.front());
        INFO(index_of_current_gene);
        for (auto parsed_it = parsed_vector.begin() + 1;
                parsed_it != parsed_vector.end();
                ++parsed_it, ++palidrome_len_it) {
            deletion_table_[index_of_current_gene][*palidrome_len_it] = std::stod(*parsed_it);
        }
    }
}

PalindromeDeletionModel::PalindromeDeletionModel(std::ifstream& in,
                                                 const IgGeneDatabase& ig_gene_database) :
        PalindromeDeletionModel(in, make_shared<const IgGeneDatabase>(ig_gene_database)) { }

std::ostream& operator<<(std::ostream& out, const PalindromeDeletionModel& model) {
    if (model.GetDeletionTable().empty())
        return out;
    out << "Gene id, Length of palindrome\n";

    for (int length : model.GetDeletionLength())
       out << length << " ";
    out << "\n";

    for (size_t i = 0; i < model.size(); ++i) {
        out << "Gene_id: " << model.GetIgGeneDatabase() -> GetByIndex(i) -> name();
        for (auto it = model.GetDeletionTable().at(i).cbegin();
                it != model.GetDeletionTable().at(i).cend();
                ++it)
            out << *it << " ";
        out << "\n";
    }
    return out;
}
/**************************************************************************************************/
