#include "recombination_model.hpp"

#include <cassert>

#include <map>
#include <string>
#include <algorithm>
#include <vector>

#include <boost/tokenizer.hpp>

using boost::tokenizer;
using boost::escaped_list_separator;
using Tokenizer = tokenizer<escaped_list_separator<char>>;

IgGeneProbabilityModel::IgGeneProbabilityModel(IgGeneProbabilityMap& ig_gene_probabilities) :
        ig_gene_probabilities_(ig_gene_probabilities) { }


double IgGeneProbabilityModel::GetProbabilityByGeneName(const string& index_string) const {
    assert(ig_gene_probabilities_.find(index_string) != ig_gene_probabilities_.end());
    return ig_gene_probabilities_.at(index_string);
}

size_t IgGeneProbabilityModel::size() const { return ig_gene_probabilities_.size(); }

IgGeneProbabilityModel::IgGeneProbabilityModel(std::ifstream& in) {
    assert(in.is_open());
    string line;
    vector<string> parsed_vector;
    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == 2);
        ig_gene_probabilities_[parsed_vector.front()] = std::stod(parsed_vector.back());
    }
}

std::ostream& operator<<(std::ostream& out, const IgGeneProbabilityModel& model) {
    for (auto& gene : model.GetGeneProbabilities())
        out << "Gene_id: " << gene.first << ", " << "Gene probability: " << gene.second << "\n";
    return out;
}


/************************************************************************************************/

const size_t NongenomicInsertionModel::alphabet_size_ = 4;

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

double NongenomicInsertionModel::GetInsertionProbabilityByLength(const unsigned int ins_length) const {
    assert(ins_length < insertion_probabilities_.size());
    return insertion_probabilities_[ins_length];
}

double NongenomicInsertionModel::GetTransitionProbability(char in, char out) const {
    string alphabet = "ACGT";
    size_t in_pos = alphabet.find(in);
    size_t out_pos = alphabet.find(out);
    assert(in_pos != string::npos);
    assert(out_pos != string::npos);
    return transition_matrix_[in_pos][out_pos];
}

std::ostream& operator<<(std::ostream& out, const NongenomicInsertionModel& model) {
    out << "Length, Probability" << "\n";
    for (auto it = model.GetInsertionProbabilities().begin();
            it != model.GetInsertionProbabilities().end(); ++it)
        out << it - model.GetInsertionProbabilities().begin() << ", " << *it << "\n";
    out << "\n";

    out << "Markov chain transition matrix" << "\n";
    out << "\tA\tC\tG\tT" << "\n";
    string alphabet = "ACGT";
    for (auto it1 = model.GetTransitionMatrix().begin();
            it1 != model.GetTransitionMatrix().end(); ++it1) {
        out << alphabet[it1 - model.GetTransitionMatrix().begin()] << ": ";
        for (auto it2 = it1->begin(); it2 != it1->end(); ++it2)
            out << *it2 << " \n"[it2 == it1->end() - 1];
    }
    return out;
}

/**************************************************************************************************/

PalindromeDeletionModel::PalindromeDeletionModel(DeletionTableMap deletion_table) :
    deletion_table_(deletion_table) { }


PalindromeDeletionModel::PalindromeDeletionModel(std::ifstream& in) {
    assert(in.is_open());
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

    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == palidrome_len_diversity + 1);
        auto palidrome_len_it = palindrome_lengths.begin();
        for (auto parsed_it = parsed_vector.begin() + 1;
                parsed_it != parsed_vector.end();
                ++parsed_it, ++palidrome_len_it) {
            deletion_table_[parsed_vector.front()][*palidrome_len_it] = std::stod(*parsed_it);
        }
    }
}

std::ostream& operator<<(std::ostream& out, const PalindromeDeletionModel& model) {
    if (model.GetDeletionTable().empty())
        return out;
    out << "Gene id, Length of palindrome\n";

    auto& gene = (*(model.GetDeletionTable().begin())).second;
    for (auto& key_value : gene)
        out << key_value.first << " ";
    out << "\n";

    for (auto& gene : model.GetDeletionTable()) {
        out << gene.first << ": ";
        for (auto& probability : gene.second)
            out << probability.second << " ";
        out << "\n";
    }
    return out;
}
/**************************************************************************************************/

HCProbabilityRecombinationModel::HCProbabilityRecombinationModel(std::ifstream& in) :
        V_gene_probability_model_(in),
        D_gene_probability_model_(in),
        J_gene_probability_model_(in),
        VD_nongenomic_insertion_model_(in),
        DJ_nongenomic_insertion_model_(in),
        V_palindrome_deletion_model_(in),
        J_palindrome_deletion_model_(in),
        DLeft_palindrome_deletion_model_(in),
        DRight_palindrome_deletion_model_(in) {
}

std::ostream& operator<<(std::ostream& out, const HCProbabilityRecombinationModel& model) {
    out << "V gene probabilities:\n";
    out << model.GetVGeneProbabilityModel();
    out << "\nD gene probabilities:\n";
    out << model.GetDGeneProbabilityModel();
    out << "\nJ gene probabilities:\n";
    out << model.GetJGeneProbabilityModel();

    out << "\nVD nongenomic insertions model:\n";
    out << model.GetVDNongenomicInsertionModel();
    out << "\nDJ nongenomic insertions model:\n";
    out << model.GetDJNongenomicInsertionModel();

    out << "\nV palindrome deletion model:\n";
    out << model.GetVPalindromeDeletionModel();
    out << "\nJ palindrome deletion model:\n";
    out << model.GetJPalindromeDeletionModel();
    out << "\nD left palindrome deletion model:\n";
    out << model.GetDLeftPalindromeDeletionModel();
    out << "\nD right palindrome deletion model:\n";
    out << model.GetDRightPalindromeDeletionModel();
    return out;
}
