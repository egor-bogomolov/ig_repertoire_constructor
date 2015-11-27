#pragma once

#include "standard.hpp"

#include <map>
#include <string>
#include <fstream>
#include <vector>

using std::string;
using std::map;
using std::vector;

class IgGeneProbabilityModel {
    using IgGeneProbabilityMap = map<string, double>;
    IgGeneProbabilityMap ig_gene_probabilities_;

 public:
    IgGeneProbabilityModel() = delete;

    explicit IgGeneProbabilityModel(IgGeneProbabilityMap&);

    IgGeneProbabilityModel(const IgGeneProbabilityModel&) = default;

    IgGeneProbabilityModel(IgGeneProbabilityModel&&) = default;

    IgGeneProbabilityModel& operator=(const IgGeneProbabilityModel&) = default;

    IgGeneProbabilityModel& operator=(IgGeneProbabilityModel&&) = default;

    virtual ~IgGeneProbabilityModel() = default;

    const IgGeneProbabilityMap& GetGeneProbabilities() const {
        return ig_gene_probabilities_;
    }

    void SetGeneProbabilities(const IgGeneProbabilityMap& ig_gene_probabilities) {
        ig_gene_probabilities_ = ig_gene_probabilities;
    }

    using citerator = IgGeneProbabilityMap::const_iterator;

    citerator cbegin() const { return ig_gene_probabilities_.cbegin(); }
    citerator cend() const { return ig_gene_probabilities_.cend(); }

    double GetProbabilityByGeneName(const string&) const;

    size_t size() const;

    explicit IgGeneProbabilityModel(std::ifstream&);
};

std::ostream& operator<<(std::ostream&, const IgGeneProbabilityModel&);

class NongenomicInsertionModel {
    static const size_t alphabet_size_;
    vector<double> insertion_probabilities_;
    using NongenomicInsertionMatrix=vector<vector<double>>;
    NongenomicInsertionMatrix transition_matrix_;

 public:
    NongenomicInsertionModel() = delete;

    NongenomicInsertionModel(vector<double>, NongenomicInsertionMatrix);

    NongenomicInsertionModel(const NongenomicInsertionModel&) = default;

    NongenomicInsertionModel(NongenomicInsertionModel&&) = default;

    NongenomicInsertionModel& operator=(const NongenomicInsertionModel&) = default;

    NongenomicInsertionModel& operator=(NongenomicInsertionModel&&) = default;

    virtual ~NongenomicInsertionModel() = default;

    const vector<double>& GetInsertionProbabilities() const { return insertion_probabilities_; }

    double GetInsertionProbabilityByLength(const unsigned int) const;

    const NongenomicInsertionMatrix& GetTransitionMatrix() const { return transition_matrix_; }

    void SetInsertionProbabilities(const vector<double>& insertion_probabilities) {
        insertion_probabilities_ = insertion_probabilities;
    }

    void SetTransitionMatrix(const NongenomicInsertionMatrix& transition_matrix) {
        transition_matrix_ = transition_matrix;
    }

    explicit NongenomicInsertionModel(std::ifstream&);

    double GetTransitionProbability(char, char) const;
};

std::ostream& operator<<(std::ostream&, const NongenomicInsertionModel&);


class PalindromeDeletionModel {
    using DeletionTableMap = map<string, map<int, double>>;
    DeletionTableMap deletion_table_;

 public:
    PalindromeDeletionModel() = delete;

    explicit PalindromeDeletionModel(DeletionTableMap deletion_table);

    PalindromeDeletionModel(const PalindromeDeletionModel&) = default;

    PalindromeDeletionModel(PalindromeDeletionModel&&) = default;

    PalindromeDeletionModel& operator=(const PalindromeDeletionModel&) = default;

    PalindromeDeletionModel& operator=(PalindromeDeletionModel&&) = default;

    virtual ~PalindromeDeletionModel() = default;

    const DeletionTableMap& GetDeletionTable() const { return deletion_table_; }

    void SetDeletionTable(const DeletionTableMap& deletion_table) {
        deletion_table_ = deletion_table;
    }

    using citerator = DeletionTableMap::const_iterator;

    citerator cbegin() const { return deletion_table_.cbegin(); }
    citerator cend() const { return deletion_table_.cend(); }

    explicit PalindromeDeletionModel(std::ifstream&);
};

std::ostream& operator<<(std::ostream&, const PalindromeDeletionModel&);

class HCProbabilityRecombinationModel {
    IgGeneProbabilityModel V_gene_probability_model_;
    IgGeneProbabilityModel D_gene_probability_model_;
    IgGeneProbabilityModel J_gene_probability_model_;

    NongenomicInsertionModel VD_nongenomic_insertion_model_;
    NongenomicInsertionModel DJ_nongenomic_insertion_model_;

    PalindromeDeletionModel V_palindrome_deletion_model_;
    PalindromeDeletionModel J_palindrome_deletion_model_;
    PalindromeDeletionModel DLeft_palindrome_deletion_model_;
    PalindromeDeletionModel DRight_palindrome_deletion_model_;

 public:
    HCProbabilityRecombinationModel() = delete;

    HCProbabilityRecombinationModel(const HCProbabilityRecombinationModel&) = default;

    HCProbabilityRecombinationModel(HCProbabilityRecombinationModel&&) = default;

    HCProbabilityRecombinationModel& operator=(const HCProbabilityRecombinationModel&) = default;

    HCProbabilityRecombinationModel& operator=(HCProbabilityRecombinationModel&&) = default;

    virtual ~HCProbabilityRecombinationModel() = default;

    const IgGeneProbabilityModel& GetVGeneProbabilityModel() const {
        return V_gene_probability_model_;
    }

    const IgGeneProbabilityModel& GetDGeneProbabilityModel() const {
        return D_gene_probability_model_;
    }

    const IgGeneProbabilityModel& GetJGeneProbabilityModel() const {
        return J_gene_probability_model_;
    }

    const NongenomicInsertionModel& GetVDNongenomicInsertionModel() const {
        return VD_nongenomic_insertion_model_;
    } 

    const NongenomicInsertionModel& GetDJNongenomicInsertionModel() const {
        return DJ_nongenomic_insertion_model_;
    } 

    const PalindromeDeletionModel& GetVPalindromeDeletionModel() const {
        return V_palindrome_deletion_model_;
    } 

    const PalindromeDeletionModel& GetJPalindromeDeletionModel() const {
        return J_palindrome_deletion_model_;
    } 

    const PalindromeDeletionModel& GetDLeftPalindromeDeletionModel() const {
        return DLeft_palindrome_deletion_model_;
    } 

    const PalindromeDeletionModel& GetDRightPalindromeDeletionModel() const {
        return DRight_palindrome_deletion_model_;
    } 

    void SetVGeneProbabilityModel(const IgGeneProbabilityModel V_gene_probability_model) {
        V_gene_probability_model_ = V_gene_probability_model;
    }

    void SetDGeneProbabilityModel(const IgGeneProbabilityModel D_gene_probability_model) {
        D_gene_probability_model_ = D_gene_probability_model;
    }

    void SetJGeneProbabilityModel(const IgGeneProbabilityModel J_gene_probability_model) {
        J_gene_probability_model_ = J_gene_probability_model;
    }

    explicit HCProbabilityRecombinationModel(std::ifstream&);
};

std::ostream& operator<<(std::ostream&, const HCProbabilityRecombinationModel&);
