#pragma once

#include "standard.hpp"
#include "../vdj_alignments/gene_database.hpp"

#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <memory>

using std::string;
using std::map;
using std::vector;
using std::shared_ptr;

using IgGeneDatabasePtrConst = shared_ptr<const IgGeneDatabase>;
using HC_GenesDatabase_PtrConst = shared_ptr<const HC_GenesDatabase>;

class IgGeneProbabilityModel {
    using IgGeneProbabilityVector = vector<double>;
    IgGeneProbabilityVector ig_gene_probabilities_;
    IgGeneDatabasePtrConst ig_gene_database_;

 public:
    IgGeneProbabilityModel() = delete;

    IgGeneProbabilityModel(const IgGeneProbabilityVector&, const IgGeneDatabasePtrConst&);

    IgGeneProbabilityModel(const IgGeneProbabilityModel&) = default;

    IgGeneProbabilityModel(IgGeneProbabilityModel&&) = default;

    IgGeneProbabilityModel& operator=(const IgGeneProbabilityModel&) = default;

    IgGeneProbabilityModel& operator=(IgGeneProbabilityModel&&) = default;

    virtual ~IgGeneProbabilityModel() = default;

    const IgGeneProbabilityVector& GetIgGeneProbabilities() const {
        return ig_gene_probabilities_;
    }

    const IgGeneDatabasePtrConst GetIgGeneDatabase() const {
        return ig_gene_database_;
    }

    void SetGeneProbabilities(const IgGeneProbabilityVector& ig_gene_probabilities) {
        ig_gene_probabilities_ = ig_gene_probabilities;
    }

    using citerator = IgGeneProbabilityVector::const_iterator;

    citerator cbegin() const { return ig_gene_probabilities_.cbegin(); }
    citerator cend() const { return ig_gene_probabilities_.cend(); }

    size_t size() const;

    IgGeneProbabilityModel(std::ifstream&, const IgGeneDatabasePtrConst&);

    IgGeneProbabilityModel(std::ifstream&, const IgGeneDatabase&);

    //double GetProbabilityByGenId(const IgGene&) const;
    double GetProbabilityByGenId(const size_t&) const;
};

std::ostream& operator<<(std::ostream&, const IgGeneProbabilityModel&);


class NongenomicInsertionModel {
    vector<double> insertion_probabilities_;
    using NongenomicInsertionMatrix = vector<vector<double>>;
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

    double GetInsertionProbabilityByLength(const long unsigned int) const;

    const NongenomicInsertionMatrix& GetTransitionMatrix() const { return transition_matrix_; }

    void SetInsertionProbabilities(const vector<double>& insertion_probabilities) {
        insertion_probabilities_ = insertion_probabilities;
    }

    void SetTransitionMatrix(const NongenomicInsertionMatrix& transition_matrix) {
        transition_matrix_ = transition_matrix;
    }

    explicit NongenomicInsertionModel(std::ifstream&);

    double GetTransitionProbability(char, char) const;

    double GetTransitionProbability(const std::pair<char, char>&) const;
};

std::ostream& operator<<(std::ostream&, const NongenomicInsertionModel&);


class PalindromeDeletionModel {
    using DeletionTableVector = vector<vector<double>>;
    DeletionTableVector deletion_table_;
    vector<int> deletion_length_;
    IgGeneDatabasePtrConst ig_gene_database_;

 public:
    PalindromeDeletionModel() = delete;

    PalindromeDeletionModel(const DeletionTableVector&,
                            const vector<int>&,
                            const IgGeneDatabasePtrConst&);

    PalindromeDeletionModel(const PalindromeDeletionModel&) = default;

    PalindromeDeletionModel(PalindromeDeletionModel&&) = default;

    PalindromeDeletionModel& operator=(const PalindromeDeletionModel&) = default;

    PalindromeDeletionModel& operator=(PalindromeDeletionModel&&) = default;

    virtual ~PalindromeDeletionModel() = default;

    const DeletionTableVector& GetDeletionTable() const { return deletion_table_; }

    void SetDeletionTable(const DeletionTableVector& deletion_table) {
        deletion_table_ = deletion_table;
    }

    const IgGeneDatabasePtrConst GetIgGeneDatabase() const {
        return ig_gene_database_;
    }

    const vector<int>& GetDeletionLength() const { return deletion_length_; }

    using citerator = DeletionTableVector::const_iterator;

    citerator cbegin() const { return deletion_table_.cbegin(); }
    citerator cend() const { return deletion_table_.cend(); }

    size_t size() const { return deletion_table_.size(); }

    PalindromeDeletionModel(std::ifstream&, const IgGeneDatabasePtrConst&);

    PalindromeDeletionModel(std::ifstream&, const IgGeneDatabase&);

    double GetDeletionProbability(const size_t& gene_id, const int deletion_length) const;
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

    HC_GenesDatabase_PtrConst HC_database;

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

    HCProbabilityRecombinationModel(std::ifstream&, const HC_GenesDatabase_PtrConst&);

    HCProbabilityRecombinationModel(std::ifstream&, const HC_GenesDatabase&);

    //double GetProbabilityByGenId(const IgGene& ig_gene) const;
    double GetProbabilityByGenId(enum IgGeneType ig_gene_type,
            const size_t& id) const;
};

std::ostream& operator<<(std::ostream&, const HCProbabilityRecombinationModel&);
