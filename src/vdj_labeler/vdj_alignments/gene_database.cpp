#include "gene_database.hpp"

#include "seqan/sequence.h"

using namespace std;

std::string IgGeneTypeToString(IgGeneType gene_type) {
    if(gene_type == variable_gene)
        return "variable";
    if(gene_type == diversity_gene)
        return "diversity";
    return "join";
}

//-------------------------------------------------------------------------------------------

std::ostream& operator<< (std::ostream &out, const IgGene &obj) {
    out << "Name: " << obj.name() << ". Seq (len " << obj.length() << "): " << obj.seq();
    return out;
}

//-------------------------------------------------------------------------------------------

void IgGeneDatabase::AddGenesFromFile(std::string filename) {
    std::vector<CharString> read_headers;
    std::vector<Dna5String> read_seqs;
    seqan::SeqFileIn seqFileIn_reads(filename.c_str());
    readRecords(read_headers, read_seqs, seqFileIn_reads);
    for(size_t i = 0; i < read_headers.size(); i++) {
        IgGenePtr ig_gene_ptr(new IgGene(gene_type_, read_headers[i], read_seqs[i]));
        ig_genes_.push_back(ig_gene_ptr);
        gene_name_map_[string(toCString(read_headers[i]))] = ig_gene_ptr;
        gene_index_map_[string(toCString(read_headers[i]))] = i;
    }
}

IgGenePtr IgGeneDatabase::GetByIndex(size_t index) const {
    assert(index < size());
    return ig_genes_[index];
}

IgGenePtr IgGeneDatabase::GetByName(std::string gene_name) const {
    assert(gene_name_map_.find(gene_name) != gene_name_map_.end());
    return gene_name_map_.at(gene_name);
}

size_t IgGeneDatabase::GetIndexByName(std::string gene_name) const {
    assert(gene_index_map_.find(gene_name) != gene_index_map_.end());
    return gene_index_map_.at(gene_name);
}

std::ostream& operator<<(std::ostream &out, const IgGeneDatabase &ig_gene_db) {
    out << "Ig genes database. Gene type: " << IgGeneTypeToString(ig_gene_db.GeneType()) << ". # records: " <<
            ig_gene_db.size() << std::endl;
    for(auto it = ig_gene_db.cbegin(); it != ig_gene_db.cend(); it++)
        out << *(*it) << std::endl;
    return out;
}

//-------------------------------------------------------------------------------------------

void HC_GenesDatabase::AddGenesFromFile(IgGeneType gene_type, std::string filename){
    if(gene_type == variable_gene)
        variable_genes_.AddGenesFromFile(filename);
    else
    if(gene_type == diversity_gene)
        diversity_genes_.AddGenesFromFile(filename);
    else
        join_genes_.AddGenesFromFile(filename);
}

size_t HC_GenesDatabase::GenesNumber(IgGeneType gene_type) const {
    if(gene_type == variable_gene)
        return variable_genes_.size();
    if(gene_type == diversity_gene)
        return diversity_genes_.size();
    return join_genes_.size();
}

IgGenePtr HC_GenesDatabase::GetByIndex(IgGeneType gene_type, size_t index) const {
    if(gene_type == variable_gene)
        return variable_genes_.GetByIndex(index);
    if(gene_type == diversity_gene)
        return diversity_genes_.GetByIndex(index);
    return join_genes_.GetByIndex(index);
}

std::ostream& operator<<(std::ostream &out, const HC_GenesDatabase& obj) {
    out << obj.VariableGenes();
    out << "--------------------" << std::endl;
    out << obj.DiversityGenes();
    out << "--------------------" << std::endl;
    out << obj.JoinGenes();
    return out;
}

//-------------------------------------------------------------------------------------------

void LC_GenesDatabase::AddGenesFromFile(IgGeneType gene_type, std::string filename){
    if(gene_type == variable_gene)
        variable_genes_.AddGenesFromFile(filename);
    else
        join_genes_.AddGenesFromFile(filename);
}

size_t LC_GenesDatabase::GenesNumber(IgGeneType gene_type) const {
    if(gene_type == variable_gene)
        return variable_genes_.size();
    return join_genes_.size();
}

IgGenePtr LC_GenesDatabase::GetByIndex(IgGeneType gene_type, size_t index) const {
    if(gene_type == variable_gene)
        return variable_genes_.GetByIndex(index);
    return join_genes_.GetByIndex(index);
}

std::ostream& operator<<(std::ostream &out, const LC_GenesDatabase& obj) {
    out << obj.VariableGenes();
    out << "--------------------" << std::endl;
    out << obj.JoinGenes();
    return out;
}