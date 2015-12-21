#include "vdj_config.hpp"
#include "../include/config_common.hpp"

void load(vdj_config::io_params::input_params::germline_genes::igh_genes &igh,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(igh.variable_genes, pt, "variable_genes");
    load(igh.diversity_genes, pt, "diversity_genes");
    load(igh.join_genes, pt, "join_genes");
}

void load(vdj_config::io_params::input_params::germline_genes::igl_genes &igl,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(igl.variable_genes, pt, "variable_genes");
    load(igl.join_genes, pt, "join_genes");
}

void load(vdj_config::io_params::input_params::germline_genes::igk_genes &igk,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(igk.variable_genes, pt, "variable_genes");
    load(igk.join_genes, pt, "join_genes");
}

void load(vdj_config::io_params::input_params::germline_genes &germlines,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(germlines.igh, pt, "igh");
    load(germlines.igl, pt, "igl");
    load(germlines.igk, pt, "igk");
}

void load(vdj_config::io_params::input_params &input, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input.input_sequences, pt, "input_sequences");
    load(input.vj_alignment_file, pt, "vj_alignment_file");
    load(input.germlines, pt, "germlines");
}

void load(vdj_config::io_params::output_params &output, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output.log_filename, pt, "log_filename");
    load(output.output_dir, pt, "output_dir");
}

void load(vdj_config::io_params &iop, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(iop.input, pt, "input_params");
    load(iop.output, pt, "output_params");
}

void load(vdj_config::run_params &rp, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(rp.max_memory, pt, "max_memory");
    load(rp.threads_count, pt, "thread_count");
}

void load(vdj_config &vdj, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(vdj.io, pt, "io_params");
    load(vdj.rp, pt, "run_params");
}

void load(vdj_config &cfg, std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}