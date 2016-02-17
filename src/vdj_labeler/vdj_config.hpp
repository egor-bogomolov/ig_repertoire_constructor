#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

struct vdj_config {
    struct io_params {
        struct input_params {
            std::string     input_sequences;
            std::string     vj_alignment_file;

            struct germline_genes {
                struct igh_genes {
                    std::string variable_genes;
                    std::string diversity_genes;
                    std::string join_genes;
                };

                struct igl_genes {
                    std::string variable_genes;
                    std::string join_genes;
                };

                struct igk_genes {
                    std::string variable_genes;
                    std::string join_genes;
                };

                igh_genes igh;
                igl_genes igl;
                igk_genes igk;
            };

            germline_genes germlines;
        };

        struct output_params {
            std::string     log_filename;
            std::string     output_dir;
        };

        input_params input;
        output_params output;
    };

    struct run_params {
        unsigned        threads_count;
        unsigned        max_memory;
    };

    io_params io;
    run_params rp;
};

void load(vdj_config &cfg, const std::string &filename);

typedef config_common::config<vdj_config> vdj_cfg;
