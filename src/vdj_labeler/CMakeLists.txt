############################################################################
# Copyright (c) 2015-2015 Saint-Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(vdj_labeler CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})

find_path(SEQAN_INCLUDE_DIR seqan/seq_io.h
  HINTS ../../ext/include ${EXT_DIR}/include)

if(SEQAN_INCLUDE_DIR)
  message("seqan headers found in ${SEQAN_INCLUDE_DIR}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})

add_executable(vdj_labeler
                fastq_read_archive.cpp
                vdj_alignments/gene_database.cpp
                vdj_alignments/alignment_structs.cpp
                vdj_alignments/vj_alignment_info.cpp
                vdj_alignments/aligners/right_v_tail_aligner.cpp
                vdj_alignments/aligners/left_j_tail_aligner.cpp
                vdj_alignments/aligners/simple_d_aligner.cpp
                vdj_alignments/vdj_hits.cpp
                vdj_alignments/vdj_hits_storage.cpp
                vdj_alignments/vdj_hits_calculators/custom_vdj_hits_calculator.cpp
                vdj_alignments/vdj_hits_calculators/info_based_vj_hits_calculator.cpp
                vdj_alignments/vdj_hits_calculators/info_based_d_hits_calculator.cpp
                recombination/cleaved_gene.cpp
                recombination/nongenomic_insertion.cpp
                recombination/recombination.cpp
                recombination_generation/gene_events_generators/shms_calculators/left_event_shms_calculator.cpp
                recombination_generation/gene_events_generators/shms_calculators/right_event_shms_calculator.cpp
                recombination_generation/gene_events_generators/shms_calculators/versatile_shms_calculator.cpp
                recombination_generation/recombination_event_storage.cpp
                recombination_generation/custom_hc_recombination_generator.cpp
                recombination_generation/gene_events_generators/v_recombination_event_generator.cpp
                recombination_generation/gene_events_generators/d_recombination_event_generator.cpp
                recombination_generation/gene_events_generators/j_recombination_event_generator.cpp
                recombination_generation/insertion_events_generators/versatile_insertion_event_generator.cpp
                recombination_calculator/hc_model_based_recombination_calculator.cpp
                model/recombination_model.cpp
                main.cpp)

target_link_libraries(vdj_labeler input yaml-cpp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(vdj_labeler PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS vdj_labeler
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/vdj_labeler"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
