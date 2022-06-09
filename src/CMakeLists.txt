set(CONFIG_SOURCES
    config/config.hpp
    config/config.cpp
    config/common.hpp
    config/common.cpp
    config/option_parser.hpp
    config/option_parser.cpp
    config/option_collation.hpp
    config/option_collation.cpp
    config/octopus_vcf.hpp
    config/octopus_vcf.cpp
)

set(EXCEPTIONS_SOURCES
    exceptions/error.hpp
    exceptions/error.cpp
    exceptions/user_error.hpp
    exceptions/program_error.hpp
    exceptions/system_error.hpp
    exceptions/missing_file_error.hpp
    exceptions/missing_file_error.cpp
    exceptions/malformed_file_error.hpp
    exceptions/malformed_file_error.cpp
    exceptions/missing_index_error.hpp
    exceptions/missing_index_error.cpp
    exceptions/unwritable_file_error.hpp
    exceptions/unwritable_file_error.cpp
    exceptions/unimplemented_feature_error.hpp
    exceptions/unimplemented_feature_error.cpp
    exceptions/file_open_error.hpp
    exceptions/file_open_error.cpp
)

set(CONCEPTS_SOURCES
    concepts/equitable.hpp
    concepts/comparable.hpp
    concepts/mappable.hpp
    concepts/mappable_range.hpp
    concepts/indexed.hpp
)

set(BASICS_SOURCES
    basics/contig_region.hpp
    basics/genomic_region.hpp
    basics/phred.hpp
    basics/cigar_string.hpp
    basics/cigar_string.cpp
    basics/aligned_read.hpp
    basics/aligned_read.cpp
    basics/mappable_reference_wrapper.hpp
    basics/ploidy_map.hpp
    basics/ploidy_map.cpp
    basics/pedigree.hpp
    basics/pedigree.cpp
    basics/trio.hpp
    basics/trio.cpp
    basics/read_pileup.hpp
    basics/read_pileup.cpp
    basics/tandem_repeat.hpp
    basics/tandem_repeat.cpp
    basics/aligned_template.hpp
    basics/aligned_template.cpp
)

set(CONTAINERS_SOURCES
    containers/mappable_flat_multi_set.hpp
    containers/mappable_flat_set.hpp
    containers/mappable_map.hpp
    containers/matrix_map.hpp
    containers/probability_matrix.hpp
    containers/mappable_block.hpp
)

set(LOGGING_SOURCES
    logging/logging.hpp
    logging/logging.cpp
    logging/progress_meter.hpp
    logging/progress_meter.cpp
    logging/error_handler.hpp
    logging/error_handler.cpp
    logging/main_logging.hpp
    logging/main_logging.cpp
)

set(IO_SOURCES
    io/reference/caching_fasta.hpp
    io/reference/caching_fasta.cpp
    io/reference/fasta.hpp
    io/reference/fasta.cpp
    io/reference/reference_genome.hpp
    io/reference/reference_genome.cpp
    io/reference/reference_reader.hpp
    io/reference/threadsafe_fasta.hpp
    io/reference/threadsafe_fasta.cpp

    io/region/region_parser.hpp
    io/region/region_parser.cpp

    io/pedigree/pedigree_reader.hpp
    io/pedigree/pedigree_reader.cpp

    io/read/htslib_sam_facade.hpp
    io/read/htslib_sam_facade.cpp
    io/read/read_manager.hpp
    io/read/read_manager.cpp
    io/read/read_reader_impl.hpp
    io/read/read_reader.hpp
    io/read/read_reader.cpp
    io/read/read_writer.hpp
    io/read/read_writer.cpp
    io/read/buffered_read_writer.hpp
    
    io/variant/htslib_bcf_facade.hpp
    io/variant/htslib_bcf_facade.cpp
    io/variant/vcf_header.hpp
    io/variant/vcf_header.cpp
    io/variant/vcf_parser.hpp
    io/variant/vcf_parser.cpp
    io/variant/vcf_reader_impl.hpp
    io/variant/vcf_reader.hpp
    io/variant/vcf_reader.cpp
    io/variant/vcf_record.hpp
    io/variant/vcf_record.cpp
    io/variant/vcf_type.hpp
    io/variant/vcf_type.cpp
    io/variant/vcf_utils.hpp
    io/variant/vcf_utils.cpp
    io/variant/vcf_writer.hpp
    io/variant/vcf_writer.cpp
    io/variant/vcf.hpp
    io/variant/vcf_spec.hpp
)

set(READPIPE_SOURCES
    readpipe/read_pipe_fwd.hpp
    readpipe/read_pipe.hpp
    readpipe/read_pipe.cpp
    readpipe/buffered_read_pipe.hpp
    readpipe/buffered_read_pipe.cpp
    
    readpipe/downsampling/downsampler.hpp
    readpipe/downsampling/downsampler.cpp
    
    readpipe/filtering/read_filter.hpp
    readpipe/filtering/read_filter.cpp
    readpipe/filtering/read_filterer.hpp
    
    readpipe/transformers/read_transform.hpp
    readpipe/transformers/read_transform.cpp
    readpipe/transformers/read_transformer.hpp
    readpipe/transformers/read_transformer.cpp
)

set(UTILS_SOURCES
    utils/append.hpp
    utils/compression.hpp
    utils/compression.cpp
    utils/hash_functions.hpp
    utils/map_utils.hpp
    utils/mappable_algorithms.hpp
    utils/maths.hpp
    utils/merge_transform.hpp
    utils/path_utils.hpp
    utils/path_utils.cpp
    utils/read_algorithms.hpp
    utils/read_stats.hpp
    utils/read_stats.cpp
    utils/sequence_utils.hpp
    utils/string_utils.hpp
    utils/string_utils.cpp
    utils/timing.hpp
    utils/type_tricks.hpp
    utils/coverage_tracker.hpp
    utils/input_reads_profiler.hpp
    utils/input_reads_profiler.cpp
    utils/kmer_mapper.hpp
    utils/kmer_mapper.cpp
    utils/memory_footprint.hpp
    utils/memory_footprint.cpp
    utils/emplace_iterator.hpp
    utils/repeat_finder.hpp
    utils/repeat_finder.cpp
    utils/genotype_reader.hpp
    utils/genotype_reader.cpp
    utils/beta_distribution.hpp
    utils/parallel_transform.hpp
    utils/thread_pool.hpp
    utils/thread_pool.cpp
    utils/concat.hpp
    utils/select_top_k.hpp
    utils/system_utils.hpp
    utils/system_utils.cpp
    utils/k_medoids.hpp
    utils/random_select.hpp
    utils/read_duplicates.hpp
    utils/read_duplicates.cpp
    utils/array_tricks.hpp
    utils/reorder.hpp
    utils/free_memory.hpp
    utils/erase_if.hpp
    utils/set_partition.hpp
)

set(CORE_SOURCES
    core/callers/caller_builder.hpp
    core/callers/caller_builder.cpp
    core/callers/caller_factory.hpp
    core/callers/caller_factory.cpp
    core/callers/caller.hpp
    core/callers/caller.cpp
    core/callers/cancer_caller.hpp
    core/callers/cancer_caller.cpp
    core/callers/individual_caller.hpp
    core/callers/individual_caller.cpp
    core/callers/population_caller.hpp
    core/callers/population_caller.cpp
    core/callers/trio_caller.hpp
    core/callers/trio_caller.cpp
    core/callers/polyclone_caller.hpp
    core/callers/polyclone_caller.cpp
    core/callers/cell_caller.hpp
    core/callers/cell_caller.cpp

    core/types/calls/call_types.hpp
    core/types/calls/call.hpp
    core/types/calls/call.cpp
    core/types/calls/germline_variant_call.hpp
    core/types/calls/germline_variant_call.cpp
    core/types/calls/reference_call.hpp
    core/types/calls/reference_call.cpp
    core/types/calls/somatic_call.hpp
    core/types/calls/somatic_call.cpp
    core/types/calls/variant_call.hpp
    core/types/calls/variant_call.cpp
    core/types/calls/denovo_call.hpp
    core/types/calls/denovo_call.cpp
    core/types/calls/denovo_reference_reversion_call.hpp
    core/types/calls/denovo_reference_reversion_call.cpp
    core/types/calls/cell_variant_call.hpp
    core/types/calls/cell_variant_call.cpp
    core/types/calls/polyclone_variant_call.hpp
    core/types/calls/polyclone_variant_call.cpp
    core/types/calls/cnv_call.hpp
    core/types/calls/cnv_call.cpp
    core/types/calls/call_wrapper.hpp
    core/types/calls/call_wrapper.cpp
    core/types/calls/call_utils.hpp
    core/types/calls/call_utils.cpp
    core/tools/vcf_header_factory.hpp
    core/tools/vcf_header_factory.cpp
    core/tools/vcf_record_factory.hpp
    core/tools/vcf_record_factory.cpp

    core/csr/facets/facet.hpp
    core/csr/facets/facet.cpp
    core/csr/facets/samples.hpp
    core/csr/facets/samples.cpp
    core/csr/facets/overlapping_reads.hpp
    core/csr/facets/overlapping_reads.cpp
    core/csr/facets/read_assignments.hpp
    core/csr/facets/read_assignments.cpp
    core/csr/facets/reference_context.hpp
    core/csr/facets/reference_context.cpp
    core/csr/facets/genotypes.hpp
    core/csr/facets/genotypes.cpp
    core/csr/facets/alleles.hpp
    core/csr/facets/alleles.cpp
    core/csr/facets/ploidies.hpp
    core/csr/facets/ploidies.cpp
    core/csr/facets/pedigree.hpp
    core/csr/facets/pedigree.cpp
    core/csr/facets/repeat_context.hpp
    core/csr/facets/repeat_context.cpp
    core/csr/facets/reads_summary.hpp
    core/csr/facets/reads_summary.cpp
    core/csr/facets/facet_factory.hpp
    core/csr/facets/facet_factory.cpp

    core/csr/filters/variant_call_filter.hpp
    core/csr/filters/variant_call_filter.cpp
    core/csr/filters/single_pass_variant_call_filter.hpp
    core/csr/filters/single_pass_variant_call_filter.cpp
    core/csr/filters/double_pass_variant_call_filter.hpp
    core/csr/filters/double_pass_variant_call_filter.cpp
    core/csr/filters/threshold_filter.hpp
    core/csr/filters/threshold_filter.cpp
    core/csr/filters/unsupervised_clustering_filter.hpp
    core/csr/filters/unsupervised_clustering_filter.cpp
    core/csr/filters/variant_call_filter_factory.hpp
    core/csr/filters/variant_call_filter_factory.cpp
    core/csr/filters/threshold_filter_factory.hpp
    core/csr/filters/threshold_filter_factory.cpp
    core/csr/filters/unsupervised_clustering_filter_factory.hpp
    core/csr/filters/unsupervised_clustering_filter_factory.cpp
    core/csr/filters/passing_filter.hpp
    core/csr/filters/passing_filter.cpp
    core/csr/filters/training_filter_factory.hpp
    core/csr/filters/training_filter_factory.cpp
    core/csr/filters/conditional_threshold_filter.hpp
    core/csr/filters/conditional_threshold_filter.cpp
    core/csr/filters/somatic_threshold_filter.hpp
    core/csr/filters/somatic_threshold_filter.cpp
    core/csr/filters/denovo_threshold_filter.hpp
    core/csr/filters/denovo_threshold_filter.cpp
    core/csr/filters/random_forest_filter.hpp
    core/csr/filters/random_forest_filter.cpp
    core/csr/filters/random_forest_filter_factory.hpp
    core/csr/filters/random_forest_filter_factory.cpp
    core/csr/filters/somatic_random_forest_filter.hpp
    core/csr/filters/somatic_random_forest_filter.cpp
    core/csr/filters/denovo_random_forest_filter.hpp
    core/csr/filters/denovo_random_forest_filter.cpp
    core/csr/filters/variant_filter_utils.hpp
    core/csr/filters/variant_filter_utils.cpp
    
    core/csr/measures/measure.hpp
    core/csr/measures/measure.cpp
    core/csr/measures/quality.hpp
    core/csr/measures/quality.cpp
    core/csr/measures/depth.hpp
    core/csr/measures/depth.cpp
    core/csr/measures/quality_by_depth.hpp
    core/csr/measures/quality_by_depth.cpp
    core/csr/measures/genotype_quality.hpp
    core/csr/measures/genotype_quality.cpp
    core/csr/measures/genotype_quality_by_depth.hpp
    core/csr/measures/genotype_quality_by_depth.cpp
    core/csr/measures/mapping_quality_zero_count.hpp
    core/csr/measures/mapping_quality_zero_count.cpp
    core/csr/measures/mean_mapping_quality.hpp
    core/csr/measures/mean_mapping_quality.cpp
    core/csr/measures/model_posterior.hpp
    core/csr/measures/model_posterior.cpp
    core/csr/measures/model_posterior_by_depth.hpp
    core/csr/measures/model_posterior_by_depth.cpp
    core/csr/measures/allele_depth.hpp
    core/csr/measures/allele_depth.cpp
    core/csr/measures/allele_frequency.hpp
    core/csr/measures/allele_frequency.cpp
    core/csr/measures/allele_frequency_bias.hpp
    core/csr/measures/allele_frequency_bias.cpp
    core/csr/measures/strand_bias.hpp
    core/csr/measures/strand_bias.cpp
    core/csr/measures/mapping_quality_divergence.hpp
    core/csr/measures/mapping_quality_divergence.cpp
    core/csr/measures/gc_content.hpp
    core/csr/measures/gc_content.cpp
    core/csr/measures/filtered_read_fraction.hpp
    core/csr/measures/filtered_read_fraction.cpp
    core/csr/measures/clipped_read_fraction.hpp
    core/csr/measures/clipped_read_fraction.cpp
    core/csr/measures/is_denovo.hpp
    core/csr/measures/is_denovo.cpp
    core/csr/measures/is_somatic.hpp
    core/csr/measures/is_somatic.cpp
    core/csr/measures/measures_fwd.hpp
    core/csr/measures/measure_factory.hpp
    core/csr/measures/measure_factory.cpp
    core/csr/measures/ambiguous_read_fraction.hpp
    core/csr/measures/ambiguous_read_fraction.cpp
    core/csr/measures/median_base_quality.hpp
    core/csr/measures/median_base_quality.cpp
    core/csr/measures/mismatch_count.hpp
    core/csr/measures/mismatch_count.cpp
    core/csr/measures/mismatch_fraction.hpp
    core/csr/measures/mismatch_fraction.cpp
    core/csr/measures/is_refcall.hpp
    core/csr/measures/is_refcall.cpp
    core/csr/measures/normal_contamination.hpp
    core/csr/measures/normal_contamination.cpp
    core/csr/measures/denovo_contamination.hpp
    core/csr/measures/denovo_contamination.cpp
    core/csr/measures/read_side_bias.hpp
    core/csr/measures/read_side_bias.cpp
    core/csr/measures/alt_allele_count.hpp
    core/csr/measures/alt_allele_count.cpp
    core/csr/measures/str_length.hpp
    core/csr/measures/str_length.cpp
    core/csr/measures/str_period.hpp
    core/csr/measures/str_period.cpp
    core/csr/measures/posterior_probability.hpp
    core/csr/measures/posterior_probability.cpp
    core/csr/measures/posterior_probability_by_depth.hpp
    core/csr/measures/posterior_probability_by_depth.cpp
    core/csr/measures/classification_confidence.hpp
    core/csr/measures/classification_confidence.cpp
    core/csr/measures/somatic_haplotype_count.hpp
    core/csr/measures/somatic_haplotype_count.cpp
    core/csr/measures/median_somatic_mapping_quality.hpp
    core/csr/measures/median_somatic_mapping_quality.cpp
    core/csr/measures/strand_disequilibrium.hpp
    core/csr/measures/strand_disequilibrium.cpp
    core/csr/measures/supplementary_fraction.hpp
    core/csr/measures/supplementary_fraction.cpp
    core/csr/measures/misaligned_read_count.hpp
    core/csr/measures/misaligned_read_count.cpp
    core/csr/measures/read_tail_bias.hpp
    core/csr/measures/read_tail_bias.cpp
    core/csr/measures/read_end_bias.hpp
    core/csr/measures/read_end_bias.cpp
    core/csr/measures/variant_length.hpp
    core/csr/measures/variant_length.cpp
    core/csr/measures/base_mismatch_count.hpp
    core/csr/measures/base_mismatch_count.cpp
    core/csr/measures/base_mismatch_fraction.hpp
    core/csr/measures/base_mismatch_fraction.cpp
    core/csr/measures/base_mismatch_quality.hpp
    core/csr/measures/base_mismatch_quality.cpp
    core/csr/measures/assigned_depth.hpp
    core/csr/measures/assigned_depth.cpp
    core/csr/measures/duplicate_concordance.hpp
    core/csr/measures/duplicate_concordance.cpp
    core/csr/measures/duplicate_allele_depth.hpp
    core/csr/measures/duplicate_allele_depth.cpp
    core/csr/measures/duplicate_allele_fraction.hpp
    core/csr/measures/duplicate_allele_fraction.cpp
    core/csr/measures/error_rate.hpp
    core/csr/measures/error_rate.cpp
    core/csr/measures/error_rate_stdev.hpp
    core/csr/measures/error_rate_stdev.cpp
    core/csr/measures/is_transversion.hpp
    core/csr/measures/is_transversion.cpp
    core/csr/measures/phase_length.hpp
    core/csr/measures/phase_length.cpp
    core/csr/measures/max_read_length.hpp
    core/csr/measures/max_read_length.cpp
    core/csr/measures/mean_likelihood.hpp
    core/csr/measures/mean_likelihood.cpp
    core/csr/measures/allele_mapping_quality.hpp
    core/csr/measures/allele_mapping_quality.cpp
    core/csr/measures/phylogeny_posterior.hpp
    core/csr/measures/phylogeny_posterior.cpp
    
    core/models/haplotype_likelihood_array.hpp
    core/models/haplotype_likelihood_array.cpp
    core/models/haplotype_likelihood_model.hpp
    core/models/haplotype_likelihood_model.cpp

    core/models/genotype/subclone_model.hpp
    core/models/genotype/subclone_model.cpp
    core/models/genotype/constant_mixture_genotype_likelihood_model.hpp
    core/models/genotype/constant_mixture_genotype_likelihood_model.cpp
    core/models/genotype/individual_model.hpp
    core/models/genotype/individual_model.cpp
    core/models/genotype/independent_population_model.hpp
    core/models/genotype/independent_population_model.cpp
    core/models/genotype/population_model.hpp
    core/models/genotype/population_model.cpp
    core/models/genotype/variational_bayes_mixture_model.hpp
    core/models/genotype/trio_model.hpp
    core/models/genotype/trio_model.cpp
    core/models/genotype/genotype_prior_model.hpp
    core/models/genotype/uniform_genotype_prior_model.hpp
    core/models/genotype/coalescent_genotype_prior_model.hpp
    core/models/genotype/cancer_genotype_prior_model.hpp
    core/models/genotype/cancer_genotype_prior_model.cpp
    core/models/genotype/population_prior_model.hpp
    core/models/genotype/uniform_population_prior_model.hpp
    core/models/genotype/coalescent_population_prior_model.hpp
    core/models/genotype/coalescent_population_prior_model.cpp
    core/models/genotype/hardy_weinberg_model.hpp
    core/models/genotype/hardy_weinberg_model.cpp
    core/models/genotype/variable_mixture_genotype_likelihood_model.hpp
    core/models/genotype/variable_mixture_genotype_likelihood_model.cpp
    core/models/genotype/variational_bayes_mixture_mixture_model.hpp
    core/models/genotype/variational_bayes_mixture_mixture_model.cpp
    core/models/genotype/single_cell_prior_model.hpp
    core/models/genotype/single_cell_prior_model.cpp
    core/models/genotype/single_cell_model.hpp
    core/models/genotype/single_cell_model.cpp
    core/models/genotype/haplogroup_prior_model.hpp
    core/models/genotype/haplogroup_prior_model.cpp

    core/models/pairhmm/pair_hmm.hpp
    core/models/pairhmm/simd_pair_hmm.hpp
    core/models/pairhmm/rolling_initializer.hpp
    core/models/pairhmm/sse2_pair_hmm_impl.hpp
    core/models/pairhmm/avx2_pair_hmm_impl.hpp
    core/models/pairhmm/avx512_pair_hmm_impl.hpp
    core/models/pairhmm/simd_pair_hmm_factory.hpp
    core/models/pairhmm/simd_pair_hmm_wrapper.hpp

    core/models/error/indel_error_model.hpp
    core/models/error/indel_error_model.cpp
    core/models/error/repeat_based_indel_error_model.hpp
    core/models/error/repeat_based_indel_error_model.cpp
    core/models/error/repeat_based_snv_error_model.hpp
    core/models/error/repeat_based_snv_error_model.cpp
    core/models/error/snv_error_model.hpp
    core/models/error/snv_error_model.cpp
    core/models/error/error_model_factory.hpp
    core/models/error/error_model_factory.cpp
    core/models/error/basic_repeat_based_indel_error_model.hpp
    core/models/error/basic_repeat_based_indel_error_model.cpp
    core/models/error/custom_repeat_based_indel_error_model.hpp
    core/models/error/custom_repeat_based_indel_error_model.cpp

    core/models/mutation/somatic_mutation_model.hpp
    core/models/mutation/somatic_mutation_model.cpp
    core/models/mutation/coalescent_model.hpp
    core/models/mutation/coalescent_model.cpp
    core/models/mutation/denovo_model.hpp
    core/models/mutation/denovo_model.cpp
    core/models/mutation/indel_mutation_model.hpp
    core/models/mutation/indel_mutation_model.cpp

    core/models/reference/individual_reference_likelihood_model.hpp
    core/models/reference/individual_reference_likelihood_model.cpp

    core/tools/coretools.hpp
    core/tools/haplotype_filter.hpp
    core/tools/haplotype_filter.cpp
    core/tools/read_assigner.hpp
    core/tools/read_assigner.cpp
    core/tools/read_realigner.hpp
    core/tools/read_realigner.cpp
    core/tools/bam_realigner.hpp
    core/tools/bam_realigner.cpp
    core/tools/indel_profiler.hpp
    core/tools/indel_profiler.cpp
    core/tools/bad_region_detector.hpp
    core/tools/bad_region_detector.cpp

    core/tools/hapgen/genome_walker.hpp
    core/tools/hapgen/genome_walker.cpp
    core/tools/hapgen/haplotype_generator.hpp
    core/tools/hapgen/haplotype_generator.cpp
    core/tools/hapgen/haplotype_tree.hpp
    core/tools/hapgen/haplotype_tree.cpp
    
    core/tools/phaser/phaser.hpp
    core/tools/phaser/phaser.cpp

    core/tools/vargen/cigar_scanner.hpp
    core/tools/vargen/cigar_scanner.cpp
    core/tools/vargen/downloader.hpp
    core/tools/vargen/downloader.cpp
    core/tools/vargen/local_reassembler.hpp
    core/tools/vargen/local_reassembler.cpp
    core/tools/vargen/randomiser.hpp
    core/tools/vargen/randomiser.cpp
    core/tools/vargen/variant_generator.hpp
    core/tools/vargen/variant_generator.cpp
    core/tools/vargen/vcf_extractor.hpp
    core/tools/vargen/vcf_extractor.cpp
    core/tools/vargen/variant_generator_builder.hpp
    core/tools/vargen/variant_generator_builder.cpp
    core/tools/vargen/active_region_generator.hpp
    core/tools/vargen/active_region_generator.cpp
    core/tools/vargen/repeat_scanner.hpp
    core/tools/vargen/repeat_scanner.cpp
    
    core/tools/vargen/utils/assembler.hpp
    core/tools/vargen/utils/assembler.cpp
    core/tools/vargen/utils/global_aligner.hpp
    core/tools/vargen/utils/global_aligner.cpp
    core/tools/vargen/utils/assembler_active_region_generator.hpp
    core/tools/vargen/utils/assembler_active_region_generator.cpp
    core/tools/vargen/utils/misaligned_reads_detector.hpp
    core/tools/vargen/utils/misaligned_reads_detector.cpp

    core/types/allele.hpp
    core/types/allele.cpp
    core/types/cancer_genotype.hpp
    core/types/cancer_genotype.cpp
    core/types/genotype.hpp
    core/types/genotype.cpp
    core/types/haplotype.hpp
    core/types/haplotype.cpp
    core/types/variant.hpp
    core/types/variant.cpp
    core/types/phylogeny.hpp
    core/types/phylogeny.cpp
    core/types/indexed_haplotype.hpp
    core/types/shared_haplotype.hpp
    core/types/partitioned_genotype.hpp

    core/calling_components.hpp
    core/calling_components.cpp

    core/octopus.hpp
    core/octopus.cpp
)

set(MISC_SOURCES
    ${octopus_SOURCE_DIR}/src/timers.hpp
    ${octopus_SOURCE_DIR}/src/timers.cpp
)

set(OCTOPUS_SOURCES
    ${CONFIG_SOURCES}
    ${EXCEPTIONS_SOURCES}
    ${CONCEPTS_SOURCES}
    ${BASICS_SOURCES}
    ${CONTAINERS_SOURCES}
    ${LOGGING_SOURCES}
    ${IO_SOURCES}
    ${READPIPE_SOURCES}
    ${UTILS_SOURCES}
    ${CORE_SOURCES}
    ${MISC_SOURCES}
)

set(INCLUDE_SOURCES
    ${octopus_SOURCE_DIR}/lib/bioio.hpp
    ${octopus_SOURCE_DIR}/lib/tandem/tandem.hpp
    ${octopus_SOURCE_DIR}/lib/ksp/custom_dijkstra_call.hpp
    ${octopus_SOURCE_DIR}/lib/ksp/yen_ksp.hpp
    ${octopus_SOURCE_DIR}/lib/ranger/Forest.h
)

set(REQUIRED_BOOST_LIBRARIES
    system
    filesystem
    program_options
    date_time
    log_setup
    log
    iostreams
    timer
    thread
)

set(WarningIgnores
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-missing-braces
    -Wno-deprecated-copy
)

set(GCCWarningIgnores
    -Wno-maybe-uninitialized
    -Wno-noexcept-type
    )

# Compile options for all builds
add_compile_options(-Wall -Wextra -Werror ${WarningIgnores})

if(CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(${GCCWarningIgnores})
endif()

# find_package(SSE)
# set(AVX512_FOUND false)
# if (AVX512F_FOUND AND AVX512BW_FOUND)
#     add_compile_options(-mavx512f -mavx512bw)
# 	set(AVX512_FOUND true)
# elseif (AVX2_FOUND)
#     add_compile_options(-mavx2)
# elseif (SSE4_1_FOUND)
#     add_compile_options(-msse4.1)
# elseif (SSSE3_FOUND)
#     add_compile_options(-msse3)
# else ()
#     add_compile_options(-msse2)
# endif()

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)

if (NOT BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(HTSlib_USE_STATIC_LIBS ON)
endif()

if (CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
set(CXX_OPTIMIZATION_FLAGS -ffast-math -mcpu=${COMPILER_ARCHITECTURE})
else()
set(CXX_OPTIMIZATION_FLAGS -ffast-math -march=${COMPILER_ARCHITECTURE})
endif()
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CXX_OPTIMIZATION_FLAGS ${CXX_OPTIMIZATION_FLAGS} -mfpmath=both)
endif()

if (BUILD_TESTING)
    # Make a library of all octopus non-main.cpp sources so can be used with tests
    add_library(Octopus ${OCTOPUS_SOURCES} ${INCLUDE_SOURCES})
    target_compile_features(Octopus PRIVATE cxx_thread_local)
    target_compile_options(Octopus PRIVATE ${CXX_OPTIMIZATION_FLAGS})
    target_include_directories(Octopus PUBLIC ${octopus_SOURCE_DIR}/lib ${octopus_SOURCE_DIR}/src)
    target_link_libraries(Octopus tandem ranger date-tz)
    target_compile_definitions(Octopus PRIVATE -DBOOST_LOG_DYN_LINK) # Required for log
    find_package (Boost 1.65 REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBRARIES} REQUIRED)
    message(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    message(STATUS "Boost libraries: " ${Boost_LIBRARIES})
    find_package(GMP REQUIRED)
    message(STATUS "GMP include dir: " ${GMP_INCLUDES})
    message(STATUS "GMP libraries: " ${GMP_LIBRARIES})
    find_package (HTSlib 1.4 REQUIRED)
    find_package(Threads REQUIRED)
    target_include_directories (Octopus PUBLIC
        ${octopus_SOURCE_DIR}/lib
        ${octopus_SOURCE_DIR}/src
        ${Boost_INCLUDE_DIR}
        ${GMP_INCLUDES}
        ${HTSlib_INCLUDE_DIRS})
    target_link_libraries (Octopus tandem ranger ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${HTSlib_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    check_ipo_supported(RESULT ipo_supported OUTPUT output)
    if(ipo_supported)
        message(STATUS "IPO is supported!")
        set_property(TARGET Octopus PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "IPO is not supported: ${output}")
    endif()
    configure_file (
        "${PROJECT_SOURCE_DIR}/src/config/system.h.in"
        "${PROJECT_BINARY_DIR}/generated/system.hpp"
    )
elseif (CMAKE_BUILD_TYPE MATCHES Debug)
    add_executable(octopus-debug main.cpp ${OCTOPUS_SOURCES} ${INCLUDE_SOURCES})
    target_compile_features(octopus-debug PRIVATE cxx_thread_local)
    target_include_directories(octopus-debug PUBLIC ${octopus_SOURCE_DIR}/lib ${octopus_SOURCE_DIR}/src)
    target_link_libraries(octopus-debug tandem ranger date-tz)
    if (NOT BUILD_SHARED_LIBS)
        message(STATUS "Linking against boost static libraries")
        set(Boost_USE_STATIC_LIBS ON)
    else()
        message(STATUS "Linking against boost dynamic libraries")
        target_compile_definitions(octopus-debug PRIVATE -DBOOST_LOG_DYN_LINK) # Required for log
    endif()
    set(Boost_NO_BOOST_CMAKE TRUE)
    find_package (Boost 1.65 REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBRARIES} REQUIRED)
    message(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    message(STATUS "Boost libraries: " ${Boost_LIBRARIES})
    find_package(GMP REQUIRED)
    message(STATUS "GMP include dir: " ${GMP_INCLUDES})
    message(STATUS "GMP libraries: " ${GMP_LIBRARIES})
    find_package (HTSlib 1.4 REQUIRED)
    find_package(Threads REQUIRED)
    target_include_directories (octopus-debug PUBLIC
        ${octopus_SOURCE_DIR}/lib
        ${octopus_SOURCE_DIR}/src
        ${Boost_INCLUDE_DIR}
        ${GMP_INCLUDES}
        ${HTSlib_INCLUDE_DIRS})
    target_link_libraries (octopus-debug tandem ranger date-tz ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${HTSlib_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    configure_file (
        "${PROJECT_SOURCE_DIR}/src/config/system.h.in"
        "${PROJECT_BINARY_DIR}/generated/system.hpp"
    )
    install(TARGETS octopus-debug DESTINATION ${CMAKE_INSTALL_PREFIX})
elseif (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    set(SanitizeFlags -fsanitize=address,undefined)
    add_executable(octopus-sanitize main.cpp ${OCTOPUS_SOURCES} ${INCLUDE_SOURCES})
    target_compile_features(octopus-sanitize PRIVATE cxx_thread_local)
    target_compile_options(octopus-sanitize PRIVATE ${SanitizeFlags} -fno-omit-frame-pointer)
    if (NOT BUILD_SHARED_LIBS)
        message(STATUS "Linking against boost static libraries")
        set(Boost_USE_STATIC_LIBS ON)
    else()
        message(STATUS "Linking against boost dynamic libraries")
        target_compile_definitions(octopus-sanitize PRIVATE -DBOOST_LOG_DYN_LINK) # Required for log
    endif()
    set(Boost_NO_BOOST_CMAKE TRUE)
    find_package (Boost 1.65 REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBRARIES} REQUIRED)
    message(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    message(STATUS "Boost libraries: " ${Boost_LIBRARIES})
    find_package(GMP REQUIRED)
    message(STATUS "GMP include dir: " ${GMP_INCLUDES})
    message(STATUS "GMP libraries: " ${GMP_LIBRARIES})
    find_package (HTSlib 1.4 REQUIRED)
    find_package(Threads REQUIRED)
    target_include_directories (octopus-sanitize PUBLIC
        ${octopus_SOURCE_DIR}/lib
        ${octopus_SOURCE_DIR}/src
        ${Boost_INCLUDE_DIR}
        ${GMP_INCLUDES}
        ${HTSlib_INCLUDE_DIRS})
    target_link_libraries (octopus-sanitize tandem ranger date-tz ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${HTSlib_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${SanitizeFlags})
    configure_file (
        "${PROJECT_SOURCE_DIR}/src/config/system.h.in"
        "${PROJECT_BINARY_DIR}/generated/system.hpp"
    )
    install(TARGETS octopus-sanitize DESTINATION ${CMAKE_INSTALL_PREFIX})
else()
    add_executable(octopus main.cpp ${OCTOPUS_SOURCES} ${INCLUDE_SOURCES})
    target_compile_features(octopus PRIVATE cxx_thread_local)
    target_compile_options(octopus PRIVATE ${CXX_OPTIMIZATION_FLAGS})
    if (BUILD_SHARED_LIBS)
        message(STATUS "Linking against boost dynamic libraries")
        target_compile_definitions(octopus PRIVATE -DBOOST_LOG_DYN_LINK) # Required for log
    else()
        message(STATUS "Linking against boost static libraries")
        set(Boost_USE_STATIC_LIBS ON)
    endif()
    set(Boost_NO_BOOST_CMAKE TRUE)
    find_package (Boost 1.65 REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBRARIES} REQUIRED)
    message(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    message(STATUS "Boost libraries: " ${Boost_LIBRARIES})
    find_package(GMP REQUIRED)
    message(STATUS "GMP include dir: " ${GMP_INCLUDES})
    message(STATUS "GMP libraries: " ${GMP_LIBRARIES})
    find_package (HTSlib 1.4 REQUIRED)
    find_package(Threads REQUIRED)
    target_include_directories (octopus PUBLIC
        ${octopus_SOURCE_DIR}/lib
        ${octopus_SOURCE_DIR}/src
        ${Boost_INCLUDE_DIR}
        ${GMP_INCLUDES}
        ${HTSlib_INCLUDE_DIRS})
    target_link_libraries (octopus tandem ranger date-tz ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${HTSlib_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    if (NOT BUILD_SHARED_LIBS)
        if (CMAKE_COMPILER_IS_GNUCXX)
            target_link_libraries(octopus m mvec -static -static-libgcc -static-libstdc++)
        endif()
    endif()
    check_ipo_supported(RESULT ipo_supported OUTPUT output)
    if(ipo_supported)
        message(STATUS "IPO is supported!")
        set_property(TARGET octopus PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        target_compile_options(octopus PRIVATE -ipo)
    else()
        message(WARNING "IPO is not supported: ${output}")
    endif()
    configure_file (
        "${PROJECT_SOURCE_DIR}/src/config/system.h.in"
        "${PROJECT_BINARY_DIR}/generated/system.hpp"
    )
    install(TARGETS octopus DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()
