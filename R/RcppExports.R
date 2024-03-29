# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

calculate_allele_spectrum_cpp <- function(input_population, markers, progress_bar) {
    .Call(`_GenomeAdmixR_calculate_allele_spectrum_cpp`, input_population, markers, progress_bar)
}

simulation_data_to_genomeadmixr_data_cpp <- function(input_population, markers) {
    .Call(`_GenomeAdmixR_simulation_data_to_genomeadmixr_data_cpp`, input_population, markers)
}

simulation_data_to_plink_cpp <- function(input_population, markers) {
    .Call(`_GenomeAdmixR_simulation_data_to_plink_cpp`, input_population, markers)
}

calculate_heterozygosity_cpp <- function(input_population, markers, progress_bar) {
    .Call(`_GenomeAdmixR_calculate_heterozygosity_cpp`, input_population, markers, progress_bar)
}

vcf_to_matrix_cpp <- function(input_mat, allele_1, allele_2) {
    .Call(`_GenomeAdmixR_vcf_to_matrix_cpp`, input_mat, allele_1, allele_2)
}

simulate_cpp <- function(input_population, select, pop_size, number_of_founders, starting_proportions, total_runtime, morgan, verbose, track_frequency, track_markers, track_junctions, multiplicative_selection, num_threads) {
    .Call(`_GenomeAdmixR_simulate_cpp`, input_population, select, pop_size, number_of_founders, starting_proportions, total_runtime, morgan, verbose, track_frequency, track_markers, track_junctions, multiplicative_selection, num_threads)
}

simulate_emp_cpp <- function(input_population, marker_positions_R, select, pop_size, total_runtime, morgan, verbose, track_frequency, track_markers_R, multiplicative_selection, mutation_rate, substitution_matrix_R, num_threads, recombination_map) {
    .Call(`_GenomeAdmixR_simulate_emp_cpp`, input_population, marker_positions_R, select, pop_size, total_runtime, morgan, verbose, track_frequency, track_markers_R, multiplicative_selection, mutation_rate, substitution_matrix_R, num_threads, recombination_map)
}

simulate_migration_cpp <- function(input_population_1, input_population_2, select, pop_size, starting_frequencies, total_runtime, morgan, verbose, track_frequency, track_markers, track_junctions, multiplicative_selection, migration_rate, num_threads) {
    .Call(`_GenomeAdmixR_simulate_migration_cpp`, input_population_1, input_population_2, select, pop_size, starting_frequencies, total_runtime, morgan, verbose, track_frequency, track_markers, track_junctions, multiplicative_selection, migration_rate, num_threads)
}

simulate_migration_emp_cpp <- function(input_population_1, input_population_2, marker_positions_R, select, pop_sizes, total_runtime, morgan, verbose, track_frequency, track_markers_R, multiplicative_selection, migration_rate, mutation_rate, substitution_matrix_R, num_threads, recombination_map) {
    .Call(`_GenomeAdmixR_simulate_migration_emp_cpp`, input_population_1, input_population_2, marker_positions_R, select, pop_sizes, total_runtime, morgan, verbose, track_frequency, track_markers_R, multiplicative_selection, migration_rate, mutation_rate, substitution_matrix_R, num_threads, recombination_map)
}

