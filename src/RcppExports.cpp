// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calculate_allele_spectrum_cpp
arma::mat calculate_allele_spectrum_cpp(Rcpp::NumericVector input_population, Rcpp::NumericVector markers, bool progress_bar);
RcppExport SEXP _GenomeAdmixR_calculate_allele_spectrum_cpp(SEXP input_populationSEXP, SEXP markersSEXP, SEXP progress_barSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population(input_populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type markers(markersSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_allele_spectrum_cpp(input_population, markers, progress_bar));
    return rcpp_result_gen;
END_RCPP
}
// calculate_heterozygosity_cpp
arma::mat calculate_heterozygosity_cpp(Rcpp::NumericVector input_population, Rcpp::NumericVector markers, bool progress_bar);
RcppExport SEXP _GenomeAdmixR_calculate_heterozygosity_cpp(SEXP input_populationSEXP, SEXP markersSEXP, SEXP progress_barSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population(input_populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type markers(markersSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_heterozygosity_cpp(input_population, markers, progress_bar));
    return rcpp_result_gen;
END_RCPP
}
// simulate_cpp
List simulate_cpp(Rcpp::NumericVector input_population, NumericMatrix select, size_t pop_size, size_t number_of_founders, Rcpp::NumericVector starting_proportions, size_t total_runtime, double morgan, bool verbose, bool track_frequency, NumericVector track_markers, bool track_junctions, bool multiplicative_selection, int seed, int num_threads);
RcppExport SEXP _GenomeAdmixR_simulate_cpp(SEXP input_populationSEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP starting_proportionsSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP verboseSEXP, SEXP track_frequencySEXP, SEXP track_markersSEXP, SEXP track_junctionsSEXP, SEXP multiplicative_selectionSEXP, SEXP seedSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population(input_populationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< size_t >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type starting_proportions(starting_proportionsSEXP);
    Rcpp::traits::input_parameter< size_t >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type track_markers(track_markersSEXP);
    Rcpp::traits::input_parameter< bool >::type track_junctions(track_junctionsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_cpp(input_population, select, pop_size, number_of_founders, starting_proportions, total_runtime, morgan, verbose, track_frequency, track_markers, track_junctions, multiplicative_selection, seed, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// simulate_emp_cpp
List simulate_emp_cpp(Rcpp::NumericMatrix input_population, Rcpp::NumericVector marker_positions_R, Rcpp::NumericMatrix select, size_t pop_size, size_t total_runtime, double morgan, bool verbose, bool track_frequency, Rcpp::NumericVector track_markers_R, bool multiplicative_selection, int seed, double mutation_rate, NumericMatrix sub_matrix);
RcppExport SEXP _GenomeAdmixR_simulate_emp_cpp(SEXP input_populationSEXP, SEXP marker_positions_RSEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP verboseSEXP, SEXP track_frequencySEXP, SEXP track_markers_RSEXP, SEXP multiplicative_selectionSEXP, SEXP seedSEXP, SEXP mutation_rateSEXP, SEXP sub_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type input_population(input_populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type marker_positions_R(marker_positions_RSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< size_t >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< size_t >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type track_markers_R(track_markers_RSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sub_matrix(sub_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_emp_cpp(input_population, marker_positions_R, select, pop_size, total_runtime, morgan, verbose, track_frequency, track_markers_R, multiplicative_selection, seed, mutation_rate, sub_matrix));
    return rcpp_result_gen;
END_RCPP
}
// simulate_migration_cpp
List simulate_migration_cpp(NumericVector input_population_1, NumericVector input_population_2, NumericMatrix select, NumericVector pop_size, NumericMatrix starting_frequencies, int total_runtime, double morgan, bool verbose, bool track_frequency, NumericVector track_markers, bool track_junctions, bool multiplicative_selection, double migration_rate, int seed);
RcppExport SEXP _GenomeAdmixR_simulate_migration_cpp(SEXP input_population_1SEXP, SEXP input_population_2SEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP starting_frequenciesSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP verboseSEXP, SEXP track_frequencySEXP, SEXP track_markersSEXP, SEXP track_junctionsSEXP, SEXP multiplicative_selectionSEXP, SEXP migration_rateSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input_population_1(input_population_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type input_population_2(input_population_2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type starting_frequencies(starting_frequenciesSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type track_markers(track_markersSEXP);
    Rcpp::traits::input_parameter< bool >::type track_junctions(track_junctionsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    Rcpp::traits::input_parameter< double >::type migration_rate(migration_rateSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_migration_cpp(input_population_1, input_population_2, select, pop_size, starting_frequencies, total_runtime, morgan, verbose, track_frequency, track_markers, track_junctions, multiplicative_selection, migration_rate, seed));
    return rcpp_result_gen;
END_RCPP
}
// simulate_migration_emp_cpp
List simulate_migration_emp_cpp(const NumericMatrix& input_population_1, const NumericMatrix& input_population_2, const NumericVector& marker_positions_R, NumericMatrix select, NumericVector pop_size, int total_runtime, double morgan, bool verbose, bool track_frequency, const NumericVector& track_markers_R, bool multiplicative_selection, double migration_rate, int seed, double mutation_rate, const NumericMatrix& substitution_matrix);
RcppExport SEXP _GenomeAdmixR_simulate_migration_emp_cpp(SEXP input_population_1SEXP, SEXP input_population_2SEXP, SEXP marker_positions_RSEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP verboseSEXP, SEXP track_frequencySEXP, SEXP track_markers_RSEXP, SEXP multiplicative_selectionSEXP, SEXP migration_rateSEXP, SEXP seedSEXP, SEXP mutation_rateSEXP, SEXP substitution_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type input_population_1(input_population_1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type input_population_2(input_population_2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type marker_positions_R(marker_positions_RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type track_markers_R(track_markers_RSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    Rcpp::traits::input_parameter< double >::type migration_rate(migration_rateSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type substitution_matrix(substitution_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_migration_emp_cpp(input_population_1, input_population_2, marker_positions_R, select, pop_size, total_runtime, morgan, verbose, track_frequency, track_markers_R, multiplicative_selection, migration_rate, seed, mutation_rate, substitution_matrix));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GenomeAdmixR_calculate_allele_spectrum_cpp", (DL_FUNC) &_GenomeAdmixR_calculate_allele_spectrum_cpp, 3},
    {"_GenomeAdmixR_calculate_heterozygosity_cpp", (DL_FUNC) &_GenomeAdmixR_calculate_heterozygosity_cpp, 3},
    {"_GenomeAdmixR_simulate_cpp", (DL_FUNC) &_GenomeAdmixR_simulate_cpp, 14},
    {"_GenomeAdmixR_simulate_emp_cpp", (DL_FUNC) &_GenomeAdmixR_simulate_emp_cpp, 13},
    {"_GenomeAdmixR_simulate_migration_cpp", (DL_FUNC) &_GenomeAdmixR_simulate_migration_cpp, 14},
    {"_GenomeAdmixR_simulate_migration_emp_cpp", (DL_FUNC) &_GenomeAdmixR_simulate_migration_emp_cpp, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_GenomeAdmixR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
