// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
// simulate_cpp
List simulate_cpp(Rcpp::NumericVector input_population, NumericMatrix select, int pop_size, int number_of_founders, Rcpp::NumericVector starting_proportions, int total_runtime, double morgan, bool progress_bar, bool track_frequency, NumericVector track_markers, bool track_junctions, bool multiplicative_selection);
RcppExport SEXP _GenomeAdmixR_simulate_cpp(SEXP input_populationSEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP starting_proportionsSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP progress_barSEXP, SEXP track_frequencySEXP, SEXP track_markersSEXP, SEXP track_junctionsSEXP, SEXP multiplicative_selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population(input_populationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type starting_proportions(starting_proportionsSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type track_markers(track_markersSEXP);
    Rcpp::traits::input_parameter< bool >::type track_junctions(track_junctionsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_cpp(input_population, select, pop_size, number_of_founders, starting_proportions, total_runtime, morgan, progress_bar, track_frequency, track_markers, track_junctions, multiplicative_selection));
    return rcpp_result_gen;
END_RCPP
}
// simulate_migration_cpp
List simulate_migration_cpp(Rcpp::NumericVector input_population_1, Rcpp::NumericVector input_population_2, NumericMatrix select, Rcpp::NumericVector pop_size, int number_of_founders, Rcpp::NumericVector starting_proportions, int total_runtime, double morgan, bool progress_bar, bool track_frequency, NumericVector track_markers, bool track_junctions, bool multiplicative_selection, double migration_rate);
RcppExport SEXP _GenomeAdmixR_simulate_migration_cpp(SEXP input_population_1SEXP, SEXP input_population_2SEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP starting_proportionsSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP progress_barSEXP, SEXP track_frequencySEXP, SEXP track_markersSEXP, SEXP track_junctionsSEXP, SEXP multiplicative_selectionSEXP, SEXP migration_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population_1(input_population_1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population_2(input_population_2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type starting_proportions(starting_proportionsSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type track_markers(track_markersSEXP);
    Rcpp::traits::input_parameter< bool >::type track_junctions(track_junctionsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    Rcpp::traits::input_parameter< double >::type migration_rate(migration_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_migration_cpp(input_population_1, input_population_2, select, pop_size, number_of_founders, starting_proportions, total_runtime, morgan, progress_bar, track_frequency, track_markers, track_junctions, multiplicative_selection, migration_rate));
    return rcpp_result_gen;
END_RCPP
}
// test_fish_functions
void test_fish_functions();
RcppExport SEXP _GenomeAdmixR_test_fish_functions() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_fish_functions();
    return R_NilValue;
END_RCPP
}
// simulate_migration_cpp
List simulate_migration_cpp(Rcpp::NumericVector input_population_1, Rcpp::NumericVector input_population_2, NumericMatrix select, Rcpp::NumericVector pop_size, int number_of_founders, Rcpp::NumericVector starting_proportions, int total_runtime, double morgan, bool progress_bar, bool track_frequency, NumericVector track_markers, bool track_junctions, bool multiplicative_selection, double migration_rate);
RcppExport SEXP _GenomeAdmixR_simulate_migration_cpp(SEXP input_population_1SEXP, SEXP input_population_2SEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP starting_proportionsSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP progress_barSEXP, SEXP track_frequencySEXP, SEXP track_markersSEXP, SEXP track_junctionsSEXP, SEXP multiplicative_selectionSEXP, SEXP migration_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population_1(input_population_1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population_2(input_population_2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type starting_proportions(starting_proportionsSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type track_markers(track_markersSEXP);
    Rcpp::traits::input_parameter< bool >::type track_junctions(track_junctionsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    Rcpp::traits::input_parameter< double >::type migration_rate(migration_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_migration_cpp(input_population_1, input_population_2, select, pop_size, number_of_founders, starting_proportions, total_runtime, morgan, progress_bar, track_frequency, track_markers, track_junctions, multiplicative_selection, migration_rate));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GenomeAdmixR_calculate_allele_spectrum_cpp", (DL_FUNC) &_GenomeAdmixR_calculate_allele_spectrum_cpp, 3},
    {"_GenomeAdmixR_simulate_cpp", (DL_FUNC) &_GenomeAdmixR_simulate_cpp, 12},
    {"_GenomeAdmixR_simulate_migration_cpp", (DL_FUNC) &_GenomeAdmixR_simulate_migration_cpp, 14},
    {"_GenomeAdmixR_test_fish_functions", (DL_FUNC) &_GenomeAdmixR_test_fish_functions, 0},
    {"_GenomeAdmixR_simulate_migration_cpp", (DL_FUNC) &_GenomeAdmixR_simulate_migration_cpp, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_GenomeAdmixR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
