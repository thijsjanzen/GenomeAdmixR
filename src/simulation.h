//
//  simulation.h
//
//  Generic Wright-Fisher simulation of one population, or of two populations
//  connected by migration.
//
//  Everything that depends on how a genome is stored (as ancestry blocks or as
//  sequence data) is handled by a Model, see models.h. A Model provides:
//
//    using individual_type = ...;
//
//    individual_type reproduce(const individual_type& first_parent,
//                              const individual_type& second_parent,
//                              rnd_t& random_generator) const;
//    bool use_selection() const;
//    double fitness(const individual_type& individual) const;
//    int number_of_alleles() const;
//    arma::mat allele_frequencies(const std::vector<individual_type>& population,
//                                 const std::vector<double>& markers,
//                                 int generation) const;
//        returns number_of_markers * number_of_alleles rows with the columns
//        time, location, allele and frequency. Markers with a negative
//        location are not tracked and give rows of zeros.
//    double count_junctions(const std::vector<individual_type>& population) const;
//    bool may_stop_when_fixed() const;
//    Rcpp::List to_list(const std::vector<individual_type>& population) const;
//

#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>

#include <RcppArmadillo.h>
#include "random_functions.h"
#include "helper_functions.h"
#include "util.h"


// Returns a copy of 'source_population' with 'target_size' individuals. If the
// sizes differ, the individuals are drawn at random (with replacement).
// Written by Claude (Anthropic).
template <class Individual>
std::vector<Individual> resize_population(
    const std::vector<Individual>& source_population,
    size_t target_size,
    rnd_t& random_generator) {

  if (source_population.size() == target_size) {
    return source_population;
  }
  if (source_population.empty()) {
    throw std::runtime_error("can not seed a population from an empty input");
  }

  const int source_size = static_cast<int>(source_population.size());
  std::vector<Individual> resized_population;
  resized_population.reserve(target_size);
  for (size_t i = 0; i < target_size; ++i) {
    const int drawn_index = random_generator.random_number(source_size);
    resized_population.push_back(source_population[drawn_index]);
  }
  return resized_population;
}


// A population is fixed when every chromosome of every individual is the same.
// Written by Claude (Anthropic).
template <class Individual>
bool is_fixed(const std::vector<Individual>& population) {
  if (population.empty()) {
    return true;
  }
  const auto& reference_chromosome = population.front().chromosome1;
  for (const auto& individual : population) {
    const bool first_matches  = individual.chromosome1 == reference_chromosome;
    const bool second_matches = individual.chromosome2 == reference_chromosome;
    if (!first_matches || !second_matches) {
      return false;
    }
  }
  return true;
}


// One population, together with the fitness of its individuals.
// Written by Claude (Anthropic), including all member functions.
template <class Model>
struct PopulationState {
  using individual_type = typename Model::individual_type;

  std::vector<individual_type> individuals;
  std::vector<double> fitness;
  double maximum_fitness = 0.0;

  void update_fitness(const Model& model) {
    if (!model.use_selection()) {
      return;
    }
    fitness.resize(individuals.size());
    for (size_t i = 0; i < individuals.size(); ++i) {
      fitness[i] = model.fitness(individuals[i]);
    }
    maximum_fitness = *std::max_element(fitness.begin(), fitness.end());
  }

  size_t number_of_fit_individuals() const {
    return std::count_if(fitness.begin(), fitness.end(),
                         [](double value) { return value > 0.0; });
  }

  // Index of a random parent: proportional to fitness under selection,
  // uniform otherwise.
  int draw_parent_index(bool use_selection, rnd_t& random_generator) const {
    if (use_selection) {
      return draw_prop_fitness(fitness, maximum_fitness, random_generator);
    }
    return random_generator.random_number(static_cast<int>(individuals.size()));
  }
};


// Draws parents for the offspring of one population. With probability
// 'migration_rate' a parent comes from the other population (if any).
// Written by Claude (Anthropic), including all member functions.
template <class Model>
class ParentSampler {
 public:
  using individual_type = typename Model::individual_type;

  ParentSampler(const PopulationState<Model>& resident_population,
                const PopulationState<Model>* migrant_population,
                double migration_rate,
                bool use_selection) :
    resident_population_(resident_population),
    migrant_population_(migrant_population),
    migration_rate_(migration_rate),
    use_selection_(use_selection) {
  }

  // Returns the parent. 'parent_id' identifies the parent across both
  // populations: migrants are numbered after the residents.
  const individual_type& draw_parent(int& parent_id,
                                     rnd_t& random_generator) const {
    if (draw_migrant(random_generator)) {
      const int migrant_index =
        migrant_population_->draw_parent_index(use_selection_, random_generator);
      parent_id = number_of_residents() + migrant_index;
      return migrant_population_->individuals[migrant_index];
    }
    const int resident_index =
      resident_population_.draw_parent_index(use_selection_, random_generator);
    parent_id = resident_index;
    return resident_population_.individuals[resident_index];
  }

  size_t number_of_possible_parents() const {
    size_t number_of_parents = resident_population_.individuals.size();
    if (migrant_population_ != nullptr && migration_rate_ > 0.0) {
      number_of_parents += migrant_population_->individuals.size();
    }
    return number_of_parents;
  }

 private:
  bool draw_migrant(rnd_t& random_generator) const {
    if (migrant_population_ == nullptr) {
      return false;  // a single population: no random number is used
    }
    return random_generator.uniform() < migration_rate_;
  }

  int number_of_residents() const {
    return static_cast<int>(resident_population_.individuals.size());
  }

  const PopulationState<Model>& resident_population_;
  const PopulationState<Model>* migrant_population_;
  double migration_rate_;
  bool use_selection_;
};


// Produces 'offspring_count' offspring, each from two different parents.
// The work is spread over 'num_threads' threads; 'seed' determines the result.
// Written by Claude (Anthropic).
template <class Model>
std::vector<typename Model::individual_type> produce_offspring(
    const Model& model,
    const ParentSampler<Model>& parent_sampler,
    size_t offspring_count,
    int num_threads,
    uint64_t seed) {

  using individual_type = typename Model::individual_type;

  if (parent_sampler.number_of_possible_parents() < 2) {
    Rcpp::stop("need at least two individuals to produce offspring");
  }

  std::vector<individual_type> offspring(offspring_count);

  auto produce_one_offspring = [&](size_t offspring_index,
                                   rnd_t& random_generator) {
    int first_parent_id = -1;
    int second_parent_id = -1;
    const individual_type* first_parent =
      &parent_sampler.draw_parent(first_parent_id, random_generator);
    const individual_type* second_parent =
      &parent_sampler.draw_parent(second_parent_id, random_generator);
    while (second_parent_id == first_parent_id) {
      second_parent = &parent_sampler.draw_parent(second_parent_id,
                                                  random_generator);
    }
    offspring[offspring_index] = model.reproduce(*first_parent,
                                                 *second_parent,
                                                 random_generator);
  };

  parallel_blocks(offspring_count, num_threads, seed, produce_one_offspring);
  return offspring;
}


// Written by Claude (Anthropic).
struct SimulationSettings {
  std::vector<size_t> population_sizes;   // one entry per population
  double migration_rate = 0.0;
  int total_runtime = 0;
  bool track_frequency = false;
  std::vector<double> tracked_markers;    // recorded every generation
  bool track_junctions = false;
  std::vector<double> summary_markers;    // recorded at the start and the end
  bool verbose = false;
  int num_threads = 1;
};


// Written by Claude (Anthropic).
template <class Model>
struct SimulationResult {
  std::vector< std::vector<typename Model::individual_type> > populations;
  arma::mat frequencies;          // tracked markers, every generation
  arma::mat initial_frequencies;
  arma::mat final_frequencies;
  std::vector<double> junctions;  // mean number of junctions per chromosome
  int generations_run = 0;
};


// Runs a simulation of one population, or of two populations connected by
// migration, for settings.total_runtime generations, or until all
// populations are fixed.
// Written by Claude (Anthropic), including all member functions.
template <class Model>
class Simulation {
 public:
  using individual_type = typename Model::individual_type;

  Simulation(const Model& model,
             const SimulationSettings& settings,
             const std::vector< std::vector<individual_type> >& starting_populations,
             rnd_t& random_generator) :
    model_(model),
    settings_(settings),
    random_generator_(random_generator) {
    check_settings(starting_populations);
    for (const auto& starting_population : starting_populations) {
      PopulationState<Model> population;
      population.individuals = starting_population;
      populations_.push_back(population);
    }
  }

  SimulationResult<Model> run() {
    update_all_fitness();
    result_.initial_frequencies = frequencies_of_all_populations(
                                      settings_.summary_markers, 0);
    allocate_frequency_table();
    print_progress_header();

    for (int generation = 0; generation < settings_.total_runtime; ++generation) {
      record_junctions();
      record_tracked_frequencies(generation);

      replace_by_next_generation();
      update_all_fitness();
      result_.generations_run = generation + 1;

      print_progress(generation);
      if (generation > 2 && can_stop_because_fixed()) {
        report_fixation();
        break;
      }
      Rcpp::checkUserInterrupt();
    }
    if (settings_.verbose) Rcpp::Rcout << "\n";

    remove_unused_frequency_rows();
    result_.final_frequencies = frequencies_of_all_populations(
                                    settings_.summary_markers,
                                    result_.generations_run);
    move_populations_to_result();
    return result_;
  }

 private:
  void check_settings(
      const std::vector< std::vector<individual_type> >& starting_populations) const {
    const size_t number_of_populations = starting_populations.size();
    if (number_of_populations < 1 || number_of_populations > 2) {
      throw std::runtime_error("need one or two populations");
    }
    if (settings_.population_sizes.size() != number_of_populations) {
      throw std::runtime_error("need one population size per population");
    }
    for (const auto& starting_population : starting_populations) {
      if (starting_population.empty()) {
        throw std::runtime_error("empty population");
      }
    }
  }

  size_t number_of_populations() const {
    return populations_.size();
  }

  // ----- selection -----

  void update_all_fitness() {
    for (auto& population : populations_) {
      population.update_fitness(model_);
    }
    check_enough_fit_individuals();
  }

  // Without migration, selection has to leave at least two possible parents,
  // otherwise drawing two different parents would never finish.
  void check_enough_fit_individuals() const {
    if (!model_.use_selection() || number_of_populations() != 1) {
      return;
    }
    const auto& population = populations_.front();
    if (population.maximum_fitness > 0.0 &&
        population.number_of_fit_individuals() < 2) {
      Rcpp::stop("fewer than two individuals with positive fitness remain");
    }
  }

  // ----- reproduction -----

  void replace_by_next_generation() {
    std::vector< std::vector<individual_type> > next_generation;
    for (size_t index = 0; index < number_of_populations(); ++index) {
      next_generation.push_back(offspring_of_population(index));
    }
    for (size_t index = 0; index < number_of_populations(); ++index) {
      populations_[index].individuals.swap(next_generation[index]);
    }
  }

  std::vector<individual_type> offspring_of_population(size_t index) {
    const PopulationState<Model>* migrant_population = nullptr;
    if (number_of_populations() == 2) {
      migrant_population = &populations_[1 - index];
    }
    ParentSampler<Model> parent_sampler(populations_[index],
                                        migrant_population,
                                        settings_.migration_rate,
                                        model_.use_selection());
    const uint64_t seed = random_generator_.next_seed();
    return produce_offspring(model_,
                             parent_sampler,
                             settings_.population_sizes[index],
                             settings_.num_threads,
                             seed);
  }

  // ----- fixation -----

  bool can_stop_because_fixed() const {
    if (!model_.may_stop_when_fixed()) {
      return false;
    }
    for (const auto& population : populations_) {
      if (!is_fixed(population.individuals)) {
        return false;
      }
    }
    return true;
  }

  void report_fixation() const {
    if (settings_.verbose) {
      Rcpp::Rcout << "\n After " << result_.generations_run
                  << " generations, the population has become completely"
                  << " homozygous and fixed\n";
    }
    R_FlushConsole();
  }

  // ----- junctions -----

  void record_junctions() {
    if (!settings_.track_junctions) {
      return;
    }
    double total_number_of_junctions = 0.0;
    double total_number_of_chromosomes = 0.0;
    for (const auto& population : populations_) {
      total_number_of_junctions += model_.count_junctions(population.individuals);
      total_number_of_chromosomes += 2.0 * population.individuals.size();
    }
    result_.junctions.push_back(total_number_of_junctions /
                                total_number_of_chromosomes);
  }

  // ----- allele frequencies -----

  arma::uword rows_per_generation() const {
    return settings_.tracked_markers.size() *
           model_.number_of_alleles() *
           number_of_populations();
  }

  arma::uword number_of_frequency_columns() const {
    // time, location, allele, frequency (and population, if there are two)
    return number_of_populations() == 1 ? 4 : 5;
  }

  void allocate_frequency_table() {
    if (!settings_.track_frequency) {
      return;
    }
    result_.frequencies.zeros(rows_per_generation() * settings_.total_runtime,
                              number_of_frequency_columns());
  }

  void record_tracked_frequencies(int generation) {
    if (!settings_.track_frequency) {
      return;
    }
    const arma::mat frequencies = frequencies_of_all_populations(
                                      settings_.tracked_markers, generation);
    const arma::uword first_row = generation * rows_per_generation();
    const arma::uword last_row = first_row + rows_per_generation() - 1;
    result_.frequencies.rows(first_row, last_row) = frequencies;
  }

  // After early fixation, the rows of the generations that were not
  // simulated are removed.
  void remove_unused_frequency_rows() {
    if (!settings_.track_frequency) {
      return;
    }
    result_.frequencies.resize(rows_per_generation() * result_.generations_run,
                               number_of_frequency_columns());
  }

  arma::mat frequencies_of_all_populations(const std::vector<double>& markers,
                                           int generation) const {
    if (number_of_populations() == 1) {
      return model_.allele_frequencies(populations_.front().individuals,
                                       markers, generation);
    }
    arma::mat all_frequencies;
    for (size_t index = 0; index < number_of_populations(); ++index) {
      const arma::mat frequencies = frequencies_with_population_column(
                                        index, markers, generation);
      all_frequencies = arma::join_cols(all_frequencies, frequencies);
    }
    return all_frequencies;
  }

  // Frequencies of one population, with the population number (1 or 2) in a
  // fifth column. Rows of untracked (negative) markers stay zero.
  arma::mat frequencies_with_population_column(size_t population_index,
                                               const std::vector<double>& markers,
                                               int generation) const {
    const arma::mat frequencies = model_.allele_frequencies(
                         populations_[population_index].individuals,
                         markers, generation);
    arma::mat output(frequencies.n_rows, 5, arma::fill::zeros);
    output.cols(0, 3) = frequencies;

    const arma::uword number_of_alleles = model_.number_of_alleles();
    const double population_number = population_index + 1.0;
    for (size_t marker_index = 0; marker_index < markers.size(); ++marker_index) {
      if (markers[marker_index] < 0) {
        continue;
      }
      const arma::uword first_row = marker_index * number_of_alleles;
      const arma::uword last_row = first_row + number_of_alleles - 1;
      output.submat(first_row, 4, last_row, 4).fill(population_number);
    }
    return output;
  }

  // ----- progress -----

  void print_progress_header() const {
    if (!settings_.verbose) {
      return;
    }
    Rcpp::Rcout << "0--------25--------50--------75--------100\n";
    Rcpp::Rcout << "*";
    R_FlushConsole();
  }

  void print_progress(int generation) const {
    const int generations_per_star = std::max(settings_.total_runtime / 20, 1);
    if (settings_.verbose && generation % generations_per_star == 0) {
      Rcpp::Rcout << "**";
    }
  }

  // ----- output -----

  void move_populations_to_result() {
    for (auto& population : populations_) {
      result_.populations.push_back(std::move(population.individuals));
    }
  }

  const Model& model_;
  const SimulationSettings& settings_;
  rnd_t& random_generator_;
  std::vector< PopulationState<Model> > populations_;
  SimulationResult<Model> result_;
};


// Written by Claude (Anthropic).
template <class Model>
SimulationResult<Model> run_simulation(
    const Model& model,
    const SimulationSettings& settings,
    const std::vector< std::vector<typename Model::individual_type> >& starting_populations,
    rnd_t& random_generator) {
  Simulation<Model> simulation(model, settings, starting_populations,
                               random_generator);
  return simulation.run();
}
