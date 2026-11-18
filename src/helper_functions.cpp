//
//  helper_functions.cpp
//
//
//  Created by Thijs Janzen on 20/11/2018.
//

#include "helper_functions.h"
#include <vector>
#include <array>
#include <string>
#include <thread>
#include <chrono>
#include <algorithm>
#include <stdexcept>

void force_output() {
//  std::this_thread::sleep_for(std::chrono::nanoseconds(100));
  std::this_thread::sleep_for(std::chrono::milliseconds(30));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}


// ancestry at position 'marker' (relative position in [0, 1])
// Ancestry at 'marker' (a relative position in [0, 1]), or -1 if the marker
// lies before the start of the chromosome.
int get_ancestry(const std::vector< junction >& chromosome,
                 double marker) {
  if (chromosome.empty() || marker < chromosome.front().pos) {
    return -1;
  }
  auto marker_is_left_of = [](double location, const junction& block_start) {
    return location < block_start.pos;
  };
  // the first junction to the right of the marker starts the next block
  auto next_block = std::upper_bound(chromosome.begin(), chromosome.end(),
                                     marker, marker_is_left_of);
  auto block_containing_marker = next_block - 1;
  return block_containing_marker->right;
}

int find_index(const std::vector<int>& v, int value) {
  for (size_t i = 0; i < v.size(); ++i) {
    if (v[i] == value) return i;
  }
  return -1;
}

void update_founder_labels(const std::vector<junction> chrom,
                           std::vector<int>& founder_labels) {
  for(auto i = chrom.begin(); i != chrom.end(); ++i) {
    if (founder_labels.empty()) {
      if((*i).right != -1) founder_labels.push_back((*i).right);
    } else {
      if(find_index(founder_labels, (*i).right) == -1) {
        if((*i).right != -1) founder_labels.push_back((*i).right);
      }
    }
  }
  return;
}

void update_anc_chrom(const std::vector<junction>& chromosome,
                      const std::vector<int>& founder_labels,
                      double marker,
                      arma::mat& allele_matrix) {
  int ancestry = get_ancestry(chromosome, marker);
  if (ancestry == -1 && !chromosome.empty()) {
    // location before the first site (sequence data): use the first site
    ancestry = chromosome.front().right;
  }
  const int founder_index = find_index(founder_labels, ancestry);
  if (founder_index < 0) {
    throw std::runtime_error("could not find ancestry at marker " +
                             std::to_string(marker));
  }
  allele_matrix(founder_index, 3)++;
}



arma::mat update_frequency_tibble(const std::vector< Fish >& v,
                                  double m,
                                  const std::vector<int>& founder_labels,
                                  int t,
                                  double morgan) {

  int num_alleles = founder_labels.size();
  arma::mat allele_matrix(num_alleles, 4);
  // initialize results
  for (int i = 0; i < num_alleles; ++i) {
    allele_matrix(i, 0) = t;
    allele_matrix(i, 1) = m * morgan;
    allele_matrix(i, 2) = founder_labels[i];
    allele_matrix(i, 3) = 0;
  }

 // for (auto it = v.begin(); it != v.end(); ++it) {
 for (int i = 0; i < v.size(); ++i) {
    update_anc_chrom(v[i].chromosome1, founder_labels,
                     m, allele_matrix);
    update_anc_chrom(v[i].chromosome2, founder_labels,
                     m, allele_matrix);
  }

  for (int i = 0; i < num_alleles; ++i) {
    allele_matrix(i, 3) *= 1.0 / (2 * v.size());
  }

  return(allele_matrix);
}

arma::mat update_all_frequencies_tibble(const std::vector< Fish >& pop,
                                        const std::vector<double>& markers,
                                        const std::vector<int>& founder_labels,
                                        int t,
                                        double morgan) {

  int number_of_alleles = founder_labels.size();
  arma::mat output(markers.size() * number_of_alleles, 4);

  for (size_t i = 0; i < markers.size(); ++i) {
     if (markers[i] < 0) continue;
     arma::mat local_mat = update_frequency_tibble(pop,
                                                   markers[i],
                                                   founder_labels,
                                                   t,
                                                   morgan);
    // now we have a (markers x alleles) x 3 tibble, e.g. [loc, anc, freq]
    // and we have to put that in the right place in the output matrix
    int start = i * number_of_alleles;
    int end = start + number_of_alleles;
    for (int j = start; j < end; ++j) {
      for (int k = 0; k < 4; ++k) {
        output(j, k) = local_mat(j - start, k);
      }
    }
  }
  return(output);
}



// total number of crossover junctions in a population. Every chromosome
// stores a start junction at 0.0 (and no end junction), so a chromosome
// carries size() - 1 crossover junctions.
double total_junctions(const std::vector< Fish >& population) {
  double number_of_junctions = 0.0;
  for (const auto& individual : population) {
    number_of_junctions += static_cast<double>(individual.chromosome1.size()) - 1.0;
    number_of_junctions += static_cast<double>(individual.chromosome2.size()) - 1.0;
  }
  return number_of_junctions;
}

int draw_prop_fitness(const std::vector<double>& fitness,
                      const double& maxFitness,
                      rnd_t& rndgen) {
  if (fitness.empty()) {
    throw std::runtime_error("Fitness is empty");
  }
  if (maxFitness <= 0.0) {
    return rndgen.random_number(fitness.size());
  }

  size_t fitness_size = fitness.size();
  double inv_fitness = 1.0 / maxFitness;
  while (true) {
    int index = rndgen.random_number(fitness_size);
    double prob = inv_fitness * fitness[index];
    if (rndgen.uniform() < prob) {
      return index;
    }
  }
}

std::vector< Fish > convert_NumericVector_to_fishVector(const NumericVector& v) {
  std::vector< Fish > output;

  Fish temp;
  int indic_chrom = 1;
  bool add_indiv = false;

  junction prev_j(-1, 0);

  for (int i = 0; i < v.size(); i += 2) {
    junction temp_j;
    temp_j.pos = v[i];
    temp_j.right = v[i+1];

    if (temp_j.pos <= prev_j.pos) {
      if(indic_chrom == 1) {
        indic_chrom = 2;
      } else {
        add_indiv = true;
      }
    }

    if (add_indiv) {
      output.push_back(temp);
      add_indiv = false;
      indic_chrom = 1;
      temp.chromosome1.clear();
      temp.chromosome2.clear();
    }

    if (indic_chrom == 1) {
      temp.chromosome1.push_back(temp_j);
    } else {
      temp.chromosome2.push_back(temp_j);
    }
    prev_j = temp_j;
  }

  // last individual is not added:
  output.push_back(temp);

  return(output);
}

List convert_to_list(const std::vector<Fish>& v) {
  int list_size = (int)v.size();
  List output(list_size);

  for (size_t i = 0; i < v.size(); ++i) {

    Fish focal = v[i];

    NumericMatrix chrom1(focal.chromosome1.size(), 2); // nrow = number of junctions, ncol = 2
    for (size_t j = 0; j < focal.chromosome1.size(); ++j) {
      chrom1(j, 0) = focal.chromosome1[j].pos;
      chrom1(j, 1) = focal.chromosome1[j].right;
    }

    NumericMatrix chrom2(focal.chromosome2.size(), 2); // nrow = number of junctions, ncol = 2
    for (size_t j = 0; j < focal.chromosome2.size(); ++j) {
      chrom2(j, 0) = focal.chromosome2[j].pos;
      chrom2(j, 1) = focal.chromosome2[j].right;
    }

    List toAdd = List::create( Named("chromosome1") = chrom1,
                               Named("chromosome2") = chrom2
    );

    output(i) = toAdd;
  }

  return output;
}

// select: one row per site; columns: location (relative, in [0, 1]),
// fitness aa, Aa, AA, ancestor A. Rows with a negative ancestor are ignored.
double calculate_fitness(const Fish& focal,
                         const std::vector< std::array<double, 5>>& select,
                         bool multiplicative_selection) {
  double fitness = multiplicative_selection ? 1.0 : 0.0;
  for (const auto& site : select) {
    const double location = site[0];
    const double selected_ancestor = site[4];
    if (selected_ancestor < 0) {
      continue;
    }
    const int ancestor = static_cast<int>(selected_ancestor);
    const bool first_copy  = get_ancestry(focal.chromosome1, location) == ancestor;
    const bool second_copy = get_ancestry(focal.chromosome2, location) == ancestor;
    const int number_of_copies = static_cast<int>(first_copy) +
                                 static_cast<int>(second_copy);
    // columns 1, 2 and 3 hold the fitness for 0, 1 and 2 copies
    const double site_fitness = site[1 + number_of_copies];
    if (multiplicative_selection) {
      fitness *= site_fitness;
    } else {
      fitness += site_fitness;
    }
  }
  return fitness;
}


int draw_random_founder(const NumericVector& v,
                        rnd_t& rndgen) {
  double r = rndgen.uniform();
  for (size_t i = 0; i < static_cast<size_t>(v.size()); ++i) {
    r -= v[i];
    if(r <= 0) {
      return(i);
    }
  }
  return(v.size() - 1);
}

// [[Rcpp::export]]
arma::mat calculate_allele_spectrum_cpp(Rcpp::NumericVector input_population,
                                        Rcpp::NumericVector markers,
                                        bool progress_bar) {
try {
  std::vector< Fish > Pop;

  Pop = convert_NumericVector_to_fishVector(input_population);
  std::vector<int> founder_labels;
  for (auto it = Pop.begin(); it != Pop.end(); ++it) {
    update_founder_labels((*it).chromosome1, founder_labels);
    update_founder_labels((*it).chromosome2, founder_labels);
  }

  // locations are relative positions on the chromosome, in [0, 1]
  arma::mat frequencies = update_all_frequencies_tibble(Pop,
                                                        Rcpp::as<std::vector<double>>(markers),
                                                        founder_labels,
                                                        0,
                                                        1.0);

  return frequencies;
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  Rcpp::stop("c++ exception (unknown reason)");
}
return arma::mat();
}

// [[Rcpp::export]]
NumericMatrix simulation_data_to_genomeadmixr_data_cpp(Rcpp::NumericVector input_population,
                                                       Rcpp::NumericVector markers) {

try {
  std::vector< Fish > Pop;
  Pop = convert_NumericVector_to_fishVector(input_population);

  NumericMatrix output(Pop.size() * 2, markers.size());

  for(size_t i = 0; i < Pop.size(); ++i) {
    int index_c1 = i * 2;
    int index_c2 = index_c1 + 1;

    for (size_t j = 0; j < static_cast<size_t>(markers.size()); ++j) {
      output(index_c1, j) = get_ancestry(Pop[i].chromosome1, markers[j]);
      output(index_c2, j) = get_ancestry(Pop[i].chromosome2, markers[j]);
    }
  }
  return(output);
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  Rcpp::stop("c++ exception (unknown reason)");
}
return NA_REAL;
}

std::string int_to_base(int a) {
  std::string  b1 = "0";
  if (a == 1) b1 = "a";
  if (a == 2) b1 = "c";
  if (a == 3) b1 = "t";
  if (a == 4) b1 = "g";
  return b1;
}

std::vector<std::string> combine_alleles(int a1, int a2) {

  std::vector<std::string> output = {"0", "0"};
  if (a1 != 0 && a2 != 0) {
    output[0] = int_to_base(a1);
    output[1] = int_to_base(a2);
  }
  return output;
}

// [[Rcpp::export]]
StringMatrix simulation_data_to_plink_cpp(Rcpp::NumericVector input_population,
                                          Rcpp::NumericVector markers) {

try {
  std::vector< Fish > Pop;
  Pop = convert_NumericVector_to_fishVector(input_population);

  StringMatrix output(Pop.size(), markers.size() * 2);
  for(size_t i = 0; i < static_cast<size_t>(Pop.size()); ++i) {
    for (size_t j = 0; j < static_cast<size_t>(markers.size()); ++j) {
      auto allele_1 = get_ancestry(Pop[i].chromosome1, markers[j]);
      auto allele_2 = get_ancestry(Pop[i].chromosome2, markers[j]);
      std::vector<std::string> entry = combine_alleles(allele_1, allele_2);

      int marker_start = j * 2;
      int second_marker = marker_start + 1;
      output(i, marker_start) = entry[0];
      output(i, second_marker) = entry[1];
    }
  }
  return(output);
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  Rcpp::stop("c++ exception (unknown reason)");
}
return NA_REAL;
}


float calc_het(const Fish& indiv, double marker) {
  int allele1 = get_ancestry(indiv.chromosome1, marker);
  int allele2 = get_ancestry(indiv.chromosome2, marker);

  if(allele1 != allele2) return 1.0f;

  return 0.f;
}


float calculate_heterozygosity_pop(const std::vector< Fish >& pop,
                                   double marker) {

  float freq_heterozygosity = 0.f;
  for(const auto& i : pop) {
    freq_heterozygosity += calc_het(i, marker);
  }
  freq_heterozygosity *= 1.0f / pop.size();
  return(freq_heterozygosity);
}

arma::mat update_heterozygosities_tibble(const std::vector< Fish >& pop,
                                         const NumericVector& markers,
                                         bool progress_bar) {

  arma::mat output(markers.size(), 2);

  int updateFreq = markers.size() / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(progress_bar) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }

  for(int i = 0; i < markers.size(); ++i) {
    output(i, 0) = markers[i];
    output(i, 1) = calculate_heterozygosity_pop(pop, markers[i]);
    if( i % updateFreq == 0 && progress_bar) {
      Rcout << "**";
    }
  }
  return(output);
}



// [[Rcpp::export]]
arma::mat calculate_heterozygosity_cpp(Rcpp::NumericVector input_population,
                                       Rcpp::NumericVector markers,
                                       bool progress_bar) {
try {
  std::vector< Fish > Pop;

  Pop = convert_NumericVector_to_fishVector(input_population);



  arma::mat heterozygosities = update_heterozygosities_tibble(Pop,
                                                              markers,
                                                              progress_bar);
  return heterozygosities;
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  Rcpp::stop("c++ exception (unknown reason)");
}
return arma::mat();
}

NumericVector scale_markers(const Rcpp::NumericVector& markers,
                            double morgan) {
  Rcpp::NumericVector outputmarkers(markers.size());

  for(int i = 0; i < markers.size(); ++i) {
    // negative values are 'no marker' placeholders and stay negative
    outputmarkers[i] = markers[i] < 0 ? markers[i] : markers[i] / morgan;
  }
  return outputmarkers;
}


//// EMP helper functions ///
//
//
//


std::vector< Fish_emp > convert_numeric_matrix_to_fish_vector(
  const Rcpp::NumericMatrix& input_population) {

  std::vector< Fish_emp > output;
  for(int i = 0; i < (input_population.nrow() - 1); i += 2) {
    Rcpp::NumericVector cc1 = input_population(i, _);
    Rcpp::NumericVector cc2 = input_population(i + 1 , _);
    std::vector<int> c1(cc1.begin(), cc1.end());
    std::vector<int> c2(cc2.begin(), cc2.end());

    output.push_back( Fish_emp(c1, c2) );
  }
  return(output);
}

double calculate_fitness(const Fish_emp& focal,
                         const std::vector<std::array<double, 5>>& select,
                         const std::vector<double>& locations,
                         bool multiplicative_selection) {

  int number_of_markers = select.size();
  std::vector<double> fitness_vec(number_of_markers,
                                  multiplicative_selection ? 1.0 : 0.0);

  for (int i = 0; i < number_of_markers; ++i) {
    auto focal_pos = select[i][0];
    auto focal_anc = select[i][4];
    if (focal_anc == -1) continue; // do not take into account
    int focal_index = find_location(locations, focal_pos);

    if (focal_index < 0 ||
        static_cast<size_t>(focal_index) >= focal.chromosome1.size() ||
        static_cast<size_t>(focal_index) >= focal.chromosome2.size()) {
      throw std::runtime_error("site under selection not found in data");
    }

    auto a1 = focal.chromosome1[focal_index];
    auto a2 = focal.chromosome2[focal_index];
    int fit_index = 1 + (a1 == focal_anc) + (a2 == focal_anc);
    fitness_vec[i] = select[i][fit_index];
 //   Rcout << a1 << " " << a2 << " " << select(i, fit_index) << "\n";
  }

  double output = multiplicative_selection ?
      std::accumulate(fitness_vec.begin(), fitness_vec.end(), 1.0,
                             std::multiplies<>()) :
      std::accumulate(fitness_vec.begin(), fitness_vec.end(), 0.0);

  return output;
}

List convert_to_list(const std::vector<Fish_emp>& v,
                     const std::vector<double>& locations) {
  int list_size = (int)v.size();
  List output(list_size);

  for(size_t i = 0; i < v.size(); ++i) {

    Fish_emp focal = v[i];

    NumericMatrix chrom1(focal.chromosome1.size(), 2); // nrow = number of junctions, ncol = 2
    for(size_t j = 0; j < focal.chromosome1.size(); ++j) {
      chrom1(j, 0) = locations[j];
      chrom1(j, 1) = focal.chromosome1[j];
    }

    NumericMatrix chrom2(focal.chromosome2.size(), 2); // nrow = number of junctions, ncol = 2
    for(size_t j = 0; j < focal.chromosome2.size(); ++j) {
      chrom2(j, 0) = locations[j];
      chrom2(j, 1) = focal.chromosome2[j];
    }

    List toAdd = List::create( Named("chromosome1") = chrom1,
                               Named("chromosome2") = chrom2
                              );
    output(i) = toAdd;
  }

  return output;
}

std::vector< std::vector<double > > update_frequency_tibble(const std::vector< Fish_emp >& pop,
                                                            size_t marker_index,
                                                            double pos,
                                                            int t) {

  int num_alleles = 5;
  std::vector< std::vector< double >> allele_matrix(num_alleles,
                                                    std::vector<double>(4, 0));

  // initialize results
  for (int i = 0; i < num_alleles; ++i) {
    allele_matrix[i][0] = t;
    allele_matrix[i][1] = pos;
    allele_matrix[i][2] = i;
    allele_matrix[i][3] = 0;
  }

  for (size_t i = 0; i < pop.size(); ++i) {
    if (marker_index >= pop[i].chromosome1.size() ||
        marker_index >= pop[i].chromosome2.size()) {
      throw std::runtime_error("marker index out of bounds");
    }
    for (int local_anc : {pop[i].chromosome1[marker_index],
                          pop[i].chromosome2[marker_index]}) {
      if (local_anc < 0 || local_anc >= num_alleles) {
        throw std::runtime_error("sequence data contains values other than "
                                 "0 (missing), 1 (a), 2 (c), 3 (t), 4 (g)");
      }
      allele_matrix[local_anc][3]++;
    }
  }

  for (size_t i = 0; i < allele_matrix.size(); ++i) {
    allele_matrix[i][3] *= 1.0 / (2 * pop.size());
  }
  return(allele_matrix);
}



int find_location(const std::vector<double>& markers,
                  double pos) {

  auto loc = std::lower_bound(markers.begin(), markers.end(), pos);

  if (loc != markers.end()) {
    if (*loc == pos) {
      return std::distance(markers.begin(), loc);
    }
  }
  return -1;
}

arma::mat update_all_frequencies_tibble(const std::vector< Fish_emp >& pop,
                                        const std::vector<double>& markers,
                                        const std::vector<double>& locations,
                                        int t,
                                        double morgan) {

  int number_of_alleles = 5; // always the same: [0, 1, 2, 3, 4]
  arma::mat output(markers.size() * number_of_alleles, 4);

  for (size_t i = 0; i < markers.size(); ++i) {
    if (markers[i] < 0) continue;  // not tracked: rows stay empty

    auto index = find_location(locations, markers[i]);

    if (index < 0) {
      throw std::runtime_error("tracked marker not found in data");
    }

    std::vector<std::vector<double >> local_mat = update_frequency_tibble(pop,
                                                  index,
                                                  markers[i],
                                                  t);
    int start = i * number_of_alleles;
    int end = start + number_of_alleles;
    for(int j = start; j < end; ++j) {
      for(int k = 0; k < 4; ++k) {
        output(j, k) = local_mat[j - start][k];
      }
    }
  }
  return(output);
}

int draw_mutated_base(int source_base,
                      const std::vector< std::vector< double >>& sub_matrix,
                      rnd_t& rndgen) {
  if (source_base == 0) // no data
    return 0;

  static int alleles[4] = {1, 2, 3, 4};

  double r = rndgen.uniform();
  for (int i = 0; i < 4; ++i) {
    r -= sub_matrix[source_base - 1][i];
    if (r <= 0.0) {
      return( alleles[i] );
    }
  }
  return alleles[3];
}

int num_mutations(int num_markers,
                  double mu,
                  rnd_t& rndgen)  {
  std::binomial_distribution<> mutate_num(num_markers, mu);
  int result = mutate_num(rndgen.rndgen_);
  return result;
}

void mutate_chrom(std::vector<int>& chrom,
                  const std::vector<std::vector<double>>& sub_matrix,
                  const double& mutation_rate,
                  rnd_t& rndgen) {

  int num_mutated_bases = num_mutations(chrom.size(),
                                        mutation_rate,
                                        rndgen);

  for (int i = 0; i < num_mutated_bases; ++i) {
    int index = rndgen.random_number(chrom.size());
    int mutated_base = draw_mutated_base(chrom[index],
                                         sub_matrix,
                                         rndgen);
    chrom[index] = mutated_base;
  }
  return;
}

void mutate(Fish_emp& indiv,
            const std::vector<std::vector<double>>& sub_matrix,
            const double& mutation_rate,
            rnd_t& rndgen) {

  mutate_chrom(indiv.chromosome1, sub_matrix, mutation_rate, rndgen);
  mutate_chrom(indiv.chromosome2, sub_matrix, mutation_rate, rndgen);

  return;
}

std::vector<int> get_alleles(int focal_genotype,
                             int allele_1,
                             int allele_2) {

  if (focal_genotype == 1) {
    return {allele_1, allele_1};
  }
  if (focal_genotype == 2) {
      return {allele_1, allele_2};
  }
  if (focal_genotype == 3) {
    return {allele_2, allele_2};
  }

  return {0, 0};
}

// [[Rcpp::export]]
NumericMatrix vcf_to_matrix_cpp(const Rcpp::NumericMatrix input_mat,
                            const NumericVector& allele_1,
                            const NumericVector& allele_2) {
  int num_indiv = input_mat.nrow();
  int num_markers = allele_1.size();

  NumericMatrix output(2 * num_indiv, num_markers);

  for (int i = 0; i < num_indiv; ++i) {
    int index_c1 = i * 2;
    int index_c2 = index_c1 + 1;

    for (int j = 0; j < num_markers; ++j) {
      std::vector< int > alleles = get_alleles(input_mat(i, j),
                                               allele_1[j],
                                               allele_2[j]);
      output(index_c1, j) = alleles[0];
      output(index_c2, j) = alleles[1];
    }
  }

  return output;
}

// Converts the R selection matrix. Returns an empty vector (no selection)
// when the matrix is the "no selection" placeholder. Locations are divided
// by 'scale' (morgan, for ancestry simulations; 1 for sequence data).
std::vector< std::array<double, 5> > convert_select_from_r(const Rcpp::NumericMatrix& select,
                                                           double scale) {
  std::vector< std::array<double, 5> > selection_sites;
  const bool is_placeholder = select.ncol() != 5;
  if (is_placeholder) {
    return selection_sites;
  }
  for (int row = 0; row < select.nrow(); ++row) {
    std::array<double, 5> site;
    for (int column = 0; column < 5; ++column) {
      site[column] = select(row, column);
    }
    site[0] /= scale;  // location
    selection_sites.push_back(site);
  }
  return selection_sites;
}
