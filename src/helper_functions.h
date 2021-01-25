//
//  helper_functions.hpp
//
//
//  Created by Thijs Janzen on 20/11/2018.
//

#ifndef helper_functions_hpp
#define helper_functions_hpp

#include <vector>
#include <RcppArmadillo.h>

#include "random_functions.h"

// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

template <typename FISH>
double calc_fitness(const FISH & focal,
                    const Rcpp::NumericMatrix& select,
                    bool multiplicative_selection) {

  int number_of_markers = select.nrow();
  std::vector<double> fitness_vec(number_of_markers);

  for (int i = 0; i < number_of_markers; ++i) {
    double focal_pos = select(i, 0);
    double focal_anc = select(i, 4);
    if (focal_anc < 0) continue; // do not take into account

    int a1 = focal.chromosome1.get_ancestry(focal_pos);
    int a2 = focal.chromosome2.get_ancestry(focal_pos);
    int fit_index = 1 + (a1 == focal_anc) + (a2 == focal_anc);
    fitness_vec[i] = select(i, fit_index);
  }

  if (!multiplicative_selection) {
    return std::accumulate(fitness_vec.begin(), fitness_vec.end(), 0.0);
  }

  return std::accumulate(fitness_vec.begin(), fitness_vec.end(), 0.0,
                         std::multiplies<>());
}

NumericVector scale_markers(const Rcpp::NumericVector& markers,
                            double morgan) {
  if (markers.size() == 1) {
    return markers;
  }

  Rcpp::NumericVector outputmarkers(markers.size());

  for(int i = 0; i < markers.size(); ++i) {
    outputmarkers[i] = markers[i] * 1.0 / morgan;
  }
  return outputmarkers;
}

template <typename CHROMO>
bool matching_chromosomes(const CHROMO& v1,
                          const CHROMO& v2)
{
  if(v1.genome.size() != v2.genome.size()) {
    return false;
  }
  for(size_t i = 0; i < v1.genome.size(); ++i) {
    if(v1.genome[i] != v2.genome[i]) {
      return false;
    }
  }
  return true;
}

template<typename FISH>
bool is_fixed(const std::vector< FISH >& v) {

  if(!matching_chromosomes(v[0].chromosome1, v[0].chromosome2)) {
    return false;
  }

  for(auto it = v.begin(); it != v.end(); ++it) {
    if(!matching_chromosomes((*it).chromosome1, v[0].chromosome1)) {
      return false;
    }
    if(!matching_chromosomes((*it).chromosome1, (*it).chromosome2)) {
      return false;
    }
  }
  return true;
}


int find_index(const std::vector<int>& v, int value) {
  for(size_t i = 0; i < v.size(); ++i) {
    if(v[i] == value) return i;
  }
  return -1;
}

template <typename CHROMO>
void update_founder_labels(const CHROMO& chrom,
                           std::vector<int>& founder_labels) {
  for(auto i = chrom.genome.begin(); i != chrom.genome.end(); ++i) {
    if(founder_labels.empty()) {
      if((*i).get_ancestry() != -1) founder_labels.push_back((*i).get_ancestry());
    } else {
      if(find_index(founder_labels, (*i).right) == -1) {
        if((*i).get_ancestry() != -1) founder_labels.push_back((*i).get_ancestry());
      }
    }
  }
  return;
}

template <typename FISH>
arma::mat update_frequency_tibble(const std::vector< FISH >& v,
                                  double m,
                                  const std::vector<int>& founder_labels,
                                  int t,
                                  double morgan) {

  int num_alleles = founder_labels.size();
  arma::mat allele_matrix(num_alleles, 4);
  // initialize results
  for(int i = 0; i < num_alleles; ++i) {
    allele_matrix(i, 0) = t;
    allele_matrix(i, 1) = m * morgan;
    allele_matrix(i, 2) = founder_labels[i];
    allele_matrix(i, 3) = 0;
  }

  for(auto it = v.begin(); it != v.end(); ++it) {

    for(auto i = ((*it).chromosome1.genome.begin()+1);
        i != (*it).chromosome1.genome.end(); ++i) {
      if((*i).pos > m) {
        auto local_anc = (*(i-1)).right;
        int index = find_index(founder_labels, local_anc);
        allele_matrix(index, 3)++;
        break;
      }
    }

    for(auto i = ((*it).chromosome2.genome.begin()+1);
        i != (*it).chromosome2.genome.end(); ++i) {
      if((*i).pos > m) {
        auto local_anc = (*(i-1)).right;
        int index = find_index(founder_labels, local_anc);
        allele_matrix(index, 3)++;
        break;
      }
    }
  }

  for(int i = 0; i < num_alleles; ++i) {
    allele_matrix(i, 3) *= 1.0 / (2 * v.size());
  }

  return(allele_matrix);
}

template <typename FISH>
arma::mat update_all_frequencies_tibble(const std::vector< FISH >& pop,
                                        const NumericVector& markers,
                                        const std::vector<int>& founder_labels,
                                        int t,
                                        double morgan) {

  int number_of_alleles = founder_labels.size();
  arma::mat output(markers.size() * number_of_alleles, 4);

  for(int i = 0; i < markers.size(); ++i) {
    arma::mat local_mat = update_frequency_tibble<FISH>(pop,
                                                        markers[i],
                                                               founder_labels,
                                                               t,
                                                               morgan);
    // now we have a (markers x alleles) x 3 tibble, e.g. [loc, anc, freq]
    // and we have to put that in the right place in the output matrix
    int start = i * number_of_alleles;
    int end = start + number_of_alleles;
    for(int j = start; j < end; ++j) {
      for(int k = 0; k < 4; ++k) {
        output(j, k) = local_mat(j - start, k);
      }
    }
  }
  return(output);
}


template <typename FISH>
arma::mat record_frequencies_pop(const std::vector< FISH >& pop,
                                 const NumericVector& markers,
                                 const std::vector<int>& founder_labels,
                                 int t,
                                 int pop_indicator,
                                 double morgan) {
  int number_of_alleles = founder_labels.size();
  arma::mat output(markers.size() * number_of_alleles, 5);

  for(int i = 0; i < markers.size(); ++i) {
    arma::mat local_mat = update_frequency_tibble<FISH>(pop,
                                                        markers[i],
                                                               founder_labels,
                                                               t,
                                                               morgan);
    // now we have a (markers x alleles) x 5 tibble, e.g. [loc, anc, freq, pop]
    // and we have to put that in the right place in the output matrix
    int start = i * number_of_alleles;
    int end = start + number_of_alleles;
    for(int j = start; j < end; ++j) {
      for(int k = 0; k < 4; ++k) {
        output(j, k) = local_mat(j - start, k);
      }
      output(j, 4) = pop_indicator;
    }
  }
  return(output);
}

template <typename FISH>
arma::mat update_all_frequencies_tibble_dual_pop(const std::vector< FISH >& pop_1,
                                                 const std::vector< FISH >& pop_2,
                                                 const NumericVector& markers,
                                                 const std::vector<int>& founder_labels,
                                                 int t,
                                                 double morgan) {
  arma::mat output_1 = record_frequencies_pop<FISH>(pop_1, markers, founder_labels, t, 1, morgan);
  arma::mat output_2 = record_frequencies_pop<FISH>(pop_2, markers, founder_labels, t, 2, morgan);

  arma::mat output = arma::join_cols(output_1, output_2);
  return(output);
}

template <typename FISH>
double calc_mean_junctions(const std::vector< FISH > & pop) {

  double mean_junctions = 0.0;
  for(auto it = pop.begin(); it != pop.end(); ++it) {
    mean_junctions += (*it).chromosome1.get_num_junctions(); // start and end don't count
    mean_junctions += (*it).chromosome2.get_num_junctions();
  }
  mean_junctions *= 1.0 / (pop.size() * 2); // diploid

  return(mean_junctions);
}

int draw_prop_fitness(const std::vector<double>& fitness,
                      double maxFitness) {

  if (maxFitness <= 0.0) {
    return random_number(fitness.size());
  }

  while(true) {
    int index = random_number(fitness.size());
    double prob = 1.0 * fitness[index] / maxFitness;
    if(uniform() < prob) {
      return index;
    }
  }
}


template <typename FISH, typename MARKER>
std::vector< FISH >
convert_NumericVector_to_fishVector(const Rcpp::NumericVector& v) {

  std::vector< FISH > output;

  FISH temp;
  int indic_chrom = 1;
  bool add_indiv = false;

  for(size_t i = 0; i < v.size(); i += 2) {
    MARKER temp_j;
    temp_j.pos = v[i];
    temp_j.right = v[i+1];

    if(indic_chrom == 1) {
      temp.chromosome1.genome.push_back(temp_j);
    } else {
      temp.chromosome2.genome.push_back(temp_j);
    }

    if(temp_j.right == -1) {
      if(indic_chrom == 1) {
        indic_chrom = 2;
      } else {
        add_indiv = true;
      }
    }

    if(add_indiv) {
      output.push_back(temp);
      add_indiv = false;
      indic_chrom = 1;
      temp.chromosome1.genome.clear();
      temp.chromosome2.genome.clear();
    }
  }

  return(output);
}


template <typename FISH>
List convert_to_list(const std::vector<FISH>& v) {
  int list_size = (int)v.size();
  List output(list_size);

  for(size_t i = 0; i < v.size(); ++i) {

    FISH focal = v[i];

    NumericMatrix chrom1(focal.chromosome1.genome.size(), 2); // nrow = number of junctions, ncol = 2
    for(size_t j = 0; j < focal.chromosome1.genome.size(); ++j) {
      chrom1(j, 0) = focal.chromosome1.genome[j].pos;
      chrom1(j, 1) = focal.chromosome1.genome[j].right;
    }

    NumericMatrix chrom2(focal.chromosome2.genome.size(), 2); // nrow = number of junctions, ncol = 2
    for(size_t j = 0; j < focal.chromosome2.genome.size(); ++j) {
      chrom2(j, 0) = focal.chromosome2.genome[j].pos;
      chrom2(j, 1) = focal.chromosome2.genome[j].right;
    }

    List toAdd = List::create( Named("chromosome1") = chrom1,
                               Named("chromosome2") = chrom2
    );

    output(i) = toAdd;
  }

  return output;
}

int draw_random_founder(const NumericVector& v) {
  double r = uniform();
  for(int i = 0; i < v.size(); ++i) {
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
                                        bool progress_bar)
{
  std::vector< Fish<chromosome_junctions> > Pop =
    convert_NumericVector_to_fishVector< Fish<chromosome_junctions>,
                                         junction> (input_population);

  std::vector<int> founder_labels;
  for(auto it = Pop.begin(); it != Pop.end(); ++it) {
    update_founder_labels<chromosome_junctions>((*it).chromosome1, founder_labels);
    update_founder_labels<chromosome_junctions>((*it).chromosome2, founder_labels);
  }
  //Rcout << "number of alleles: " << founder_labels.size() << "\n";

  double morgan = markers[markers.size() - 1];
  markers = scale_markers(markers, morgan); // make sure they are in [0, 1];

  arma::mat frequencies = update_all_frequencies_tibble<Fish<chromosome_junctions>>(Pop,
                                                                                    markers,
                                                                                    founder_labels,
                                                                                    0,
                                                                                    morgan);

  return frequencies;
}

template <typename FISH>
float calc_het(const FISH& indiv, float marker) {
  int allele1 = indiv.chromosome1.get_ancestry(marker);
  int allele2 = indiv.chromosome2.get_ancestry(marker);
  if(allele1 != allele2) return 1.0f;

  return 0.f;
}


template <typename FISH>
float calculate_heterozygosity_pop(const std::vector< FISH >& pop,
                                   float marker) {

  float freq_heterozygosity = 0.f;
  for(const auto& i : pop) {
    freq_heterozygosity += calc_het<FISH>(i, marker);
  }
  freq_heterozygosity *= 1.0f / pop.size();
  return(freq_heterozygosity);
}

template <typename FISH>
arma::mat update_heterozygosities_tibble(const std::vector< FISH >& pop,
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
    output(i, 1) = calculate_heterozygosity_pop<FISH>(pop, markers[i]);
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
  std::vector< Fish<chromosome_junctions> > Pop;

  Pop = convert_NumericVector_to_fishVector<Fish<chromosome_junctions>,
                                            junction>(input_population);

  arma::mat heterozygosities = update_heterozygosities_tibble<Fish<chromosome_junctions>>(Pop,
                                                                                          markers,
                                                                                          progress_bar);
  return heterozygosities;
}

#endif /* helper_functions_hpp */
