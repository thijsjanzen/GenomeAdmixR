//
//  helper_functions.cpp
//
//
//  Created by Thijs Janzen on 20/11/2018.
//

#include "helper_functions.h"
#include <vector>
#include <string>
#include <thread>
#include <chrono>

#include <chrono>
#include <thread>

void force_output() {
//  std::this_thread::sleep_for(std::chrono::nanoseconds(100));
  std::this_thread::sleep_for(std::chrono::milliseconds(300));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}


bool is_fixed(const std::vector< Fish >& v) {

  if (!matching_chromosomes(v[0].chromosome1, v[0].chromosome2)) {
    return false;
  }

  for(auto it = v.begin(); it != v.end(); ++it) {
    if (!matching_chromosomes((*it).chromosome1, v[0].chromosome1)) {
      return false;
    }
    if (!matching_chromosomes((*it).chromosome1, (*it).chromosome2)) {
      return false;
    }
  }
  return true;
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

void update_anc_chrom(const std::vector<junction>& chrom,
                      const std::vector<int>& founder_labels,
                      double m,
                      arma::mat& allele_matrix) {
  if (chrom.size() == 1) {
    if (m >= chrom[0].pos) {
      int local_anc = chrom[0].right;
      int index = find_index(founder_labels, local_anc);
      allele_matrix(index, 3)++;
      return;
    }
  }

  if (m > chrom.back().pos) {
    int local_anc = chrom.back().right;
    int index = find_index(founder_labels, local_anc);
    allele_matrix(index, 3)++;
    return;
  }


  for(auto i = chrom.begin(); i != chrom.end(); ++i) {
    if (i->pos == m) {
      int local_anc = i->right;
      int index = find_index(founder_labels, local_anc);
      allele_matrix(index, 3)++;
      return;
    }

    if (i->pos > m) {
      if (i != chrom.begin()) {

        int local_anc = (*(i-1)).right;
        int index = find_index(founder_labels, local_anc);
        allele_matrix(index, 3)++;
        return;
      }
    }
  }

  Rcpp::Rcout << m << "\n"; force_output();
  Rcpp::stop("ERROR DID NOT FIND MARKER");
  return;
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
                                        const NumericVector& markers,
                                        const std::vector<int>& founder_labels,
                                        int t,
                                        double morgan) {

  int number_of_alleles = founder_labels.size();
  arma::mat output(markers.size() * number_of_alleles, 4);

  for (int i = 0; i < markers.size(); ++i) {
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



arma::mat record_frequencies_pop(const std::vector< Fish >& pop,
                                 const NumericVector& markers,
                                 const std::vector<int>& founder_labels,
                                 int t,
                                 int pop_indicator,
                                 double morgan) {
  int number_of_alleles = founder_labels.size();
  arma::mat output(markers.size() * number_of_alleles, 5);
  if (markers.size() < 1) {
    Rcout << "markers empty\n"; force_output();
    return(output);
  }

  for(int i = 0; i < markers.size(); ++i) {
    if (markers[i] < 0) continue;
    arma::mat local_mat = update_frequency_tibble(pop,
                                                  markers[i],
                                                         founder_labels,
                                                         t,
                                                         morgan);
    // now we have a (markers x alleles) x 5 tibble, e.g. [loc, anc, freq, pop]
    // and we have to put that in the right place in the output matrix
    int start = i * number_of_alleles;
    int end = start + number_of_alleles;
    for (int j = start; j < end; ++j) {
      for (int k = 0; k < 4; ++k) {
        output(j, k) = local_mat(j - start, k);
      }
      output(j, 4) = pop_indicator;
    }
  }
  return(output);
}

arma::mat update_all_frequencies_tibble_dual_pop(const std::vector< Fish >& pop_1,
                                                 const std::vector< Fish >& pop_2,
                                                 const NumericVector& markers,
                                                 const std::vector<int>& founder_labels,
                                                 int t,
                                                 double morgan) {
  arma::mat output_1 = record_frequencies_pop(pop_1, markers, founder_labels, t, 1, morgan);
  arma::mat output_2 = record_frequencies_pop(pop_2, markers, founder_labels, t, 2, morgan);

  arma::mat output = arma::join_cols(output_1, output_2);
  return(output);
}

double calc_mean_junctions(const std::vector< Fish> & pop) {

  double mean_junctions = 0.0;
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    mean_junctions += (*it).chromosome1.size() - 2; // start and end don't count
    mean_junctions += (*it).chromosome2.size() - 2;
  }
  mean_junctions *= 1.0 / (pop.size() * 2); // diploid

  return(mean_junctions);
}

int draw_prop_fitness(const std::vector<double>& fitness,
                      const double& maxFitness,
                      rnd_t& rndgen) {

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

  int num_indiv = 0;

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

      num_indiv++;
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
//  Rcout << "converting to list\n"; force_output();
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

double calculate_fitness(const Fish& focal,
                         const NumericMatrix& select,
                         bool multiplicative_selection) {

  int number_of_markers = select.nrow();
  std::vector< int > num_alleles(number_of_markers, 0);

  int focal_marker = 0;
  double pos = select(focal_marker, 0);
  double anc = select(focal_marker, 4);
  // loc aa  Aa  AA ancestor
  //  0  1   2  3  4

  for(auto it = (focal.chromosome1.begin()+1); it != focal.chromosome1.end(); ++it) {
    if ((*it).pos > pos) {
      if ((*(it-1)).right == anc) num_alleles[focal_marker]++;
      focal_marker++;
      if(focal_marker >= number_of_markers) {
        break;
      }
      pos = select(focal_marker, 0);
      anc = select(focal_marker, 4);
    }
    if (anc < 0) break; // these entries are only for tracking alleles over time, not for selection calculation
  }

  focal_marker = 0;
  pos = select(focal_marker, 0);
  anc = select(focal_marker, 4);

  for(auto it = (focal.chromosome2.begin()+1); it != focal.chromosome2.end(); ++it) {
    if ((*it).pos > pos) {
      if ((*(it-1)).right == anc) num_alleles[focal_marker]++;
      focal_marker++;
      if (focal_marker >= number_of_markers) {
        break;
      }
      pos = select(focal_marker, 0);
      anc = select(focal_marker, 4);
    }
    if(anc < 0) break; // these entries are only for tracking alleles over time, not for selection calculation
  }

  double fitness = 0.0;
  if (multiplicative_selection) fitness = 1.0;
  for (size_t i = 0; i < num_alleles.size(); ++i) {
    if(select(i, 4) < 0) break; // these entries are only for tracking alleles over time, not for selection calculation

    int fitness_index = 1 + num_alleles[i];
    if (multiplicative_selection) {
      fitness *= select(i, fitness_index);
    } else {
      fitness += select(i, fitness_index);
    }
  }

  return(fitness);
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

  double morgan = markers[markers.size() - 1];
  markers = scale_markers(markers, morgan); // make sure they are in [0, 1];

  arma::mat frequencies = update_all_frequencies_tibble(Pop,
                                                        markers,
                                                        founder_labels,
                                                        0,
                                                        morgan);

  return frequencies;
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return arma::mat();
}

int get_ancestry(const std::vector< junction >& chrom,
                 float marker) {

  if (marker < chrom.front().pos) {
    return(-1);
  }

  if (marker > chrom.back().pos) {
    return(chrom.back().right);
  }

  for(auto i = chrom.begin(); i != chrom.end(); ++i) {
    if ((*i).pos == marker) {
      return (*i).right;
    }

    if ((*i).pos > marker) {
      return((i-1)->right);
    }
  }
  return chrom.back().right;
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
  ::Rf_error("c++ exception (unknown reason)");
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
      //   Rcout << i << " " << j << "\n"; force_output();
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
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}


float calc_het(const Fish& indiv, float marker) {
  int allele1 = get_ancestry(indiv.chromosome1, marker);
  int allele2 = get_ancestry(indiv.chromosome2, marker);

  if(allele1 != allele2) return 1.0f;

  return 0.f;
}


float calculate_heterozygosity_pop(const std::vector< Fish >& pop,
                                   float marker) {

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
  ::Rf_error("c++ exception (unknown reason)");
}
return arma::mat();
}

NumericVector scale_markers(const Rcpp::NumericVector& markers,
                            double morgan) {
  if (markers.size() == 1) {
    return markers;
  }

  Rcpp::NumericVector outputmarkers(markers.size());

 // double max = *std::max_element(markers.begin(), markers.end());

  for(int i = 0; i < markers.size(); ++i) {
    outputmarkers[i] = markers[i] * 1.0 / morgan;
  }
  return outputmarkers;
}

std::vector<double> scale_markers(const std::vector<double>& markers,
                                  double morgan) {
  if (markers.size() == 1) {
    return markers;
  }

  std::vector<double> outputmarkers(markers.size());

 // double max = *std::max_element(markers.begin(), markers.end());

  for (size_t i = 0; i < markers.size(); ++i) {
    outputmarkers[i] = markers[i] * 1.0 / morgan;
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
  for(int i = 0; i < (input_population.nrow() - 1); i+= 2) {

    Rcpp::NumericVector cc1 = input_population(i, _);
    Rcpp::NumericVector cc2 = input_population(i + 1 , _);
    std::vector<int> c1(cc1.begin(), cc1.end());
    std::vector<int> c2(cc2.begin(), cc2.end());
    output.emplace_back( Fish_emp(c1, c2));
  }
  return(output);
}

void update_founder_labels(const std::vector<int>& chrom,
                           std::vector<int>& founder_labels) {
  for(auto i = chrom.begin(); i != chrom.end(); ++i) {
    if(founder_labels.empty()) {
       founder_labels.push_back((*i));
    } else {
      if(find_index(founder_labels, (*i)) == -1) {
        founder_labels.push_back((*i));
      }
    }
  }
  return;
}



double calculate_fitness(const Fish_emp& focal,
                         const NumericMatrix& select,
                         const std::vector<double>& locations,
                         bool multiplicative_selection) {

  int number_of_markers = select.nrow();
  std::vector<double> fitness_vec(number_of_markers);

  for (int i = 0; i < number_of_markers; ++i) {
    auto focal_pos = select(i, 0);
    auto focal_anc = select(i, 4);
    if (focal_anc == -1) continue; // do not take into account
    int focal_index = find_location(locations, focal_pos);

    auto a1 = focal.chromosome1[focal_index];
    auto a2 = focal.chromosome2[focal_index];
    int fit_index = 1 + (a1 == focal_anc) + (a2 == focal_anc);
    fitness_vec[i] = select(i, fit_index);
 //   Rcout << a1 << " " << a2 << " " << select(i, fit_index) << "\n";
  }

  if (!multiplicative_selection) {
    return std::accumulate(fitness_vec.begin(), fitness_vec.end(), 0.0);
  }

  return std::accumulate(fitness_vec.begin(), fitness_vec.end(), 1.0,
                         std::multiplies<>());



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

//  Rcout << "start update_frequency_tibble\n"; force_output();
  // initialize results
  for (int i = 0; i < num_alleles; ++i) {
    allele_matrix[i][0] = t;
    allele_matrix[i][1] = pos;
    allele_matrix[i][2] = i;
    allele_matrix[i][3] = 0;
  }

  for (size_t i = 0; i < pop.size(); ++i) {
    if (marker_index >= pop[i].chromosome1.size()) {
      Rcpp::stop("marker index out of bounds");
    }
    if (marker_index >= pop[i].chromosome2.size()) {
      Rcpp::stop("marker index out of bounds");
    }
    size_t local_anc1 = pop[i].chromosome1[marker_index];
    allele_matrix[local_anc1][3]++;
    size_t local_anc2 = pop[i].chromosome2[marker_index];
    allele_matrix[local_anc2][3]++;
  }

 // Rcout << "starting normalization\n"; force_output();
  for (size_t i = 0; i < allele_matrix.size(); ++i) {
    allele_matrix[i][3] *= 1.0 / (2 * pop.size());
  }
  return(allele_matrix);
}



int find_location(const std::vector<double>& markers,
                  double pos) {

 // Rcout << "this is find location with: " << pos << "\n"; force_output();

  //auto loc = std::find(markers.begin(), markers.end(), pos);
  auto loc = std::lower_bound(markers.begin(), markers.end(), pos);

  if (loc != markers.end()) {
    if (*loc == pos) {
      return std::distance(markers.begin(), loc);
    }
  }

  /*for (int i = 0; i < markers.size(); ++i) {
    if (markers[i] == pos) {
      return i;
    }
  }*/
  return - 1;

//  return 0;
}

arma::mat update_all_frequencies_tibble(const std::vector< Fish_emp >& pop,
                                        const std::vector<double>& markers,
                                        const std::vector<double>& locations,
                                        int t,
                                        double morgan) {

  int number_of_alleles = 5; // always the same: [0, 1, 2, 3, 4]
  arma::mat output(markers.size() * number_of_alleles, 4);

  for (size_t i = 0; i < markers.size(); ++i) {

    auto index = find_location(locations, markers[i]);

    if (index < 0) {
      Rcout << i << " " << index << "\n"; force_output();
      Rcpp::stop("out of bounds in update_all_frequencies_tibble\n");
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

bool is_fixed(const std::vector< Fish_emp >& v) {

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

bool matching_chromosomes(const std::vector< junction >& v1,
                          const std::vector< junction >& v2) {
  if(v1.size() != v2.size()) {
    return false;
  }
  for(size_t i = 0; i < v1.size(); ++i) {
    if(v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

bool matching_chromosomes(const std::vector< int >& v1,
                          const std::vector< int >& v2) {
  if(v1.size() != v2.size()) {
    return false;
  }
  for(size_t i = 0; i < v1.size(); ++i) {
    if(v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

int count_num_j(const std::vector< int >& chrom) {
  int num_j = 0;
  for (size_t i = 1; i < chrom.size(); ++i) {
    if (chrom[i] != chrom[i - 1]) num_j++;
  }
  return(num_j);
}


double number_of_junctions(const std::vector< Fish_emp>& pop) {
  double num_j = 0.0;
  for(const auto& i : pop) {
    num_j += count_num_j(i.chromosome1);
    num_j += count_num_j(i.chromosome2);
  }
  num_j *= 1.0 / (2 * pop.size());
  return(num_j);
}

arma::mat record_frequencies_pop(const std::vector< Fish_emp >& pop,
                                 const std::vector<double>& markers,
                                 const std::vector<double>& locations,
                                 int t,
                                 int pop_indicator,
                                 double morgan) {
  int number_of_alleles = 5;
  arma::mat output(markers.size() * number_of_alleles, 5);
  if (markers.size() < 1) {
    Rcout << "markers empty\n"; force_output();
    return(output);
  }

  for(size_t i = 0; i < markers.size(); ++i) {
   if (markers[i] < 0) continue;

   int index = find_location(locations, markers[i]);
   std::vector< std::vector< double> > local_mat = update_frequency_tibble(pop,
                                                                           index,
                                                            markers[i],
                                                            t);
    // now we have a (markers x alleles) x 5 tibble, e.g. [loc, anc, freq, pop]
    // and we have to put that in the right place in the output matrix
    int start = i * number_of_alleles;
    int end = start + number_of_alleles;
    for(int j = start; j < end; ++j) {
      for(int k = 0; k < 4; ++k) {
        output(j, k) = local_mat[j - start][k];
      }
      output(j, 4) = pop_indicator;
    }
  }
  return(output);
}

arma::mat update_all_frequencies_tibble_dual_pop(const std::vector< Fish_emp >& pop_1,
                                                 const std::vector< Fish_emp >& pop_2,
                                                 const std::vector<double>& markers,
                                                 const std::vector<double>& locations,
                                                 int t,
                                                 double morgan) {
  arma::mat output_1 = record_frequencies_pop(pop_1, markers, locations, t, 1, morgan);
  arma::mat output_2 = record_frequencies_pop(pop_2, markers, locations, t, 2, morgan);

  arma::mat output = arma::join_cols(output_1, output_2);
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
                             int allele_2,
                             rnd_t& rndgen) {

  if (focal_genotype == 1) {
    return {allele_1, allele_1};
  }
  if (focal_genotype == 2) {
   // if (rndgen.random_number(2) == 0) {
      return {allele_1, allele_2};
  //  } else {
    //  return {allele_2, allele_1};
  //  }
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
try {
  int num_indiv = input_mat.nrow();
  int num_markers = allele_1.size();

  NumericMatrix output(2 * num_indiv, num_markers);

  rnd_t rndgen;

  for (int i = 0; i < num_indiv; ++i) {
    int index_c1 = i * 2;
    int index_c2 = i + 1;

    for (int j = 0; j < num_markers; ++j) {
      std::vector< int > alleles = get_alleles(input_mat(i, j),
                                               allele_1[j],
                                               allele_2[j],
                                               rndgen);
      output(index_c1, j) = alleles[0];
      output(index_c2, j) = alleles[1];
    }
  }

  return output;
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}
