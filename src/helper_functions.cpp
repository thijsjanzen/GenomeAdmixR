//
//  helper_functions.cpp
//  
//
//  Created by Thijs Janzen on 20/11/2018.
//

#include "helper_functions.h"
#include <vector>

bool verify_individual_cpp(const Fish& Nemo) {
    for(int i = 0; i < Nemo.chromosome1.size(); ++i) {
        if(Nemo.chromosome1[i].right >  1000 |
           Nemo.chromosome1[i].right < -1000) {
            return false;
        }
    }

    for(int i = 0; i < Nemo.chromosome2.size(); ++i) {
        if(Nemo.chromosome2[i].right >  1000 |
           Nemo.chromosome2[i].right < -1000) {
            return false;
        }
    }

    return true;
}


bool verify_pop_cpp(const std::vector< Fish >& pop) {
    for(auto it = pop.begin(); it != pop.end(); ++it) {
        if(!verify_individual_cpp((*it))) {
            return false;
        }
    }

    return true;
}

bool matching_chromosomes(const std::vector< junction >& v1,
                          const std::vector< junction >& v2)
{
    if(v1.size() != v2.size()) {
        return false;
    }
    for(int i = 0; i < v1.size(); ++i) {
        if(v1[i] != v2[i]) {
            return false;
        }
    }
    return true;
}

bool is_fixed(const std::vector< Fish >& v) {

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

void update_founder_labels(const std::vector<junction> chrom,
                           std::vector<int>& founder_labels) {
    for(auto i = chrom.begin(); i != chrom.end(); ++i) {
        if(founder_labels.empty()) {
            founder_labels.push_back((*i).right);
        } else {
            if(std::find(founder_labels.begin(), founder_labels.end(), (*i).right) == founder_labels.end()) {
                founder_labels.push_back((*i).right);
            }
        }
    }
    return;
}

NumericVector update_frequency(const std::vector< Fish >& v,
                               double m,
                               int num_alleles) {

    NumericVector freq(num_alleles, 0.0);

    for(auto it = v.begin(); it != v.end(); ++it) {
        for(auto i = ((*it).chromosome1.begin()+1); i != (*it).chromosome1.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                if(index >= num_alleles || index < 0) {
                    Rcout << "ERROR!!\n";
                    Rcout << "trying to access NumericVector freq outside bounds\n";
                }
                freq(index)++;
                break;
            }
        }

        for(auto i = ((*it).chromosome2.begin()+1); i != (*it).chromosome2.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                if(index >= num_alleles || index < 0) {
                    Rcout << "ERROR!!\n";
                    Rcout << "trying to access NumericVector freq outside bounds\n";
                }
                freq(index)++;
                break;
            }
        }
    }

    for(int i = 0; i < freq.size(); ++i) {
        freq(i) = freq(i) * 1.0 / (2*v.size());
    }

    return(freq);
}

int find_index(const std::vector<int>& v, int value) {
    for(int i = 0; i < v.size(); ++i) {
        if(v[i] == value) return i;
    }
    Rcout << "ERROR! Could not find ancestry label, returning -1, expect out of range error soon\n";
    return -1;
}


arma::mat update_frequency_tibble(const std::vector< Fish >& v,
                                  double m,
                                  const std::vector<int>& founder_labels) {

    int num_alleles = founder_labels.size();
    arma::mat allele_matrix(num_alleles,3);
    for(int i = 0; i < num_alleles; ++i) {
        allele_matrix(i, 0) = m;
        allele_matrix(i, 1) = founder_labels[i];
        allele_matrix(i, 2) = 0;
    }

   for(auto it = v.begin(); it != v.end(); ++it) {
       for(auto i = ((*it).chromosome1.begin()+1); i != (*it).chromosome1.end(); ++i) {
           if((*i).pos > m) {
               int local_anc = (*(i-1)).right;
               int index = find_index(founder_labels, local_anc);
               allele_matrix(index, 2)++;
               break;
           }
       }

       for(auto i = ((*it).chromosome2.begin()+1); i != (*it).chromosome2.end(); ++i) {
           if((*i).pos > m) {
               int local_anc = (*(i-1)).right;
               int index = find_index(founder_labels, local_anc);
               allele_matrix(index, 2)++;
               break;
           }
       }
   }

   for(int i = 0; i < num_alleles; ++i) {
       allele_matrix(i, 2) *= 1.0 / (2 * v.size());
   }

   return(allele_matrix);
}


arma::mat update_all_frequencies_tibble(const std::vector< Fish >& pop,
                                        const NumericVector& markers,
                                        const std::vector<int>& founder_labels) {

    //Rcout << "this is update_all_frequencies_tibble\n";
    int number_of_alleles = founder_labels.size();
    Rcout << number_of_alleles << "\n";
    arma::mat output(markers.size() * number_of_alleles, 3);

    for(int i = 0; i < markers.size(); ++i) {
        //Rcout << "collect local_mat\n";
        arma::mat local_mat = update_frequency_tibble(pop,
                                                 markers[i],
                                                 founder_labels);
       // now we have a (markers x alleles) x 3 tibble, e.g. [loc, anc, freq]
       // and we have to put that in the right place in the output matrix
        //Rcout << "now we feed local mat to output:\n";
        int start = i * number_of_alleles;
        int end = start + number_of_alleles;
        for(int j = start; j < end; ++j) {
           for(int k = 0; k < 3; ++k) {
              // Rcout << j << "\t" << k << "\t" << j - start << "\n";
               output(j, k) = local_mat(j - start, k);
           }
        }
    }
   return(output);
}

arma::mat update_all_frequencies(const std::vector< Fish >& pop,
                                 const NumericVector& markers,
                                 int number_of_alleles) {

    arma::mat output(markers.size(), number_of_alleles);

    for(int i = 0; i < markers.size(); ++i) {
        NumericVector v = update_frequency(pop,
                                           markers[i],
                                           number_of_alleles);
        for(int j = 0; j < v.size(); ++j) {
            output(i, j) = v(j);
        }
    }
    return(output);
}

double calc_mean_junctions(const std::vector< Fish> & pop) {

    double mean_junctions = 0.0;
    for(auto it = pop.begin(); it != pop.end(); ++it) {
        mean_junctions += (*it).chromosome1.size() - 2; // start and end don't count
        mean_junctions += (*it).chromosome2.size() - 2;
    }
    mean_junctions *= 1.0 / (pop.size() * 2); // diploid

    return(mean_junctions);
}

int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness) {

    if(maxFitness <= 0.0) {
        Rcout << "maxFitness = " << maxFitness << "\n";
        Rcpp::stop("Cannot draw fitness if maxFitness <= 0");
        return(-1);
    }

    if(maxFitness > 10000.0) {
        Rcout << "maxFitness = " << maxFitness << "\n";
        Rcpp::stop("It appears maxfitness has encountered a memory access violation\n");
        return(-1);
    }

    for(int i = 0; i < 1e6; ++i) {
        int index = random_number(fitness.size());
        double prob = 1.0 * fitness[index] / maxFitness;
        if(uniform() < prob) {
            return index;
        }
    }
    Rcout << maxFitness << "\n";
    Rcpp::stop("ERROR!Couldn't pick proportional to fitness");
    return -1;
}




std::vector< Fish > convert_NumericVector_to_fishVector(const NumericVector v) {
    std::vector< Fish > output;

    Fish temp;
    int indic_chrom = 1;
    bool add_indiv = false;

    for(int i = 0; i < v.size(); i += 2) {
        junction temp_j;
        temp_j.pos = v[i];
        temp_j.right = v[i+1];

        if(indic_chrom == 1) {
            temp.chromosome1.push_back(temp_j);
        } else {
            temp.chromosome2.push_back(temp_j);
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
            temp.chromosome1.clear();
            temp.chromosome2.clear();
        }
    }

    return(output);
}

List convert_to_list(const std::vector<Fish>& v) {
    int list_size = (int)v.size();
    List output(list_size);

    for(int i = 0; i < v.size(); ++i) {

        Fish focal = v[i];

        NumericMatrix chrom1(focal.chromosome1.size(), 2); // nrow = number of junctions, ncol = 2
        for(int j = 0; j < focal.chromosome1.size(); ++j) {
            chrom1(j, 0) = focal.chromosome1[j].pos;
            chrom1(j, 1) = focal.chromosome1[j].right;
        }

        NumericMatrix chrom2(focal.chromosome2.size(), 2); // nrow = number of junctions, ncol = 2
        for(int j = 0; j < focal.chromosome2.size(); ++j) {
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
        if((*it).pos > pos) {
            if((*(it-1)).right == anc) num_alleles[focal_marker]++;
            focal_marker++;
            if(focal_marker >= number_of_markers) {
                break;
            }
            pos = select(focal_marker, 0);
            anc = select(focal_marker, 4);
        }
        if(anc < 0) break; // these entries are only for tracking alleles over time, not for selection calculation
    }

    focal_marker = 0;
    pos = select(focal_marker, 0);
    anc = select(focal_marker, 4);

    for(auto it = (focal.chromosome2.begin()+1); it != focal.chromosome2.end(); ++it) {
        if((*it).pos > pos) {
            if((*(it-1)).right == anc) num_alleles[focal_marker]++;
            focal_marker++;
            if(focal_marker >= number_of_markers) {
                break;
            }
            pos = select(focal_marker, 0);
            anc = select(focal_marker, 4);
        }
        if(anc < 0) break; // these entries are only for tracking alleles over time, not for selection calculation
    }

    double fitness = 0.0;
    if (multiplicative_selection) fitness = 1.0;
    for(int i = 0; i < num_alleles.size(); ++i) {
        if(select(i, 4) < 0) break; // these entries are only for tracking alleles over time, not for selection calculation

        int fitness_index = 1 + num_alleles[i];
        if(multiplicative_selection) {
            fitness *= select(i, fitness_index);
        } else {
            fitness += select(i, fitness_index);
        }
    }

    return(fitness);
}

int draw_random_founder(const std::vector<double>& v) {
    double r = uniform();
    for(int i = 0; i < v.size(); ++i) {
        r -= v[i];
        if(r <= 0) {
            return(i);
        }
    }
    return(v.back());
}
