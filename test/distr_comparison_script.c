//
//  distr_script.c
//  
//
//  Created by Jordan Eizenga on 1/17/16.
//
//

#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "nanopore_hdp.h"

#define CYTOSINE_MODS "CHM"

typedef enum DistributionMetric {
    SYMMETRIC_KL_DIVERGENCE,
    HELLINGER,
    L2,
    SHANNON_JENSEN
} DistributionMetric;

void print_help() {
    fprintf(stderr, "Usage: distr_comparison_script nanopore_hdp_file output_file (options)\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-h print this message and exit.\n");
    fprintf(stderr, "\t-d distance metric (kl,hellinger,l2,shannon) [default: hellinger]\n");
    fprintf(stderr, "\t-n nucleotides to compare [default: CHM]\n");
    fprintf(stderr, "\t-c canonical comparison nucleotide (will appear in output) [default: first listed in -n]\n");
    fprintf(stderr, "\t-m max number of nucleotide combos to compare [default: 1]\n");
}

char* kmer_from_index(int index, char* alphabet, int alphabet_size, int kmer_length) {
    char* kmer = (char*) malloc(sizeof(char) * kmer_length);
    int remainder = index;
    for (int i = kmer_length - 1; i >= 0; i--) {
        kmer[i] = alphabet[remainder % alphabet_size];
        remainder /= alphabet_size;
    }
    return kmer;
}

//int power(int n, int k) {
//    int res = 1;
//    for (int i = 0; i < k; i++) {
//        res *= n;
//    }
//    return res;
//}

int choose(int n, int r) {
    int lim = n - r;
    if (lim < r) {
        int tmp = r;
        r = lim;
        lim = tmp;
    }
    
    int res = 1;
    for (int val = n; val > lim; val--) {
        res *= val;
    }
    for (int val = 2; val <= r; val++) {
        res /= val;
    }
    return res;
}

bool* choice_index_to_binary_vector(int index, int n, int r) {
    index++;
    bool* binary_vector = (bool*) malloc(sizeof(bool) * n);
    int ncr;
    
    for (int i = 0; i < n ; i++) {
        if (n - i == r) {
            for (; i < n; i++) {
                binary_vector[i] = true;
            }
            break;
        }
        ncr = choose((n - i - 1), r);
        if (index > ncr) {
            binary_vector[i] = true;
            r--;
            index -= ncr;
            
        }
        else {
            binary_vector[i] = false;
        }
    }
    
    return binary_vector;
}


void count_occurrences(char* kmer, int kmer_length, char* count_set, int set_size,
                       int* count_out, int** positions_out) {
    
    int count = 0;
    for (int i = 0; i < kmer_length; i++) {
        for (int j = 0; j < set_size; j++) {
            if (kmer[i] == count_set[j]) {
                count++;
                break;
            }
        }
    }
    *count_out = count;
    
    int* positions = (int*) malloc(sizeof(int) * count);
    *positions_out = positions;
    
    if (count == 0) {
        return;
    }
    
    int k = 0;
    for (int i = 0; i < kmer_length; i++) {
        for (int j = 0; j < set_size; j++) {
            if (kmer[i] == count_set[j]) {
                positions[k] = i;
                k++;
                break;
            }
        }
    }
}


double mean_pairwise_distance(char* restricted_kmer, int kmer_length, int* positions, int num_positions,
                              char* compare_set, int set_size, NanoporeDistributionMetricMemo* metric_memo) {
    
    double sum_distance = 0.0;
    int num_nuclotide_combos = power(set_size, num_positions);
    int num_pairs = choose(num_nuclotide_combos, 2);
    
    char* query_kmer_1 = (char*) malloc(sizeof(char) * kmer_length);
    char* query_kmer_2 = (char*) malloc(sizeof(char) * kmer_length);
    for (int i = 0; i < kmer_length; i++) {
        query_kmer_1[i] = restricted_kmer[i];
        query_kmer_2[i] = restricted_kmer[i];
    }
    
    char* substitutions_1;
    char* substitutions_2;
    for (int i = 1; i < num_nuclotide_combos; i++) {
        substitutions_1 = kmer_from_index(i, compare_set, set_size, num_positions);
        
        for (int k = 0; k < num_positions; k++) {
            query_kmer_1[positions[k]] = substitutions_1[k];
        }
        
        for (int j = 0; j < i; j++) {
            substitutions_2 = kmer_from_index(j, compare_set, set_size, num_positions);
            
            for (int k = 0; k < num_positions; k++) {
                query_kmer_2[positions[k]] = substitutions_2[k];
            }
            
            sum_distance += get_kmer_distr_distance(metric_memo, query_kmer_1, query_kmer_2);
            
            free(substitutions_2);
        }
        free(substitutions_1);
    }
    
    free(query_kmer_1);
    free(query_kmer_2);
    
    return sum_distance / num_pairs;
}

char* get_labeled_kmer(char* restricted_kmer, int kmer_length, int* positions, int num_positions) {
    
    char* labeled_kmer = (char*) malloc(sizeof(char) * (kmer_length + num_positions + 1));
    
    int j = 0;
    for (int i = 0; i < kmer_length; i++) {
        labeled_kmer[i + j] = restricted_kmer[i];
        if (i == positions[j]) {
            j++;
            labeled_kmer[i + j] = '*';
        }
    }
    labeled_kmer[kmer_length + num_positions] = '\0';
    return labeled_kmer;
}

int main(int argc, char** argv) {
    
    DistributionMetric metric = HELLINGER;
    char* comparison_nucleotides = CYTOSINE_MODS;
    char root_comparison_nucleotide = '\0';
    int max_combo_size = 1;
    
    int flag;
    while ((flag = getopt(argc, argv, "hdncm")) != -1) {
        switch (flag) {
            case 'h':
                print_help();
                exit(EXIT_SUCCESS);
                
            case 'd':
                if (strcmp(optarg, "kl") == 0) {
                    metric = SYMMETRIC_KL_DIVERGENCE;
                }
                else if (strcmp(optarg, "hellinger") == 0) {
                    metric = HELLINGER;
                }
                else if (strcmp(optarg, "shannon") == 0) {
                    metric = SHANNON_JENSEN;
                }
                else if (strcmp(optarg, "l2") == 0) {
                    metric = L2;
                }
                else {
                    fprintf(stderr, "ERROR: Unrecognized distribution metric '%s'.\n", optarg);
                    fprintf(stderr, "Accepted metrics are: kl, hellinger, l2, shannon.\n");
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 'n':
                comparison_nucleotides = optarg;
                break;
                
            case 'c':
                sscanf(optarg, "%c", &root_comparison_nucleotide);
                break;
                
            case 'm':
                sscanf(optarg, "%d", &max_combo_size);
                break;
                
            default:
                fprintf(stderr, "ERROR: Unrecognized option '-%c'.\n", flag);
                print_help();
                exit(EXIT_FAILURE);
                break;
        }
    }
    
    // get hdp and pull its attributes
    NanoporeHDP* nhdp = deserialize_nhdp(argv[optind]);
    char* alphabet = get_nanopore_hdp_alphabet(nhdp);
    int alphabet_size = get_nanopore_hdp_alphabet_size(nhdp);
    int kmer_length = get_nanopore_hdp_kmer_length(nhdp);
    
    int num_comparison_nucleotides = strlen(comparison_nucleotides);
    
    // validate comparison nucleotides belong to same alphabet
    for (int i = 0; i < num_comparison_nucleotides; i++) {
        for (int j = 0; j < alphabet_size; j++) {
            if (comparison_nucleotides[i] == alphabet[j]) {
                break;
            }
            else if (j == alphabet_size - 1) {
                fprintf(stderr, "ERROR: Comparison nucleotide '%c' not included in alphabet '%s'.\n",
                        comparison_nucleotides[j], alphabet);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if (root_comparison_nucleotide == '\0') {
        // default to first nucleotide from comparison group
        root_comparison_nucleotide = comparison_nucleotides[0];
    }
    else {
        // else verify that provided nucleotide is in comparison group
        bool in_comparison_nucleotides = false;
        for (int i = 0; i < num_comparison_nucleotides; i++) {
            if (root_comparison_nucleotide == comparison_nucleotides[i]) {
                in_comparison_nucleotides = true;
                break;
            }
        }
        if (!in_comparison_nucleotides) {
            fprintf(stderr, "ERROR: root comparison nucleotide '%c' not in comparison group '%s'.\n",
                    root_comparison_nucleotide, comparison_nucleotides);
        }
    }
    
    FILE* output_file = fopen(argv[optind + 1], "w");
    
    NanoporeDistributionMetricMemo* metric_memo;
    char* metric_col_name;
    
    switch (metric) {
        case SYMMETRIC_KL_DIVERGENCE:
            metric_memo = new_nhdp_kl_divergence_memo(nhdp);
            metric_col_name = "mean_symmetric_KL_divergence";
            break;
        case HELLINGER:
            metric_memo = new_nhdp_hellinger_distance_memo(nhdp);
            metric_col_name = "mean_Hellinger_distance";
            break;
        case L2:
            metric_memo = new_nhdp_l2_distance_memo(nhdp);
            metric_col_name = "mean_L2_distance";
            break;
        case SHANNON_JENSEN:
            metric_memo = new_nhdp_shannon_jensen_distance_memo(nhdp);
            metric_col_name = "mean_Shannon-Jensen_distance";
            break;
            
        default:
            fprintf(stderr, "ERROR: could not create distance metric memo.\n");
            exit(EXIT_FAILURE);
            break;
    }
    
    // create restricted alphabet with only one of the comparison nucleotides
    int restricted_alphabet_size = alphabet_size - num_comparison_nucleotides + 1;
    char* restricted_alphabet = (char*) malloc(sizeof(char) * (restricted_alphabet_size + 1));
    
    int k = 0;
    bool in_restricted;
    for (int i = 0; i < alphabet_size; i++) {
        in_restricted = true;
        for (int j = 0; j < num_comparison_nucleotides; j++) {
            if (alphabet[i] == comparison_nucleotides[j] &&
                comparison_nucleotides[j] != root_comparison_nucleotide) {
                in_restricted = false;
                break;
            }
        }
        if (in_restricted) {
            restricted_alphabet[k] = alphabet[i];
            k++;
        }
    }
    
    // add header
    fprintf(output_file, "#k-mer\t%s\n", metric_col_name);
    
    int num_restricted_kmers = power(restricted_alphabet_size, kmer_length);
    
    char* labeled_kmer;
    double mean_dist;
    int* filtered_positions;
    bool* which_positions;
    
    for (int kmer_num = 0; kmer_num < num_restricted_kmers; kmer_num++) {
        char* kmer = kmer_from_index(kmer_num, restricted_alphabet, restricted_alphabet_size, kmer_length);
        
        int count;
        int* positions;
        
        count_occurrences(kmer, kmer_length, comparison_nucleotides, num_comparison_nucleotides, &count, &positions);
        
        // only consider kmers that contain comparison set
        if (count == 0) {
            continue;
        }
        // not enough varying sites to need to break them up
        else if (count <= max_combo_size) {
            labeled_kmer = get_labeled_kmer(kmer, kmer_length, positions, count);
            mean_dist = mean_pairwise_distance(kmer, kmer_length, positions, count, comparison_nucleotides,
                                               num_comparison_nucleotides, metric_memo);
            
            fprintf(output_file, "%s\t%lf\n", labeled_kmer, mean_dist);
            
            free(labeled_kmer);
        }
        // break up into all combos less than max size
        else {
            int num_combos = choose(count, max_combo_size);
            filtered_positions = (int*) malloc(sizeof(int) * max_combo_size);
            
            for (int i = 0; i < num_combos; i++) {
                // get which sites should be included for this combo index
                which_positions = choice_index_to_binary_vector(i, count, max_combo_size);
                
                // limit position array to these sites
                int j = 0, k = 0;
                while (j < max_combo_size) {
                    if (which_positions[k]) {
                        filtered_positions[j] = positions[k];
                        j++;
                    }
                    k++;
                }
                free(which_positions);
                
                labeled_kmer = get_labeled_kmer(kmer, kmer_length, filtered_positions, max_combo_size);
                mean_dist = mean_pairwise_distance(kmer, kmer_length, filtered_positions, max_combo_size, comparison_nucleotides,
                                                   num_comparison_nucleotides, metric_memo);
                
                fprintf(output_file, "%s\t%lf\n", labeled_kmer, mean_dist);
                free(labeled_kmer);
            }
            
            free(filtered_positions);
        }
        
        free(kmer);
        free(positions);
    }
    
}
