//
//  nanopore_hdp_tests.c
//  
//
//  Created by Jordan Eizenga on 1/12/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "CuTest.h"
#include "nanopore_hdp.h"

void test_first_kmer_index(CuTest* ct) {
    int64_t kmer_id = 0;
    int64_t length = 5;
    int64_t alphabet_size = 4;
    int64_t* kmer = get_word(kmer_id, alphabet_size, length);
    
    for (int64_t i = 0; i < length; i++) {
        CuAssertIntEquals(ct, 0, kmer[i]);
    }
    
    free(kmer);
}

void test_second_kmer_index(CuTest* ct) {
    int64_t kmer_id = 1;
    int64_t length = 5;
    int64_t alphabet_size = 4;
    int64_t* kmer = get_word(kmer_id, alphabet_size, length);
    
    for (int i = 0; i < length - 1; i++) {
        CuAssertIntEquals(ct, 0, kmer[i]);
    }
    CuAssertIntEquals(ct, 1, kmer[length - 1]);
    
    free(kmer);
}

void test_sixth_kmer_index(CuTest* ct) {
    int64_t kmer_id = 6;
    int64_t length = 5;
    int64_t alphabet_size = 4;
    int64_t* kmer = get_word(kmer_id, alphabet_size, length);
    
    for (int64_t i = 0; i < length - 2; i++) {
        CuAssertIntEquals(ct, 0, kmer[i]);
    }
    CuAssertIntEquals(ct, 1, kmer[length - 2]);
    CuAssertIntEquals(ct, 2, kmer[length - 1]);
    
    free(kmer);
}

void test_multiset_creation(CuTest* ct) {
    int64_t length = 6;
    int64_t alphabet_size = 4;
    int64_t* multiset_1 = get_word_multiset(1, alphabet_size, length);
    int64_t* multiset_2 = get_word_multiset(4, alphabet_size, length);
    int64_t* multiset_3 = get_word_multiset(16, alphabet_size, length);
    
    for (int64_t i = 0; i < length; i++) {
        CuAssertIntEquals(ct, multiset_1[i], multiset_2[i]);
        CuAssertIntEquals(ct, multiset_2[i], multiset_3[i]);
    }
    
    free(multiset_1);
    free(multiset_2);
    free(multiset_3);
}

void test_word_id_to_multiset_id(CuTest* ct) {
    int64_t length = 8;
    int64_t alphabet_size = 4;
    
    CuAssertIntEquals(ct, word_id_to_multiset_id(0, alphabet_size, length), 0);
    CuAssertIntEquals(ct, word_id_to_multiset_id(1, alphabet_size, length), 1);
    CuAssertIntEquals(ct, word_id_to_multiset_id(2, alphabet_size, length), 2);
    CuAssertIntEquals(ct, word_id_to_multiset_id(3, alphabet_size, length), 3);
    CuAssertIntEquals(ct, word_id_to_multiset_id(4, alphabet_size, length), 1);
    CuAssertIntEquals(ct, word_id_to_multiset_id(5, alphabet_size, length), 4);
    CuAssertIntEquals(ct, word_id_to_multiset_id(6, alphabet_size, length), 5);
    CuAssertIntEquals(ct, word_id_to_multiset_id(7, alphabet_size, length), 6);
    CuAssertIntEquals(ct, word_id_to_multiset_id(8, alphabet_size, length), 2);
    CuAssertIntEquals(ct, word_id_to_multiset_id(10, alphabet_size, length), 7);
    CuAssertIntEquals(ct, word_id_to_multiset_id(11, alphabet_size, length), 8);
    CuAssertIntEquals(ct, word_id_to_multiset_id(12, alphabet_size, length), 3);
    CuAssertIntEquals(ct, word_id_to_multiset_id(13, alphabet_size, length), 6);
    CuAssertIntEquals(ct, word_id_to_multiset_id(14, alphabet_size, length), 8);
    CuAssertIntEquals(ct, word_id_to_multiset_id(15, alphabet_size, length), 9);
    CuAssertIntEquals(ct, word_id_to_multiset_id(16, alphabet_size, length), 1);
    
}

void test_kmer_id(CuTest* ct) {
    CuAssertIntEquals(ct, kmer_id("AAAC", "ACGT", 4, 4), 1);
    CuAssertIntEquals(ct, kmer_id("AAAT", "ACGT", 4, 4), 3);
    CuAssertIntEquals(ct, kmer_id("AAAT", "ACT", 3, 4), 2);
    CuAssertIntEquals(ct, kmer_id("GGGG", "ABCDEFG", 7, 4), power(7, 4) - 1);
    CuAssertIntEquals(ct, standard_kmer_id("AACAA", 5), 16);
}

void create_flat_model(CuTest* ct) {
    printf("flat\n");
    NanoporeHDP* nhdp = flat_hdp_model("ACGT", 4, 6, 1.0, 1.0, -1.0 ,1.0, 10,
                                       "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_model.model");
    update_nhdp_from_alignment(nhdp, "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_alignment.tsv",
                               false);
    
    CuAssert(ct, "did not complete nhdp construction\n", true);
    destroy_nanopore_hdp(nhdp);
}

void create_multiset_model(CuTest* ct) {
    printf("multiset\n");
    NanoporeHDP* nhdp = multiset_hdp_model("ACGT", 4, 6, 1.0, 1.0, 1.0, -1.0 ,1.0, 10,
                                       "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_model.model");
    update_nhdp_from_alignment(nhdp, "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_alignment.tsv",
                               false);
    
    CuAssert(ct, "did not complete nhdp construction\n", true);
}

void create_purine_model(CuTest* ct) {
    printf("purine\n");
    NanoporeHDP* nhdp = purine_composition_hdp_model("AG", 2, "CT", 2, 6, 1.0, 1.0, 1.0, -1.0 ,1.0, 10,
                                                     "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_model.model");
    update_nhdp_from_alignment(nhdp, "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_alignment.tsv",
                               false);
    
    CuAssert(ct, "did not complete nhdp construction\n", true);
    destroy_nanopore_hdp(nhdp);
}

void create_middle_nts_model(CuTest* ct) {
    printf("middle\n");
    NanoporeHDP* nhdp = middle_2_nts_hdp_model("ACGT", 4, 6, 1.0, 1.0, 1.0, -1.0 ,1.0, 10,
                                               "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_model.model");
    update_nhdp_from_alignment(nhdp, "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_alignment.tsv",
                               false);
    
    CuAssert(ct, "did not complete nhdp construction\n", true);
    destroy_nanopore_hdp(nhdp);
}


CuSuite* get_suite() {
    
    CuSuite* suite = CuSuiteNew();
    
    SUITE_ADD_TEST(suite, test_first_kmer_index);
    SUITE_ADD_TEST(suite, test_second_kmer_index);
    SUITE_ADD_TEST(suite, test_sixth_kmer_index);
    SUITE_ADD_TEST(suite, test_multiset_creation);
    SUITE_ADD_TEST(suite, test_word_id_to_multiset_id);
    SUITE_ADD_TEST(suite, test_kmer_id);
    SUITE_ADD_TEST(suite, create_flat_model);
    SUITE_ADD_TEST(suite, create_multiset_model);
    SUITE_ADD_TEST(suite, create_purine_model);
    SUITE_ADD_TEST(suite, create_middle_nts_model);
    
    return suite;
}

int main(void) {
    CuSuite* suite = get_suite();
    
    CuString* output = CuStringNew();
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
}








