#include <iostream>

#include "simd_pair_hmm.hpp"

using namespace std;
using namespace octopus::hmm::simd;

const int MAXLEN = 256;

struct TestCase {
    string name;
    char truth[MAXLEN];
    char target[MAXLEN];
    int8_t quals[MAXLEN];
    int8_t gap_open[MAXLEN];
    short gap_extend;
    short nuc_prior;
    int first_pos;
    char aln1[2 * MAXLEN];
    char aln2[2 * MAXLEN];
    int score;
    int bandwidth = 8;
};

TestCase case1 = {
    "case 1",
    "ACGTACGTACGTACGAAAA",
    "AAAA",
    {40,40,40,40,-1},
    {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,-1},
    1,
    4,
    15,
    "AAAA",
    "AAAA",
    0
};

TestCase case2 = {
    "case 2",
    "ACGTACGTACGTACGAATA",
    "AAAA",
    {40,40,40,40,-1},
    {90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,-1},
    1,
    4,
    15,
    "AATA",
    "AAAA",
    40
};

TestCase case3 = {
    "case 3",
    "ACGTACGTACGTACGAAGC",
    "CGGC",
    {40,40,40,40,-1},
    {90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,70,90,90,90,-1},
    1,
    4,
    13,
    "CGAAGC",
    "CG--GC",
    71
};

TestCase case4 = {
    "case 4",
    "CGAAGCACGTACGTACGTA",
    "CGGC",
    {40,40,40,40,-1},
    {90,90,70,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,-1},
    1,
    4,
    0,
    "CGAAGC",
    "CG--GC",
    71
};

TestCase testcases[] = {case1, case2, case3, case4};



int ilen( int8_t* a ) {
    int i=0;
    while (a[i] != -1) i++;
    return i;
}

int runTest( TestCase& testcase ) {
    
    size_t truth_len = std::strlen( testcase.truth );
    size_t target_len = std::strlen( testcase.target );
    char aln1[ 2*truth_len + 1 ];
    char aln2[ 2*truth_len + 1 ];
    int first_pos;
    int fail = 0;

    cout << "Test " << testcase.name;
    
    if (truth_len != target_len + 2 * testcase.bandwidth - 1) {
        cout << endl << "Truth and target lengths are " << truth_len << " " << target_len << " and do not differ by " <<
            2 * testcase.bandwidth - 1 << endl;
        return 1;
    }

    if (testcase.quals[ target_len ] != -1) {
        cout << endl << "Qualities has wrong length " << ilen(testcase.quals) << " expected " << target_len << endl;
        return 1;
    }

    if (testcase.gap_open[ truth_len ] != -1) {
        cout << endl << "Gap openings has wrong length " << ilen(testcase.gap_open) << " expected " << truth_len << endl;
        return 1;
    }

    int score = align( testcase.truth,
                       testcase.target,
                       testcase.quals,
                       truth_len,
                       target_len,
                       testcase.gap_open,
                       testcase.gap_extend,
                       testcase.nuc_prior,
                       first_pos,
                       aln1,
                       aln2 );

    int i = 0;
    for (; aln1[i] && testcase.aln1[i]; i++) {
        if (aln1[i] != testcase.aln1[i] || aln2[i] != testcase.aln2[i]) 
            fail = 1;
    }
    if (aln1[i] != testcase.aln1[i]) 
        fail = 1;
    
    if (fail) {
        cout << endl << "Alignments mismatch.  Expected" << endl << testcase.aln1 << endl << testcase.aln2 << endl;
        cout << "Found" << endl << aln1 << endl << aln2;
    }

    if ( first_pos != testcase.first_pos ) {
        cout << endl << "First pos expected " << testcase.first_pos << " got " << first_pos;
        fail = 1;
    }

    if ( score != testcase.score ) {
        cout << endl << "Score expected " << testcase.score << " got " << score << endl;
        fail = 1;
    }

    if (fail) {
        cout << endl;
    } else {
        cout << " Success" << endl;
    }
    return (fail);
}



int main()
{
    for (auto& tc : testcases) {
        runTest( tc );
    }
}
