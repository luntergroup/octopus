/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Dec 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
 ******************************************************************************************************************/

#ifndef align_hpp
#define align_hpp

#include <cstdint>

int fastAlignmentRoutine(const char* seq1, const char* seq2, const char* qual2, int len1, int len2,
                         int gapextend, int nucprior, const std::int8_t* localgapopen);

int fastAlignmentRoutine(const char* seq1, const char* seq2, const char* qual2, int len1, int len2,
                         int gapextend, int nucprior, const std::int8_t* localgapopen,
                         char* aln1, char* aln2, int* firstpos);

int calculateFlankScore(int hapLen, int leftHapFlank, int rightHapFlank, const char* quals,
                        const std::int8_t* localGapOpen, int gapExtend, int nucprior,
                        int firstpos, const char* aln1, const char* aln2);

#endif
