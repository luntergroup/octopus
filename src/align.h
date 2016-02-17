#ifndef ALIGN_H
#define ALIGN_H

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009, Nov 2014
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/
int fastAlignmentRoutine(const char* seq1, const char* seq2, const char* qual2, int len1, int len2, int gapextend, int nucprior,
			 const char* localgapopen, char* aln1, char* aln2, int* firstpos);

int calculateFlankScore(int hapLen, int hapFlank, const char* quals, const char* localGapOpen, int gapExtend, int nucprior,
			int firstpos, const char* aln1, const char* aln2);

#endif
