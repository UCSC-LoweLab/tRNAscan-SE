/* eufindtRNA  - Eukaryotic tRNA finder
 *
 * pavesi.h - functions for finding transcriptional control regions 
 *
 * C implementation of algorithm described by Pavesi, Conterio, 
 * Bolchi, Dieci, & Ottonello in NAR 22:1247-56 (94)
 * "Identification of new eukaryotic tRNA genes in genomic DNA
 * databases by a multistep weight matix analysis of transcriptional
 * control regions"
 *
 * To be used in tRNAscan-SE package to increase sensitivity by
 * complementing tRNAscan 1.3 first-pass scan
 *
 * by Todd MJ Lowe    4/8/96
 *
 * Uses Sean Eddy's function library for biological sequence analysis
 * (Squid v1.5g)
 *
 */

#include "squid.h"
#include "eufind_const.h"

void Init_tRNA(TRNA_TYPE *tRNA);

int IntEncodeSeq (char *intseq, char *seq, int seqlen);

int GetBbox (float *score, int *seqidx, char *iseq, int seqlen, int strand, int verbose);

float  Get_ABdist_weight(int ABdist);

int GetSecABox(TRNA_TYPE *tRNA, char *seq);

void GetBestABox (TRNA_TYPE *tRNA, char *seq, char *iseq, int seqlen, int strand, int verbose, int Max_AB_dist, int prev_Abox_st);

int GetBestTrxTerm (TRNA_TYPE *tRNA, char *seq, int seqlen, float TermPenalty);

void Get_IsoType (TRNA_TYPE *tRNA);

void Get_anticodon (TRNA_TYPE *tRNA, char *seq);

void Get_tRNA_stats (TRNA_TYPE *tRNA, char *seq, int seqlen, int strand);

void Save_tRNA (TRNA_TYPE *tRNA, SQINFO *sqinfo, char *seq, int strand, int ShowScores, long int sqoffset);

int tRNAOverlap (TRNA_TYPE *tRNA1, TRNA_TYPE *tRNA2, int strand);
