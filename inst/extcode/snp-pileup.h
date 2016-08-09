#ifndef _MAIN_H_
#define _MAIN_H_

#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h" 
#include "htslib/kstring.h"
#include "htslib/khash_str2int.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"

// bam flags
#define READ_PAIRED 0x1
#define READ_MAPPED_PROPER_PAIR 0x2
#define READ_UNMAPPED 0x4
#define MATE_UNMAPPED 0x8
#define READ_REVERSE_STRAND 0x10
#define MATE_REVERSE_STRAND 0x20
#define READ_FIRST_IN_PAIR 0x40
#define READ_SECOND_IN_PAIR 0x80
#define READ_SECONDARY_ALIGNMENT 0x100
#define READ_FAILED_QUALITY_CHECKS 0x200
#define READ_PCR_OPTICAL_DUPLICATE 0x400
#define READ_SUPPLEMENTARY 0x800
#define READ_PRIMARY 0x900

using namespace std;

struct arguments {
	vector<char*> args;
	bool count_orphans;
	bool gzipped;
	BGZF* gzippedPointer;
	bool ignore_overlaps;
	int min_base_quality;
	int min_map_quality;
	vector<int> min_read_counts;
	int max_depth;
	bool progress;
	int pseudo_snps;
	bool verbose;
	void (*outFunc)(arguments, string, FILE *);
};

struct file_info {
	int refs;
	int alts;
	int errors;
	int deletions;
};

typedef struct { 
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag, all; 
    int rflag_require, rflag_filter; 
    int openQ, extQ, tandemQ, min_support; // for indels 
    double min_frac; // for indels 
    char *reg, *pl_list, *fai_fname, *output_fname; 
    faidx_t *fai; 
    void *bed, *rghash; 
    int argc; 
    char **argv;
} mplp_conf_t; 


#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}
typedef struct { 
    char *ref[2]; 
    int ref_id[2]; 
    int ref_len[2]; 
} mplp_ref_t; 


typedef struct { 
    samFile *fp; 
	struct arguments *conf; 
    bam_hdr_t *h;
    mplp_ref_t *ref;
} mplp_aux_t; 


#endif
