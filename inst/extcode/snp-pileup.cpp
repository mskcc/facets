#include "snp-pileup.h"

#include <argp.h>

static char args_doc[] = "<vcf file> <output file> <sequence files...>";

static struct argp_option options[] = {
	{"count-orphans", 'A', 0, 0, "Do not discard anomalous read pairs."},
	{"gzip", 'g', 0, 0, "Compresses the output file with BGZF."},
	{"ignore-overlaps", 'x', 0, 0, "Disable read-pair overlap detection."},
	{"min-base-quality", 'Q', "QUALITY", 0, "Sets the minimum threshold for base quality. Default is 0."},
	{"min-map-quality", 'q', "QUALITY", 0, "Sets the minimum threshold for mapping quality. Default is 0."},
	{"min-read-counts", 'r', "READS", 0, "Comma separated list of minimum read counts for a position to be output. Default is 0."},
	{"max-depth", 'd',"DEPTH", 0, "Sets the maximum depth. Default is 4000."},
	{"progress", 'p', 0, 0, "Show a progress bar. WARNING: requires additional time to calculate number of SNPs, and will take longer than normal."},
	{"pseudo-snps", 'P', "MULTIPLE", 0, "Every MULTIPLE positions, if there is no SNP, insert a blank record with the total count at the position."},
	{"verbose", 'v', 0, 0, "Show detailed messages."},
	{ 0 }
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
	struct arguments *arguments = (struct arguments*) state->input;

	switch (key) {
	case 'A':
		arguments->count_orphans = true;
		break;

	case 'd':
		arguments->max_depth = atoi(arg);
		break;

	case 'g':
		arguments->gzipped = true;
		break;

	case 'P':
		arguments->pseudo_snps = atoi(arg);
		break;

	case 'p':
		arguments->progress = true;
		break;

	case 'r':
		char * token;
		token = strtok(arg, ",");
		while (token) {
			arguments->min_read_counts.push_back(atoi(token));
			token = strtok(NULL, ",");
		}
		break;

	case 'Q':
		arguments->min_base_quality = atoi(arg);
		break;

	case 'q':
		arguments->min_map_quality = atoi(arg);
		break;

	case 'x':
		arguments->ignore_overlaps = true;
		break;

	case 'v':
		arguments->verbose = true;
		break;

	case ARGP_KEY_ARG:
		arguments->args.push_back(arg);
		break;

	case ARGP_KEY_END:
		if (state->arg_num < 3) {
			// Not enough arguments.
			argp_usage (state);
		}
		for (int i = 0; i < state->arg_num - 2; i++) {
			arguments->min_read_counts.push_back(0);
		}
		break;

	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp = { options, parse_opt, args_doc, 0 };

static int mplp_func(void *data, bam1_t *b) {
	// it seems to me that this function is run once every read, and lets you control if the read should be skipped or not
	// here is where you would put quality checks and things like that

//	printf("hello\n");
	char *ref;
	mplp_aux_t *ma = (mplp_aux_t*)data;
	int ret, skip = 0, ref_len;
	do {
		int has_ref;
		ret = sam_read1(ma->fp, ma->h, b);
		if (ret < 0) break;

		if (b->core.flag & READ_FAILED_QUALITY_CHECKS || b->core.flag & READ_SECONDARY_ALIGNMENT || b->core.flag & READ_PCR_OPTICAL_DUPLICATE || b->core.flag & READ_UNMAPPED) {
			skip = 1;
			continue;
		}

		if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
			skip = 1;
			continue;
		}

		/*if (ma->conf->bed && ma->conf->all == 0) { // test overlap
			skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
			if (skip) continue;
		}*/

		skip = 0;
		/*if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
		if (has_ref && ma->conf->capQ_thres > 10) {
			int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
			if (q < 0) skip = 1;
			else if (b->core.qual > q) b->core.qual = q;
			} */
		if (b->core.qual < ma->conf->min_map_quality) {
			skip = 1;
		} else if (b->core.qual == 0) {
			skip = 1;
		} else if (!ma->conf->count_orphans && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) {
			skip = 1;
		}
	} while (skip);
	return ret;
}

static int vcf_chr_to_bam(char * chromosome, char ** bam_chrs, int32_t n_targets) {
	// try to remove chr prefix if it exists
	if (!strncmp(chromosome, "chr", 3)) {
		chromosome += 3;
	}
	bool bam_chr_prefix = (!strncmp(bam_chrs[0], "chr", 3));
	for (int i = 0; i < n_targets; i++) {
		if (!bam_chr_prefix) {
			if(!strcmp(bam_chrs[i], chromosome)) {
				return i;
			}
		} else {
			if(!strcmp(bam_chrs[i] + 3, chromosome)) {
				return i;
			}
		}
	}
	return -1;
}

uint64_t get_snp_count(char * file) {
	uint64_t count = 0;
	// load vcf file
	bcf_srs_t * vcfReader = bcf_sr_init();
	if (!bcf_sr_add_reader(vcfReader, file)) {
		printf("Failed to read VCF file: %s\n", bcf_sr_strerror(vcfReader->errnum));
		return 1;
	}
	while ( bcf_sr_next_line(vcfReader))
	{
		count++;
	}

	bcf_sr_destroy(vcfReader);
	return count;
}

void print_output(arguments arguments, string str, FILE * fp) {
	fputs(str.c_str(), fp);
}

void gzip_output(arguments arguments, string str, FILE * fp) {
    if (bgzf_write(arguments.gzippedPointer, str.c_str(), str.length()) < 0) {
		printf("failed to write to file, terminating.\n");
		exit(1);
	}
}

inline bool ends_with(std::string const & value, std::string const & ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

int program_main(arguments arguments) {
	clock_t start = clock();

	int i = 0;
	int n = arguments.args.size() - 2; // n is the number of files. currently hardcoded at 2
	bam_mplp_t iter;
	const bam_pileup1_t ** plp = (const bam_pileup1_t **) calloc(n, sizeof(bam_pileup1_t*));
	int * n_plp = (int*) calloc(n, sizeof(int));
	char *ref;
	int ref_len;
	mplp_ref_t mp_ref = MPLP_REF_INIT;
	struct arguments * conf = &arguments;
	hts_verbose = 1;

	// load vcf file
	bcf_srs_t * vcfReader = bcf_sr_init();
	if (!bcf_sr_add_reader(vcfReader, arguments.args[0])) {
		printf("Failed to read VCF file: %s\n", bcf_sr_strerror(vcfReader->errnum));
		return 1;
	}
	bcf_hdr_t * vcfHdr = vcfReader->readers[0].header;

	uint64_t count = 0;
	if (arguments.progress) {
		printf("Calculating SNP count...");
		cout.flush();
		count = get_snp_count(arguments.args[0]);
		printf("done.\n");
		cout.flush();
	}

	// construct data to pass to pileup engine
	mplp_aux_t **data;
	data = (mplp_aux_t **) calloc(n, sizeof(mplp_aux_t*)); // allocate memory for data
	for (i = 0; i < n; ++i) {
		data[i] = (mplp_aux_t*) calloc(1, sizeof(mplp_aux_t));

		// open file
		hFILE *hfp = hopen(arguments.args[i + 2], "r");
		htsFormat fmt;

		if (!hfp) {
			printf("Failed to read sequence file %s (%s).\n", arguments.args[i + 2], strerror(errno));
			return 1;
		}

		hts_detect_format(hfp, &fmt);

		if (arguments.verbose) {
			printf("Detected format for file %d: %s\n", i + 1, hts_format_description(&fmt));
		}

		data[i]->fp = hts_hopen(hfp, arguments.args[i + 2], "rb");
		if (!data[i]->fp) {
			// couldn't open file, laugh at user
			printf("Couldn't open sequence file %s! (%s)\n", arguments.args[i + 2], strerror(errno));
			return 1;
		}

		hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0); // not sure what this does but samtools does it so I guess it's important

		data[i]->conf = conf;
		data[i]->h = sam_hdr_read(data[i]->fp);
		data[i]->ref = &mp_ref;
	}

	// read first header
	bam_hdr_t * hdr;
	hdr = data[0]->h;

	// start pileup engine
	iter = bam_mplp_init(n, mplp_func, (void**)data);
	if (!arguments.ignore_overlaps) {
		bam_mplp_init_overlaps(iter);
	}
	int max_depth = arguments.max_depth;
	bam_mplp_set_maxcnt(iter, max_depth);

	if (arguments.verbose) {
		printf("Max per-file depth set to %d.\n", max_depth);
	}

	string fname = string(arguments.args[1]);
	if (arguments.gzipped) {
		if (!ends_with(fname, ".gz")) {
			fname += ".gz";
		}
		arguments.outFunc = gzip_output;
	}
	// check if output exists
	FILE * test_output = fopen(fname.c_str(), "r");
	if (test_output) {
		printf("Output file %s already exists!\n", arguments.args[1]);
		fclose(test_output);
		return 1;
	}
	// DON'T CLOSE test_output HERE because if you are here, test_output is null
	FILE * output_file = NULL;
	if (arguments.gzipped) {
		arguments.gzippedPointer = bgzf_open(fname.c_str(), "w+");
	} else {
		// open output file
		output_file = fopen(fname.c_str(), "w+");
		if (!output_file) {
			printf("Failed to open output file for writing: %s\n", strerror(errno));
			return 1;
		}
	}

	ostringstream output;
	// output header to file
	output << "Chromosome,Position,Ref,Alt";
	for (i = 0; i < n; ++i) {
		output << ",File";
		output << (i + 1);
		output << "R,File";
		output << (i + 1);
		output << "A,File";
		output << (i + 1);
		output << "E,File";
		output << (i + 1);
		output << "D";
	}
	output << "\n";
	(arguments.outFunc)(arguments, output.str(), output_file);

	// go through it
	int ret;
	int tid, pos, vcf_tid;
	vector<file_info> f_info = vector<file_info>();
	bool first = true;
	uint64_t current_count = 0;
	float last_progress = 100.0;
	bool have_snpped = false;
	while ( bcf_sr_next_line(vcfReader) ) {
		bcf1_t *line = vcfReader->readers[0].buffer[0];
		vcf_tid = vcf_chr_to_bam((char*) vcfHdr->id[BCF_DT_CTG][line->rid].key, hdr->target_name, hdr->n_targets);
		current_count++;

//		printf("hi\n");
		float progress = ((float)current_count/count);
		if (arguments.progress && (int)(last_progress*100) != (int)(progress*100)) {
			int barWidth = 70;

			cout << "[";
			int pos = barWidth * progress;
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos) cout << "=";
				else if (i == pos) cout << ">";
				else cout << " ";
			}
			cout << "] " << int(progress * 100.0) << " %\r";
			cout.flush();

			last_progress = progress;
		}
		if (line->n_allele != 2 || strlen(line->d.allele[0]) != 1 || strlen(line->d.allele[1]) != 1) {
			continue;
		}
		//printf("vcf: chrom %s, pos %d, id %s, ref %s, alt %s\n", vcfHdr->id[BCF_DT_CTG][line->rid].key, line->pos + 1, line->d.id, line->d.allele[0], line->d.allele[1]);
		if (!first && tid > vcf_tid) {
			// keep going...
			continue;
		}
		if (!first && tid == vcf_tid && line->pos < pos) {
			//printf("line->pos < pos, skipping\n");
			continue;
		}
		while ( (ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
			if (first) {
				first = false;
			}
			//printf("position %d, chr %s, n_plp %d, file %d\n", pos + 1, hdr->target_name[tid], n_plp[i], i);
			have_snpped = false;
			if (vcf_tid > tid) {
				// snp's tid comes after where I am at
				// skip ahead
				continue;
			}
			if (tid == vcf_tid && pos == line->pos) {
				// we found what we're looking for!
				// output it to the file
				bool is_not_zero = false;
				bool fails_min = false;
				for (i = 0; i < n; ++i) {
					// this is once for each file
					file_info this_file = file_info();
					this_file.refs = 0;
					this_file.alts = 0;
					this_file.errors = 0;
					this_file.deletions = 0;
					if (n_plp[i] >= arguments.min_read_counts[i]) {
						for (int j = 0; j < n_plp[i]; ++j) {
							const bam_pileup1_t *p = plp[i] + j;
							int c = p->qpos < p->b->core.l_qseq
								? bam_get_qual(p->b)[p->qpos]
								: 0;
							if (c == 0) {
								continue; // no
							}
							if (c < arguments.min_base_quality) {
								continue; // skip anything with quality below threshold
							}
							// this is once for each nucleotide
							if (p->is_del) {
								this_file.deletions++;
							} else if (seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)] == line->d.allele[0][0]) {
								this_file.refs++;
							} else if (seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)] == line->d.allele[1][0]) {
								this_file.alts++;
							} else {
								this_file.errors++;
							}
						}
					} else {
						fails_min = true;
					}
					if (this_file.refs != 0 || this_file.alts != 0 || this_file.errors != 0) {
						is_not_zero = true;
					}
					f_info.push_back(this_file);
				}
				if (is_not_zero && !fails_min) {
					output.str("");
					output.clear();
					output << hdr->target_name[tid];
					output << ",";
					output << (line->pos + 1);
					output << ",";
					output << line->d.allele[0][0];
					output << ",";
					output << line->d.allele[1][0];
					for (i = 0; i < n; i++) {
						output << ",";
						output << f_info[i].refs;
						output << ",";
						output << f_info[i].alts;
						output << ",";
						output << f_info[i].errors;
						output << ",";
						output << f_info[i].deletions;
					}
					output << "\n";
					(arguments.outFunc)(arguments, output.str(), output_file);
				}
				f_info.clear();
				have_snpped = true;
				break;
			}
			if (arguments.pseudo_snps && !have_snpped && (((pos + 1) % arguments.pseudo_snps) == 0)) {
				bool is_not_zero = false;
				bool fails_min = false;
				for (i = 0; i < n; i++) {
					file_info this_file = file_info();
					this_file.refs = 0;
					this_file.alts = 0;
					this_file.errors = 0;
					this_file.deletions = 0;
					if (n_plp[i] >= arguments.min_read_counts[i]) {
						for (int j = 0; j < n_plp[i]; ++j) {
							const bam_pileup1_t *p = plp[i] + j;
							int c = p->qpos < p->b->core.l_qseq
								? bam_get_qual(p->b)[p->qpos]
								: 0;
							if (c == 0) {
								continue; // no
							}
							if (c < arguments.min_base_quality) {
								continue; // skip anything with quality below threshold
							}
							this_file.refs++;
						}
					} else {
						fails_min = true;
					}
					if (this_file.refs != 0) {
						is_not_zero = true;
					}
					f_info.push_back(this_file);
				}
				// add a pseudo snp!
				if (is_not_zero && !fails_min) {
					output.str("");
					output.clear();
					output << hdr->target_name[tid];
					output << ",";
					output << (pos + 1);
					output << ",.,.";
					for (i = 0; i < n; i++) {
						output << ",";
						output << f_info[i].refs;
						output << ",0,0,0";
					}
					output << "\n";
					(arguments.outFunc)(arguments, output.str(), output_file);
				}
				f_info.clear();
			}
			if (tid > vcf_tid || (tid == vcf_tid && pos > line->pos)) {
				// we missed it, go to the next one
				//printf("i'm at %d, vcf is at %d\n", tid, vcf_tid);
				break;
			}
		}
	}

	double duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("Finished in %f seconds.\n", duration);

	if (arguments.gzipped) {
		bgzf_close(arguments.gzippedPointer);
	} else {
		fclose(output_file);
	}

	return 0; // added return code for successful completion 
}

int main(int argc, char ** argv) {
	struct arguments arguments;

	arguments.args = vector<char *>();
	arguments.count_orphans = false;
	arguments.gzipped = false;
	arguments.ignore_overlaps = false;
	arguments.min_base_quality = 0;
	arguments.min_map_quality = 0;
	arguments.min_read_counts = vector<int>();
	arguments.max_depth = 4000;
	arguments.outFunc = print_output;
	arguments.progress = false;
	arguments.pseudo_snps = 0;
	arguments.verbose = false;

	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	return program_main(arguments);
}
