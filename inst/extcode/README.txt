snp-pileup
----------

This application will, given a VCF file containing SNP locations, output for
each SNP the counts of the reference nucleotide, alternative nucleotide,
errors, and deletions. These counts can then be used in facets.

Was developed on a linux machine and tested with htslib v1.3.1

Installation
------------
First, HTSlib must be installed on your system. To do that, download it from
http://www.htslib.org/download/ and follow the "Building and installing" 
instructions on that page. If installed systemwide (in /usr/local/lib) using
"make install" ensure that the libraries are available with the command 
"sudo ldconfig" (only needs to be run once).

This code can be compiled using 

     g++ -std=c++11 snp-pileup.cpp -lhts -o snp-pileup

when htslib is available systemwide, or

     g++ -std=c++11 -I/path/htslib/include snp-pileup.cpp \\
     -L/path/htslib/lib -lhts -Wl,-rpath=/path/htslib/lib -o snp-pileup 

when it is installed locally and path is the location where it is available.


Usage
-----
     snp-pileup <vcf file> <output file> <sequence files...>

Usage of snp-pileup requires a VCF file and one (or multiple) sequence files
containing DNA. The sequence files should be in the BAM format, and both the
VCF and all sequence files must be sorted. A suitable option for VCF is one
from NCBI. (WARNING: the link below points to the VCF for the current version
of genome build)

 ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz

Use a VCF file that is consistent with the genome build used for aligning the
sequencing data. [Available snp/genome build versions are in directories
human_9606_b149_GRCh37p13, human_9606_b150_GRCh37p13, human_9606_b149_GRCh38p7,
human_9606_b150_GRCh38p7.]


Parameters

Here is a list of parameters snp-pileup accepts and information about what
they do. Some of them, such as -q, -Q, -A, and -x, are the same as their
equivalent in samtools mpileup, and are used the same way.

  -A, --count-orphans        Do not discard anomalous read pairs.
  -d, --max-depth=DEPTH      Sets the maximum depth. Default is 4000.
  -g, --gzip                 Compresses the output file with BGZF.
  -p, --progress             Show a progress bar. WARNING: requires additional
                             time to calculate number of SNPs, and will take
                             longer than normal.
  -P, --pseudo-snps=MULTIPLE Every MULTIPLE positions, if there is no SNP,
                             insert a blank record with the total count at the
                             position.
  -q, --min-map-quality=QUALITY   Sets the minimum threshold for mapping
                             quality. Default is 0.
  -Q, --min-base-quality=QUALITY   Sets the minimum threshold for base quality.
                             Default is 0.
  -r, --min-read-counts=READS   Comma separated list of minimum read counts for
                             a position to be output. Default is 0.
  -v, --verbose              Show detailed messages.
  -x, --ignore-overlaps      Disable read-pair overlap detection.


Unlike mpileup which splits the max-depth (-d option) value specified evenly
among the bam files, it is applied to each file individually.

The pseudo-snps (-P option) specification gives the depths for non-polymorphic
loci along a grid of positions. These are useful for estimating total copy
number in regions that are sparse in snps.

The min-read-counts (-r option) is used to output only loci that are captured
at sufficient depth. This helps to limit the output file size by eliminating
off-target reads. It is used to specify lower bound for the normal sample. For
instance if the normal is sequenced at 50x depth we use -r20,0 where 20 is the
minimum depth for the normal bam and 0 for the tumor bam.

This program outputs a csv file and -g option produces a gzip compressed one.

You can view this list at any time by using --help.

Limitations
-----------
SNPs where there are multiple nucleotides changing will be ignored, and all
minimum thresholds (except for the minimum read count) apply equally to all
files—there is no way to set them on a per-file basis.
