//Shared information for linked optimisation
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

using namespace std;

struct run_params {
	int seed; //Random seed.  Get from command line
	double c; //Overdispersion parameter
	int n_haps; //Number of haplotypes to include in model.  Default is 1.
    int n_samples; //Number of samples in multi-locus trajectory file.
	int hap_its; //Number of sets of haplotypes to consider in optimisation.  Default is 1000.
	int hap_term; //Number of consecutive changes to a haplotype not giving a previously observed set of haplotypes before terminating
	int freq_its; //Number of iterations to calculate over in frequency optimisation. Default is 10000.
	int freq_term; //Number of consecutive iterations after which to terminate frequency optimisation given no improvements found.  Default is 1000.
	int freq_rep; //Number of independent frequency optimisations to run for each set of full haplotypes.  Default is 1.
	int run_time; //Allows the code to be run for a specific length of time.  Default is 1000000 seconds.  Setting this fixes hap_its to be large.
	int read_hap; //Read in a fixed set of haplotypes rather than finding an optimal set
    int read_prev; //Read in a previous set of haplotypes as part of a starting point for a new optimisation
	int print_sample; //Flag for resampling code. Outputs resample data
	int sample_reps; //Number of samples to conduct
	int suppress_output; //Anti-verbose output
	int exp_like; //Use explicit likelihood function
	int growth; //Within-host viral growth rate for compound likelihood method
	int max_n; //Maximum bottleneck size considered by inference method
	double extinct; //Cutoff frequency for regarding a haplotype to be extinct
	int find_max; //Find maximum likelihood Nt in bottleneck inference
	int print_like; //Print likelihoods for multiple Nt into file
	double null_freq; //Used in likelihood to replace zeros
	double delta_bic; //Change in BIC required to infer an additional haplotype during reconstruction
	const char* err_flag; //Name of input file containing multi-locus polymorphisms
	const char* traj_in; //Name of input file containing multi-locus polymorphisms
	const char* traj_out; //Name of output file to contain multi-locus polymorphisms
	const char* var_name; //Prefix for variant data
	const char* like_name; //Prefix for variant data
	const char* fullhap_in; //Name of input file containing haplotypes
    const char* prev_in; //Name of input file containing previous haplotypes
	int verb; //Verbose output
	int err; //Error flag
};


struct mhap {
	int nvar;
	string st;
	vector<char> seq;
	vector<int> loci;
	vector<int> n_loci;
	vector<int> times;
	vector<int> obs;
	vector<int> match;
	int tot; //Number of observations of this haplotype
	int set; //Index by loci covered
};

struct haplo {
	string st;
	vector<char> seq;
};

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

void GetOptions (run_params& p, int argc, const char **argv);

int ReconstructHaplotypes (run_params p, double& best_bic, vector<haplo>& haplotypes, vector<double>& pre_freqs, vector<double>& post_freqs);
int ReconstructHaplotypesMulti (run_params p, double& best_bic, vector<haplo>& haplotypes, vector< vector<double> >& multi_freqs);


//get_files
void GetVariantData (run_params p, vector<mhap>& haps);
void GetFullHaplotypes (run_params p, vector<haplo>& full_haps);
