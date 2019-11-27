#include "shared_haps.h"
#include "sampling.h"
#include "sample_haplotypes.h"
#include "reconstruct_haplotypes.h"
#include "io.h"

//rec data;
int SampleData (run_params p, vector< vector< vector<mhap> > >& all_hap_data_sets) {

	double C=p.c;
	
	//Initialise random number generator
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, p.seed);
	
	
	//Read in multi-locus trajectories file
	vector<mhap> hap_data;
	GetVariantData(p,hap_data);
	//Calculate full set of loci
	vector<int> loci;
	CalculateFullLocusSet(hap_data,loci);
	sort(loci.begin(),loci.end());

	//Update data with n_loci values: As 0, 1, 2, 3, rather than positions in genome
	CalculateNLociData(hap_data,loci);
	
	//Update data with set i.e. which reads cover which loci
	int n_sets;
	CalculateSetData(n_sets,hap_data);
	
	//Calculate N values - number of samples of each multi-locus type
	vector<int> Npre;
	vector<int> Npost;
	CalculateObsTotalN (n_sets,Npre,Npost,hap_data);

	//Split partial haplotype data according to which loci are spanned by each read
	vector< vector<mhap> > hap_data_sets;
	SplitHaplotypeData (n_sets,hap_data,hap_data_sets);
	
	//Next step: Read in inferred haplotypes and their frequencies
	vector<haplo> full_haps; //Structure for full haplotypes
	vector<double> freq_pre;
	vector<double> freq_post;
	GetFullHaplotypesFreq(p,full_haps,freq_pre,freq_post);
	//Remove XXXX haplotype
	freq_pre.pop_back();
	freq_post.pop_back();
	full_haps.pop_back();

	//Generate a resampled population given the read structure full haplotype reconstruction.
	//Output as a new Multi_locus_trajectories type file
	for (int i=0;i<p.sample_reps;i++) {
		GenerateSample (p,C,Npre,Npost,freq_pre,freq_post,full_haps,hap_data_sets,all_hap_data_sets,rgen);
	}

	return 0;

	
}
