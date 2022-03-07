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
void GetOptions (run_params& p, int argc, const char **argv);
void GetVariantData (run_params p, vector<mhap>& haps);
void GetVariantDataMulti (run_params p, vector<mhap>& haps);
void GetFullHaplotypes (run_params p, vector<haplo>& full_haps);
void GetFullHaplotypesFreq (run_params p, vector<haplo>& full_haps, vector<double>& freq_pre, vector<double>& freq_post);
void PrintResampleMulti (vector< vector< vector<mhap> > >& all_hap_data_sets);
void PrintSample (int i, vector<mhap>& hap_data);
void PrintVariances (run_params p, vector<double>& pre_var, vector<double>& post_var);
void GetLikelihoods (run_params p, vector< vector<double> >& likelihoods);
void GetVarianceMatrices(run_params p, int dim, vector< vector<double> >& var_pre_provis, vector< vector<double> >& var_post_provis);
void PrintReducedData(run_params p, int dim, vector<double> freq_pre, vector<double> freq_post, gsl_matrix *var_pre, gsl_matrix *var_post);
void OutputReconstructHaplotypes (vector<haplo> haplotypes,vector<double> pre_freqs,vector<double> post_freqs,double best_bic);
void OutputReconstructHaplotypesMulti (run_params p, vector<haplo> haplotypes, vector< vector<double> >& multi_freqs, double best_bic);

