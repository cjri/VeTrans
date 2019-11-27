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

//utilities
void GetOptions (run_params& p, int argc, const char **argv);
void CalculateFullLocusSet (vector<mhap>& hap_data, vector<int>& loci);
void CalculateNLociData (vector<mhap>& hap_data, vector<int>& loci);
void CalculateSetData (int& n_sets, vector<mhap>& hap_data);
void ConstructHaplotypeSet (int n_haps, vector<int>& loci, vector<mhap>& hap_data, vector<haplo>& full_haps, gsl_rng *rgen);
void ConstructHaplotypeDeNovo (haplo& hap, vector<int>& loci, vector<mhap>& hap_data, gsl_rng *rgen);
void ConstructHaplotypeX (haplo& hap, vector<int>& loci);
void CheckUniqueHaplotypes (int& check, vector<mhap>& hap_data, vector<haplo>& full_haps, gsl_rng *rgen);
int ChangeHaplotypeRandom (run_params p, vector<haplo>& full_haps, vector<haplo>& haplotypes, vector< vector<mhap> >& hap_data_sets, vector< vector<haplo> >& full_haps_archive, gsl_rng *rgen);
void CheckUniqueHaplotypesSplit (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, gsl_rng *rgen);
int CheckHaplotypesArchive (vector<haplo>& full_haps, vector< vector<haplo> >& full_haps_archive);
void SplitHaplotypeData (int n_sets, vector<mhap>& hap_data, vector< vector<mhap> >& hap_data_sets);
void CheckXHap (vector<int>& loci, vector<haplo>& full_haps);
void MatchFullHapsSplit (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps);
void CalculateObsTotalN (int n_sets, vector<int>& Npre, vector<int>& Npost, vector<mhap>& hap_data);
void ConstructObservationsX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector<int> >& obs_freqs_pre, vector< vector<int> >& obs_freqs_post);
void FindLogFact(vector<double>& fact_store, vector<int>& Npre, vector<int>& Npost);
void SetInitialFreqs (vector<double>& hap_freqs, vector<haplo> full_haps, gsl_rng *rgen);
void NormaliseFreqs (vector<double>& hap_freqs);
void CheckNoiseFreq (vector<double>& hap_freqs);
void ChangeFreq (double changex, vector<haplo>& full_haps, vector<double>& hap_freqs, gsl_rng *rgen);
void OptimiseFrequencies (run_params p, double& lL_pre, double& lL_post, double C, vector<int>& Npre, vector<int>& Npost, vector<haplo>& full_haps, vector<double>& hap_freqs_pre, vector<double>& hap_freqs_post, vector< vector<int> >& obs_freqs_pre, vector< vector<int> >& obs_freqs_post, vector< vector<mhap> >& hap_data_sets, vector<double>& fact_store, gsl_rng *rgen);
void MissingFreqRedist(vector<double>& inf_freqs);
void FindBestRep (double& lL_tot, vector<double>& hap_freqs_pre, vector<double>& hap_freqs_post, vector<double>& logs_store, vector< vector<double> >& pre_store, vector< vector<double> >& post_store);
void OptimiseFrequenciesFixedFullHaps (run_params p, double C, vector<double>& pre_freqs, vector<double>& post_freqs, vector<haplo>& full_haps, vector<int>& Npre, vector<int>& Npost, vector<double>& fact_store, vector< vector<mhap> >& hap_data_sets, gsl_rng *rgen);
int CheckTime (run_params p, time_t timer_s);
double DirichletMultiCalc (run_params p, int N, double c, vector<int>& obs, vector<double>& inf, vector<double>& fact_store);
double FindBIC(double L, vector<haplo>& full_haps);
void PrintInference (vector<haplo>& full_haps, vector<double>& hap_freqs_pre, vector<double>& hap_freqs_post);
void FinalOutput (run_params p, double lL_hapbest, vector<haplo>& haplotypes, vector<double>& hap_freqs_pre_best, vector<double>& hap_freqs_post_best, vector< vector<haplo> >& full_haps_store);
int ReconstructHaplotypesFixedFromInternal (run_params& p, double& best_bic, vector<haplo>& haplotypes, vector<double>& pre_freqs, vector<double>& post_freqs, vector<mhap>& hap_data, vector<haplo>& full_haps);

