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

void CalculateObsTotalNMulti (run_params& p, int n_sets, vector< vector<int> >& Nmulti, vector<mhap>& hap_data);
void ConstructObservationsXMulti (run_params& p, vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector< vector<int> > >& obs_freqs_multi);
void FindLogFactMulti(vector<double>& fact_store, vector< vector<int> >& Nmulti);
void OptimiseFrequenciesMulti (run_params p, vector<double>& lL, double C, vector< vector<int> >& Nmulti, vector<haplo>& full_haps, vector< vector<double> >& hap_freqs_multi, vector< vector< vector<int> > >& obs_freqs_multi, vector< vector<mhap> >& hap_data_sets, vector<double>& fact_store, gsl_rng *rgen);
void FindBestRepMulti (run_params& p, double& lL_tot, vector< vector<double> >& hap_freqs_multi, vector<double>& logs_store, vector< vector< vector<double> > >& multi_store);
void OptimiseFrequenciesFixedFullHapsMulti (run_params p, double C, double& best_bic, vector< vector<double> >& multi_freqs, vector<haplo>& full_haps, vector< vector<int> >& Nmulti, vector<double>& fact_store, vector< vector<mhap> >& hap_data_sets, gsl_rng *rgen);
void PrintInferenceMulti (run_params& p, vector<haplo>& full_haps, vector< vector<double> >& hap_freqs_multi);
void FinalOutputMulti (run_params p, double lL_hapbest, vector<haplo>& haplotypes, vector< vector<double> >& hap_freqs_best_multi, vector< vector<haplo> >& full_haps_store);

