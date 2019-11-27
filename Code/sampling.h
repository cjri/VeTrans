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
void GenerateSample (run_params p, double C, vector<int>& Npre, vector<int>& Npost, vector<double> freq_pre, vector<double> freq_post, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets, vector< vector< vector<mhap> > >& all_hap_data_sets, gsl_rng *rgen);
vector<int> DirMultSampling(int N, vector<double> &freqs, double C, const gsl_rng *r);
vector<int> multinomialSampling(int N, vector<double> p, const gsl_rng *r);
void FlattenData (int i, vector< vector< vector<mhap> > >& all_hap_data_sets, vector<mhap>& hap_data);
void CalculateVarianceMatrices(vector< vector<double> >& all_pre_freqs, vector< vector<double> >& all_post_freqs, vector<double>& pre_var, vector<double>& post_var);





