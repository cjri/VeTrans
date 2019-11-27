//Explicit likelihood code
void SetupExplicit1 (run_params& p, vector<int>& Npre, vector<int>& Npost, vector< vector<mhap> >& hap_data_sets);
void SetupExplicit2 (run_params p, vector<haplo> full_haps, vector<int>& Npost, vector< vector<mhap> >& hap_data_sets, vector< vector<int> >& obs_freqs_post, vector<double>& fact_store);
void MatchFullHapsSplitX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps);
void MatchFullHapsSplitNoX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps);
void ConstructObservationsX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector<int> >& obs_freqs_post);
void ConstructObservationsNoX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector<int> >& obs_freqs_post);
void FindLogFact(vector<int> Npost, vector<double>& fact_store);
double CalculateExplicitLikelihood (run_params p, int Nt, double C, vector<int>& Npost, vector<double>& freq_pre, vector<double>& freq_post, vector<double>& fact_store, vector< vector<int> >& list, vector< vector<int> >& obs_freqs_post, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets);
void GetMultinomialOutcomes (int Nt, vector<double>& freq_pre, vector< vector<int> >& list);
int cprint (int n, int r, int pos, int sum, vector<int>& choose, vector< vector<int> >& list);
void FindMnomProbs (int Nt, vector< vector<int> >& list, vector<double>& freq_pre, vector<double>& mnom_probs, vector<double>& fact_store);
void GetInferredHaplotypeFrequencies(int i, int Nt, vector<double> freq_post, vector<double>& freq_post_prob, vector< vector<int> >& list);
void GetScaledPartialHaplotypeLikelihood (run_params p, double C, double& L, vector<int>& Npost, vector< vector<int> >& obs_freqs_post, vector<double>& freq_post_prob, vector<double>& inf_freqs_post, vector<double>& fact_store, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets);

