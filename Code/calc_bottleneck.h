//Bottleneck calculation
int CalcBottleneck(run_params p);
void RemoveNotSeenPreTransmission (run_params p, int& dim, vector<double>& freq_pre, vector<double>& freq_post, vector< vector<double> >& var_pre_provis, vector< vector<double> >& var_post_provis, vector<haplo>& full_haps);
void DefineVariances(int dim, vector< vector<double> > var_pre_provis, vector< vector<double> > var_post_provis, gsl_matrix* var_pre, gsl_matrix* var_post);
void NormaliseFreqs(vector<double> &freqs);
vector<vector<int> > findSwaps(std::vector<int>& extinctions, int dimensions);
void swapMatrix(gsl_matrix *mat, vector<vector<int> > swaps);
void swapVector(vector<double>& v, vector<vector<int> >& swaps);



