//Compound likelihood functions
void constructMatrixM(vector<double> x, gsl_matrix *M);
double computeLikelihoodCont(vector<double> x, vector<double> mean, gsl_matrix *var);
double logMultivariateNormalPDF(vector<double> &x, vector<double> &mu, gsl_matrix *sigma);
vector<double> subtractVectors(vector<double> &a, vector<double> &b);
double FindCompoundLikelihood (run_params p, int Nt, int dim, vector<double>& freq_pre, vector<double>& freq_post, gsl_matrix *var_pre, gsl_matrix *var_post);
void FindCombinedLikelihood (vector< vector<double> >& likelihoods);
void FindMaxLikelihood (run_params p, int dim, int& maxN, double& maxL,double C, vector<int>& Npost, vector<double>& freq_pre, vector<double>& freq_post, vector<double>& fact_store, vector< vector<int> >& list, vector< vector<int> >& obs_freqs_post, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets, gsl_matrix *var_pre, gsl_matrix *var_post);






