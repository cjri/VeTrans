#include "shared_haps.h"
#include "explicit.h"
#include "compound.h"
#include "calc_bottleneck.h"
#include "reconstruct_haplotypes.h"
#include "io.h"

int CalcBottleneck(run_params p) {
	double C=p.c;
	
	//Set up for explicit likelihood code
	vector< vector<mhap> > hap_data_sets;
	vector<int> Npre;
	vector<int> Npost;
	if (p.exp_like==1) {
        if (p.verb==1) {
            cout << "Explicit method\n";
        }
		SetupExplicit1(p,Npre,Npost,hap_data_sets);
	}
	
	//Read in haplotype frequency data
	if (p.verb==1) {
		cout << "Reading in frequencies from " << p.fullhap_in << "\n";
	}
	vector<haplo> full_haps; //Structure for full haplotypes
	vector<double> freq_pre;
	vector<double> freq_post;
	GetFullHaplotypesFreq(p,full_haps,freq_pre,freq_post);
	if (full_haps.size()==0) {
		cout << "Error: No data in haplotype frequencies\n";
		return 0;
	}
	int dim=freq_pre.size();
	
	//Read in provisional variance matrices
	vector< vector<double> > var_pre_provis;
	vector< vector<double> > var_post_provis;
	GetVarianceMatrices(p,dim,var_pre_provis,var_post_provis);
	
	//Remove haplotypes that are only seen post-transmission.  These are assumed to result from de novo mutation
	RemoveNotSeenPreTransmission (p,dim,freq_pre,freq_post,var_pre_provis,var_post_provis,full_haps);
	
	//Setup for explicit likelihood calculation.
	vector< vector<int> > obs_freqs_post;
	vector<double> fact_store;
	if (p.exp_like==1) {
		SetupExplicit2 (p,full_haps,Npost,hap_data_sets,obs_freqs_post,fact_store);
	}

	gsl_matrix *var_pre = gsl_matrix_alloc(dim, dim);
	gsl_matrix *var_post = gsl_matrix_alloc(dim, dim);
	DefineVariances(dim,var_pre_provis,var_post_provis,var_pre,var_post);

	//Normalise frequencies to add to 1
	NormaliseFreqs(freq_pre);
	NormaliseFreqs(freq_post);

	if (p.verb==1) {
		PrintReducedData(p,dim,freq_pre,freq_post,var_pre,var_post);
	}
	
	//Calculate the maximum likelihood
	double maxL=-1e80;
	int maxN=-1;
	int chk=0;
	if (p.print_like==1) {
		ofstream like_file;
		if (p.exp_like==1) {
			like_file.open("Likelihoods_exp.out");
		} else {
			like_file.open("Likelihoods.out");
		}
		for (int Nt=1;Nt<=p.max_n;Nt++) {
            if (p.verb==1) {
                cout << "Nt " << Nt;
            }
			if (p.exp_like==1) {
				vector< vector<int> > list;
				if (p.max_n==21) {
					cout << "Warning: Explicit code may take a long time for Nt>20.";
				}
				//Get explicit likelihood
				double lL=CalculateExplicitLikelihood (p,Nt,C,Npost,freq_pre,freq_post,fact_store,list,obs_freqs_post,full_haps,hap_data_sets);
				//cout << "Nt " << Nt << " " << log(lL) << "\n";
				if (isinf(log(lL))) {
					like_file << Nt << " " << -1e80 << "\n";
                    if (p.verb==1) {
                        cout << " " << -1e80 << "\n";
                    }

				} else {
					like_file << Nt << " " << log(lL) << "\n";
                    if (p.verb==1) {
                        cout << " " << log(lL) << "\n";
                    }

				}
			} else {
				if (full_haps.size()<=1) {
					if (p.verb==1&&chk==0) {
						chk=1;
						cout << "Error: Require at least two haplotypes for compound likelihood to work\n";
					}
					like_file << Nt << " " << 0 << "\n";
                    if (p.verb==1) {
                        cout << " " << 0 << "\n";
                    }

				} else {
					//Get approximate likelihood
					double L=FindCompoundLikelihood (p,Nt,dim,freq_pre,freq_post,var_pre,var_post);
					like_file << Nt << " " << L << "\n";
                    if (p.verb==1) {
                        cout << " " << L << "\n";
                    }

				}
			}

		}
	}
	
	if (p.find_max==1) {
		vector< vector<int> > list;
		//Find maximum likelihood
		FindMaxLikelihood (p,dim,maxN,maxL,C,Npost,freq_pre,freq_post,fact_store,list,obs_freqs_post,full_haps,hap_data_sets,var_pre,var_post);
		
		if (p.exp_like) {
		cout << "Maximum likelihood bottleneck size was " << maxN << " with log likelihood of " << log(maxL) << "\n";
		} else {
			cout << "Maximum likelihood bottleneck size was " << maxN << " with log likelihood of " << maxL << "\n";
		}
	}
	return 0;
}

void RemoveNotSeenPreTransmission (run_params p, int& dim, vector<double>& freq_pre, vector<double>& freq_post, vector< vector<double> >& var_pre_provis, vector< vector<double> >& var_post_provis, vector<haplo>& full_haps) {
	vector< vector<double> > var_pre_red;
	vector< vector<double> > var_post_red;
	p.verb=0;
	int ns=0;
	vector<int> never_seen;
	for (int i=dim-1;i>=0;i--) {
		if (freq_pre[i]<p.extinct) {
			never_seen.push_back(i);
			ns++;
		}
	}
	if (p.verb==1) {
		cout << "Number of haplotypes not seen before transmission = " << ns << "\n";
	}
	if (ns==0) {
		var_pre_red=var_pre_provis;
		var_post_red=var_post_provis;
	} else {
		//Remove unseen frequencies from haplotype frequency vectors and haplotype list
		for (int i=0;i<never_seen.size();i++) {
			freq_pre.erase(freq_pre.begin()+never_seen[i]);
			freq_post.erase(freq_post.begin()+never_seen[i]);
			full_haps.erase(full_haps.begin()+never_seen[i]);
		}
		//Sort out variance matrices
		for (int i=0;i<dim;i++) {
			vector<double> vpB;
			vector<double> vpA;
			int inc=1;
			for (int k=0;k<never_seen.size();k++) {
				if (never_seen[k]==i) {
					inc=0;
				}
			}
			if (inc==1) {
				for (int j=0;j<dim;j++) {
					inc=1;
					for (int k=0;k<never_seen.size();k++) {
						if (never_seen[k]==j) {
							inc=0;
						}
					}
					if (inc==1) {
						vpB.push_back(var_pre_provis[i][j]);
						vpA.push_back(var_post_provis[i][j]);
					}
				}
				var_pre_red.push_back(vpB);
				var_post_red.push_back(vpA);
			}
		}
		
		if (p.verb==1) {
			cout << "Reduced variance : Pre transmission:\n";
			for(int i=0;i<var_pre_red.size();i++) {
				for(int j=0;j<var_pre_red[0].size();j++) {
					cout << var_pre_red[i][j] << " ";
				}
				cout << "\n";
			}
			cout << "Reduced variance : Post transmission:\n";
			for(int i=0;i<var_post_red.size();i++) {
				for(int j=0;j<var_post_red[0].size();j++) {
					cout << var_post_red[i][j] << " ";
				}
				cout << "\n";
			}
		}
		
		dim=freq_pre.size();
		if (p.verb==1) {
			cout << "Reduced dimension to " << dim << "\n";
		}
	}
	var_pre_provis=var_pre_red;
	var_post_provis=var_post_red;
}

void DefineVariances(int dim, vector< vector<double> > var_pre_provis, vector< vector<double> > var_post_provis, gsl_matrix* var_pre, gsl_matrix* var_post) {
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) {
			gsl_matrix_set(var_pre,i,j,var_pre_provis[i][j]);
			gsl_matrix_set(var_post,i,j,var_post_provis[i][j]);
		}
	}
}

vector<vector<int> > findSwaps(std::vector<int>& extinctions, int dimensions) {
	vector<vector<int> > swaps;
	//Go through extinctions one at a time (largest first), then swap with furhest to the right non-extinction
	//Example: All sites = {0,1,2,3,4,5,6}, extinction sites = {5,3,1}.
	//Then first swap 5 with 6, so get order {0,1,2,3,4,6,5} and swap pair {5,6}
	//Then swap 3 with 5, i.e. get {0,1,2,6,4,3,5} and swap pairs {5,6}, {3,5}
	//Finally swap 1 with 4, i.e. get {0,4,2,6,1,3,5} and swap pairs {5,6}, {3,5}, {1,4}.
	//We could further swap {0,4,2,6} to get them in order, but it doesn't matter much/is time consuming to do additional swaps
	for(unsigned int i=0; i<extinctions.size(); i++) { //Loop over extinction entries
		
		//Check if extinction is too far to the left, if so, move it to right position.
		//E.g. if extinction[i] = 5 and dimensions = 7 and (i+1) = 1, then 5 < 7-1 is true, so make swap with site 7-1=6
		//Similarly, if extinction[i] = 3 and dimensions = 7 and (i+1) = 2, then 3 < 7-2 is true, so make swap with position 7-2=5
		if(extinctions[i] < (int)(dimensions-(i+1))) {
			
			int swapPartner = dimensions - (i+1);
			vector<int> swap;
			swap.push_back(extinctions[i]);
			swap.push_back(swapPartner);
			swaps.push_back(swap);
		}
	}
	
	return swaps;
	
}


//Perform a set of swaps (e.g. swaps = {{5,6}, {3,5}, {1,4}}) for a matrix mat using gsl methods
void swapMatrix(gsl_matrix *mat, vector<vector<int> > swaps) {
	for(unsigned int i=0; i<swaps.size(); i++) {
		int swap1 = swaps[i][0]; //First entry for swapping, e.g. 5
		int swap2 = swaps[i][1]; //Second enttry for swappging, e.g. 6
		gsl_matrix_swap_rows(mat, swap1, swap2);
		gsl_matrix_swap_columns(mat, swap1, swap2);
	}
}


//Perform a set of swaps (e.g. swaps = {{5,6}, {3,5}, {1,4}}) for a double vector (e.g. {a,b,c,d,e,f,g} becomes {a,e,c,g,b,d,f}
void swapVector(vector<double>& v, vector<vector<int> >& swaps) {
	//	cout << "Swapping vector: "; printVector(v);
	//	cout << "Swap vector is: ";
	//	for(unsigned int i=0; i<swaps.size(); i++) {
	//
	//		for(unsigned int j=0; j<swaps[i].size(); j++) {
	//			cout << swaps[i][j] << " ";
	//		}
	//		cout << "\n";
	//	}
	for(unsigned int i=0; i<swaps.size(); i++) {
		int swap1 = swaps[i][0]; //First entry for swapping, e.g. 5
		int swap2 = swaps[i][1]; //Second enttry for swappging, e.g. 6
		iter_swap(v.begin() + swap1, v.begin() + swap2); //Defined in algorithm
		//		cout << "Vector after swap " << i << ": "; printVector(v);
	}
}
