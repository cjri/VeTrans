#include "shared_haps.h"
#include "reconstruct_haplotypes.h"
#include "explicit.h"
#include <iostream>
#include <string>
#include <sstream>

//Code to run the explicit likelihood
void SetupExplicit1 (run_params& p, vector<int>& Npre, vector<int>& Npost, vector< vector<mhap> >& hap_data_sets) {
    if (p.verb==1) {
        cout << "Set up for explicit calcluation\n";
    }
	p.null_freq=p.extinct/10;
	//Need the multi-locus data to run the explicit likelihood calculation
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
	CalculateObsTotalN (n_sets,Npre,Npost,hap_data);
	//Split partial haplotype data according to which loci are spanned by each read
	SplitHaplotypeData (n_sets,hap_data,hap_data_sets);
}

void SetupExplicit2 (run_params p, vector<haplo> full_haps, vector<int>& Npost, vector< vector<mhap> >& hap_data_sets, vector< vector<int> >& obs_freqs_post, vector<double>& fact_store) {
	if (full_haps[full_haps.size()-1].st.find("X")==0) {
	//	cout << "Match X\n";
		MatchFullHapsSplitX(hap_data_sets,full_haps);
		ConstructObservationsX (hap_data_sets,full_haps,obs_freqs_post);
	} else {
	//	cout << "Match NoX \n";
		MatchFullHapsSplitNoX(hap_data_sets,full_haps);
		ConstructObservationsNoX (hap_data_sets,full_haps,obs_freqs_post);
	}
	FindLogFact(Npost,fact_store);
	if (p.verb==1) {
		cout << "Data with match information\n";
		for (int i=0;i<hap_data_sets.size();i++) {
			for (int j=0;j<hap_data_sets[i].size();j++) {
				cout << hap_data_sets[i][j].n_loci.size() << " ";
				for (int k=0;k<hap_data_sets[i][j].n_loci.size();k++) {
					cout << hap_data_sets[i][j].n_loci[k] << " ";
				}
				cout << hap_data_sets[i][j].st << " ";
				cout << hap_data_sets[i][j].times.size() << " ";
				for (int k=0;k<hap_data_sets[i][j].times.size();k++) {
					cout << hap_data_sets[i][j].times[k] << " " << hap_data_sets[i][j].obs[k] << " ";
				}
				for (int k=0;k<hap_data_sets[i][j].match.size();k++) {
					cout << hap_data_sets[i][j].match[k] << " ";
				}
				cout << "\n";

			}
		}
	}
}

void MatchFullHapsSplitX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps) {
	//Identify which partial haplotype observations come from which full haplotypes
	for (int i=0;i<hap_data_sets.size();i++) {
		for (int j=0;j<hap_data_sets[i].size();j++) {
			//Hap_data j
			hap_data_sets[i][j].match.clear();
			for (int k=0;k<full_haps.size()-1;k++) {
				//Full haplotype k
				int match=1;
				for (int l=0;l<hap_data_sets[i][j].n_loci.size();l++) {
					if (hap_data_sets[i][j].st[l]!=full_haps[k].seq[hap_data_sets[i][j].n_loci[l]]) {
						match=0;
					}
				}
				if (match==1) {
					hap_data_sets[i][j].match.push_back(k);
				}
				
			}
			if (hap_data_sets[i][j].match.size()==0) { //No matches: Match to X
				hap_data_sets[i][j].match.push_back(full_haps.size()-1);
			}
		}
	}
}

void MatchFullHapsSplitNoX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps) {
	for (int i=0;i<hap_data_sets.size();i++) {
		for (int j=0;j<hap_data_sets[i].size();j++) {
			//Hap_data j
			hap_data_sets[i][j].match.clear();
			for (int k=0;k<full_haps.size();k++) {
				//Full haplotype k
				int match=1;
				for (int l=0;l<hap_data_sets[i][j].n_loci.size();l++) {
					if (hap_data_sets[i][j].st[l]!=full_haps[k].seq[hap_data_sets[i][j].n_loci[l]]) {
						match=0;
					}
				}
				if (match==1) {
					hap_data_sets[i][j].match.push_back(k);
				}
				
			}
		}
	}
}

void ConstructObservationsX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector<int> >& obs_freqs_post) {
	//Calculate sets of numbers of observed partial haplotypes
	for (int i=0;i<hap_data_sets.size();i++) {
		vector<int> opp;
		int oxp=0;
//		cout << "i= " << i << " " << hap_data_sets[i].size() << "\n";
		for (int j=0;j<hap_data_sets[i].size();j++) {  //Construct set of all matching observations + sum of all non-observed observations
			if (hap_data_sets[i][j].match[0]<full_haps.size()-1) {
				opp.push_back(hap_data_sets[i][j].obs[1]);
			} else {
				oxp=oxp+hap_data_sets[i][j].obs[1];
			}
		}
		if (opp.size()==0) {
			opp.push_back(0);
		}
//		cout << oxp << "\n";
		opp.push_back(oxp);
		obs_freqs_post.push_back(opp);
	/*	for (int j=0;j<opp.size();j++) {
			cout << opp[j] << " ";
		}
		cout << "\n";*/
	}
}


void ConstructObservationsNoX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector<int> >& obs_freqs_post) {
	for (int i=0;i<hap_data_sets.size();i++) {
		vector<int> opp;
		for (int j=0;j<hap_data_sets[i].size();j++) {  //Construct set of all matching observations
			if (hap_data_sets[i][j].match.size()>0) {
				if (hap_data_sets[i][j].match[0]<full_haps.size()) {
					opp.push_back(hap_data_sets[i][j].obs[1]);
				}
			}
		}
		if (opp.size()==0) {
			opp.push_back(0);
		}
		obs_freqs_post.push_back(opp);
	}
}

void FindLogFact(vector<int> Npost, vector<double>& fact_store){
	int max=0;
	for (int i=0;i<Npost.size();i++) {
		if (Npost[i]>max) {
			max=Npost[i];
		}
	}
	int Nt=max+10;
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=Nt+1;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
	}
}

double CalculateExplicitLikelihood (run_params p, int Nt, double C, vector<int>& Npost, vector<double>& freq_pre, vector<double>& freq_post, vector<double>& fact_store, vector< vector<int> >& list, vector< vector<int> >& obs_freqs_post, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets) {
	GetMultinomialOutcomes (Nt,freq_pre,list);
	//Probability of each outcome in the list
	vector<double> mnom_probs;
	FindMnomProbs(Nt,list,freq_pre,mnom_probs,fact_store);
	//Evaluate likeihood
	double lL=0;
	vector<double> inf_freqs_post;
	vector<double> freq_post_prob;
	for (int i=0;i<mnom_probs.size();i++) { //Set of multinomial outcomes
		//cout << i << " " << mnom_probs.size() << "\n";
		GetInferredHaplotypeFrequencies(i,Nt,freq_post,freq_post_prob,list); //From multinomial output
		//Get the multinomial outcome here.
		double L=0;
		//Calculate likelihood for this multinomial outcome
		GetScaledPartialHaplotypeLikelihood (p,C,L,Npost,obs_freqs_post,freq_post_prob,inf_freqs_post,fact_store,full_haps, hap_data_sets);
		//cout << "i " << i << " " << mnom_probs[i] << " L " << L << " e^L " << exp(L) << "\n";
		L=exp(L)*mnom_probs[i];
		lL=lL+L;
	}
	return lL;
}

void GetMultinomialOutcomes (int Nt, vector<double>& freq_pre, vector< vector<int> >& list) {
//	cout << "GetMultinomial\n";
	list.clear();
	vector<int> choose;
	for (int i=0;i<freq_pre.size();i++) {
		choose.push_back(0);
	}
	cprint(Nt,freq_pre.size(),0,0,choose,list);
}

//Generate list of multinomial possibilities
int cprint (int n, int r, int pos, int sum, vector<int>& choose, vector< vector<int> >& list) {
	if (pos==r-1) {
		choose[choose.size()-1]=n-sum;
		list.push_back(choose);
		return 0;
	} else {
		int csum=sum;
		for (int p=0;p<=n-csum;p++) {
			choose[pos]=p;
			sum=csum+p;
			cprint(n,r,pos+1,sum,choose,list);
		}
	}
	return 0;
}

void FindMnomProbs (int Nt, vector< vector<int> >& list, vector<double>& freq_pre, vector<double>& mnom_probs, vector<double>& fact_store) {
	mnom_probs.clear();
	//cout << "Size is " << list.size() << " " << fact_store.size() << "\n";
	for (int i=0;i<list.size();i++) {
		//Find log multinomial coefficient
		double L=fact_store[Nt];
		for (int j=0;j<list[0].size();j++) {
			L=L-fact_store[list[i][j]];
		}
		//Find probabilities
		for (int j=0;j<list[0].size();j++) {
			L=L+(list[i][j]*log(freq_pre[j]));
		}
		L=exp(L);
		mnom_probs.push_back(L);
	}
}

void GetInferredHaplotypeFrequencies(int i, int Nt, vector<double> freq_post, vector<double>& freq_post_prob, vector< vector<int> >& list) {
	freq_post_prob.clear();
	for (int j=0;j<freq_post.size();j++) {
		double q=(list[i][j]+0.)/(Nt+0.);
		freq_post_prob.push_back(q);
	}
}

void GetScaledPartialHaplotypeLikelihood (run_params p, double C, double& L, vector<int>& Npost, vector< vector<int> >& obs_freqs_post, vector<double>& freq_post_prob, vector<double>& inf_freqs_post, vector<double>& fact_store, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets) {
	//cout << "Start GetScale\n";
	for (int j=0;j<Npost.size();j++) { //Sets of partial haplotypes
		inf_freqs_post.clear();
		for (int k=0;k<obs_freqs_post[j].size();k++) {
			inf_freqs_post.push_back(0);
		}
		int xmatch=0;
		int indexk=0;
		if (full_haps[full_haps.size()-1].st.find("X")==0) { //List of haplotypes includes an XXXX type
			//Code for when there is an X type
			for (int k=0;k<hap_data_sets[j].size();k++) { //Partial haplotypes in subset
				if (hap_data_sets[j][k].match.size()>0) {
					if (hap_data_sets[j][k].match[0]==full_haps.size()-1) { //Corresponds to unobserved haplotype namely X
						if (xmatch==0) { //Only need to do this once
							inf_freqs_post[inf_freqs_post.size()-1]=freq_post_prob[hap_data_sets[j][k].match[0]];
							xmatch=1;
						}
					} else { //Corresponds to observed haplotype - find sum of frequencies of matching haplotypes
						for (int l=0;l<hap_data_sets[j][k].match.size();l++) {
							inf_freqs_post[indexk]=inf_freqs_post[indexk]+freq_post_prob[hap_data_sets[j][k].match[l]];
						}
						
						indexk++;
					}
				}
			}
		} else {
			//Code for when there is no X type
			for (int k=0;k<hap_data_sets[j].size();k++) { //Partial haplotypes in subset
				//Corresponds to observed haplotype - find sum of frequencies of matching haplotypes
				if (hap_data_sets[j][k].match.size()>0) {
					for (int l=0;l<hap_data_sets[j][k].match.size();l++) {
						inf_freqs_post[indexk]=inf_freqs_post[indexk]+freq_post_prob[hap_data_sets[j][k].match[l]];
					}
					indexk++;
				}
			}
		}
		if (inf_freqs_post.size()>1) { //Don't bother calculating likelihoods where there is only one observed haplotype
			//For example, if the two haplotypes are GA and GC, don't assess single-locus data from the first locus
			//This can in theory be done, but can cause problems with the overall likelihood value being large
			//In practice we only care about relative likelihoods
			int rdepth=0;
			//N.B. Here Npost might not total the matching observations - calculate read depth for relevant samples
			for (int k=0;k<obs_freqs_post[j].size();k++) {
				rdepth=rdepth+obs_freqs_post[j][k];
			}
			L=L+DirichletMultiCalc(p,rdepth,C,obs_freqs_post[j],inf_freqs_post,fact_store);
		}
	}
}

