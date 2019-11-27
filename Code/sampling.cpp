#include "shared_haps.h"
#include "sampling.h"
#include <iostream>
#include <string>
#include <sstream>

void GenerateSample (run_params p, double C, vector<int>& Npre, vector<int>& Npost, vector<double> freq_pre, vector<double> freq_post, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets, vector< vector< vector<mhap> > >& all_hap_data_sets, gsl_rng *rgen) {
	ofstream mhap_file;
	if (p.print_sample==1) {
		mhap_file.open(p.traj_out);
	}
	vector< vector<mhap> > new_hap_data;
	for (int i=0;i<Npre.size();i++) {
		vector<int> sample_pre = DirMultSampling(Npre[i],freq_pre,C,rgen);
		vector<int> sample_post = DirMultSampling(Npost[i],freq_post,C,rgen);
		vector<int> sites=hap_data_sets[i][0].n_loci;
		//Combine partial haplotypes into single strings
		vector<string> p_haps;
		for (int j=0;j<full_haps.size();j++) {
			string str="";
			for (int k=0;k<sites.size();k++) {
				str=str+full_haps[j].st[sites[k]];
			}
			p_haps.push_back(str);
		}
		//Find unique haplotypes - numerical values that increase given a not-previously-seen partial haplotype
		vector<int> hap_uniq;
		int max=0; //Increment to get number of unique haplotypes
		for (int i=0;i<p_haps.size();i++) {
			hap_uniq.push_back(max);
			for (int j=0;j<i;j++) {
				if (p_haps[j]==p_haps[i]) {
					hap_uniq[i]=hap_uniq[j];
					max--;
					break;
				}
			}
			max++;
		}
		//Find accompanying unique haplotype sequences
		int index=0;
		vector<string> phap_uniq;
		for (int j=0;j<hap_uniq.size();j++) {
			if (hap_uniq[j]==index) {
				phap_uniq.push_back(p_haps[j]);
				index++;
			}
		}
		//Construct sums of numbers of observations of each partial haplotype.
		//Note here that sample_pre and sample_post are the sampled full haplotype frequencies
		vector<int> phap_sample_pre;
		vector<int> phap_sample_post;
		for(int i=0;i<=max;i++) {
			phap_sample_pre.push_back(0);
			phap_sample_post.push_back(0);
		}
		for (int j=0;j<hap_uniq.size();j++) {
			phap_sample_pre[hap_uniq[j]]=phap_sample_pre[hap_uniq[j]]+sample_pre[j];
			phap_sample_post[hap_uniq[j]]=phap_sample_post[hap_uniq[j]]+sample_post[j];
		}
		//Generate multi-locus output here.
		vector<mhap> new_mh;
		for (int j=0;j<phap_uniq.size();j++) {
			mhap mh;

			if (p.print_sample==1) {
				mhap_file << sites.size() << " ";
			}
			for (int k=0;k<sites.size();k++) {
				if (p.print_sample==1) {
					mhap_file << sites[k] << " ";
				}
				mh.loci.push_back(sites[k]);
			}
			mh.st=phap_uniq[j];
			mh.times.push_back(0);
			mh.times.push_back(1);
			mh.obs.push_back(phap_sample_pre[j]);
			mh.obs.push_back(phap_sample_post[j]);
			if (p.print_sample==1) {
				mhap_file << phap_uniq[j] << " 2 0 " << phap_sample_pre[j] << " 1 " << phap_sample_post[j] << "\n";
			}
			new_mh.push_back(mh);
		}
		new_hap_data.push_back(new_mh);
	}
	all_hap_data_sets.push_back(new_hap_data);
}


vector<int> DirMultSampling(int N, vector<double> &freqs, double C, const gsl_rng *r) {
	
	//Define the prior, alpha, for the Dirichlet distribution
	double alpha[freqs.size()];
	for(unsigned int i=0;i<freqs.size();i++) {
		alpha[i] = C*freqs[i];
	}
	
	//Sample frequencies p from the Dirichlet distribution
	double p[freqs.size()]; //Placeholder
	gsl_ran_dirichlet(r, freqs.size(), alpha, p);
	
	//Convert p to a vector
	vector<double> pVec;
	for(unsigned int i=0;i<freqs.size();i++) {
		pVec.push_back(p[i]);
	}
	
	//Return multinomial sample based on Dirichlet prior
	return multinomialSampling(N,pVec,r);
}

vector<int> multinomialSampling(int N, vector<double> p, const gsl_rng *r) {
	
	size_t dim = p.size();
	unsigned int n[dim];
	double* pPointer = &p[0];
	
	gsl_ran_multinomial(r,dim,N,pPointer,n);
	
	vector<int> result(n, n + sizeof n / sizeof n[0]);
	
	return(result);
	
}

void FlattenData (int i, vector< vector< vector<mhap> > >& all_hap_data_sets, vector<mhap>& hap_data) {
	vector< vector<mhap> > hap_data_sets=all_hap_data_sets[i];
	for (int j=0;j<hap_data_sets.size();j++) {
		for (int k=0;k<hap_data_sets[j].size();k++) {
			hap_data.push_back(hap_data_sets[j][k]);
		}
	}
}

void CalculateVarianceMatrices(vector< vector<double> >& all_pre_freqs, vector< vector<double> >& all_post_freqs, vector<double>& pre_var, vector<double>& post_var) {
	//N.B. These are diagonal matrices by approximation
	//Calculate mean frequencies from data.  Only use these to get the sample variance
	vector<double> freq_pre;
	vector<double> freq_post;
	for (int i=0;i<all_pre_freqs[0].size();i++) {
		freq_pre.push_back(0);
		freq_post.push_back(0);
	}
	for (int i=0;i<all_pre_freqs.size();i++) {
		for (int j=0;j<all_pre_freqs[0].size();j++) {
			freq_pre[j]=freq_pre[j]+all_pre_freqs[i][j];
			freq_post[j]=freq_post[j]+all_post_freqs[i][j];
		}
	}
	for (int i=0;i<all_pre_freqs[0].size();i++) {
		freq_pre[i]=freq_pre[i]/all_pre_freqs.size();
		freq_post[i]=freq_post[i]/all_post_freqs.size();
	}

	//Initialise the variances
	for (int i=0;i<all_pre_freqs[0].size();i++) {
		pre_var.push_back(0);
		post_var.push_back(0);
	}
	//Calculate variance terms
	for (int i=0;i<all_pre_freqs.size();i++) {
		for (int j=0;j<all_pre_freqs[0].size();j++) {
			double diff=all_pre_freqs[i][j]-freq_pre[j];
			diff=pow(diff,2);
			pre_var[j]=pre_var[j]+diff;
		}
	}
	for (int i=0;i<all_post_freqs.size();i++) {
		for (int j=0;j<all_post_freqs[0].size();j++) {
			double diff=all_post_freqs[i][j]-freq_post[j];
			diff=pow(diff,2);
			post_var[j]=post_var[j]+diff;
		}
	}
	for (int j=0;j<all_post_freqs[0].size();j++) {
		pre_var[j]=pre_var[j]/(all_post_freqs.size()-1);
		post_var[j]=post_var[j]/(all_post_freqs.size()-1);
	}
}
