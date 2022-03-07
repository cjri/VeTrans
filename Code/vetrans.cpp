//Applies linear model to unobserved haplotypes.  Retains loci for which there is no variation.

//Efficiency: Around 72% of time in Dirichlet Multinomial likelihood calculation.  Not too shabby

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <time.h>

using namespace std;

#include "shared_haps.h"
#include "reconstruct_haplotypes.h"
#include "sample_haplotypes.h"
#include "sampling.h"
#include "combine_likelihoods.h"
#include "calc_bottleneck.h"
#include "io.h"

int main(int argc, const char **argv) {

	if (argc==1) {
		cout << "Error: A method must be speficied here.  Methods are:\n";
		cout << "  reconstruct_haplotypes: Does haplotype reconstruction\n";
		return 0;
	}
	
	run_params p;
	GetOptions(p,argc,argv);
	if (p.err==1) {
		cout << "Error in input flag " << p.err_flag << " is not understood\n";
		return 0;
	}

	
	string method=argv[1];
	
	//Structures for recovering reconstruction data
	vector<haplo> haplotypes;
	vector<double> pre_freqs;
	vector<double> post_freqs;
    vector< vector<double> > multi_freqs;
	double best_bic;
	
	if (method.compare("reconstruct_haplotypes")==0) {
		//May want to: 1.  Call haplotypes and print.  2.  Call haplotypes and return
		ReconstructHaplotypes (p,best_bic,haplotypes,pre_freqs,post_freqs);

		//Output from ReconstructHaplotypes
		OutputReconstructHaplotypes(haplotypes,pre_freqs,post_freqs,best_bic);

    } else if (method.compare("reconstruct_haplotypes_multi")==0) {
        cout << "Here now\n";
        
        //Designed for more general haplotype reconstruction
        ReconstructHaplotypesMulti (p,best_bic,haplotypes,multi_freqs);

        OutputReconstructHaplotypesMulti(p,haplotypes,multi_freqs,best_bic);
        
	} else if (method.compare("find_haplotypes")==0) {
		vector<haplo> haplotypes_store;
		vector<double> pre_freqs_store;
		vector<double> post_freqs_store;

		p.n_haps=1;
		double opt_bic=1e10;
		int done=0;
		while (done==0) {
			ReconstructHaplotypes (p,best_bic,haplotypes,pre_freqs,post_freqs);
			cout << "# of haplotypes = " << p.n_haps << "; BIC = " << best_bic << "\n";
			if (best_bic<opt_bic-p.delta_bic) {
				opt_bic=best_bic;
				p.n_haps++;
				haplotypes_store=haplotypes;
				pre_freqs_store=pre_freqs;
				post_freqs_store=post_freqs;
			} else {
				cout << "No significant BIC improvement: End calculation\n";
				ofstream out_file;
				out_file.open("Inferred_haplotypes.out");
				//out_file << "BIC " << opt_bic << "\n";
				for (int i=0;i<haplotypes_store.size();i++) {
					for (int j=0;j<haplotypes_store[i].seq.size();j++) {
						out_file << haplotypes_store[i].seq[j];
					}
					out_file << " " << pre_freqs_store[i] << " " << post_freqs_store[i] << "\n";
				}
				done=1;
			}
		}
        
    } else if (method.compare("find_haplotypes_multi")==0) {
        vector<haplo> haplotypes_store;
        vector< vector<double> > multi_freqs_store;

        p.n_haps=1;
        double opt_bic=1e10;
        int done=0;
        while (done==0) {
            ReconstructHaplotypesMulti (p,best_bic,haplotypes,multi_freqs);
            cout << "# of haplotypes = " << p.n_haps << "; BIC = " << best_bic << "\n";
            if (best_bic<opt_bic-p.delta_bic) {
                opt_bic=best_bic;
                p.n_haps++;
                haplotypes_store=haplotypes;
                multi_freqs_store=multi_freqs;
            } else {
                cout << "No significant BIC improvement: End calculation\n";
                ofstream out_file;
                out_file.open("Inferred_haplotypes.out");
                //out_file << "BIC " << opt_bic << "\n";
                for (int i=0;i<haplotypes_store.size();i++) {
                    for (int j=0;j<haplotypes_store[i].seq.size();j++) {
                        out_file << haplotypes_store[i].seq[j];
                    }
                    for (int s=0;s<p.n_samples;s++) {
                        out_file << " " << multi_freqs[s][i];
                    }
                }
                done=1;
            }
        }

	
	} else if (method.compare("resample")==0) {
		//Randomly resamples Multi_locus_trajectories.out using the inferred haplotypes.  Output to Multi_locus_trajectories_sample.out
		p.sample_reps=1;
		p.print_sample=1;
		vector< vector< vector<mhap> > > all_hap_data_sets;
		SampleData(p,all_hap_data_sets);
	
	} else if (method.compare("resample_multi")==0) {
		//Multiple resampling.  Output to Multi_locus_trajectories_sample?.out
		//--reps flag says how many resamples
		vector< vector< vector<mhap> > > all_hap_data_sets;
		p.print_sample=0;
		SampleData(p,all_hap_data_sets);
		PrintResampleMulti(all_hap_data_sets);
		
	} else if (method.compare("find_variance")==0) {
		vector< vector< vector<mhap> > > all_hap_data_sets;
		p.print_sample=0;
		SampleData(p,all_hap_data_sets);
		vector< vector<double> > all_pre_freqs;
		vector< vector<double> > all_post_freqs;
		vector<haplo> full_haps; //Structure for full haplotypes
		//Read in haplotype set
		GetFullHaplotypes(p,full_haps);
		for (int i=0;i<all_hap_data_sets.size();i++) {
			vector<mhap> hap_data;
			FlattenData (i,all_hap_data_sets,hap_data);
			vector<haplo> haplotypes;
			vector<double> pre_freqs;
			vector<double> post_freqs;
			ReconstructHaplotypesFixedFromInternal (p,best_bic,haplotypes,pre_freqs,post_freqs,hap_data,full_haps);
			all_pre_freqs.push_back(pre_freqs);
			all_post_freqs.push_back(post_freqs);
		}
		vector<double> pre_var;
		vector<double> post_var;
		CalculateVarianceMatrices(all_pre_freqs,all_post_freqs,pre_var,post_var);
		//Print variance matrices.  Name can be specified by --var_name option
		PrintVariances(p,pre_var,post_var);
		
	} else if (method.compare("calc_bottleneck")==0) {
		CalcBottleneck(p);
		
	} else if (method.compare("combine_likelihoods")==0) {
		CombineLikelihoods(p);
		
	} else {
		cout << "Method not recognised\n";
	}
	
	return 0;
}
