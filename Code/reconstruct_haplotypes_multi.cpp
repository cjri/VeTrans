#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <time.h>

using namespace std;

#include "shared_haps.h"
#include "io.h"
#include "reconstruct_haplotypes.h"
#include "reconstruct_haplotypes_multi.h"

int ReconstructHaplotypesMulti (run_params p, double& best_bic, vector<haplo>& haplotypes, vector< vector<double> >& multi_freqs)  {
	if (p.verb==1) {
		cout << "Running code to reconstruct haplotypes:\n";
		cout << " Read multi-locus trajectory data from " << p.traj_in << "\n";
		cout << "   Change with --traj_in\n";
		cout << " Conduct reconstruction with " << p.n_haps << " haplotypes\n";
		cout << "   Change with --n_haps\n";
		if (p.read_hap==1) {
			cout << " Read in haplotype data from " << p.fullhap_in << "\n";
			cout << "   Change with --fullhap_in\n";
		}
		cout << "\n";
	}
	
	double C=p.c;
	time_t timer_s;
	time(&timer_s);

	//Initialise random number generator
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, p.seed);
	
    cout << "Get variant data\n";
	vector<mhap> hap_data;
	GetVariantDataMulti(p,hap_data);
	
    
    cout << "Calculate loci\n";
	//Calculate full set of loci
	vector<int> loci;
	CalculateFullLocusSet(hap_data,loci);
	sort(loci.begin(),loci.end());

	//Update data with n_loci values: Loci numbers become 0, 1, 2, 3, etc.
	CalculateNLociData(hap_data,loci);
	
	if (p.verb==1) {
		cout << "Loci\n";
		for (int i=0;i<loci.size();i++) {
			cout << loci[i] << "\n";
		}
	}
	
	//Update data with set e.g. which cover which loci
	int n_sets;
	CalculateSetData(n_sets,hap_data);
	
	if (p.verb==1) {
		cout << "Sets\n";
		for (int i=0;i<hap_data.size();i++) {
			cout << i << " " << hap_data[i].set << "\n";
		}
	}
	
	//Calculate N values
    vector< vector<int> > Nmulti;
    for (int s=0;s<p.n_samples;s++) {
        vector<int> temp;
        Nmulti.push_back(temp);
    }
	CalculateObsTotalNMulti (p,n_sets,Nmulti,hap_data);

    /*cout << "Nmulti" << "\n";
    for (int i=0;i<Nmulti.size();i++) {
        for (int j=0;j<Nmulti[i].size();j++) {
            cout << Nmulti[i][j] << " ";
        }
        cout << "\n";
    }*/
    
	vector<double> fact_store;
	FindLogFactMulti(fact_store,Nmulti);

	//Split partial haplotype data according to which loci are spanned by each read
	vector< vector<mhap> > hap_data_sets;
	SplitHaplotypeData (n_sets,hap_data,hap_data_sets);
	
	vector<haplo> full_haps; //Structure for full haplotypes
	
	//Break point here - to simply read in a set of haplotypes and do the frequency optimisation for these haplotypes
	if (p.read_hap==1) {
		//Read in haplotype set
		GetFullHaplotypes(p,full_haps);

		//Add an X haplotype if none such exists
		CheckXHap (loci,full_haps);

		//Optimise full haplotype frequencies
		OptimiseFrequenciesFixedFullHapsMulti (p,C,multi_freqs,full_haps,Nmulti,fact_store,hap_data_sets,rgen);
		haplotypes=full_haps;
		return 0;
	}
	
	//Generate initial set of full haplotypes - begin full optimisation process
	ConstructHaplotypeSet(p.n_haps,loci,hap_data,full_haps,rgen);
	//Check all haplotypes are different
	int check=0;
	CheckUniqueHaplotypes(check,hap_data,full_haps,rgen);
	if (check>=1000) {
		cout << "Can't create sufficient haplotypes.  Terminating...\n";
		return 0;
	}
	if (p.verb==1) {
		cout << "Full haps\n";
		for (int i=0;i<full_haps.size();i++) {
			for (int j=0;j<full_haps[i].seq.size();j++) {
				cout << full_haps[i].seq[j];
			}
			cout << "\n";
		}
	}
	
	//From here: Haplotype specific after initial generation
	vector< vector<haplo> > full_haps_store; //Keep best reconstruction with likelihoods
	vector< vector<haplo> > full_haps_archive; //Keep all reconstructions - don't do more than once
	vector<double> lL_store;
	
	int firsthapset=1;
    vector<double> lL;
    for (int s=0;s<p.n_samples;s++) {
        lL.push_back(-1e10);
    }
	double lL_tot=-1e10; //Total likelihood for current set of full haplotypes
	double lL_hapbest=-1e10; //Maximum likelihood for a set of full haplotypes
    vector< vector<double> > hap_freqs_multi;
    vector< vector<double> > hap_freqs_best_multi;
    
	double run_time=0;
    vector< vector< vector<double> > > multi_store;
	vector<double> logs_store;
	int ct=1;
	
	for (int hapset=0;hapset<p.hap_its;hapset++) { //Loop over sets of full haplotypes
        //cout << "Hapset " << hapset << "\n";
		ct=CheckTime (p,timer_s);
		if (ct==0) {
			return 0;
		}
		if (firsthapset==0) {
			//Full set of attempted haplotype sets
			full_haps_archive.push_back(full_haps);
			
			//Evaluate full haplotype set
			if (lL_tot>lL_hapbest) {
				full_haps_store.push_back(full_haps);
				lL_store.push_back(lL_tot);
				lL_hapbest=lL_tot;
                hap_freqs_best_multi=hap_freqs_multi;
                multi_freqs=hap_freqs_multi;
				//cout << "Better likelihood: " << lL_tot << "\n";
				double bic=FindBIC(lL_tot,full_haps);
				best_bic=bic;
				if (p.verb==1) {
					cout << "Better BIC " << bic << "\n";
                    PrintInferenceMulti (p,full_haps,hap_freqs_multi);
				}
			} else {
				if (run_time>p.run_time/10||hapset>p.hap_its/10) {
					full_haps=full_haps_store[full_haps_store.size()-1];
				}
			}
			
			//Change full haplotype set
			int changes=ChangeHaplotypeRandom(p,full_haps,haplotypes,hap_data_sets,full_haps_archive,rgen);
			//cout << "Changes required " << changes << "\n";
			if (changes==100) {
				//cout << "End\n";
				//Final output
				FinalOutputMulti (p,lL_hapbest,haplotypes,hap_freqs_best_multi,full_haps_store);
				return 0;
			}
		}
		firsthapset=0;
	
		//Calculate matching between partial haplotype data and full haplotypes
		MatchFullHapsSplit(hap_data_sets,full_haps);
	
		//Clear storage for replicate frequency inferences
        multi_store.clear();
        vector< vector<double> > temp;
        for (int s=0;s<p.n_samples;s++) {
            multi_store.push_back(temp);
        }
        logs_store.clear();
       // cout << "Start reps\n";
		for (int rep=0;rep<p.freq_rep;rep++) { //Independent frequency starting points
            //cout << "Rep " << rep << "\n";
			//Set initial haplotype frequencies
            hap_freqs_multi.clear();

            for (int s=0;s<p.n_samples;s++) {
                vector<double> freqs;
                SetInitialFreqs(freqs,full_haps,rgen);
                hap_freqs_multi.push_back(freqs);
            }

            //Check initial frequencies - last cannot be >1% if it is an X-type haplotype
            for (int s=0;s<p.n_samples;s++) {
                CheckNoiseFreq(hap_freqs_multi[s]);
            }
            
			//Construct observations X.  Here the last obs combines all of the 'noise' observations
            vector< vector< vector<int> > > obs_freqs_multi;
            for (int s=0;s<p.n_samples;s++) {
                vector< vector<int> > ofreq;
                obs_freqs_multi.push_back(ofreq);
            }
            ConstructObservationsXMulti (p,hap_data_sets,full_haps,obs_freqs_multi);

			//Optimise frequencies to fit to data
            for (int s=0;s<p.n_samples;s++) {
                lL[s]=-1e10;
            }
	
            OptimiseFrequenciesMulti(p,lL,C,Nmulti,full_haps,hap_freqs_multi,obs_freqs_multi,hap_data_sets,fact_store,rgen);
            //Here in the Multi conversion.
            lL_tot=0;
            for (int s=0;s<p.n_samples;s++) {
                lL_tot=lL_tot+lL[s];
            }
            logs_store.push_back(lL_tot);
            for (int s=0;s<p.n_samples;s++) {
                multi_store[s].push_back(hap_freqs_multi[s]);
            }
		}
		//Find best replicate
        FindBestRepMulti (p,lL_tot,hap_freqs_multi,logs_store,multi_store);
	}
	FinalOutputMulti (p,lL_hapbest,haplotypes,hap_freqs_best_multi,full_haps_store);
    for (int s=0;s<p.n_samples;s++) {
        multi_freqs[s]=hap_freqs_best_multi[s];
    }

	return 0;
}

void CalculateObsTotalNMulti (run_params& p, int n_sets, vector< vector<int> >& Nmulti, vector<mhap>& hap_data) {
    for (int set=0;set<n_sets;set++) {
        vector<int> n;
        for (int s=0;s<p.n_samples;s++) {
            n.push_back(0);
        }
        for (int i=0;i<hap_data.size();i++) {
            if (hap_data[i].set==set) {
                for (int s=0;s<p.n_samples;s++) {
                    n[s]=n[s]+hap_data[i].obs[s];
                }
            }
        }
        for (int s=0;s<p.n_samples;s++) {
            Nmulti[s].push_back(n[s]);
        }
    }
}


void ConstructObservationsXMulti (run_params& p, vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector< vector<int> > >& obs_freqs_multi) {
    //cout << "Construct ObsXMulti\n";
    for (int i=0;i<hap_data_sets.size();i++) {
        //Observations corresponding to known haplotypes
        vector< vector<int> > opmult;
        vector<int> op;
        for (int s=0;s<p.n_samples;s++) {
            opmult.push_back(op);
        }
        //Observations corresponding to X haplotype
        vector<int> oxmult;
        for (int s=0;s<p.n_samples;s++) {
            oxmult.push_back(0);
        }
        for (int j=0;j<hap_data_sets[i].size();j++) {  //Construct set of all matching observations + sum of all non-observed observations
            if (hap_data_sets[i][j].match[0]<full_haps.size()-1) {
                for (int s=0;s<p.n_samples;s++) {
                    opmult[s].push_back(hap_data_sets[i][j].obs[s]);
                }
            } else {
                for (int s=0;s<p.n_samples;s++) {
                    oxmult[s]=oxmult[s]+hap_data_sets[i][j].obs[s];
                }
            }
        }
        for (int s=0;s<p.n_samples;s++) {
            if (opmult[s].size()==0) {
                opmult[s].push_back(0);
            }
        }
        for (int s=0;s<p.n_samples;s++) {
            opmult[s].push_back(oxmult[s]);
            obs_freqs_multi[s].push_back(opmult[s]);
        }
    }
}

void FindLogFactMulti(vector<double>& fact_store, vector< vector<int> >& Nmulti){
    int maxn=0;
    for (int i=0;i<Nmulti.size();i++) {
        for (int j=0;j<Nmulti[i].size();j++) {
            if (Nmulti[i][j]>maxn) {
                maxn=Nmulti[i][j]+100;
            }
        }
    }
    
    double logN=0;
    fact_store.push_back(0);
    for (int i=1;i<=maxn;i++) {
        logN=logN+log(i);
        fact_store.push_back(logN);
    }
}

void OptimiseFrequenciesMulti (run_params p, vector<double>& lL, double C, vector< vector<int> >& Nmulti, vector<haplo>& full_haps, vector< vector<double> >& hap_freqs_multi, vector< vector< vector<int> > >& obs_freqs_multi, vector< vector<mhap> >& hap_data_sets, vector<double>& fact_store, gsl_rng *rgen) {
    //cout << "OptimiseFrequencesMulti\n";
    //Optimise haplotype frequencies for the data given a set of haplotypes
    double changex=0.01;
    vector<double> bestL;
    for (int s=0;s<p.n_samples;s++) {
        bestL.push_back(-1e10);
    }
    vector< vector<double> > inf_freqs_multi;
    for (int s=0;s<p.n_samples;s++) {
        vector<double> temp;
        inf_freqs_multi.push_back(temp);
    }
    vector< vector<double> > hap_freqs_multi_store;
    for (int s=0;s<p.n_samples;s++) {
        vector<double> temp;
        hap_freqs_multi_store.push_back(temp);
    }

    double L=0;
    int first=1;
    int fails=0;
    for (int it=0;it<p.freq_its;it++) {
        //cout << "Iteration " << it << "\n";
        if (first==0) {
            for (int s=0;s<p.n_samples;s++) {
                if (lL[s]>bestL[s]) {
                    hap_freqs_multi_store[s]=hap_freqs_multi[s];
                    bestL[s]=lL[s];
                    fails=0;
                } else {
                    hap_freqs_multi[s]=hap_freqs_multi_store[s];
                    fails++;
                }
            }
        }
        if (fails==p.freq_term) {
            for (int s=0;s<p.n_samples;s++) {
                lL[s]=bestL[s];
            }
            break;
        }
        first=0;
        
        //Change frequencies
        for (int s=0;s<p.n_samples;s++) {
            ChangeFreq (changex,full_haps,hap_freqs_multi[s],rgen);
        }

        //Calculate likelihood - multiple samples
        for (int s=0;s<p.n_samples;s++) {
            lL[s]=0;
        }
        
        for (int i=0;i<hap_data_sets.size();i++) { //Sets of partial haplotypes
            //cout << "Set " << i << "\n";
            for (int s=0;s<p.n_samples;s++) {
                inf_freqs_multi[s].clear();
            }
            for (int s=0;s<p.n_samples;s++) {
                for (int k=0;k<obs_freqs_multi[s][i].size();k++) { //Inference has same number of elements as observations
                    inf_freqs_multi[s].push_back(0);
                }
            }
            int xmatch=0;
            int indexj=0;
            for (int j=0;j<hap_data_sets[i].size();j++) { //Partial haplotypes in subset
                //cout << "PH " << j << "\n";
                if (hap_data_sets[i][j].match[0]==full_haps.size()-1) { //Corresponds to unobserved haplotype namely X
                    if (xmatch==0) { //Only need to do this once
                        for (int s=0;s<p.n_samples;s++) {
                            inf_freqs_multi[s][inf_freqs_multi[s].size()-1]=hap_freqs_multi[s][hap_data_sets[i][j].match[0]];
                        }
                        xmatch=1;
                    }
                } else { //Corresponds to observed haplotype - find sum of frequencies of matching haplotypes
                    for (int k=0;k<hap_data_sets[i][j].match.size();k++) {
                        for (int s=0;s<p.n_samples;s++) {
                            inf_freqs_multi[s][indexj]=inf_freqs_multi[s][indexj]+hap_freqs_multi[s][hap_data_sets[i][j].match[k]];
                        }
                    }
                    indexj++;
                }
            }
            
            //NB: Sometimes there is no match to a full haplotype.  In this case the inferred frequencies will not sum to 1.  Deal with this by putting them in the noise frequency, then correcting if this is >1%.  This seems to work out as a proportionate penalty term though perhaps there is a better solution.
            for (int s=0;s<p.n_samples;s++) {
                MissingFreqRedist(inf_freqs_multi[s]);
            }
            
            //Calculate actual likelihood here...
            for (int s=0;s<p.n_samples;s++) {
                /*cout << "N " << Nmulti[s][i] << "\n";
                cout << "obs_freqs_multi\n";
                for (int k=0;k<obs_freqs_multi[s][i].size();k++){
                    cout << obs_freqs_multi[s][i][k] << " ";
                }
                cout << "\n";
                cout << "inf_freqs_multi\n";
                for (int k=0;k<inf_freqs_multi[s].size();k++){
                    cout << inf_freqs_multi[s][k] << " ";
                }
                cout << "\n";*/
                L=DirichletMultiCalc(p,Nmulti[s][i],C,obs_freqs_multi[s][i],inf_freqs_multi[s],fact_store);
                //cout << "L " << L << "\n";
                lL[s]=lL[s]+L;
            }
            /*for (int s=0;s<p.n_samples;s++) {
                cout << "s " << s << " " << lL[s] << "\n";
            }*/
        }
    }
    for (int s=0;s<p.n_samples;s++) {
        lL[s]=bestL[s];
    }
}

void FindBestRepMulti (run_params& p, double& lL_tot, vector< vector<double> >& hap_freqs_multi, vector<double>& logs_store, vector< vector< vector<double> > >& multi_store) {
    double maxL=-1e10;
    int index=-1;
    for (int i=0;i<logs_store.size();i++) {
        if (logs_store[i]>maxL) {
            index=i;
            maxL=logs_store[i];
        }
    }
    lL_tot=logs_store[index];
    for (int s=0;s<p.n_samples;s++) {
        hap_freqs_multi[s]=multi_store[s][index];
    }
}

void OptimiseFrequenciesFixedFullHapsMulti (run_params p, double C, vector< vector<double> >& multi_freqs, vector<haplo>& full_haps, vector< vector<int> >& Nmulti, vector<double>& fact_store, vector< vector<mhap> >& hap_data_sets, gsl_rng *rgen) {
    //Full optimisation for case in which we read in a fixed set of haplotypes
    MatchFullHapsSplit(hap_data_sets,full_haps);
    vector< vector< vector<double> > > multi_store;
    for (int s=0;s<p.n_samples;s++) {
        vector< vector<double> > temp;
        multi_store.push_back(temp);
    }
    vector<double> logs_store;
    double lL_tot=-1e10;
    vector< vector<double> > hap_freqs_multi;
    for (int rep=0;rep<p.freq_rep;rep++) { //Independent frequency starting points
        //Set initial frequencies
        hap_freqs_multi.clear();
        for (int s=0;s<p.n_samples;s++) {
            vector<double> freqs;
            SetInitialFreqs(freqs,full_haps,rgen);
            hap_freqs_multi.push_back(freqs);
        }
        //Check initial frequencies - last cannot be >1% if it is an X-type haplotype
        for (int s=0;s<p.n_samples;s++) {
            CheckNoiseFreq(hap_freqs_multi[s]);
        }

        
        //Match observations to full haplotype frequencies.  Here the last obs combines all of the 'noise' observations
        vector< vector< vector<int> > > obs_freqs_multi;
        ConstructObservationsXMulti (p,hap_data_sets,full_haps,obs_freqs_multi);

        //Optimise frequencies to fit to data
        vector<double> lL;
        for (int s=0;s<p.n_samples;s++) {
            lL.push_back(-1e10);
        }
        OptimiseFrequenciesMulti(p,lL,C,Nmulti,full_haps,hap_freqs_multi,obs_freqs_multi,hap_data_sets,fact_store,rgen);
        
        lL_tot=0;
        for (int s=0;s<p.n_samples;s++) {
            lL_tot=lL_tot+lL[s];
        }
        logs_store.push_back(lL_tot);
        for (int s=0;s<p.n_samples;s++) {
            multi_store[s].push_back(hap_freqs_multi[s]);
        }
    }
    
    //Find best replicate
    FindBestRepMulti (p,lL_tot,hap_freqs_multi,logs_store,multi_store);
    double bic=FindBIC(lL_tot,full_haps);
    if (p.verb==1) {
        cout << "Final BIC " << bic << "\n";
    }
    //Here in the Multi conversion.

    if (p.suppress_output==0) {
        PrintInferenceMulti(p,full_haps,hap_freqs_multi);
    }
    for (int s=0;s<p.n_samples;s++) {
        multi_freqs[s]=hap_freqs_multi[s];
    }
}

void PrintInferenceMulti (run_params& p, vector<haplo>& full_haps, vector< vector<double> >& hap_freqs_multi) {
    for (int i=0;i<full_haps.size();i++) {
        for (int j=0;j<full_haps[i].seq.size();j++) {
            cout << full_haps[i].seq[j];
        }
        for (int s=0;s<p.n_samples;s++) {
            cout << " " << hap_freqs_multi[s][i];
        }
        cout << "\n";
    }
}

void FinalOutputMulti (run_params p, double lL_hapbest, vector<haplo>& haplotypes, vector< vector<double> >& hap_freqs_best_multi, vector< vector<haplo> >& full_haps_store) {
    vector<haplo> full_haps=full_haps_store[full_haps_store.size()-1];
    vector< vector<double> > hap_freqs_multi=hap_freqs_best_multi;
    double bic=FindBIC(lL_hapbest,full_haps);
    ofstream out_file;
    string file="Inference_"+to_string(p.n_haps)+"_"+to_string(p.seed)+".out";
    out_file.open(file.c_str());
    out_file << "Final BIC " << bic << "\n";
    haplotypes=full_haps;
    for (int i=0;i<full_haps.size();i++) {
        for (int j=0;j<full_haps[i].seq.size();j++) {
            out_file << full_haps[i].seq[j];
            cout << full_haps[i].seq[j];
        }
        for (int s=0;s<p.n_samples;s++) {
            out_file << " " << hap_freqs_multi[s][i];
            cout << " " << hap_freqs_multi[s][i];
        }
        out_file << "\n";
        cout << "\n";
    }
}


