#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <time.h>

using namespace std;

#include "shared_haps.h"
#include "io.h"
#include "reconstruct_haplotypes.h"

int ReconstructHaplotypes (run_params p, double& best_bic, vector<haplo>& haplotypes, vector<double>& pre_freqs, vector<double>& post_freqs)  {
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
	
	vector<mhap> hap_data;
	GetVariantData(p,hap_data);
	
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
	vector<int> Npre;
	vector<int> Npost;
	CalculateObsTotalN (n_sets,Npre,Npost,hap_data);

	vector<double> fact_store;
	FindLogFact(fact_store,Npre,Npost);

	//Split partial haplotype data according to which loci are spanned by each read
	vector< vector<mhap> > hap_data_sets;
	SplitHaplotypeData (n_sets,hap_data,hap_data_sets);
	
	vector<haplo> full_haps; //Structure for full haplotypes
	
	//Break point here - to read in set haplotypes and do the frequency optimisation...
	if (p.read_hap==1) {
		//Read in haplotype set
		GetFullHaplotypes(p,full_haps);

		//Add an X haplotype if none such exists
		CheckXHap (loci,full_haps);

		//Optimise full haplotype frequencies
		OptimiseFrequenciesFixedFullHaps (p,C,pre_freqs,post_freqs,full_haps,Npre,Npost,fact_store,hap_data_sets,rgen);
		haplotypes=full_haps;
		return 0;
	}
	
	//Generate initial set of full haplotypes
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
	double lL_pre=-1e10; //Post-transmission likelihood
	double lL_post=-1e10; //Pre-transmission likelihood
	double lL_tot=-1e10; //Total likelihood for current set of full haplotypes
	double lL_hapbest=-1e10; //Maximum likelihood for a set of full haplotypes
	vector<double> hap_freqs_pre;
	vector<double> hap_freqs_post;
	vector<double> hap_freqs_pre_best;
	vector<double> hap_freqs_post_best;

	double run_time=0;
	vector< vector<double> > pre_store;
	vector< vector<double> > post_store;
	vector<double> logs_store;
	int ct=1;
	
	for (int hapset=0;hapset<p.hap_its;hapset++) { //Loop over sets of full haplotypes
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
				hap_freqs_pre_best=hap_freqs_pre;
				hap_freqs_post_best=hap_freqs_post;
				pre_freqs=hap_freqs_pre;
				post_freqs=hap_freqs_post;
				//cout << "Better likelihood: " << lL_tot << "\n";
				double bic=FindBIC(lL_tot,full_haps);
				best_bic=bic;
				if (p.verb==1) {
					cout << "Better BIC " << bic << "\n";
					PrintInference(full_haps,hap_freqs_pre,hap_freqs_post);
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
				FinalOutput (p,lL_hapbest,haplotypes,hap_freqs_pre_best,hap_freqs_post_best,full_haps_store);
				return 0;
			}
		}
		firsthapset=0;
	
		//Calculate matching between partial haplotype data and full haplotypes
		MatchFullHapsSplit(hap_data_sets,full_haps);
	
		//Clear storage for replicate frequency inferences
		pre_store.clear();
		post_store.clear();
		logs_store.clear();
		
		for (int rep=0;rep<p.freq_rep;rep++) { //Independent frequency starting points
		
			//Set initial haplotype frequencies
			hap_freqs_pre.clear();
			hap_freqs_post.clear();
			SetInitialFreqs(hap_freqs_pre,full_haps,rgen);
			SetInitialFreqs(hap_freqs_post,full_haps,rgen);
			//Check frequencies - last cannot be >1%
			CheckNoiseFreq(hap_freqs_pre);
			CheckNoiseFreq(hap_freqs_post);

			//Construct observations X.  Here the last obs combines all of the 'noise' observations
			vector< vector<int> > obs_freqs_pre;
			vector< vector<int> > obs_freqs_post;
			ConstructObservationsX (hap_data_sets,full_haps,obs_freqs_pre,obs_freqs_post);
	
			//Optimise frequencies to fit to data
			lL_pre=-1e10;
			lL_post=-1e10;
	
			OptimiseFrequencies (p,lL_pre,lL_post,C,Npre,Npost,full_haps,hap_freqs_pre,hap_freqs_post,obs_freqs_pre,obs_freqs_post,hap_data_sets,fact_store,rgen);
	
			lL_tot=lL_pre+lL_post;
		
			logs_store.push_back(lL_tot);
			pre_store.push_back(hap_freqs_pre);
			post_store.push_back(hap_freqs_post);
			
		}
		//Find best replicate
		FindBestRep (lL_tot,hap_freqs_pre,hap_freqs_post,logs_store,pre_store,post_store);
	}
	FinalOutput (p,lL_hapbest,haplotypes,hap_freqs_pre_best,hap_freqs_post_best,full_haps_store);
	pre_freqs=hap_freqs_pre_best;
	post_freqs=hap_freqs_post_best;

	return 0;
}

void CalculateFullLocusSet (vector<mhap>& hap_data, vector<int>& loci) {
	for (int i=0;i<hap_data.size();i++) {
		for (int j=0;j<hap_data[i].loci.size();j++) {
			int match=0;
			for (int k=0;k<loci.size();k++) {
				if (loci[k]==hap_data[i].loci[j]) {
					match=1;
				}
			}
			if (match==0) {
				loci.push_back(hap_data[i].loci[j]);
			}
		}
	}
}

void CalculateNLociData (vector<mhap>& hap_data, vector<int>& loci) { //Loci renumbered to 0, 1, 2, etc.
	for (int i=0;i<hap_data.size();i++) {
		for (int j=0;j<hap_data[i].loci.size();j++) {
			int match=-1;
			for (int k=0;k<loci.size();k++) {
				if (hap_data[i].loci[j]==loci[k]) {
					match=k;
				}
			}
			hap_data[i].n_loci.push_back(match);
		}
	}
}

void CalculateSetData (int& n_sets, vector<mhap>& hap_data) {  //Split ML reads by loci covered
	int index=0;
	n_sets=1;
	hap_data[0].set=0;
	if (hap_data.size()>1) {
		for (int i=1;i<hap_data.size();i++) {
			if (hap_data[i].loci!=hap_data[i-1].loci) {
				index++;
				n_sets++;
			}
			hap_data[i].set=index;
		}
	}
}

void ConstructHaplotypeSet (int n_haps, vector<int>& loci, vector<mhap>& hap_data, vector<haplo>& full_haps, gsl_rng *rgen) {
	for (int i=0;i<n_haps;i++) {
		haplo hap;
		ConstructHaplotypeDeNovo (hap,loci,hap_data,rgen);
		full_haps.push_back(hap);
	}
	haplo hapx;
	ConstructHaplotypeX(hapx,loci);
	full_haps.push_back(hapx);
}


void ConstructHaplotypeDeNovo (haplo& hap, vector<int>& loci, vector<mhap>& hap_data, gsl_rng *rgen) {
	//Start with all X then randomly replace with haplotypes from data
	for (int i=0;i<loci.size();i++) {
		hap.seq.push_back('X');
	}
	int nx=loci.size();
	
	while (nx>0) {
		//Pick a haplotype and add it
		int h=(gsl_rng_uniform(rgen)*hap_data.size());
		for (int j=0;j<hap_data[h].n_loci.size();j++) {
			hap.seq[hap_data[h].n_loci[j]]=hap_data[h].st[j];
		}
		nx=0;
		for (int i=0;i<loci.size();i++) {
			//cout << hap.seq[i];
			if (hap.seq[i]=='X') {
				nx++;
			}
		}
	}
}

void ConstructHaplotypeX (haplo& hap, vector<int>& loci) {
	for (int i=0;i<loci.size();i++) {
		hap.seq.push_back('X');
	}
}

void CheckUniqueHaplotypes (int& check, vector<mhap>& hap_data, vector<haplo>& full_haps, gsl_rng *rgen) {
	//Ensure that all of the haplotypes are different from one another
	check=0;
	for (int i=1;i<full_haps.size();i++) {
		int uniq=0;
		while (uniq==0&&check<1001) {
			check++;
			uniq=1;
			for (int j=0;j<i;j++) {
				if (full_haps[j].seq==full_haps[i].seq) {
					uniq=0;
				}
			}
			if (uniq==0) {
				//Make random change to i
				int h=(gsl_rng_uniform(rgen)*hap_data.size());
				for (int j=0;j<hap_data[h].n_loci.size();j++) {
					full_haps[i].seq[hap_data[h].n_loci[j]]=hap_data[h].st[j];
				}
			}
		}
	}
}



int ChangeHaplotypeRandom (run_params p, vector<haplo>& full_haps, vector<haplo>& haplotypes, vector< vector<mhap> >& hap_data_sets, vector< vector<haplo> >& full_haps_archive, gsl_rng *rgen) {
	//Uniform random change to a single haplotype
	vector<haplo> full_haps_old = full_haps;
	int seen=1;
	int count=0;
	while (seen==1) {
		//cout << "Change haplotype\n";
		int samehaps=1;
		while (samehaps==1) {
			int samehap=1;
			while (samehap==1) {
				//Choose a random haplotype
				int r=(gsl_rng_uniform(rgen)*(full_haps.size()-1)); //Last haplotype is the X set
				vector<char> orig=full_haps[r].seq;
				//Make a random change - overwrite some loci with a partial haplotype
				int s=(gsl_rng_uniform(rgen)*hap_data_sets.size()); //Random set
				int h=(gsl_rng_uniform(rgen)*hap_data_sets[s].size()); //Random partial haplotype in set.  N.B. Could skew random choice here?
				for (int j=0;j<hap_data_sets[s][h].n_loci.size();j++) {
					full_haps[r].seq[hap_data_sets[s][h].n_loci[j]]=hap_data_sets[s][h].st[j];
				}
				if (full_haps[r].seq!=orig) {
					samehap=0;
				}
			}
			
			//Check haplotypes are unique in new model
			CheckUniqueHaplotypesSplit (hap_data_sets,full_haps,rgen);
			
			//Check haplotype set against archive of haplotype sets
			seen=CheckHaplotypesArchive(full_haps,full_haps_archive);
			count++;
			if (count==p.hap_term) {
				if (p.verb==1) {
					cout << "No further haplotypes found after " << p.hap_term << " iterations\n";
				}
				haplotypes=full_haps;

				return count;
			}
			for (int i=0;i<full_haps.size();i++) {
				if (full_haps[i].seq!=full_haps_old[i].seq) {
					samehaps=0;
				}
			}
		}
		
	}
	return count;
}


void CheckUniqueHaplotypesSplit (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, gsl_rng *rgen){
	//Ensure that all of the haplotypes are different from one another
	for (int i=1;i<full_haps.size();i++) {
		int uniq=0;
		while (uniq==0) {
			uniq=1;
			for (int j=0;j<i;j++) {
				if (full_haps[j].seq==full_haps[i].seq) {
					uniq=0;
				}
			}
			if (uniq==0) {
				//Make random change to i
				int s=(gsl_rng_uniform(rgen)*hap_data_sets.size()); //Random set
				int h=(gsl_rng_uniform(rgen)*hap_data_sets[s].size()); //Random partial haplotype in set.
				for (int j=0;j<hap_data_sets[s][h].n_loci.size();j++) {
					full_haps[i].seq[hap_data_sets[s][h].n_loci[j]]=hap_data_sets[s][h].st[j];
				}
			}
		}
	}
}

int CheckHaplotypesArchive (vector<haplo>& full_haps, vector< vector<haplo> >& full_haps_archive) {
	int seen=0; //Flag - has this combination of haplotypes ever been seen before?
	for (int i=0;i<full_haps_archive.size();i++) {
		int match_tot=0;
		for (int j=0;j<full_haps.size()-1;j++) {
			int match=0;
			for (int k=0;k<full_haps_archive[i].size();k++) {
				if (full_haps[j].seq==full_haps_archive[i][k].seq) {
					match=1;  //This haplotype exists in the previous set
					break;
				}
			}
			if (match==1) {
				match_tot++;
			} else {
				break;
			}
		}
		if (match_tot==full_haps.size()-1) {
			//This set of haplotypes is identical to a previous set of haplotypes
			seen=1;
			break;
		}
	}
	return seen;
}


void SplitHaplotypeData (int n_sets, vector<mhap>& hap_data, vector< vector<mhap> >& hap_data_sets) {
	for (int s=0;s<n_sets;s++) {
		vector<mhap> hap_temp;
		for (int i=0;i<hap_data.size();i++) {
			if (hap_data[i].set==s) {
				hap_temp.push_back(hap_data[i]);
			}
		}
		hap_data_sets.push_back(hap_temp);
	}
}

void CheckXHap (vector<int>& loci, vector<haplo>& full_haps) {
	if (full_haps[full_haps.size()-1].st.find("X")!=0) {
		haplo hapx;
		ConstructHaplotypeX(hapx,loci);
		full_haps.push_back(hapx);
	}
}

void MatchFullHapsSplit (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps) {
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


void CalculateObsTotalN (int n_sets, vector<int>& Npre, vector<int>& Npost, vector<mhap>& hap_data) {
	for (int set=0;set<n_sets;set++) {
		int n1=0;
		int n2=0;
		
		for (int i=0;i<hap_data.size();i++) {
			if (hap_data[i].set==set) {
				n1=n1+hap_data[i].obs[0];
				n2=n2+hap_data[i].obs[1];
			}
		}
		Npre.push_back(n1);
		Npost.push_back(n2);
	}
}


void ConstructObservationsX (vector< vector<mhap> >& hap_data_sets, vector<haplo>& full_haps, vector< vector<int> >& obs_freqs_pre, vector< vector<int> >& obs_freqs_post) {
	for (int i=0;i<hap_data_sets.size();i++) {
		vector<int> opr;
		vector<int> opp;
		int oxr=0; //Observations corresponding to X haplotype
		int oxp=0;
		for (int j=0;j<hap_data_sets[i].size();j++) {  //Construct set of all matching observations + sum of all non-observed observations
			if (hap_data_sets[i][j].match[0]<full_haps.size()-1) {
				opr.push_back(hap_data_sets[i][j].obs[0]);
				opp.push_back(hap_data_sets[i][j].obs[1]);
			} else {
				oxr=oxr+hap_data_sets[i][j].obs[0];
				oxp=oxp+hap_data_sets[i][j].obs[1];
			}
		}
		if (opr.size()==0) {
			opr.push_back(0);
		}
		if (opp.size()==0) {
			opp.push_back(0);
		}
		opr.push_back(oxr);
		opp.push_back(oxp);
		obs_freqs_pre.push_back(opr);
		obs_freqs_post.push_back(opp);
	}
	
}

void FindLogFact(vector<double>& fact_store, vector<int>& Npre, vector<int>& Npost){
	int maxn=0;
	for (int i=0;i<Npre.size();i++) {
		if (Npre[i]>maxn) {
			maxn=Npre[i]+100;
		}
	}
	for (int i=0;i<Npost.size();i++) {
		if (Npost[i]>maxn) {
			maxn=Npost[i]+100;
		}
	}
	
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=maxn;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
		//cout << "fact_store "<<i<<" "<<gsl_vector_get(fact_store,i)<<"\n";
	}
}


void SetInitialFreqs (vector<double>& hap_freqs, vector<haplo> full_haps, gsl_rng *rgen) {
	for (int i=0;i<full_haps.size();i++) {
		double x=gsl_rng_uniform(rgen);
		hap_freqs.push_back(x);
	}
	NormaliseFreqs(hap_freqs);
}


void NormaliseFreqs (vector<double>& hap_freqs) {
	double tot=0;
	for (int i=0;i<hap_freqs.size();i++) {
		tot=tot+hap_freqs[i];
		if (hap_freqs[i]<1e-10) {
			hap_freqs[i]=1e-10;
		}
	}
	for (int i=0;i<hap_freqs.size();i++) {
		hap_freqs[i]=hap_freqs[i]/tot;
	}
}

void CheckNoiseFreq (vector<double>& hap_freqs) {
	if (hap_freqs[hap_freqs.size()-1]>0.01) {
		double remainder=1-hap_freqs[hap_freqs.size()-1];
		hap_freqs[hap_freqs.size()-1]=0.01;
		for (int i=0;i<hap_freqs.size()-1;i++) {
			hap_freqs[i]=hap_freqs[i]*(0.99/remainder);
		}
	}
	NormaliseFreqs(hap_freqs);
}


void ChangeFreq (double changex, vector<haplo>& full_haps, vector<double>& hap_freqs, gsl_rng *rgen) {
	vector<double> hap_old=hap_freqs;
	int j=0;
	while (hap_freqs==hap_old) { //Ensure change in frequencies
		j=floor(gsl_rng_uniform(rgen)*(hap_freqs.size()));
		hap_freqs[j]=hap_freqs[j]+(gsl_rng_uniform(rgen)*changex)-(changex/2);
		NormaliseFreqs(hap_freqs);
		CheckNoiseFreq(hap_freqs);
	}
}

void OptimiseFrequencies (run_params p, double& lL_pre, double& lL_post, double C, vector<int>& Npre, vector<int>& Npost, vector<haplo>& full_haps, vector<double>& hap_freqs_pre, vector<double>& hap_freqs_post, vector< vector<int> >& obs_freqs_pre, vector< vector<int> >& obs_freqs_post, vector< vector<mhap> >& hap_data_sets, vector<double>& fact_store, gsl_rng *rgen) {
	//Optimise haplotype frequencies for the data given a set of haplotypes
	double changex=0.01;
	double bestL_pre=-1e10;
	double bestL_post=-1e10;
	vector<double> inf_freqs_pre;
	vector<double> inf_freqs_post;
	vector<double> hap_freqs_pre_store=hap_freqs_pre;
	vector<double> hap_freqs_post_store=hap_freqs_post;
	double L=0;
	int first=1;
	int fails=0;
	for (int it=0;it<p.freq_its;it++) {
		if (first==0) {
			if (lL_pre>bestL_pre) { //Better pre-transmission frequencies
				hap_freqs_pre_store=hap_freqs_pre;
				bestL_pre=lL_pre;
				fails=0;
			} else {
				hap_freqs_pre=hap_freqs_pre_store;
				fails++;
			}
			if (lL_post>bestL_post) { //Better post-transmission frequencies
				hap_freqs_post_store=hap_freqs_post;
				bestL_post=lL_post;
				fails=0;
			} else {
				hap_freqs_post=hap_freqs_post_store;
				fails++;
			}
		}
		if (fails==p.freq_term) {
			lL_pre=bestL_pre;
			lL_post=bestL_post;
			break;
		}
		first=0;
		
		//Change frequencies
		ChangeFreq (changex,full_haps,hap_freqs_pre,rgen);
		ChangeFreq (changex,full_haps,hap_freqs_post,rgen);
		
		//Calculate likelihood - before and after
		lL_pre=0;
		lL_post=0;
		for (int i=0;i<hap_data_sets.size();i++) { //Sets of partial haplotypes
			inf_freqs_pre.clear();
			inf_freqs_post.clear();
			for (int k=0;k<obs_freqs_pre[i].size();k++) { //Inference has same number of elements as observations
				inf_freqs_pre.push_back(0);
			}
			for (int k=0;k<obs_freqs_post[i].size();k++) {
				inf_freqs_post.push_back(0);
			}
			int xmatch=0;
			int indexj=0;
			for (int j=0;j<hap_data_sets[i].size();j++) { //Partial haplotypes in subset
				if (hap_data_sets[i][j].match[0]==full_haps.size()-1) { //Corresponds to unobserved haplotype namely X
					if (xmatch==0) { //Only need to do this once
						inf_freqs_pre[inf_freqs_pre.size()-1]=hap_freqs_pre[hap_data_sets[i][j].match[0]];
						inf_freqs_post[inf_freqs_pre.size()-1]=hap_freqs_post[hap_data_sets[i][j].match[0]];
						xmatch=1;
					}
				} else { //Corresponds to observed haplotype - find sum of frequencies of matching haplotypes
					for (int k=0;k<hap_data_sets[i][j].match.size();k++) {
						inf_freqs_pre[indexj]=inf_freqs_pre[indexj]+hap_freqs_pre[hap_data_sets[i][j].match[k]];
						inf_freqs_post[indexj]=inf_freqs_post[indexj]+hap_freqs_post[hap_data_sets[i][j].match[k]];
					}
					indexj++;
				}
			}
			
			//NB: Sometimes there is no match to a full haplotype.  In this case the inferred frequencies will not sum to 1.  Deal with this by putting them in the noise frequency, then correcting if this is >1%.  This seems to work out as a proportionate penalty term though perhaps there is a better solution.
			
			MissingFreqRedist(inf_freqs_pre);
			MissingFreqRedist(inf_freqs_post);
			
			L=DirichletMultiCalc(p,Npre[i],C,obs_freqs_pre[i],inf_freqs_pre,fact_store);
			lL_pre=lL_pre+L;
			L=DirichletMultiCalc(p,Npost[i],C,obs_freqs_post[i],inf_freqs_post,fact_store);
			lL_post=lL_post+L;
		}
	}
	lL_pre=bestL_pre;
	lL_post=bestL_post;
}

void MissingFreqRedist(vector<double>& inf_freqs) {
	//Move missing frequencies to noise term
	double tot=0;
	for (int j=0;j<inf_freqs.size();j++) {
		tot=tot+inf_freqs[j];
	}
	tot=1-tot;
	inf_freqs[inf_freqs.size()-1]=inf_freqs[inf_freqs.size()-1]+tot;
	//Now do correction: Noise term can't be more than 1%.
	if (inf_freqs[inf_freqs.size()-1]>0.01) {
		//Distribute frequencies evenly between other sites if so
		double diff=inf_freqs[inf_freqs.size()-1]-0.01;
		for (int j=0;j<inf_freqs.size()-1;j++) {
			inf_freqs[j]=inf_freqs[j]+(diff/(inf_freqs.size()-1));
		}
		inf_freqs[inf_freqs.size()-1]=0.01;
	}
}

void FindBestRep (double& lL_tot, vector<double>& hap_freqs_pre, vector<double>& hap_freqs_post, vector<double>& logs_store, vector< vector<double> >& pre_store, vector< vector<double> >& post_store) {
	double maxL=-1e10;
	int index=-1;
	for (int i=0;i<logs_store.size();i++) {
		if (logs_store[i]>maxL) {
			index=i;
			maxL=logs_store[i];
		}
	}
	lL_tot=logs_store[index];
	hap_freqs_pre=pre_store[index];
	hap_freqs_post=post_store[index];
}

void OptimiseFrequenciesFixedFullHaps (run_params p, double C, vector<double>& pre_freqs, vector<double>& post_freqs, vector<haplo>& full_haps, vector<int>& Npre, vector<int>& Npost, vector<double>& fact_store, vector< vector<mhap> >& hap_data_sets, gsl_rng *rgen) {
	//Full optimisation for case in which we read in a fixed set of haplotypes
	MatchFullHapsSplit(hap_data_sets,full_haps);
	vector< vector<double> > pre_store;
	vector< vector<double> > post_store;
	vector<double> logs_store;
	double lL_tot=-1e10;
	vector<double> hap_freqs_pre;
	vector<double> hap_freqs_post;
	for (int rep=0;rep<p.freq_rep;rep++) { //Independent frequency starting points
		//Set initial frequencies
		hap_freqs_pre.clear();
		hap_freqs_post.clear();
		SetInitialFreqs(hap_freqs_pre,full_haps,rgen);
		SetInitialFreqs(hap_freqs_post,full_haps,rgen);
		//Check initial frequencies - last cannot be >1% if it is an X-type haplotype
		CheckNoiseFreq(hap_freqs_pre);
		CheckNoiseFreq(hap_freqs_post);
		//Match observations to full haplotype frequencies.  Here the last obs combines all of the 'noise' observations
		vector< vector<int> > obs_freqs_pre;
		vector< vector<int> > obs_freqs_post;
		ConstructObservationsX (hap_data_sets,full_haps,obs_freqs_pre,obs_freqs_post);
		//Optimise frequencies to fit to data
		double lL_pre=-1e10;
		double lL_post=-1e10;
		OptimiseFrequencies (p,lL_pre,lL_post,C,Npre,Npost,full_haps,hap_freqs_pre,hap_freqs_post,obs_freqs_pre,obs_freqs_post,hap_data_sets,fact_store,rgen);
		lL_tot=lL_pre+lL_post;
		logs_store.push_back(lL_tot);
		pre_store.push_back(hap_freqs_pre);
		post_store.push_back(hap_freqs_post);
	}
	//Find best replicate
	FindBestRep (lL_tot,hap_freqs_pre,hap_freqs_post,logs_store,pre_store,post_store);
	double bic=FindBIC(lL_tot,full_haps);
	if (p.verb==1) {
		cout << "Final BIC " << bic << "\n";
	}
	if (p.suppress_output==0) {
		PrintInference(full_haps,hap_freqs_pre,hap_freqs_post);
	}
	pre_freqs=hap_freqs_pre;
	post_freqs=hap_freqs_post;
}

int CheckTime (run_params p, time_t timer_s) {
	//Used for case in which we impose a time limit on the optimisation
	int cont=1;
	time_t timer_n;
	time(&timer_n);
	double run_time=difftime(timer_n,timer_s);
	if (run_time>p.run_time) {
		cout << "Terminated after " << p.run_time << " seconds\n";
		cont=0;
	}
	return cont;
}


double DirichletMultiCalc (run_params p, int N, double c, vector<int>& obs, vector<double>& inf, vector<double>& fact_store) {
	double bin=0;
	//Correct inference if required: No zeros or negative frequencies allowed
	for (int i=0;i<inf.size();i++) {
		if (inf[i]<p.null_freq) {
			inf[i]=p.null_freq;
		}
	}
	if (p.null_freq>1e-20) {
		double tot=0;
		for (int i=0;i<inf.size();i++) {
			tot=tot+inf[i];
		}
		for (int i=0;i<inf.size();i++) {
			inf[i]=inf[i]/tot;
		}
	}
	/*cout << "Fact_store size: " << fact_store.size() << "\n";
	cout << "N : " << N << "\n";*/
	/*cout << "Obs\n";
	 for (unsigned int i=0;i<obs.size();i++) {
	 cout << obs[i] << " ";
	 }
	 cout << "\n";
	 cout << "Inf\n";
	 for (unsigned int i=0;i<inf.size();i++) {
	 cout << inf[i] << " ";
	 }
	 cout << "\n";
	 if (obs.size()!=inf.size()) {
	 cout << "Size match error " << obs.size() << " " << inf.size() << "\n";
	 }*/
	if (N>0) {
		bin=fact_store[N];
		for (unsigned int i=0;i<obs.size();i++) {
			bin=bin-fact_store[obs[i]];
		}
		vector<double> alpha;
		for (unsigned int i=0;i<obs.size();i++) {
			alpha.push_back(c*inf[i]);
		}
		double a=0;
		for (unsigned int i=0;i<alpha.size();i++) {
			a=a+alpha[i];
			bin=bin-gsl_sf_lngamma(alpha[i]);
		}
		bin=bin+gsl_sf_lngamma(a);
		a=0;
		for (unsigned int i=0;i<alpha.size();i++) {
			double b=alpha[i]+obs[i];
			a=a+b;
			bin=bin+gsl_sf_lngamma(b);
		}
		bin=bin-gsl_sf_lngamma(a);
	}
//		cout << "L " << bin << "\n";
//		cout << "\n";
	return(bin);
}

double FindBIC(double L, vector<haplo>& full_haps) {
	int params=full_haps.size()*2;
	double bic=-2*L + log(params);
	return bic;
}

void PrintInference (vector<haplo>& full_haps, vector<double>& hap_freqs_pre, vector<double>& hap_freqs_post) {
	for (int i=0;i<full_haps.size();i++) {
		for (int j=0;j<full_haps[i].seq.size();j++) {
			cout << full_haps[i].seq[j];
		}
		cout << " " << hap_freqs_pre[i] << " " << hap_freqs_post[i] << "\n";
	}
}

void FinalOutput (run_params p, double lL_hapbest, vector<haplo>& haplotypes, vector<double>& hap_freqs_pre_best, vector<double>& hap_freqs_post_best, vector< vector<haplo> >& full_haps_store) {
	vector<haplo> full_haps=full_haps_store[full_haps_store.size()-1];
	vector<double> hap_freqs_pre=hap_freqs_pre_best;
	vector<double> hap_freqs_post=hap_freqs_post_best;
	double bic=FindBIC(lL_hapbest,full_haps);
	ofstream out_file;
	string file="Inference_"+to_string(p.n_haps)+"_"+to_string(p.seed)+".out";
	out_file.open(file.c_str());
	out_file << "Final BIC " << bic << "\n";
	haplotypes=full_haps;
	for (int i=0;i<full_haps.size();i++) {
		for (int j=0;j<full_haps[i].seq.size();j++) {
			out_file << full_haps[i].seq[j];
		}
		out_file << " " << hap_freqs_pre[i] << " " << hap_freqs_post[i] << "\n";
	}
	
}

int ReconstructHaplotypesFixedFromInternal (run_params& p, double& best_bic, vector<haplo>& haplotypes, vector<double>& pre_freqs, vector<double>& post_freqs, vector<mhap>& hap_data, vector<haplo>& full_haps)  {
	//This is a cut-down version of the main reconstruction code, used for reading in data from an internal source rather than an external file.  It implements only the fixed-haplotype inference.  Used in the variance calculation
	double C=p.c;
	
	//Initialise random number generator
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, p.seed);
	
	//Calculate full set of loci
	vector<int> loci;
	CalculateFullLocusSet(hap_data,loci);
	sort(loci.begin(),loci.end());
	
	//Update data with n_loci values: Loci numbers become 0, 1, 2, 3, etc.
	CalculateNLociData(hap_data,loci);
	
	//Update data with set e.g. which cover which loci
	int n_sets;
	CalculateSetData(n_sets,hap_data);
	
	//Calculate N values
	vector<int> Npre;
	vector<int> Npost;
	CalculateObsTotalN (n_sets,Npre,Npost,hap_data);
	
	vector<double> fact_store;
	FindLogFact(fact_store,Npre,Npost);
	
	//Split partial haplotype data according to which loci are spanned by each read
	vector< vector<mhap> > hap_data_sets;
	SplitHaplotypeData (n_sets,hap_data,hap_data_sets);
	
	//Add an X haplotype if none such exists
	CheckXHap (loci,full_haps);
	p.suppress_output=1;
	//Optimise full haplotype frequencies
	OptimiseFrequenciesFixedFullHaps (p,C,pre_freqs,post_freqs,full_haps,Npre,Npost,fact_store,hap_data_sets,rgen);
	haplotypes=full_haps;
	return 0;
}
