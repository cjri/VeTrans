#include "shared_haps.h"
#include "io.h"
#include "sampling.h"
#include <sstream>

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.seed=(int) time(NULL);
	p.c=200;
	p.n_haps=1;
    p.n_samples=2;
	p.hap_its=1000;
	p.hap_term=100;
	p.freq_its=10000;
	p.freq_term=1000;
	p.freq_rep=1;
	p.run_time=1000000;
	p.read_hap=0;
	p.err=0;
	p.sample_reps=100;
	p.suppress_output=0;
	p.exp_like=0;
	p.growth=22;
	p.max_n=1000;
	p.delta_bic=10;
	p.extinct=0.005;
	p.find_max=0;
	p.print_like=1;
	p.null_freq=1e-20;
	p.var_name="Variance";
	p.traj_out="Multi_locus_trajectories_sample.out";
	p.traj_in="Multi_locus_trajectories.out";
	p.like_name="Likelihoods.out";
	p.fullhap_in="Inferred_haplotypes.out";
	p.verb=0;
	int x=2;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--n_haps")==0) {
			x++;
			p.n_haps=atoi(argv[x]);
        } else if (p_switch.compare("--n_samples")==0) {
            x++;
            p.n_samples=atoi(argv[x]);
        } else if (p_switch.compare("--hap_its")==0) {
			x++;
			p.n_haps=atoi(argv[x]);
		} else if (p_switch.compare("--hap_term")==0) {
			x++;
			p.hap_term=atoi(argv[x]);
		} else if (p_switch.compare("--freq_its")==0) {
			x++;
			p.freq_its=atoi(argv[x]);
		} else if (p_switch.compare("--freq_term")==0) {
			x++;
			p.freq_term=atoi(argv[x]);
		} else if (p_switch.compare("--freq_rep")==0) {
			x++;
			p.freq_rep=atoi(argv[x]);
		} else if (p_switch.compare("--verb")==0) {
			x++;
			p.verb=atoi(argv[x]);
		} else if (p_switch.compare("--explicit")==0) {
			x++;
			p.exp_like=atoi(argv[x]);
            p.like_name="Likelihoods_exp.out";
		} else if (p_switch.compare("--reps")==0) {
			x++;
			p.sample_reps=atoi(argv[x]);
		} else if (p_switch.compare("--growth")==0) {
			x++;
			p.growth=atoi(argv[x]);
		} else if (p_switch.compare("--extinct")==0) {
			x++;
			p.extinct=atof(argv[x]);
		} else if (p_switch.compare("--max_n")==0) {
			x++;
			p.max_n=atoi(argv[x]);
        } else if (p_switch.compare("--seed")==0) {
            x++;
            p.seed=atoi(argv[x]);
		} else if (p_switch.compare("--find_max")==0) {
			x++;
			p.find_max=atoi(argv[x]);
		} else if (p_switch.compare("--print_like")==0) {
			x++;
			p.print_like=atoi(argv[x]);
		} else if (p_switch.compare("--read_hap")==0) {
			x++;
			p.read_hap=atoi(argv[x]);
			if (p.run_time!=1000000) {
				cout << "Warning: Setting of time limit disabled for frequency-only inference\n";
				p.run_time=1000000;
			}
		} else if (p_switch.compare("--haps_in")==0) {
			x++;
			p.fullhap_in=argv[x];
		} else if (p_switch.compare("--traj_in")==0) {
			x++;
			p.traj_in=argv[x];
		} else if (p_switch.compare("--traj_out")==0) {
			x++;
			p.traj_out=argv[x];
		} else if (p_switch.compare("--var_name")==0) {
			x++;
			p.var_name=argv[x];
		} else if (p_switch.compare("--like_name")==0) {
			x++;
			p.like_name=argv[x];
		} else if (p_switch.compare("--run_time")==0) {
			x++;
			if (p.read_hap==1) {
				cout << "Warning: Setting of time limit disabled for frequency-only inference\n";
			} else {
				p.run_time=atoi(argv[x]);
				cout << "Set time limit to " << p.run_time << " seconds\n";
				p.hap_its=1000000;
			}
		} else if (p_switch.compare("--C")==0) {
			x++;
			p.c=atof(argv[x]);
		} else {
			p.err=1;
			p.err_flag=argv[x];
		}
		x++;
	}
}


void GetVariantData (run_params p, vector<mhap>& haps) {
	string file=p.traj_in;
	//string file="../Data/"+gene+"/Multi_locus_trajectories.out";
	if (p.verb==1) {
		cout << "Read partial haplotype data from " << p.traj_in << "\n";
	}
	ifstream in_file;
	in_file.open(file.c_str());
	int n;
	string s;
	for (int i=0;i<1000000;i++) {
		mhap m;
		if (!(in_file >> n)) break;
		if (p.verb==1) {
			cout << n << " ";
		}
		m.nvar=n;
		int k=n;
		for (int j=0;j<k;j++) {
			if (!(in_file >> n)) break;
			m.loci.push_back(n);
			if (p.verb==1) {
				cout << n << " ";
			}
		}
		if (!(in_file >> s)) break;
		m.st=s;
		if (p.verb==1) {
			cout << s << " ";
		}
		if (!(in_file >> n)) break;
		if (p.verb==1) {
			cout << n << " ";
		}

		int t=n;
		int tot=0;
		for (int j=0;j<t;j++) {
			//Time of sample
			if (!(in_file >> n)) break;
			m.times.push_back(n);
			if (p.verb==1) {
				cout << n << " ";
			}

			//Observations in sample
			if (!(in_file >> n)) break;
			m.obs.push_back(n);
			if (p.verb==1) {
				cout << n << " ";
			}

			tot=tot+n;
		}
		if (p.verb==1) {
			cout << "\n";
		}

		m.tot=tot;
		haps.push_back(m);
	}
}

void GetVariantDataMulti (run_params p, vector<mhap>& haps) {
    string file=p.traj_in;
    //string file="../Data/"+gene+"/Multi_locus_trajectories.out";
    if (p.verb==1) {
        cout << "Read partial haplotype data from " << p.traj_in << "\n";
    }
    ifstream in_file;
    in_file.open(file.c_str());
    int n;
    string s;
    for (int i=0;i<1000000;i++) {
        mhap m;
        if (!(in_file >> n)) break;
        if (p.verb==1) {
            cout << n << " ";
        }
        m.nvar=n;
        int k=n;
        for (int j=0;j<k;j++) {
            if (!(in_file >> n)) break;
            m.loci.push_back(n);
            if (p.verb==1) {
                cout << n << " ";
            }
        }
        if (!(in_file >> s)) break;
        m.st=s;
        if (p.verb==1) {
            cout << s << " ";
        }
        if (!(in_file >> n)) break;
        if (p.verb==1) {
            cout << n << " ";
        }

        int t=n;
        //cout << "T is " << t << "\n";
        int tot=0;
        for (int j=0;j<t;j++) {
            //Time of sample
            if (!(in_file >> n)) break;
            m.times.push_back(n);
            if (p.verb==1) {
                cout << "T " << n << " ";
            }

            //Observations in sample
            if (!(in_file >> n)) break;
            m.obs.push_back(n);
            if (p.verb==1) {
                cout << "N " << n << " ";
            }

            tot=tot+n;
        }
        if (p.verb==1) {
            cout << "\n";
        }

        m.tot=tot;
        haps.push_back(m);
    }
}


void GetFullHaplotypes (run_params p, vector<haplo>& full_haps) {
	ifstream hap_file;
	hap_file.open(p.fullhap_in);
	string fullhap;
	while (hap_file >> fullhap) {
		//cout << fullhap << "\n";
		haplo h;
		h.st=fullhap;
		for (int i=0;i<fullhap.length();i++) {
			h.seq.push_back(fullhap[i]);
		}
		full_haps.push_back(h);
		hap_file.ignore(numeric_limits<streamsize>::max(), '\n');
	}
}

void GetFullHaplotypesFreq (run_params p, vector<haplo>& full_haps, vector<double>& freq_pre, vector<double>& freq_post) {
	ifstream hap_file;
	//cout << "Hap in " << p.hap_in << "\n";
	hap_file.open(p.fullhap_in);
    if (p.verb==1) {
        cout << "Read from file " << p.fullhap_in << "\n";
    }
	string fullhap;
	double freq1;
	double freq2;
	for (int i=0;i<1000000;i++) {
		if (!(hap_file >> fullhap)) break;
		haplo h;
		h.st=fullhap;
        if (p.verb==1) {
            cout << h.st << " ";
        }
		for (int i=0;i<fullhap.length();i++) {
			h.seq.push_back(fullhap[i]);
		}
		full_haps.push_back(h);
		if (!(hap_file >> freq1)) break;
		if (!(hap_file >> freq2)) break;
		freq_pre.push_back(freq1);
		freq_post.push_back(freq2);
        if (p.verb==1) {
            cout << freq1 << " " << freq2 << "\n";
        }

	}
	if (p.verb==1) {
		for (int i=0;i<full_haps.size();i++) {
			cout << full_haps[i].st << " " << freq_pre[i] << " " << freq_post[i] << "\n";
		}
	}
}

void PrintResampleMulti (vector< vector< vector<mhap> > >& all_hap_data_sets) {
	for (int i=0;i<all_hap_data_sets.size();i++) {
		vector<mhap> hap_data;
		FlattenData(i,all_hap_data_sets,hap_data);
		PrintSample(i,hap_data);
	}
}

void PrintSample (int i, vector<mhap>& hap_data) {
	ofstream mhap_file;
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Multi_locus_trajectories_sample"+temp+".out";
	mhap_file.open(name.c_str());
	for (int j=0;j<hap_data.size();j++) {
		mhap_file << hap_data[j].loci.size() << " ";
		for (int k=0;k<hap_data[j].loci.size();k++) {
			mhap_file << hap_data[j].loci[k] << " ";
		}
		mhap_file << hap_data[j].st << " ";
		mhap_file << hap_data[j].times.size() << " ";
		for (int k=0;k<hap_data[j].times.size();k++) {
			mhap_file << hap_data[j].times[k] << " " << hap_data[j].obs[k] << " ";
		}
		mhap_file << "\n";
	}
	mhap_file.close();
}

void PrintVariances (run_params p, vector<double>& pre_var, vector<double>& post_var) {
	ofstream pre_var_file;
	ofstream post_var_file;
	string name_pre=string(p.var_name) + "_pre.out";
	string name_post=string(p.var_name) + "_post.out";
	pre_var_file.open(name_pre.c_str());
	post_var_file.open(name_post.c_str());
	for (int i=0;i<pre_var.size();i++) {
		for (int j=0;j<pre_var.size();j++) {
			if (i==j) {
				pre_var_file << pre_var[i] << " ";
			} else {
				pre_var_file << " 0 ";
			}
		}
		pre_var_file << "\n";
	}
	for (int i=0;i<post_var.size();i++) {
		for (int j=0;j<post_var.size();j++) {
			if (i==j) {
				post_var_file << post_var[i] << " ";
			} else {
				post_var_file << " 0 ";
			}
		}
		post_var_file << "\n";
	}
}

void GetLikelihoods (run_params p, vector< vector<double> >& likelihoods) {
	ifstream gene_file;
	gene_file.open("Gene_list.in");
	vector<string> genes;
	string s;
	for (int i=0;i<1000000;i++) {
		if (!(gene_file >> s)) break;
		genes.push_back(s);
	}
	ifstream like_file;
	for (int i=0;i<genes.size();i++) {
		string name=genes[i]+"/"+string(p.like_name);
        if (p.verb==1) {
            cout << "Reading file " << name << "\n";
        }
 		like_file.open(name.c_str());
		double x;
		int n;
		vector<double> like;
		for (int j=0;j<1000000;j++) {
			if (!(like_file >> n)) break;
			if (!(like_file >> x)) break;
			like.push_back(x);
		}
		likelihoods.push_back(like);
		like_file.close();
	}
}

void GetVarianceMatrices(run_params p, int dim, vector< vector<double> >& var_pre_provis, vector< vector<double> >& var_post_provis) {
	ifstream vpre_file;
	ifstream vpost_file;
	string name_pre=string(p.var_name) + "_pre.out";
	string name_post=string(p.var_name) + "_post.out";
	vpre_file.open(name_pre.c_str());
	vpost_file.open(name_post.c_str());
	double x;
	for (int i=0;i<dim;i++) {
		vector<double> vp;
		for (int j=0;j<dim;j++) {
			vpre_file >> x;
			vp.push_back(x);
		}
		var_pre_provis.push_back(vp);
	}
	for (int i=0;i<dim;i++) {
		vector<double> vp;
		for (int j=0;j<dim;j++) {
			vpost_file >> x;
			vp.push_back(x);
		}
		var_post_provis.push_back(vp);
	}
	vpre_file.close();
	vpost_file.close();
	if (p.verb==1) {
		cout << "Read in pre-transmission variance matrix:\n";
		for (int i=0;i<dim;i++) {
			for (int j=0;j<dim;j++) {
				cout << var_pre_provis[i][j] << " ";
			}
			cout << "\n";
		}
		cout << "Read in post-transmission variance matrix:\n";
		for (int i=0;i<dim;i++) {
			for (int j=0;j<dim;j++) {
				cout << var_post_provis[i][j] << " ";
			}
			cout << "\n";
		}
	}
}

void PrintReducedData(run_params p, int dim, vector<double> freq_pre, vector<double> freq_post, gsl_matrix *var_pre, gsl_matrix *var_post) {
	cout << "Reduced frequencies\n";
	for(int i=0;i<freq_pre.size();i++) {
		cout << freq_pre[i] << "      " << freq_post[i] << "\n";
	}
	if (p.exp_like==1) {
		cout << "Reduced variance Pre transmission:\n";
		for(int i=0;i<dim;i++) {
			for(int j=0;j<dim;j++) {
				cout << gsl_matrix_get(var_pre,i,j) << " ";
			}
			cout << "\n";
		}
		cout << "Reduced variance Post transmission:\n";
		for(int i=0;i<dim;i++) {
			for(int j=0;j<dim;j++) {
				cout << gsl_matrix_get(var_post,i,j) << " ";
			}
			cout << "\n";
		}
	}
	
}

void OutputReconstructHaplotypes (vector<haplo> haplotypes,vector<double> pre_freqs,vector<double> post_freqs,double best_bic) {
	cout << "Result of calculation\n";
	cout << "BIC " << best_bic << "\n";
	for (int i=0;i<haplotypes.size();i++) {
		for (int j=0;j<haplotypes[i].seq.size();j++) {
			cout << haplotypes[i].seq[j];
		}
		cout << " " << pre_freqs[i] << " " << post_freqs[i] << "\n";
	}
}

void OutputReconstructHaplotypesMulti (run_params p, vector<haplo> haplotypes, vector< vector<double> >& multi_freqs, double best_bic) {
    cout << "Result of calculation\n";
    cout << "BIC " << best_bic << "\n";
    for (int i=0;i<haplotypes.size();i++) {
        for (int j=0;j<haplotypes[i].seq.size();j++) {
            cout << haplotypes[i].seq[j];
        }
        for (int s=0;s<p.n_samples;s++) {
            cout << " " << multi_freqs[s][i];
        }
        cout << "\n";
    }
}

