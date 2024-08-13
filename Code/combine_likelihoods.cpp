#include "shared_haps.h"
#include "combine_likelihoods.h"
#include "io.h"

void CombineLikelihoods(run_params p) {
	if (p.verb==1) {
		cout << "Combine likelihoods from different genes...\n";
	}
	vector< vector<double> > likelihoods;
	GetLikelihoods(p,likelihoods);
	FindCombinedLikelihood(p,likelihoods);
}


void FindCombinedLikelihood (run_params p, vector< vector<double> >& likelihoods) {
	ofstream l_file;
	if (p.exp_like==1) {
		l_file.open("Combined_likelihoods_exp.out");
	} else {
		l_file.open("Combined_likelihoods_cpd.out");
	}
	double maxL=-1e80;
	int maxN=-1;
    cout << likelihoods.size() << "\n" << likelihoods[0].size() << "\n";
	for (int i=0;i<likelihoods[0].size();i++) {
		double tot=0;
		for (int j=0;j<likelihoods.size();j++) {
			tot=tot+likelihoods[j][i];
		}
		l_file << i+1 << " " << tot << "\n";
		if (tot>maxL) {
			maxL=tot;
			maxN=i+1;
		}
	}
	cout << "Maximum likelihood bottleneck = " << maxN << "\n";
}
