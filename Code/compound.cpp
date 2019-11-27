//Functions for compound likelihood
#include "shared_haps.h"
#include "compound.h"
#include "explicit.h"
#include "calc_bottleneck.h"

void constructMatrixM(vector<double> x, gsl_matrix *M) {
	
	//Create subset of vec
	int dim = (int) x.size();
	
	//Hardcode the very simple cases as it gives a considerable speed-up. See mMatrixChecker.cpp.
	if(dim==1) { //Very simple case, if x={x1} then M = {{ x1-x1^2 }}, hardcode for efficiency. About 1.5 times as fast
		
		double value = x[0] - x[0]*x[0];
		gsl_matrix *mat = gsl_matrix_alloc(dim,dim);
		gsl_matrix_set(mat,0,0,value);
		gsl_matrix_memcpy(M,mat); //Copy mat into M, entry by entry
		gsl_matrix_free(mat); //Clean up
		return;
		
	} else if(dim==2) { //Also simple case. Here M[{x1,x2}] = {{x1 - x1^2, -x1 x2}, {-x1 x2, x2 - x2^2}}. About 1.75 times as fast.
		
		gsl_matrix *mat = gsl_matrix_alloc(dim,dim);
		gsl_matrix_set(mat,0,0,x[0]-x[0]*x[0]);
		gsl_matrix_set(mat,1,0,-x[0]*x[1]);
		gsl_matrix_set(mat,0,1,-x[0]*x[1]);
		gsl_matrix_set(mat,1,1,x[1]-x[1]*x[1]);
		gsl_matrix_memcpy(M,mat); //Copy mat into M, entry by entry
		gsl_matrix_free(mat); //Clean up
		return;
		
	}
	
	//Create diagonal matrix with elements of x
	gsl_matrix *diagX = gsl_matrix_alloc(dim,dim);
	for(int j=0; j<dim; j++) {
		for(int k=0; k<dim; k++) {
			
			if(j==k) {
				gsl_matrix_set(diagX,j,k,x[j]);
			} else {
				gsl_matrix_set(diagX,j,k,0);
			}
		}
	}
	
	/*
	 / Create outer product of x*x^transpose.
	 / This method does not exist in CBLAS library, but can be obtained otherwise.
	 / In particular, Outer(a,b) = Matrix(columns=a) * DiagonalMatrix(b).
	 */
	//Create matrices
	gsl_matrix *A = gsl_matrix_alloc(dim,dim);
	for(int j=0; j<dim; j++) {
		for(int k=0; k<dim; k++) {
			gsl_matrix_set(A,j,k,x[j]);
		}
	}
	
	//Multiply them and store in XX
	gsl_matrix *XX=gsl_matrix_alloc(dim,dim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
					1.0, A, diagX,
					0.0, XX);
	
	//Create M=diag(x) -xxp^t
	gsl_matrix_sub(diagX,XX); //diagX=diagX-XX, i.e. store result in diagPb
	
	//Clean up
	gsl_matrix_memcpy(M,diagX); //Copy to M, entry by entry
	gsl_matrix_free(diagX);
	gsl_matrix_free(A);
	gsl_matrix_free(XX);
}


//Use the fully continuous method for computing a likelihood
double computeLikelihoodCont(vector<double> x, vector<double> mean, gsl_matrix *var) {
	
	//Computing likelihoods depend on dimensionality
	//	double logLikelihood = -numeric_limits<double>::max(); //Set to neg infinity
	double logLikelihood = -1e80;
	double pi=3.14159265358979323846264338327950288;
	
	int dim = (int) mean.size();
	if(dim > 2) { //Dimension larger than 2
		// Reduce dimensionality by 1 to ensure non-degeneracy
		mean.pop_back(); //WLOG remove last haplotype
		x.pop_back(); //WLOG remove last haplotype
		gsl_matrix_view varReducedView = gsl_matrix_submatrix(var, 0,0, dim-1, dim-1); //Remove last row and column
		gsl_matrix* varReduced = gsl_matrix_alloc(dim-1, dim-1);
		gsl_matrix_memcpy(varReduced, &(varReducedView.matrix)); //Copy the reduced matrix view into a "normal" matrix varReduced
		
		//Compute likelihood
		logLikelihood = logMultivariateNormalPDF(x,mean,varReduced);
		
	} else { //i.e. dim==2 or dim==1 (same computation)
		
		
		//In the case where there is a single haplotype in a set (i.e. dim=1),
		//this can be thought of as a case with a single observed haplotype
		//and one (or several) unobserved haplotypes.
		//If all loci in the haplotype set are variants then there MUST
		//exist at least one other haplotype containing the unobserved alleles,
		//however, this (or these) haplotype(s) are not observed.
		//As such, WLOG we can consider this a case of dim=2, which gets reduced to a one dimensional system
		//under reduction, i.e. similar to the dim==2 case.
		double meanSingle = mean[0];
		double varSingle = gsl_matrix_get(var,0,0);
		double stDevSingle = sqrt(varSingle);
		
		logLikelihood = -(x[0]-meanSingle)*(x[0]-meanSingle)/(2*varSingle) - log(stDevSingle*sqrt(2*pi)); //log L
		
	}
	
	return logLikelihood;
}

//Compute log likelihood for multivariate normal distribution in x with mean mu and variance sigma
double logMultivariateNormalPDF(vector<double> &x, vector<double> &mu, gsl_matrix *sigma) {
	double pi=3.14159265358979323846264338327950288;
	
	//Set up permutation and matrix copy
	int dim= (int) sigma->size1;
	vector<double> xMinusMu = subtractVectors(x,mu);
	gsl_permutation *p = gsl_permutation_calloc(dim);
	int signum;
	gsl_matrix * tmp_ptr = gsl_matrix_calloc(dim,dim);
	gsl_matrix_memcpy(tmp_ptr,sigma); //Copy elements of sigma into tmp_ptr
	
	//Get LU decomposition
	gsl_linalg_LU_decomp(tmp_ptr,p,&signum);
	
	//Get determinant
	double det = gsl_linalg_LU_det(tmp_ptr, signum);
	//	cout <<  "Determinant: " << det << "\n";
	//	if(det==0) { det=1e-10; }
	
	//Get inverse
	gsl_matrix *sigmaInv = gsl_matrix_alloc(dim,dim);
	gsl_set_error_handler_off();
	int status = gsl_linalg_LU_invert(tmp_ptr,p,sigmaInv);
	if (status) {
		
		//Clean up, then return negative infinity
		gsl_matrix_free(sigmaInv);
		gsl_matrix_free(tmp_ptr);
		gsl_permutation_free(p);
		
		cout << "Matrix not positive definite. Returning probability of neg infinity.\n";
		//	cout << "Determinant should be 0 in this case. Determinant was: " << det << "\n";
		double L=-1e80;
		return L;
	}
	
	//double logPrefactor= -0.5*log(pow(2*M_PI,dim)*det);
	double logPrefactor= -0.5*dim*log(2*pi) -0.5*log(det);
	
	//Convert xMinusMu to a gsl_vector
	gsl_vector *xMinusMuGSL = gsl_vector_alloc(dim);
	for(unsigned int i=0;i<xMinusMu.size();i++) {
		gsl_vector_set(xMinusMuGSL,i,xMinusMu[i]);
	}
	
	//Perform matrix*vector multiplication
	gsl_vector *sigmaInvXminusMu = gsl_vector_alloc(dim);
	gsl_blas_dgemv(CblasNoTrans,1.0,sigmaInv,xMinusMuGSL,0,sigmaInvXminusMu);
	
	//Perform vector*vector multiplication
	double dotProd;
	gsl_blas_ddot(xMinusMuGSL,sigmaInvXminusMu,&dotProd);
	
	//Clean up
	gsl_vector_free(xMinusMuGSL);
	gsl_vector_free(sigmaInvXminusMu);
	gsl_matrix_free(sigmaInv);
	gsl_matrix_free(tmp_ptr);
	gsl_permutation_free(p);
	
	return logPrefactor -0.5*dotProd;
}

vector<double> subtractVectors(vector<double> &a, vector<double> &b) {
	
	vector<double> result;
	//Assume equal length
	for(unsigned int i=0; i<a.size(); i++) {
		result.push_back(a[i]-b[i]);
	}
	return result;
}

double FindCompoundLikelihood (run_params p, int Nt, int dim, vector<double>& freq_pre, vector<double>& freq_post, gsl_matrix *var_pre, gsl_matrix *var_post) {
	//Create factors gamma and delta
	double Ng = Nt*p.growth; //22 fold growth
	double NtNg = Nt*Ng;
	double gamma = (Nt+Ng-1)/NtNg;
	double delta = (NtNg-Nt-Ng+1)/NtNg;
	
	//Create matrix gamma*M(q_post)  where M is the outer product
	gsl_matrix * gammaMqStarB = gsl_matrix_alloc(dim,dim);
	constructMatrixM(freq_post, gammaMqStarB);
	gsl_matrix_scale(gammaMqStarB, gamma);
	
	//Create delta*var_post
	gsl_matrix * deltaVarB = gsl_matrix_alloc(dim,dim);
	gsl_matrix_memcpy(deltaVarB, var_post);
	gsl_matrix_scale(deltaVarB, delta);
	
	//construct the variance of A: Goes into the var matrix
	gsl_matrix * var = gsl_matrix_alloc(dim,dim);
	gsl_matrix_memcpy(var, var_pre);  //Copy var_pre into var
	gsl_matrix_add(var, gammaMqStarB);  //Add gammaM(q_post)
	gsl_matrix_add(var, deltaVarB);  //Add delta*var_post
	
	if (p.verb==1) {
		cout << "Compound variance:\n";
		for(int i=0;i<dim;i++) {
			for(int j=0;j<dim;j++) {
				cout << gsl_matrix_get(var,i,j) << " ";
			}
			cout << "\n";
		}
	}
	
	//Free up unused variables
	gsl_matrix_free(gammaMqStarB);
	gsl_matrix_free(deltaVarB);
	
	//Identify haplotype extinction events: Not seen post-transmission
	int n = 0; //Number of extinctions
	vector<int> extinctions;
	for(int i=dim-1;i>=0;i--) { //Get extinction entries in descending order, e.g. {5,3,1}
		if(freq_post[i] <= p.extinct) {  //Below threshold qualifies as extinction
			//cout << i << " " << freq_post[i] << "\n";
			extinctions.push_back(i);
			n++;
		}
	}
	if (p.verb==1) {
		cout << "Dimension " << dim << "\n";
		cout << "Number of extinctions = " << n << "\n";
	}
	if (Nt<dim-n) { //Bottleneck is less than the number of surviving haplotypes
		double L=-1e80;
		return L;
	} else if (n == 0) {  //No extinction events
		//Compute likelihood the fully continuous (multivariate normal) way
		double L = computeLikelihoodCont(freq_post,freq_pre,var);
		return L;
	} else { //Some extinction - use compound likelihood method
		//Start by finding swaps
		vector<vector<int> > swaps = findSwaps(extinctions, dim);
		if (p.verb==1) {
			cout << "Swaps:\n";
			for (int i=0;i<swaps.size();i++) {
				cout << swaps[i][0] << " " << swaps[i][1] << "\n";
			}
		}
		//Then update x, mean and variances: N.B. This code alters the original variables
		swapVector(freq_post, swaps);
		swapVector(freq_pre, swaps);
		swapMatrix(var, swaps);
		if (p.verb==1) {
			cout << "Swapped frequencies\n";
			for(int i=0;i<freq_pre.size();i++) {
				cout << freq_pre[i] << "      " << freq_post[i] << "\n";
			}
			
			cout << "Swapped compound variance:\n";
			for(int i=0;i<dim;i++) {
				for(int j=0;j<dim;j++) {
					cout << gsl_matrix_get(var,i,j) << " ";
				}
				cout << "\n";
			}
		}
		
		// We now want to compute the likelihood L = P(surviving haplotypes | extinct haplotypes) * P(extinct haplotypes).
		// We note that the above is written in probability space, in log likelihood space the terms are summed, not multiplied.
		
		//First we create the conditional distribution for the surviving haplotypes conditional on the extinct haplotypes, see https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions.
		//Here we have a system of the form:
		// Dimensionality of dim = freq_post.size()
		// Number of extinctions n = numExtinctions
		// Inferred frequencies freq_post = q*A = [q1, q2] with q1 being of dimension (dim-n) x 1 and q2 being of dimensopm n x 1
		// Mean = [m1, m2] specified by freq_pre with m1 being of dimension (dim-n) x 1 and m2 being of dimension n x 1
		// var = {{v11, v12}, {v21, v22}} with v11 being of dimension (dim-n) x (dim-n), v12 with dimension (dim-n) x n, v21 with dimension n x (dim-n) and v22 with dimension n x n
		// Extinction outcome a = x[dim-n, dim] (i.e. basically {0,0,0,0...,0} for n entries. Note, this is not entirely true, a some entries can be e.g. 10^-8 and still count as extinction)
		//
		// Then the conditional distribution has mean
		// meanCond = m1 + v12 . (v22)^-1 . (a - m2)   (here we use '.' to indicate matrix multiplication)
		//
		// And it has variance
		// varCond = v11 - v12 . (v22)^-1 . v21
		
		//First create m1, m2, and a:
		vector<double> m1(freq_pre.begin(), freq_pre.begin() + (dim-n));
		vector<double> m2(freq_pre.begin() + (dim-n), freq_pre.end());
		vector<double> a(freq_post.begin() + (dim-n), freq_post.end());
		
		//Then create v11: (which we will call 'varCond' - the reason for this will become clear below)
		gsl_matrix_view v11View = gsl_matrix_submatrix(var, 0,0, (dim-n), (dim-n)); //Get sub matrix with upper-left element at (0,0), numRows=numColumns = (dim-n)
		gsl_matrix * varCond = gsl_matrix_alloc((dim-n),(dim-n));
		gsl_matrix_memcpy(varCond, &(v11View.matrix)); //Copy the v11 view into a "normal" matrix v11
		
		//Then create v12:
		gsl_matrix_view v12View = gsl_matrix_submatrix(var, 0,(dim-n), (dim-n), n); //Get submatrix with upper-left element at (0,(dim-n)), numRows=(dim-n), numCols=n
		gsl_matrix * v12 = gsl_matrix_alloc((dim-n),n);
		gsl_matrix_memcpy(v12, &(v12View.matrix)); //Copy the v12 view into a "normal" matrix v12
		
		//Then create v21:
		gsl_matrix_view v21View = gsl_matrix_submatrix(var, (dim-n), 0, n, (dim-n)); //Get submatrix with upper-left element at ((dim-n),0), numRows=n, numCols=(dim-n)
		gsl_matrix * v21 = gsl_matrix_alloc(n,(dim-n));
		gsl_matrix_memcpy(v21, &(v21View.matrix)); //Copy the v21 view into a "normal" matrix v21
		
		//Then create v22:
		gsl_matrix_view v22View = gsl_matrix_submatrix(var, (dim-n), (dim-n), n, n); //Get submatrix with upper-left element at ((dim-n),(dim-n)), numRows=numColumns = n
		gsl_matrix * v22 = gsl_matrix_alloc(n,n);
		gsl_matrix_memcpy(v22, &(v22View.matrix)); //Copy the v22 view into a "normal" matrix v22
		
		if (p.verb==1) {
			cout << "V11\n";
			for (int i=0;i<varCond->size1;i++) {
				for (int j=0;j<varCond->size2;j++) {
					cout << gsl_matrix_get(varCond,i,j) << " ";
				}
				cout << "\n";
			}
			cout << "V12\n";
			for (int i=0;i<v12->size1;i++) {
				for (int j=0;j<v12->size2;j++) {
					cout << gsl_matrix_get(v12,i,j) << " ";
				}
				cout << "\n";
			}
			cout << "V21\n";
			for (int i=0;i<v21->size1;i++) {
				for (int j=0;j<v21->size2;j++) {
					cout << gsl_matrix_get(v21,i,j) << " ";
				}
				cout << "\n";
			}
			cout << "V22\n";
			for (int i=0;i<v22->size1;i++) {
				for (int j=0;j<v22->size2;j++) {
					cout << gsl_matrix_get(v22,i,j) << " ";
				}
				cout << "\n";
			}
		}
		
		/////// Now we create the conditional mean
		//First create a gsl vector representing (a-m2)
		gsl_vector *aMinusM2 = gsl_vector_alloc(n);
		for(int i=0;i<n;i++) {
			gsl_vector_set(aMinusM2,i,a[i] - m2[i]);
		}
		
		//Then we want to create the inverse of v22:
		//First set up permutation and signum
		gsl_permutation *pm = gsl_permutation_alloc(n);
		int signum;
		
		//Then compute the LU decomposition
		gsl_linalg_LU_decomp(v22,pm,&signum); //Note that v22 is changed in this process!
		
		//Finally create inverse of v22
		gsl_matrix *v22Inv = gsl_matrix_alloc(n,n);
		gsl_set_error_handler_off(); //Stops the default error handler from exiting the program - instead we deal with the error ourselves
		int status = gsl_linalg_LU_invert(v22,pm,v22Inv); //Status of 0 equal success
		if (status) {
			
			//Error occurred, so clean up and exit
			cout << "V22 matrix not positive definite, so cannot invert. Exiting.\n";
			exit(1);
		}
		
		//Next we do the multiplications, starting with v12 . (v22)^1
		gsl_matrix* v12v22Inv = gsl_matrix_alloc(dim-n,n);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, v12, v22Inv, 0.0, v12v22Inv);
		
		
		//Followed by the matrix-vector multiplication (v12.(v22)^1) . (a-m2)
		//Actually, we can be a bit clever. The matrix-vector multiplication in gsl is as follows: y=alpha*A*x +beta*y.
		//In other words, if we pre-define a vector y, we can do both matrix-vector multiplication of A times x, but also add y.
		//This is exactly what we need to compute meanCond = m1 + v12 . (v22)^-1 . (a - m2), i.e. we can let y=m1, A=(v12.(v22)^-1) and x=(a-m2).
		gsl_vector *y = gsl_vector_alloc(dim-n);
		for(int i=0; i<(dim-n); i++) {
			gsl_vector_set(y,i,m1[i]);
		}
		gsl_blas_dgemv(CblasNoTrans, 1.0, v12v22Inv, aMinusM2, 1.0, y); //Result is now stored in y
		
		//Get conditional mean as 'normal' C++ vector
		vector<double> meanCond;
		for(int i=0; i<(dim-n); i++) {
			meanCond.push_back(gsl_vector_get(y,i));
		}
		
		//////// Now we create the conditional variance:  varCond = v11 - v12 . (v22)^-1 . v21
		//Here we can be a bit clever. The matrix-matrix multiplication in gsl is as follows: C=alpha*A*B+beta*C.
		//In other words, if we pre-define a matrix C we can do all the computation in one go. To this end, we
		//associate C=v11, alpha=-1, A=v12v22Inv, B=v21, beta=1.0. To avoid renaming matrices (through wasteful copying) we
		//will refer to C=v11 as 'varCond' (we already did this at the defintion of v11 above).
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,v12v22Inv,v21, 1.0, varCond);
		
		//Create conditional observations, x (i.e. q*A subsetted to match the mean)
		vector<double> xCond(freq_post.begin(), freq_post.begin()+(dim-n));
		
		//Compute the continuous likelihood
		double logLcont = computeLikelihoodCont(xCond, meanCond, varCond);
		
		//Compute discrete likelihood
		//Noticing that x_i=0 for all the extinct haplotypes, the binomial/multinomial likelihood simplifies to
		//P(q_1,q_2,...,q_n) = (1-q_1-q_2-...-q_n)^(Nt) for n extinct haplotypes. As such
		//L(q_1,q_2,...,q_n) = Nt log(1-q_1-q_2-..._q_n) (log likelihood)
		double sumFreqs = 0;
		for(int i = (dim-n); i<dim; i++) {
			sumFreqs += freq_pre[i];
		}
		double logLdisc = Nt*log(1-sumFreqs);
		double lL=logLdisc+logLcont;
		if (p.verb==1) {
			cout << "Nt: " << Nt << "\tlogLDisc: " << logLdisc << "\tlogLcont: " << logLcont << "\n";
		}
		return lL;
	}
	
}

void FindMaxLikelihood (run_params p, int dim, int& maxN, double& maxL,double C, vector<int>& Npost, vector<double>& freq_pre, vector<double>& freq_post, vector<double>& fact_store, vector< vector<int> >& list, vector< vector<int> >& obs_freqs_post, vector<haplo>& full_haps, vector< vector<mhap> >& hap_data_sets, gsl_matrix *var_pre, gsl_matrix *var_post) {
	int max=min(20,p.max_n);
	for (int Nt=1;Nt<=max;Nt++) {
		double L=0;
		if (p.exp_like==1) {
			L=CalculateExplicitLikelihood (p,Nt,C,Npost,freq_pre,freq_post,fact_store,list,obs_freqs_post,full_haps,hap_data_sets);
		} else {
			L=FindCompoundLikelihood (p,Nt,dim,freq_pre,freq_post,var_pre,var_post);
		}
		if (L>maxL) {
			maxL=L;
			maxN=Nt;
		}
	}
	if (maxN>1) {
		max=min(100,p.max_n);
		for (int Nt=21;Nt<=max;Nt++) {
			double L=0;
			if (p.exp_like==1) {
				L=CalculateExplicitLikelihood (p,Nt,C,Npost,freq_pre,freq_post,fact_store,list,obs_freqs_post,full_haps,hap_data_sets);
			} else {
				L=FindCompoundLikelihood (p,Nt,dim,freq_pre,freq_post,var_pre,var_post);
			}
			if (L>maxL) {
				maxL=L;
				maxN=Nt;
			}
		}
	}
	if (maxN>21) {
		max=min(1000,p.max_n);
		for (int Nt=101;Nt<=max;Nt++) {
			double L=0;
			if (p.exp_like==1) {
				L=CalculateExplicitLikelihood (p,Nt,C,Npost,freq_pre,freq_post,fact_store,list,obs_freqs_post,full_haps,hap_data_sets);
			} else {
				L=FindCompoundLikelihood (p,Nt,dim,freq_pre,freq_post,var_pre,var_post);
			}
			if (L>maxL) {
				maxL=L;
				maxN=Nt;
			}
		}
	}
}

