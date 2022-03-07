CC	      = g++
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ 
VTR		= vetrans.o reconstruct_haplotypes.o reconstruct_haplotypes_multi.o io.o sampling.o sample_haplotypes.o combine_likelihoods.o explicit.o compound.o calc_bottleneck.o

vtr: $(VTR)
	$(CC) $(CC_FLAGS) $(VTR) -o VeTrans $(LD_FLAGS)
vetrans.o: vetrans.cpp
	$(CC) $(CC_FLAGS) -c vetrans.cpp
reconstruct_haplotypes.o: reconstruct_haplotypes.cpp
	$(CC) $(CC_FLAGS) -c reconstruct_haplotypes.cpp
reconstruct_haplotypes_multi.o: reconstruct_haplotypes_multi.cpp
	$(CC) $(CC_FLAGS) -c reconstruct_haplotypes_multi.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
sampling.o: sampling.cpp
	$(CC) $(CC_FLAGS) -c sampling.cpp
sample_haplotypes.o: sample_haplotypes.cpp
	$(CC) $(CC_FLAGS) -c sample_haplotypes.cpp
explicit.o: explicit.cpp
	$(CC) $(CC_FLAGS) -c explicit.cpp
compound.o: compound.cpp
	$(CC) $(CC_FLAGS) -c compound.cpp
calc_bottleneck.o: calc_bottleneck.cpp
	$(CC) $(CC_FLAGS) -c calc_bottleneck.cpp
combine_likelihoods.o: combine_likelihoods.cpp
	$(CC) $(CC_FLAGS) -c combine_likelihoods.cpp


