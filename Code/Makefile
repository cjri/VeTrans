CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I  /usr/local/include/ -std=c++0x
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lm
VTR		= vetrans reconstruct_haplotypes io sampling sample_haplotypes combine_likelihoods explicit compound calc_bottleneck
OBJ = $(addsuffix .o, $(VTR))
TARGETS	= VeTrans

$(OBJ): %.o: %.cpp
	$(CC) $(CC_FLAGS) -c $<
$(TARGETS): $(OBJ)
	$(CC) $(CC_FLAGS) $(OBJ) -o $(TARGETS) $(LD_FLAGS)


