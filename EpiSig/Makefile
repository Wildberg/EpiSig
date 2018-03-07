CC = gcc
CXX = g++
MPICXX = mpic++
CFLAGS = -O2 -Wall -Werror -Wno-unused-but-set-variable -Wno-unused-result
LFLAGS = -lm -pthread
MPIFLAGS = -lm -DCHROMASIG_MPI=1


OBJ_MT = ConfUtils.o DataUtils.o FileUtils.o utils.o genome_size.o pnorm.o thread_random.o pre_EpiSig.o filter.o filter_ad.o EpiSig_v1.0.o
OBJ_MPI = EpiSig_v1.0.MPI.o
OBJ_BI = global_data.o fast_pmf_cycle.o
OBJ_pre_EpiSig = ConfUtils.o FileUtils.o utils.o genome_size.o pre_EpiSig.o
#OBJ_filter = ConfUtils.o FileUtils.o DataUtils.o utils.o pnorm.o genome_size.o filter.o
OBJ_filter_ad = utils.o DataUtils.o FileUtils.o genome_size.o filter_ad.o
OBJ_EpiSig = ConfUtils.o FileUtils.o utils.o fast_pmf_cycle.o DataUtils.o thread_random.o global_data.o genome_size.o pnorm.o filter.o EpiSig_v1.0.o 
OBJ_EpiSig_MPI = ConfUtils.o FileUtils.o utils.o fast_pmf_cycle.o DataUtils.o thread_random.o global_data.o pnorm.o genome_size.o filter.o EpiSig_v1.0.MPI.o

MT = pre_EpiSig filter_ad EpiSig_v1.0
MPI = EpiSig_v1.0.MPI 

all: $(MT)
serial: $(MT)
mpi: $(MPI)

pre_EpiSig: $(OBJ_pre_EpiSig)
	$(CXX) $(CFLAGS) $(LFLAGS) $(OBJ_pre_EpiSig) -o pre_EpiSig
#filter: $(OBJ_filter)
#	$(CXX) $(CFLAGS) $(LFLAGS) $(OBJ_filter) -o filter
filter_ad: $(OBJ_filter_ad)
	$(CXX) $(CFLAGS) $(LFLAGS) $(OBJ_filter_ad) -o filter_ad
EpiSig_v1.0: $(OBJ_EpiSig)
	$(CXX) $(CFLAGS) $(LFLAGS) $(OBJ_EpiSig) -o EpiSig_v1.0
EpiSig_v1.0.MPI: $(OBJ_EpiSig_MPI)
	$(MPICXX) $(CFLAGS) $(MPIFLAGS) $(OBJ_EpiSig_MPI) -o EpiSig_v1.0.MPI

$(OBJ_MT) :
%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(LFLAGS) $< -o $@
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJ_MPI) : %.o: %.cpp
	$(MPICXX) -c $(CFLAGS) $(MPIFLAGS) $< -o $@

$(OBJ_BI) : %.o: %.cpp
ifeq ($(MAKECMDGOALS),mpi)
	$(MPICXX) -c $(CFLAGS) $(MPIFLAGS) $< -o $@
else
	$(CXX) -c $(CFLAGS) $(LFLAGS) $< -o $@
endif

.PHONY : clean clean_mpi
clean: 
	 rm -f *.o
	 rm -f $(MT) $(MPI)
clean_mpi:
	rm -f $(OBJ_BI)
