CXX ?= c++
GXX ?= g++
CXXFLAGS = -std=c++11 -fPIC -DNDEBUG -O3  -I./include 
LDFLAGS = -L./lib
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
INCLUDES = 

all: exactalign count_kmers_kmc

exactalign.o: exactalign.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@ 

exactalign: exactalign.o
	${CXX}  $^ ${LDFLAGS} -fPIC -o $@ ${LIBS}

	#${GXX}  -std=c++11 -DNDEBUG -O3  $^ -I${PREFIX}/include -L${PREFIX}/lib -fPIC -o program -lsdsl -ldivsufsort -ldivsufsort64
	#${GXX}    ${LDFLAGS} ${CXXFLAGS} $^ -o program -lsdsl -ldivsufsort -ldivsufsort64

	# g++ -std=c++11 -DNDEBUG -O3 -I$(CPATH) $^ -o program -L$(LIBRARY_PATH) -lsdsl -ldivsufsort -ldivsufsort64


count_kmers_kmc.o: count_kmers_kmc.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -I./KMC -I./KMC/kmc_counter $< -o $@ 

count_kmers_kmc: count_kmers_kmc.o
	${CXX}  $^ ${LDFLAGS} -I./KMC -I./KMC/kmc_counter ./KMC/kmc_api/*.o  -fPIC -o $@ ${LIBS}


#count_kmers: count_kmers_kmc.cpp 
# ${GXX} -std=c++11 -Wall -I./KMC -I./KMC/kmc_counter ./KMC/kmc_api/*.o count_kmers_kmc.cpp -o count_kmers_kmc


clean:
	rm -f exactalign exactalign.o count_kmers_kmc count_kmers_kmc.o
