GXX ?= g++
CXXFLAGS ?= -std=c++11 -DNDEBUG -O3  -I/home/pzm11/include 
LDFLAGS ?= -L/home/pzm11/lib
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
INCLUDES = 

all: program 
#all: program count_kmers


.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

program: program.o
	${GXX}  $^ ${LDFLAGS} -fPIC -o $@ ${LIBS}

	#${GXX}  -std=c++11 -DNDEBUG -O3  $^ -I${PREFIX}/include -L${PREFIX}/lib -fPIC -o program -lsdsl -ldivsufsort -ldivsufsort64
	#${GXX}    ${LDFLAGS} ${CXXFLAGS} $^ -o program -lsdsl -ldivsufsort -ldivsufsort64

	# g++ -std=c++11 -DNDEBUG -O3 -I$(CPATH) $^ -o program -L$(LIBRARY_PATH) -lsdsl -ldivsufsort -ldivsufsort64

count_kmers: count_kmers_kmc.cpp 
	${GXX} -std=c++11 -Wall -I./KMC -I./KMC/kmc_counter ./KMC/kmc_api/*.o count_kmers_kmc.cpp -o count_kmers_kmc

findCommonKmers: findCommonKmers.cpp
	${GXX} -O3 $^ -o findCommonKmers


clean:
	rm -f findCommonKmers findCommonKmers.o
	rm -f program program.o
