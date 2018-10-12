all: program count_kmers

program: program.cpp
	${GXX}  -std=c++11 -DNDEBUG -O3  $^ -I${PREFIX}/include -I${PREFIX}/lib -o program -lsdsl -ldivsufsort -ldivsufsort64

# g++ -std=c++11 -DNDEBUG -O3 -I$(CPATH) $^ -o program -L$(LIBRARY_PATH) -lsdsl -ldivsufsort -ldivsufsort64

count_kmers: count_kmers_kmc.cpp 
	${GXX} -std=c++11 -Wall -I./KMC -I./KMC/kmc_counter ./KMC/kmc_api/*.o count_kmers_kmc.cpp -o count_kmers_kmc

findCommonKmers: findCommonKmers.cpp
	${GXX} -O3 $^ -o findCommonKmers


clean:
	rm -f findCommonKmers
	rm -f program
