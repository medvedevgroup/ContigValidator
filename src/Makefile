CXX ?= c++
GXX ?= g++
CXXFLAGS = -std=c++11 -fPIC -DNDEBUG -O3  -I./include 
LDFLAGS = -L./lib
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
INCLUDES = 

all: exactalign 

exactalign.o: exactalign.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@ 

exactalign: exactalign.o
	${CXX}  $^ ${LDFLAGS} -fPIC -o $@ ${LIBS}


#${GXX}  -std=c++11 -DNDEBUG -O3  $^ -I${PREFIX}/include -L${PREFIX}/lib -fPIC -o program -lsdsl -ldivsufsort -ldivsufsort64
#${GXX}    ${LDFLAGS} ${CXXFLAGS} $^ -o program -lsdsl -ldivsufsort -ldivsufsort64
# g++ -std=c++11 -DNDEBUG -O3 -I$(CPATH) $^ -o program -L$(LIBRARY_PATH) -lsdsl -ldivsufsort -ldivsufsort64


clean:
	rm -f exactalign exactalign.o count_kmers_kmc count_kmers_kmc.o



