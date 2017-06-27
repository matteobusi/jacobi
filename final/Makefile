CC=icc

src = src/main.cpp $(wildcard src/solvers/*.cpp)

CFLAGS = -O3 -I/home/spm1611/fastflow -std=c++11 -DNO_DEFAULT_MAPPINGS
VECREPORT = -qopt-report-phase=vec -qopt-report=2
LDFLAGS = -lpthread


jacobix: 
	$(CC) $(src) $(CFLAGS) -o bin/$@ $(LDFLAGS) $(VECREPORT)

jacobim:
	$(CC) $(src) -mmic $(CFLAGS) -o bin/$@ $(LDFLAGS) $(VECREPORT)

offload: bin/jacobim
	scp bin/jacobim mic1:jacobim

.PHONY: clean
clean:
	rm -f bin/* results/* graphs/*
