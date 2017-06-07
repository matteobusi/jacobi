CC=g++

src = main.cpp $(wildcard solvers/*.cpp)

CFLAGS = -O3 -I./fastflow -std=c++11 -DNO_DEFAULT_MAPPINGS
LDFLAGS = -lpthread


jacobix:
	$(CC) $(src) $(CFLAGS) -o $@ $(LDFLAGS)

jacobim:
	$(CC) $(src) -mmic $(CFLAGS) -o $@ $(LDFLAGS)

offload: $(obj)
	scp jacobim mic1:jacobim

.PHONY: clean
clean:
	rm -f jacobim jacobix *~
