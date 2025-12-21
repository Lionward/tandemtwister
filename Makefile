# Help message
help:
	@echo "Available targets:" 
	@echo "  make			 - Build the TandemTwist binary"
	@echo "  make clean		 - Remove generated files"
	@echo "  make test		 - Run tests (if applicable)"

# Check if CONDA_PREFIX is set
.DEFAULT_GOAL := default
ifndef CONDA_PREFIX
	HTSLIB_CHECK := $(shell if [ -e "/usr/local/lib/libhts.so" ]; then echo 1; else echo 0; fi)
ifeq ($(HTSLIB_CHECK), 1)
	# Use system htslib
CXXFLAGS=-std=c++20 -Wall -Wextra -Werror -O3 -g -mtune=native -I/usr/local/include 
LDFLAGS=-lhts -Wl,-rpath=/usr/local/lib -L/usr/local/lib
else
	$(error CONDA_PREFIX not set and htslib library not found in system)
endif
else
# Check if htslib library exists in the Conda environment
HTSLIB_CHECK := $(shell if [ -e "$(CONDA_PREFIX)/lib/libhts.so" ]; then echo 1; else echo 0; fi)
ifeq ($(HTSLIB_CHECK), 0)
	$(error htslib library not found in CONDA_PREFIX)
endif
# Use Conda htslib
CXXFLAGS=-std=c++23  -funroll-loops -ftree-vectorize -fopenmp  -Wall -Wextra -Werror -O3 -g -mtune=native -Wno-unknown-pragmas -I$(CONDA_PREFIX)/include 
LDFLAGS=-lhts -lfmt -Wl,-rpath=$(CONDA_PREFIX)/lib -L$(CONDA_PREFIX)/lib 
endif

CXX=g++
.PHONY: default clean install

# Source files
SRCS := src/assemblyInput.cpp src/LongReadInput.cpp src/parseCommandLine.cpp src/ReadPolisher.cpp src/ReadHaplotyper.cpp src/MotifAlignerCounter.cpp src/helper.cpp src/TandemTwisterConstructor.cpp src/creatVcfOutput.cpp src/main.cpp 

# Object files
OBJS := $(SRCS:.cpp=.o)


# Default target
default: tandemtwister

# Rule to compile each source file into an object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@



tandemtwister: $(OBJS)
	$(CXX) $(CXXFLAGS) -o tandemtwister $^ $(LDFLAGS)

clean:
	rm -f tandemtwister $(OBJS) $(DEPS)

install: tandemtwister
	@echo "Installing tandemtwister binary to /usr/local/bin"
	@sudo install -m 755 tandemtwister /usr/local/bin

.PHONY: test

test : tandemtwister
	@echo "Running tests"
	@./tandemtwister --germline -b test/example.bam -r test/ref_test.fa  -o test/test_sample -m test/example_tr_regions.bed -sn test -rt CCS -t 1  -s 0
