# Makefile

# Compiler and compiler options
CPP = g++
OPT = -Wall -O2

# Executable name
EXE = d2ts

# Source files
SRCS = d2ts.cpp tabu_search.cpp reads.cpp main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Build rule for object files
%.o: %.cpp
	$(CPP) -c $< $(OPT) -o $@

# Default target to build the executable
all: $(EXE)

$(EXE): $(OBJS)
	$(CPP) -o $@ $(OBJS) $(OPT)

# Clean rule to remove object files and the executable
clean:
	rm -f $(EXE) $(OBJS) *~
