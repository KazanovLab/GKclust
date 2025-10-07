# Compiler
CXX = g++
LIBS = -lz

# Source files
SRC = $(wildcard src/*.cpp)

# Output executable
OUT = gkclust

# Default rule: compile all .cpp files into an executable
$(OUT): $(SRC)
	$(CXX) -o $(OUT) $(SRC) $(LIBS)

# Clean rule: remove executable
clean:
	rm -f $(OUT)
