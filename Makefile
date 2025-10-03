# Compiler
CXX = g++

# Source files
SRC = $(wildcard src/*.cpp)

# Output executable
OUT = gkclust

# Default rule: compile all .cpp files into an executable
$(OUT): $(SRC)
	$(CXX) -o $(OUT) $(SRC)

# Clean rule: remove executable
clean:
	rm -f $(OUT)
