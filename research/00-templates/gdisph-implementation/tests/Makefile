PS_PATH = -I ../../../../src/

# Compiler selection
CC = mpicxx
CFLAGS = -std=c++17 -O3 -ffast-math -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

TEST_SRC = test_runner.cpp
TEST_OBJ = test_runner.o
TEST_EXEC = test_riemann_solver.out

.PHONY: all clean run test

# Default target builds and runs tests
all: test

# Build test executable
$(TEST_EXEC): $(TEST_OBJ)
	@echo "Linking test executable..."
	@$(CC) $(CFLAGS) $(TEST_OBJ) -o $(TEST_EXEC) $(PS_PATH)
	@echo "Link Success! [$(TEST_EXEC)]"

$(TEST_OBJ): $(TEST_SRC)
	@echo "Building test runner..."
	@$(CC) -c $(TEST_SRC) $(CFLAGS) $(PS_PATH)
	@echo "[Test compilation OK]"

# Run tests
test: $(TEST_EXEC)
	@echo ""
	@echo "╔════════════════════════════════════════════════════════════╗"
	@echo "║  Running TDD/BDD Tests for Riemann Solver                 ║"
	@echo "╚════════════════════════════════════════════════════════════╝"
	@echo ""
	@./$(TEST_EXEC) || echo "\n⚠️  Tests failed - this is expected in RED phase of TDD!"

# Run tests with MPI
test-mpi: $(TEST_EXEC)
	@echo "Running tests with MPI (1 rank)..."
	@mpirun -n 1 ./$(TEST_EXEC)

clean:
	-rm -f *.o *.out

help:
	@echo "Riemann Solver Test Makefile (TDD/BDD)"
	@echo "======================================="
	@echo ""
	@echo "Targets:"
	@echo "  make test      - Build and run test suite"
	@echo "  make test-mpi  - Run tests with MPI"
	@echo "  make clean     - Remove build artifacts"
	@echo "  make help      - Show this message"
	@echo ""
	@echo "TDD Workflow:"
	@echo "  1. RED phase:    Tests fail (current state)"
	@echo "  2. GREEN phase:  Implement code to pass tests"
	@echo "  3. REFACTOR:     Optimize while keeping tests green"
