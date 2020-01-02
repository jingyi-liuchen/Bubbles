CXX = mpiicpc
LIBS = ../eigen-git-mirror/Eigen/

CXXFLAGS =  -g -std=c++11 
CXXFLAGS += -I $(LIBS)

OBJDIR=obj
SRCDIR=src

SRC := $(wildcard $(SRCDIR)/*.cpp)
OBJ := $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))
DEP := $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.d, $(SRC))

all: $(OBJDIR)/bubble_ana

$(OBJDIR)/bubble_ana: $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

-include $(DEP)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

.PHONY : clean
clean:
	rm $(OBJ) $(DEP) $(OBJDIR)/bubble_ana

.PHONY:  test
test:
	@printf "===RUNNING TESTS===\n"
	@chmod u+rx ./tests/test.sh
	@cd ./tests/; ./test.sh
	@printf "\n"
