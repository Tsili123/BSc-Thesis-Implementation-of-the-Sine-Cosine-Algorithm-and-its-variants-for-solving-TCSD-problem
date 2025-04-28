
CXX = g++

CXXFLAGS = -Wall -I. -std=c++11	

PROGRAM = sca_prog

OBJS = main.o simulator.o population_criterion.o populations.o functions.o

$(PROGRAM): clean $(OBJS)
	$(CXX) $(OBJS) -o $(PROGRAM)

run: $(PROGRAM)		# Run with default arguments
	./$(PROGRAM) config.txt 30 6 20 10.0 30.00

clean: 
	rm -f $(PROGRAM) $(OBJS)