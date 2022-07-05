CXX = g++
MAIN = main
MAIN_OBJECTS = main.o graph.o parse.o stat.o basicstat.o connectstat.o biconnected.o eigencentrality.o closenesscentrality.o betweennesscentrality.o countstat.o makeplot.o html.o
TEST = test
TEST_OBJECTS = test.o graph.o parse.o diameter.o
PYPLOT_DIR = pyplot
CXXFLAGS = -std=c++17 -Wall -O0 -g
LDFLAGS =

all: $(MAIN) $(TEST) $(PYPLOT_DIR)

$(MAIN): $(MAIN_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(TEST): $(TEST_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(PYPLOT_DIR):
	mkdir -p pyplot

clean:
	rm -rf $(MAIN) $(MAIN_OBJECTS) $(TEST) $(TEST_OBJECTS) *.html pyplot __pycache__
