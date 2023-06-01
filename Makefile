CXX := g++
MAIN := main
MAIN_OBJECTS := main.o graph.o parse.o stat.o basicstat.o connectstat.o biconnected.o eigencentrality.o closenesscentrality.o betweennesscentrality.o kclist.o nucleusdecomposition.o countstat.o makeplot.o html.o statanalyzer.o
TEST := test
TEST_OBJECTS := test.o graph.o parse.o diameter.o
PYPLOT_DIR := pyplot
CXXFLAGS := -std=c++17 -Wall -O3
LDFLAGS :=
MAIN_DEPS := $(MAIN_OBJECTS:.o=.d)
TEST_DEPS := $(TEST_OBJECTS:.o=.d)
CPPFLAGS := -MMD -MP

all: $(MAIN) $(TEST) $(PYPLOT_DIR)

$(MAIN): $(MAIN_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(TEST): $(TEST_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(PYPLOT_DIR):
	mkdir -p pyplot

clean:
	rm -rf $(MAIN) $(MAIN_OBJECTS) $(TEST) $(TEST_OBJECTS) $(MAIN_DEPS) $(TEST_DEPS) *.html pyplot __pycache__

-include $(MAIN_DEPS)
-include $(TEST_DEPS)
