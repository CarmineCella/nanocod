PRODUCT := nanocod

CXX := g++
LINKER := g++
INCDIRS := -I.
CXXFLAGS := -std=c++11 -Wall -Wextra

SRCFILES := $(wildcard *.cpp) 
OBJFILES := $(patsubst %.cpp,%.o,$(SRCFILES))
DEPFILES := $(patsubst %.cpp,%.d,$(SRCFILES))

$(PRODUCT): $(OBJFILES)
	$(LINKER) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c $< -o $@

%.d: %.cpp
	$(CXX) $(INCDIRS) -MM $< > $@

-include $(DEPFILES)


