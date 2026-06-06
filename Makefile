CXX = mpicxx
CXXFLAGS = -Wall -O3 -fopenmp 
LDFLAGS = -fopenmp 
LDLIBS = -lmuparser

# name of the executable
TARGET = laplace

SRCS = src/main.cpp 
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)