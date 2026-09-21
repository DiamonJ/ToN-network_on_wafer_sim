# dumpi2ccdg — DUMPI trace → CCDG converter
# Requires libundumpi headers/libs from the sst-dumpi build (see ../README_REPRODUCE.md).
# Override via make DUMPI_PREFIX=/path/to/install if not using ./install.
DUMPI_PREFIX ?= $(abspath ../install)
CXX      := g++
CXXFLAGS := -std=c++11 -O2 -Wall -Wextra
INCLUDES := -I$(DUMPI_PREFIX)/include
LDFLAGS  := -L$(DUMPI_PREFIX)/lib
LIBS     := -lundumpi -lpthread

TARGET   := dumpi2ccdg
SRCS     := dumpi2ccdg.cpp
OBJS     := $(SRCS:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)
