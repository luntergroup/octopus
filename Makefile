CPP      = g++
CPPFLAGS = g -Wall -std=c++14 -O3

SOURCES = main.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

TARGET  = octopus

LDFLAGS =
LDLIBS  =

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CPP) $(CPPFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJECTS) $(LDLIBS)

depend: .depend

.depend: $(SOURCES)
    rm -f ./.depend
    $(CPP) $(CPPFLAGS) -MM $^>>./.depend;

clean:
    rm -f $(OBJECTS)

dist-clean: clean
    rm -f *~ .depend

include .depend