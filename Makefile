CC     = g++
CFLAGS = g -Wall -std=c++14 -O3

SOURCES = main.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

TARGET  = octopus

LDFLAGS =

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c

clean:
	$(RM) $(TARGET)
	