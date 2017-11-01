
CC = g++
CPPFLAGS = -std=c++11 -fopenmp
OPT = -o3

SRC = main.cpp
OBJ = main.o

all: main
	@echo "-------------------------"
	@echo "Compilation Successful"
	@echo "-------------------------"

main: main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp
	$(CC) $(CPPFLAGS) main.o -o main 

clean:
	rm -f *.o main
