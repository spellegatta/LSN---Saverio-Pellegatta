CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o random.o functions.o
	$(CC) random.o main.o functions.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.x seed.out
