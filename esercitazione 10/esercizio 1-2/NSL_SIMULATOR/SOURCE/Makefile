CC =  mpic++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: simulator.exe
simulator.exe : population.o NSL_SIMULATOR.o journey.o random.o
	$(CC) population.o journey.o random.o NSL_SIMULATOR.o -o simulator.exe $(AFLAGS)
NSL_SIMULATOR.o : NSL_SIMULATOR.cpp population.h journey.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)
population.o : population.cpp population.h random.h journey.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
journey.o : journey.cpp journey.h
	$(CC) -c journey.cpp -o journey.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*

