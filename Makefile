
CC=g++

CFLAGS=-std=c++11 -O3

objects=LinearCDSfold

linearcdsfold:

		$(CC) $(CFLAGS) src/main.cpp -o LinearCDSfold 

clean:
	-rm $(objects)
