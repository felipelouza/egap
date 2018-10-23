CC = g++
CPPFLAGS = -Wall -ansi -pedantic -g -O2

BCR_BWTCollection_obs = BCR_BWTCollection.o BWTCollection.o BCRexternalBWT.o Tools.o Sorting.o TransposeFasta.o Timer.o
BCR_BWTCollection: $(BCR_BWTCollection_obs)
	$(CC) -o BCR_BWTCollection $(BCR_BWTCollection_obs)

clean:
	rm -f core *.o *~ BCR_BWTCollection

depend:
	$(CC) -MM *.cpp *.c > dependencies.mk

include dependencies.mk
