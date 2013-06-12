# Makefile May 2013 by Akash Khare

CC = g++

CFLAGS = -O -pg -g 

RM = rm -f 

OBJS = ngila_wrapper.o


NgilaWrapper: $(OBJS)
	$(CC) $(OBJS) -o NgilaWrapper 

NgilaWrapper.o: ngila_wrapper.cpp
	$(CC) $(CFLAGS) ngila_wrapper.cpp

clean:
	$(RM) *.o NgilaWrapper


