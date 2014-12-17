all: filevercmp

clean:
	rm -f filevercmp
	rm -f filevercmp.o

.PHONY: all clean

filevercmp.o: filevercmp.c main.c filevercmp.h
	gcc -c filevercmp.c

filevercmp: filevercmp.o
	gcc -o filevercmp main.c filevercmp.o
