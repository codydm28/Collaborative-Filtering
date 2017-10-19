pa2: CSR_STRUCTURE.o main.o
	g++ CSR_STRUCTURE.o main.o -o pa2
main.o: CSR_STRUCTURE.o main.cpp
	g++ -c main.cpp
CSR_STRUCTURE.o: CSR_STRUCTURE.cpp CSR_STRUCTURE.h
	g++ -c CSR_STRUCTURE.cpp
clean:
	rm -f *.o pa2.out