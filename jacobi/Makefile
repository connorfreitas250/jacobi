C = gcc
CFLAGS = -g -Wall -pthread -O3
DEPS = jacobi.h
FILES = jacobi.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< 

jacobi: $(FILES) 
	$(CC) $(CFLAGS) -o $@ $^  

clean:
	-rm jacobi $(FILES)
