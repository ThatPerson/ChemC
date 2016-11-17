CC=gcc
CFLAGS="-Wall"

debug:clean
	$(CC) $(CFLAGS) -g -lm -o chemsi main.c
stable:clean
	$(CC) $(CFLAGS) -o chemsi main.c -lm
clean:
	rm -vfr *~ chemsi
