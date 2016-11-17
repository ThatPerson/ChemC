CC=gcc
CFLAGS="-Wall"

debug:clean
	$(CC) $(CFLAGS) -g -o chemsi main.c
stable:clean
	$(CC) $(CFLAGS) -o chemsi main.c
clean:
	rm -vfr *~ chemsi
