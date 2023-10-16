BINS = newton
CFLAGS = -O2 -lpthread -march=native -lm -g

.PHONY : all
all : $(BINS) 

newton : newton.c
	gcc $(CFLAGS) -o $@ $<


.PHONY : clean
clean :
	rm -rf $(BINS)