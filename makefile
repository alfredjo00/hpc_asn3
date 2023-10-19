BINS = newton
CFLAGS = -O2 -g -march=native -lm -lpthread

.PHONY : all
all : $(BINS) 

newton : newton.c
	gcc -o $@ $< $(CFLAGS)


.PHONY : clean
clean :
	rm -rf $(BINS)
	rm *.ppm
	rm *.png