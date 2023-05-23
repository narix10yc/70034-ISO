CC = gcc

VECFLAGS = 

EXTRA_CFLAGS = 
CFLAGS = -O1
LDFLAGS = -fopenmp

OBJDIR = obj

SRC = timing.c no_sse.c scalar.c sse.c avx.c avx512.c

OBJ = $(addprefix $(OBJDIR)/,$(SRC:.c=.o))

program : main.c $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lm

$(OBJDIR)/timing.o : timing.c 
	$(CC) -O1 -funroll-loops -c $^ -o $@ 

$(OBJDIR)/no_sse.o: no_sse.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -mno-sse -c $< -o $@

$(OBJDIR)/scalar.o: scalar.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -msse -c $< -o $@

$(OBJDIR)/sse.o: sse.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -msse -c $< -o $@

$(OBJDIR)/avx.o: avx.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -mavx -c $< -o $@

$(OBJDIR)/avx512.o: avx512.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -mavx512f -c $< -o $@

.PHONY : clean clean_playground

clean :
	rm $(OBJ) program