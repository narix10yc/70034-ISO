CC = gcc


CFLAGS = 
LDFLAGS = -fopenmp

OBJDIR = obj

SRC = timing.c complex_mul.c

OBJ = $(addprefix $(OBJDIR)/,$(SRC:.c=.o))

program_scalar : program_scalar.c $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lm


$(OBJDIR)/timing.o : timing.c 
	$(CC) -O1 -funroll-loops -c $^ -o $@ 


$(OBJDIR)/%.o: %.c
	$(CC) -O1 -march=native $(CFLAGS) -c $< -o $@



playground : playground.c $(OBJDIR)/timing.o 
	$(CC) $(LDFLAGS) playground.c $(OBJDIR)/timing.o -o $@ -lm



.PHONY : clean

clean :
	rm $(OBJ) program_scalar


