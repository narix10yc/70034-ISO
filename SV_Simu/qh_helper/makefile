CC = gcc

VECFLAGS = 

CFLAGS = 
LDFLAGS = -fopenmp

OBJDIR = obj

SRC = timing.c statevector.c qh_gates.c qh_type.c

OBJ = $(addprefix $(OBJDIR)/,$(SRC:.c=.o))

program : program.c $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lm

$(OBJDIR)/timing.o : timing.c 
	$(CC) -O1 -funroll-loops -c $^ -o $@ 

$(OBJDIR)/%.o: %.c
	$(CC) -O1 $(CFLAGS) -c $< -o $@



playground : playground.c $(OBJDIR)/timing.o 
	$(CC) $(LDFLAGS) playground.c $(OBJDIR)/timing.o -o $@ -lm



.PHONY : clean clean_playground

clean :
	rm $(OBJ) program

clean_playground :
	rm playground