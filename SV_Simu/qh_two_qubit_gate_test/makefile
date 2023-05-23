CC = gcc


CFLAGS = 
LDFLAGS = -fopenmp

OBJDIR = obj

SRC = timing.c statevector.c qh_type.c two_qubit_explicit.c two_qubit_general.c two_qubit_quest.c

OBJ = $(addprefix $(OBJDIR)/,$(SRC:.c=.o))

program : program.c $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lm

$(OBJDIR)/timing.o : timing.c 
	$(CC) -O1 -funroll-loops -c $^ -o $@ 

$(OBJDIR)/qh_type.o : qh_type.c 
	$(CC) -c $^ -o $@ 


$(OBJDIR)/%.o: %.c
	$(CC) -O1 -fopenmp $(CFLAGS) -c $< -o $@



playground : playground.c $(OBJDIR)/timing.o 
	$(CC) $(LDFLAGS) playground.c $(OBJDIR)/timing.o -o $@ -lm



.PHONY : clean

clean :
	rm $(OBJ) program


