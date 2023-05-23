CC = gcc


CFLAGS = 
LDFLAGS = -fopenmp

OBJDIR = obj

SRC = timing.c single_qubit.c qh_type.c

OBJ = $(addprefix $(OBJDIR)/,$(SRC:.c=.o))

program : program.c $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ -lm


$(OBJDIR)/timing.o : timing.c 
	$(CC) -O1 -funroll-loops -c $^ -o $@ 

$(OBJDIR)/qh_type.o : qh_type.c 
	$(CC) -c $^ -o $@ 

$(OBJDIR)/single_qubit.o : single_qubit.c
	$(CC) -O1 -fopenmp $(CFLAGS) -c $^ -o $@


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@



playground : playground.c $(OBJDIR)/timing.o 
	$(CC) $(LDFLAGS) playground.c $(OBJDIR)/timing.o -o $@ -lm



.PHONY : clean

clean :
	rm $(OBJ) program


