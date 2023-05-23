
asm : main.s
	as -o main.o main.s
	ld -o main main.o

main : main.c
	gcc -o main main.c

.PHONY : clean 

clean :
	rm main.o main