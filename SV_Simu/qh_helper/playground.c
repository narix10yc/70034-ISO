#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct S1
{
    double data[4];
} S1;

typedef struct S2
{
    double *data;
} S2;



int main(int argc, char *argv[]) {
    
    char s[3][3] = { "AB", "CD", "EF" };    

    auto p = s;

    printf("%p\n", s);
    printf("%p\n", s+1);
    printf("%p\n", &s + 1);

    printf("%p\n", p);
    printf("%p\n", p+1);
    printf("%p\n", &s + 1);
    return 0;
}