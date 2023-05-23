#include <stdio.h>

int main() {
    int num1 = 10;
    int num2 = -10;

    printf("%d in binary: ", num1);
    for (int i = 8 * sizeof(int) - 1; i >= 0; i--) {
        if (i % 8 == 7)
            printf(" ");
        printf("%d", (num1 >> i) & 1);
    }
    printf("\n");

    printf("%d in binary: ", num2);
    for (int i = 8 * sizeof(int) - 1; i >= 0; i--) {
        if (i % 8 == 7)
            printf(" ");
        printf("%d", (num2 >> i) & 1);
    }
    printf("\n");

    printf("%u\n", ~0);

    return 0;
}
