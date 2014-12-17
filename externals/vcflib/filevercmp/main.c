#include "filevercmp.h"
#include <stdio.h>


int main(int argc, char** argv) {
    if (argc != 3) {
        printf("usage: %s [a] [b]\n", argv[0]);
        printf("shows version-string comparison of strings a and b\n", argv[0]);
        printf("for instance, chr1 < chr10, 1 < 10, abca < bcac\n");
        return 1;
    }
    int c = filevercmp(argv[1], argv[2]);
    printf("%s < %s = %i\n", argv[1], argv[2], c);
    return 0;
}
