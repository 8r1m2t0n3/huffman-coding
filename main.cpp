#include <iostream>

#include "Huffman.h"

void show_help() {
    std::cout << "Usage:\n\t huffman mode(0 - encode, 1 - decode) inputfile outputfile" << std::endl;
    exit(1);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        show_help();
    }
    Huffman h(argv[2], argv[3]);
    if (*argv[1] == '0') {
        h.encode();
    }
    else if (*argv[1] == '1') {
        h.decode();
    }
    else {
        show_help();
    }
    return 0;
}
