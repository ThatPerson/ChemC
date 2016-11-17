#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Define constants.
#define BOLTZMANN 1.23 * pow(10, -23) // Boltzmann constants
#define AVOGADRO 6.02*pow(10, 23) // Avogadro's number
#define PLANCK 6.626 * pow(10, -34) // Planck's Constant
#define MP_K 1.2025*pow(10, -4) // Melting point prediction constants.
#define MP_D 3.45*pow(10, -11)
#define RYDBERG 13.6057

// Haxxy stuff
#define arrlen(x)  (sizeof(x) / sizeof((x)[0]))


#include "structs.c"

struct Element periodic_table[500];
int pt_length = 0;

#include "element.c"
#include "compounds.c"



void read_periodic_table(char * filename) {
    FILE* reading;
    int current = 0;
    reading = fopen(filename, "r");
    
    while (!feof(reading)) {
        fscanf(reading, "%s %s %d %f %d %f %f", periodic_table[current].name, periodic_table[current].small, &periodic_table[current].position, &periodic_table[current].molar, &periodic_table[current].atomic_number, &periodic_table[current].electronegativity, &periodic_table[current].atomic_radius);
        periodic_table[current].num_electrons = periodic_table[current].atomic_number;
        current++;
    }
    pt_length = current - 1;
    
    fclose(reading);
    return;
}

int length(int q[50]) {
    int i = 0;
    while (q[i] != 0) {
        if (q[i] == 0) {
            return i -1;
        }
        i++;
    }
    return i-1;
}

int initialise(char * datafile) {
    int i;
    read_periodic_table(datafile);
    for (i = 0; i < pt_length; i++) {
        get_shells(&periodic_table[i]);
    }
    return 1;
}

int main(int argc, char **argv) {
    
    /*Testing code. Just so I don't have to pass args*/
    initialise("data.dat");
    /* end */
    
    
    int i;
    for (i = 0; i < argc - 1; i++) {
        if (strcmp(argv[i], "-pt") == 0) {
            printf("Reading %s\n", argv[i+1]);
            read_periodic_table(argv[i+1]);
            initialise(argv[i+1]);
        } else if (strcmp(argv[i], "-s") == 0) {
            printf("Script in %s\n", argv[i+1]);
        } else if (strcmp(argv[i], "-e") == 0) {
            print_element(&periodic_table[find_element(argv[i+1])]);
        }
    }
    
    struct Compound NaCl;
    NaCl.constituents[0] = periodic_table[find_element("Na")];
    NaCl.constituents[1] = periodic_table[find_element("Cl")];
    NaCl.num_constituents = 2;
    printf("%f///%f\n", predict_melting_point(&NaCl, 1), predict_melting_point(&NaCl, 0));
    
    return 0;
}
