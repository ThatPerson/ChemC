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
#include "reaction.c"


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
    for (i = 0; i < argc-1; i++) {
        if (strcmp(argv[i], "-pt") == 0) {
            printf("Reading %s\n", argv[i+1]);
            read_periodic_table(argv[i+1]);
            initialise(argv[i+1]);
        } else if (strcmp(argv[i], "-s") == 0) {
            printf("Script in %s\n", argv[i+1]);
        } else if (strcmp(argv[i], "-e") == 0) {
            print_element(&periodic_table[find_element(argv[i+1])]);
        } else if (strcmp(argv[i], "-l") == 0) {
            if (i < argc-2)
                printf("%0.3f\n", periodic_table[find_element(argv[i+1])].atomic_radius + periodic_table[find_element(argv[i+2])].atomic_radius);
            else
                printf("%0.3f\n", 2*periodic_table[find_element(argv[i+1])].atomic_radius);
        } else if (strcmp(argv[i], "-m") == 0) {
            if (i < argc-2) {
                float m[2] = {periodic_table[find_element(argv[i+1])].molar, periodic_table[find_element(argv[i+2])].molar};
                printf("%0.3f\n", (m[0]*m[1])/(m[0]+m[1]));
              } else {
                float m = periodic_table[find_element(argv[i+1])].molar;
                printf("%0.3f\n", m/2);
              }
        } else if (strcmp(argv[i], "-c") == 0) {
            struct Compound NaCl;
            strcpy(NaCl.name, argv[i+1]);
            find_constituents(&NaCl);

            printf("%s %f\n", find_name(&NaCl), compound_molarity(&NaCl));

            strcpy(NaCl.name, find_name(&NaCl));
            struct Bond bonds[50];

            int q = predict_bonding(bonds,50, &NaCl);

            char bond_type[5][10] = {"single", "single", "double", "triple", "quadruple"};
            int c = q;
            if (c < 1)
                c = 1;

            for (i = 0; i < c; i++) {
                if (strcmp(bonds[i].atoms[0].name, "") != 0) {
                    printf("%s is involved in a %s bond with %s\n", bonds[i].atoms[0].name, bond_type[bonds[i].num_bonds], bonds[i].atoms[1].name);
                }
            }
            float sd[50];
            int ls = predict_ir(sd, bonds, q);
            int pqw;
            for (pqw = 0; pqw < ls; pqw++) {
                if (sd[pqw] == 0) {
                    break;
                }
                printf("%0.2f\n", sd[pqw]);
            }

        } else if (strcmp(argv[i], "-r") == 0) {
            struct Mixture reactants, products;
            int lk;
            for (lk = i+1; lk<argc; lk++) {
                printf("%d\n", lk-i-1);
                strcpy(reactants.compounds[lk-i-1].name, argv[lk]);
                find_constituents(&reactants.compounds[lk-i-1]);
                reactants.num_compounds++;
            }
            predict_reaction(&reactants, &products);
            printf("%d\n", products.num_compounds);
            printf("%d\n", products.compounds[0].num_constituents);
            printf("%s\n", products.compounds[1].constituents[0].small);
            for (lk = 0; lk < products.num_compounds; lk++) {
                printf("%s\n", find_name(&products.compounds[lk]));
            }
        }
    }


    return 0;
}
