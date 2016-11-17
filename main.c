#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Define constants.
#define BOLTZMANN 1.23 * pow(10, -23); // Boltzmann constants
#define AVOGADRO 6.02*pow(10, 23); // Avogadro's number
#define PLANCK 6.626 * pow(10, -34); // Planck's Constant
#define MP_K 1.2025*pow(10, -4); // Melting point prediction constants.
#define MP_D 3.45*pow(10, -11);
#define RYDBERG 13.6057

struct Element {
    char name[100];
    char small[5];
    int position;
    float molar;
    int atomic_number;
    int num_electrons;
    float electronegativity;
    float atomic_radius;
    int valency;
    int shells[6][6][6];
};

struct Element periodic_table[500];
int pt_length = 0;

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

void get_shells(struct Element * p) {
    /**
     * Purpose:     Predicts which orbitals electrons fill at lowest energy state. 
     *              Does not take into account penetration of orbitals - so it may
     *              not fill properly.
     * Arguments:   Pointer to element struct.
     * Method:      Very ugly algorithm. Loops over all electrons present working
     *              up the shells array in the element, filling from lowest (ie 1s
     *              to highest (ie to 2s, 2p, 3s, 3p etc).
     * Issues:      Does not take into account penetration - so fills 3d before 4s.
     *              Additionally, does not take into account individual energy 
     *              levels - so Cu and Cr are not done properly. Gives reasonable
     *              approximation for < 20 electrons.
     */
    
    int n = 0, l =0, m = 0, a_n = p->num_electrons, i, q;
    while (a_n > 0) {
        q = 0;
        if (n > 7)
            break;
        if (p->shells[n][l][m + l] >= 2) {
            if (m < l) {
                m++;
            } else {
                if (l < n && l < 5) {
                    l++;
                    m = -l;
                } else {
                    n++;
                    l = 0;
                    m = 0;
                }
            }
        } else if (p->shells[n][l][m+l] == 1) {
            for (i = m + l; i < 2*l + 1; i++) {
                if (p->shells[n][l][i] == 0 && q == 0) {
                    p->shells[n][l][i]++;
                    a_n--;
                    q = 1;
                    break;
                }
            }
            if (q == 0) {
                p->shells[n][l][m+l]++;
                a_n--;
            }
        } else {
            p->shells[n][l][m+l]++;
            a_n--;
        }
    }
    
    return;
}

int length(int q[50]) {
    int i = 0;
    while (q[i] != 0) {
        printf("%d\n", i);
        if (q[i] == 0) {
            return i -1;
        }
        i++;
    }
    return i-1;
}

float predict_electron_energy(struct Element * p, int n, int l, int m, int x, int legacy) {
    /**
     * Purpose:     Predict the energy of an electron in n l m.
     * Arguments:   nlm. x - amount of electron shielding. p - pointer to atom.
     * Method:      1.
     *                  Uses Rydberg equation, taking x as the shielding.
     *              2. 
     *                  Same as 1 but uses Slater's rules to predict the shielding.
     *              Legacy switches mode - 1 works better for prediction, but 2 is more accurate.
     */
    float energy;
    if (legacy == 1) {
        int Z = p->atomic_number - x;
        energy = -RYDBERG * (pow(Z, 2) / pow(n, 2));
    } else {
        int S = 0;
        int i, j, k, s_n = n - 1, s_l = l, s_m = m;
        int p[3] = {0, 0, 0};
        for (i = 0; i < s_l; i++) {
            for (j = 0; j < length(p.shells[s_n][i]); j++) {
                p[0] += p->shells[s_n][i][j];
                p[0]--;
            }
        }
        
    }
        
    
    
    return energy;
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
    read_periodic_table("data.dat");
    /* end */
    
    
    int i;
    for (i = 0; i < argc - 1; i++) {
        if (strcmp(argv[i], "-pt") == 0) {
            printf("Reading %s\n", argv[i+1]);
            read_periodic_table(argv[i+1]);
            initialise(argv[i+1]);
        } else if (strcmp(argv[i], "-s") == 0) {
            printf("Script in %s\n", argv[i+1]);
        }
    }
    printf("%s\n", periodic_table[2].name);

    int s[5] = {2, 3, 1, 1, 1};
    printf("%d\n", length(s));
    
    return 0;
}
