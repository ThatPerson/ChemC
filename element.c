int group(struct Element *p) {
    int n;
    for (n = 0; n < 5; n++) {
        if (p->shells[n][0][0] == 0) {
            break;
        }
    }
    return n;
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
    float Z;
    if (legacy == 1) {
        Z = p->atomic_number - x;

    } else {
        float S = 0;
        int i, j, k, s_n = n - 1;

        int ps[3] = {0, 0, 0};
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {

                ps[0] += p->shells[s_n][i][j];
            }

        }
        ps[0] = ps[0] - 1;
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                if (s_n - 1 >= 0)
                    ps[1] += p->shells[s_n-1][i][j];
            }
        }
        for (i = 0; i < s_n-1; i++) {
            for (j = 0; j < 5; j++) {
                for (k = 0; k <= 5; k++) {
                    ps[2] += p->shells[i][j][k];
                }
            }
        }
        S = 0.35 * ps[0] + 0.85*ps[1] + 1*ps[2];
        Z = p->atomic_number - S;
    }

    energy = -RYDBERG * (pow(Z, 2) / pow(n, 2));

    return energy;
}

float predict_valence_electron_energy(struct Element *p) {
    /**
     * Purpose:     Find the energy of the outermost electron.
     * Arguments:   p - pointer to atom
     * Method:      Loop over all orbitals until the outermost one is found. Calculate the energy of it.
     */
    int c, w, n, l, m, c_n=0, c_l=0, c_m=0;
    float energy;
    for (n = 0; n < 5; n++) {
        for (l = 0; l < 5; l++) {
            w = 0;
            for (m = 0; m < 5; m++) {
                if (p->shells[n][l][m] != 0) {
                    c_n = n;
                    c_l = l;
                    c_m = m;
                    w += p->shells[n][l][m];
                }
            }
            c += w;
        }
    }
    energy = predict_electron_energy(p, c_n+1, c_l, c_m-c_l, c, 0);
    return energy;
}

int get_valence_electrons(struct Element *p) {
    int n,l,m, outershell = 0;
    for (n = 0; n < 5; n++) {
        if (p->shells[n][0][0] == 0) {
            break;
        }
    }
    n--;
    for (l = 0; l < 5; l++) {
        for (m = 0; m < 5; m++) {
            outershell += p->shells[n][l][m];
        }
    }
    return outershell;
}


int atom_valency(struct Element *p) {
    int valency = 1, outershell = get_valence_electrons(p);
    if (outershell < 8) {
        if (outershell < 4)
            valency = outershell;
        else
            valency = 8-outershell;
    }
    return valency;
}


void print_element(struct Element *p) {
    printf("Name:\t\t\t%s\nSmall:\t\t\t%s\nGroup:\t\t\t%d\nAtomic Number:\t\t%d\nMolar Mass:\t\t%0.2f\nAtomic Radius:\t\t%0.2f\nElectronegativity:\t%0.2f\nEnergy:\t\t\t%0.2f eV\nElectrons;  \nValency:\t\t\t%d\n", p->name, p->small, group(p), p->atomic_number, p->molar, p->atomic_radius, p->electronegativity, predict_valence_electron_energy(p), atom_valency(p));
    int c, w, n, l, m;
    char ns[4] = {'s', 'p', 'd', 'f'};
    char ls[4][5][6] = {   {"","","","",""},
                            {"x", "y", "z", "", ""},
                            {"xy", "yz", "xz", "z2", "x2-y2"},
                            {"","","","",""}};
    for (n = 0; n < 5; n++) {
        for (l = 0; l < 5; l++) {
            w = 0;
            for (m = 0; m < 5; m++) {
                if (p->shells[n][l][m] != 0) {
                    printf("\t%d%c(%s): %d (%0.2f eV)\n", n+1, ns[l], ls[l][m], p->shells[n][l][m], predict_electron_energy(p, n+1, l, m-l, c, 0));
                    w += p->shells[n][l][m];
                }
            }
            c += w;
        }
    }
    return;
}

int find_element(char * str) {
    int i;
    for (i = 0; i < pt_length; i++) {
        if (strcmp(str, periodic_table[i].small) == 0 || strcmp(str, periodic_table[i].name) == 0) {
            return i;
        }
    }
    return -1;
}


int is_element_cation(struct Element *p) {
    int outershell = get_valence_electrons(p);
    if (outershell < 8) {
        if (outershell < 4)
            return 1;
        else
            return 0;
    }
    return -1;
}
