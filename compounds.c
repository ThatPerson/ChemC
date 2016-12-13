float predict_melting_point(struct Compound * q, int algorithm) {
    float melting_point = 0;
    if (algorithm == 0) {
        int i;
        int Z[3] = {0, 0, 0}, N[3] = {0, 0, 0};
        for (i = 0; i < q->num_constituents; i++) {
            Z[1+is_element_cation(&q->constituents[i])] = q->constituents[i].atomic_number;
            N[1+is_element_cation(&q->constituents[i])] = group(&q->constituents[i]);
        }
        if (N[1] == 0 || N[2] == 0)
            return -1;
        melting_point = fabs(RYDBERG * (pow(Z[2]/N[2], 2) - pow(Z[1]/N[1], 2)));
        melting_point = -0.101334 * melting_point + 1109.81;
    } else {
        int n[3] = {0, 0, 0}, c[3] = {0, 0, 0}, i;
        float d[3] = {0, 0, 0};
        for (i = 0; i < q->num_constituents; i++) {
            int tmp = 1 + is_element_cation(&q->constituents[i]);
            n[tmp]++;
            c[tmp] += atom_valency(&q->constituents[i]);
            d[tmp] += q->constituents[i].atomic_radius;
        }
        int num = 2;
        double X = 0;
        double cat = 1.3 - (0.3 * c[2]);
        double dist = (d[1] + d[2]) * pow(10, -12);
        if (dist != 0 && n[1] != 0) {
            X = ((MP_K*num*cat*c[1])/(dist*n[1]) * (1-(MP_D/dist)));
            melting_point = 0.00148848 * X + 1.0007;
        } else {
            return -1;
        }
    }
    return melting_point;
}

/* Find constituents function - for now just don't */

int find_constituents(struct Compound *q) {
    /**
     * Sometimes in life it's best just not to ask.
     */
    int multiplicative_factor = 1;
    int i, qs = 0;
    if (48 <= q->name[0] && q->name[0] <= 57) {
        qs = 1;
        multiplicative_factor = (q->name[0] - 48);
    }

    char constituents[50][3];
    int current = 0,r;
    char buffer[3] = {' ',' ',' '};
    int current_buf = 0;
    int multiplt = 0;
    for (i = qs; i < strlen(q->name); i++) {
        if (48 <= q->name[i] && q->name[i] <= 57) {
            multiplt = (10 * multiplt) + (q->name[i] - 48);
        }
        if (current_buf > 2) {
            current_buf = 0;
        }
        if (65 <= q->name[i] && q->name[i] <= 90) {
            if (strlen(buffer) > 0) {
                if (multiplt == 0)
                    multiplt = 1;
                for (r = 0; r < multiplt; r++) {
                    strcpy(constituents[current], buffer);
                    current++;
                }
                multiplt = 0;
            }
            current_buf = 0;
            buffer[0] = 0;
            buffer[1] = 0;
            buffer[2] = 0;
            buffer[current_buf] = q->name[i];
            current_buf++;
        } else if (97 <= q->name[i] && q->name[i] <= 122) {
            buffer[current_buf] = q->name[i];
            current_buf++;
        }

    }

    if (multiplt == 0)
        multiplt = 1;
    for (r = 0; r < multiplt; r++) {
        strcpy(constituents[current], buffer);
        current++;
    }
    int n = 0;
    for (r = 0; r < multiplicative_factor; r++) {
        for (i = 1; i < current; i++) {
            q->constituents[n] = periodic_table[find_element(constituents[i])];
            n++;
        }
    }
    q->num_constituents = n;
    return n;
}

char * find_name(struct Compound * q) {
    char chemicals[50][3];
    int number_of[50];
    static char output[500] = "";
    int i, n, p, c = 0;
    for (i = 0; i < 50; i++) {
        number_of[i] = 0;
    }
    for (i = 0; i < q->num_constituents; i++) {
        n = 0;
        for (p = 0; p < 50; p++) {
            if (strcmp(chemicals[p], q->constituents[i].small) == 0) {
                n = 1;
                number_of[p]++;
            }
        }
        if (n == 0) {
            strcpy(chemicals[c], q->constituents[i].small);
            number_of[c] = 1;
            c++; // ay
        }
    }

    strcpy(output, "");
    for (i = 0; i < sizeof(chemicals)/sizeof(chemicals[0]); i++) {
        if (number_of[i] > 1)
            sprintf(output, "%s%s%d", output, chemicals[i], number_of[i]);
        else if (number_of[i] == 1)
            sprintf(output, "%s%s", output, chemicals[i]);
    }
    return output;
}

int get_electronegative(struct Compound *q, int type) {
    float current_lowest = 0;
    int current_atom = 0, i;
    if (q->num_constituents > 0) {
        current_lowest = q->constituents[0].electronegativity;
        current_atom = 0;
    } else {
        return -1; // No atoms
    }
    for (i = 0; i < q->num_constituents; i++) {
        if (type == 1) {
            if (q->constituents[i].electronegativity > current_lowest && q->constituents[i].present == 1) {
                current_lowest = q->constituents[i].electronegativity;
                current_atom = i;
            }
        } else {
            if (q->constituents[i].electronegativity < current_lowest && q->constituents[i].present == 1) {
                current_lowest = q->constituents[i].electronegativity;
                current_atom = i;
            }
        }
    }
    return current_atom;
}

int get_electronegative_m(struct Mixture *q, int type) {
    float current_lowest = 0;
    int current_atom = 0, i;
    if (q->num_constituents > 0) {
        current_lowest = q->constituents[0].electronegativity;
        current_atom = 0;
    } else {
        return -1; // No atoms
    }
    for (i = 0; i < q->num_constituents; i++) {
      VERBOSE("%d %s    %f  %d\n", q->constituents[i].present, q->constituents[i].name, q->constituents[i].electronegativity, i);
        if (type == 1) {

            if (q->constituents[i].electronegativity > current_lowest && q->constituents[i].present == 1) {
                current_lowest = q->constituents[i].electronegativity;
                current_atom = i;
            }
        } else {
            if (q->constituents[i].electronegativity < current_lowest && q->constituents[i].present == 1) {
                current_lowest = q->constituents[i].electronegativity;
                current_atom = i;
            }
        }
    }
    return current_atom;
}

int predict_bonding(struct Bond * bonds, int length, struct Compound * q) {
    int c_left = q->num_constituents;
    int i;
    for (i = 0; i < q->num_constituents; i++) {
        q->constituents[i].valency = atom_valency(&q->constituents[i]);
        q->constituents[i].present = 1;
    }
    int current_bond = 0, current_atom;
    int last_atom = get_electronegative(q, 0), p = 0;
    for (i = 0; i < c_left; i++) {
        if (strcmp(q->constituents[i].small, q->constituents[last_atom].small) == 0 && p == 0) {
            q->constituents[i].present = 0;
            p = 1;
        } else {
            q->constituents[i].present = 1;
        }
    }
    c_left--;
    int n = 0;
    while (c_left > 0) {
        n = (n == 0) ? 1 : 0;
        current_atom = get_electronegative(q, n);
        for (i = 0; i < q->num_constituents; i++) {
            if (strcmp(q->constituents[i].small, q->constituents[current_atom].small) == 0 && q->constituents[i].present == 1) {
                if (q->constituents[last_atom].valency >= q->constituents[current_atom].valency) {
                    q->constituents[last_atom].valency -= q->constituents[current_atom].valency;
                    bonds[current_bond].atoms[0] = q->constituents[last_atom];
                    bonds[current_bond].atoms[1] = q->constituents[current_atom];
                    bonds[current_bond].num_bonds = q->constituents[current_atom].valency;
                    current_bond++;
                    q->constituents[current_atom].valency = 0;
                    q->constituents[i].present = 0;
                    c_left--;
                } else if (q->constituents[last_atom].valency < q->constituents[current_atom].valency) {
                    bonds[current_bond].atoms[0] = q->constituents[last_atom];
                    bonds[current_bond].atoms[1] = q->constituents[current_atom];
                    bonds[current_bond].num_bonds = q->constituents[current_atom].valency - q->constituents[last_atom].valency;
                    current_bond++;
                    q->constituents[i].present = 0;
                    c_left--;
                    last_atom = current_atom;
                }
            }
        }
    }
    return current_bond - 1;
}

int predict_ir(float * ir, struct Bond * bonds, int num_bonds) {
    float reduced_mass, tmp;
    int i, p, j, current_id;
    for (i = 0; i < 50; i++)
        ir[i] = 0;
    for (i = 0; i < num_bonds; i++) {
        reduced_mass = (bonds[i].atoms[0].molar * bonds[i].atoms[1].molar) / (bonds[i].atoms[0].molar + bonds[i].atoms[1].molar);
        tmp = 4.12*sqrt((5*pow(10, 5) * bonds[i].num_bonds)/reduced_mass);
        j = 0;
        for (p = 0; p < current_id; p++) {
            if (ir[p] == tmp) {
                j = 1;
            }
        }
        if (j == 0) {
            ir[current_id] = tmp;
            current_id++;
        }
    }
    return current_id;
}


float compound_molarity(struct Compound *q) {
    int i;
    float molar = 0;
    for (i = 0; i < q->num_constituents; i++) {
        molar += q->constituents[i].molar;
    }
    return molar;
}



/*
Char	Dec	Hex	Oct
A	65	41	101
Z	90	5A	132
13	0D	15
10	0A	12
a	97	61	141
z	122	7A	172
13	0D	15
10	0A	12
0	48	30	60
9	57	39	71
*/
