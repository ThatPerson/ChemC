float predict_melting_point(struct Compound * q, int algorithm) {
    float melting_point = 0;
    if (algorithm == 0) {
        int i;
        int Z[3] = {0, 0, 0}, N[2] = {0, 0};
        for (i = 0; i < arrlen(q->constituents); i++) {
            Z[1+is_element_cation(&q->constituents[i])] = q->constituents[i].atomic_number;
            N[1+is_element_cation(&q->constituents[i])] = group(&q->constituents[i]);
        }
        if (N[1] == 0 || N[2] == 0)
            return -1;
        melting_point = abs(RYDBERG * (pow(Z[2]/N[2], 2) - pow(Z[1]/N[1], 2)));
        melting_point = -0.101334 * melting_point + 1109.81;
    } else {
        int n[3] = {0, 0, 0}, c[3] = {0, 0, 0}, i;
        float d[3] = {0, 0, 0};
        for (i = 0; i < arrlen(q->constituents); i++) {
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

