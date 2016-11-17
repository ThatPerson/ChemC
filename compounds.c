float predict_melting_point(struct Compound * q, int algorithm) {
    if (algorithm == 0) {
        int i;
        int Z[3] = {0, 0, 0}, N[2] = {0, 0};
        for (i = 0; i < arrlen(q->constituents); i++) {
            Z[1+is_element_cation(q->constituents[i])] = q->constituents[i]->atomic_number;
            N[1+is_element_cation(q->constituents[i])] = group(q->constituents[i]);
        }
        if (N[1] == 0 || N[2] == 0)
            return -1;
        float melting_point = abs(RYDBERG * (pow(Z[2]/N[2], 2) - pow(Z[1]/N[1], 2)));
        melting_point = -0.101334 * melting_point + 1109.81;
        return melting_point;
    } else {
    }
}
