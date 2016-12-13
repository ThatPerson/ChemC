int predict_reaction(struct Mixture * p, struct Mixture * q) {
    int i, s;
    // Get all components in initial reaction mixture.
    p->num_constituents = 0;
    for (i = 0; i < p->num_compounds; i++) {
        if (p->compounds[i].num_constituents == 0) {
            find_constituents(&p->compounds[i]);
            for (s = 0; s < p->compounds[i].num_constituents; s++) {
                p->constituents[s] = p->compounds[i].constituents[s]; // Get copying function sorted. (actually pointers so shouldn't matter.
                p->constituents[s].valency = atom_valency(&p->constituents[s]);
                p->constituents[s].present = 1;
                p->num_constituents++;
            }
        }
    }

    // loop from 1 to the number of constituent atoms.
    int c_left = q->num_constituents; // Same scheme as predict bonding to make it easier to deal with. Essentially the algorithm is identical.
    int current_out = 0;
    int last_atom = get_electronegative_m(p, 0);
    int current_element_in_compound = 0;
    int current_element_in_mixture = 0;
    for (i = 0; i < c_left; i++) {
        if (strcmp(p->constituents[i].small, p->constituents[last_atom].small) == 0) {
            p->constituents[i].present = 0;
            q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[i]; /* Need a function to copy the struct */
            q->constituents[current_element_in_mixture] = p->constituents[i];
            current_element_in_mixture++;
            current_element_in_compound++;
            break;
        }
    }
    c_left--;
    int qs=0;
    int current_atom, n = 0;
    while (c_left > 0) {
        qs = 0;
        n = (n==0)?1:0;
        current_atom = get_electronegative_m(p, n);
        for (i = 0; i < p->num_constituents; i++) {
            if (strcmp(p->constituents[i].small, p->constituents[current_atom].small) == 0 && p->constituents[i].present == 1) {
                qs = 1;
                if (p->constituents[last_atom].valency >= p->constituents[current_atom].valency) {
                    p->constituents[last_atom].valency -= p->constituents[current_atom].valency;
                    q->constituents[current_element_in_mixture] = p->constituents[current_atom];
                    q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[current_atom];
                    current_element_in_compound++;
                    current_element_in_mixture++;
                    p->constituents[current_out].valency = 0;
                    p->constituents[current_out].present = 0;
                } else {
                    if (p->constituents[last_atom].valency == 0) {
                        current_out++;
                        last_atom = current_atom;
                        current_element_in_compound = 0;
                    } else {
                        p->constituents[current_atom].valency -= p->constituents[last_atom].valency;
                        p->constituents[last_atom].valency = 0;
                        p->constituents[last_atom].present = 0;
                        q->constituents[current_element_in_mixture] = p->constituents[current_atom];
                        q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[current_atom];
                        last_atom = current_atom;
                        c_left--;
                    }
                }

            }
            if (qs == 1)
                break;
        }
    }

    return 1;
}
