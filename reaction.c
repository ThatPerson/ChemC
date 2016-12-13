int predict_reaction(struct Mixture * p, struct Mixture * q) {
    int i, s;
    // Get all components in initial reaction mixture.
    p->num_constituents = 0;
    for (i = 0; i < p->num_compounds; i++) {
        find_constituents(&p->compounds[i]);

        for (s = 0; s < p->compounds[i].num_constituents; s++) {
            p->constituents[p->num_constituents] = p->compounds[i].constituents[s]; // Get copying function sorted. (actually pointers so shouldn't matter.
            p->constituents[p->num_constituents].valency = atom_valency(&p->constituents[s]);
            p->constituents[p->num_constituents].present = 1;
            p->num_constituents++;

        }
    }
    // loop from 1 to the number of constituent atoms.
  //   // Same scheme as predict bonding to make it easier to deal with. Essentially the algorithm is identical.
    int current_out = 0;
    int last_atom = get_electronegative_m(p, 0);
    int c_left = p->num_constituents;
    int current_element_in_compound = 0;
    int current_element_in_mixture = 0;
    for (i = 0; i < c_left; i++) {
        printf("%s ::::: %s\n", p->constituents[i].small, p->constituents[last_atom].small);
        if (strcmp(p->constituents[i].small, p->constituents[last_atom].small) == 0) {
            p->constituents[i].present = 0;
            q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[i];
            q->constituents[current_element_in_mixture] = p->constituents[i];
            current_element_in_mixture++;
            current_element_in_compound++;
            break;
        }
    }
    printf("%s %d\n", q->compounds[0].constituents[0].name, current_element_in_compound);
    c_left--;
    int qs=0;
    int current_atom, n = 0;
    int runs = 0;
    printf("====================== INTO LOOP %d======================\n", c_left);
    while (c_left > 0) {

        printf("NUM%d\n", c_left);
        qs = 0;
        n = (n==0)?1:0;
        current_atom = get_electronegative_m(p, n);
        for (i = 0; i < p->num_constituents; i++) {
            if (strcmp(p->constituents[i].small, p->constituents[current_atom].small) == 0 && p->constituents[i].present == 1) {
                qs = 1;
                c_left--;
                if (p->constituents[last_atom].valency >= p->constituents[current_atom].valency) {
                    p->constituents[last_atom].valency -= p->constituents[current_atom].valency;
                    q->constituents[current_element_in_mixture] = p->constituents[current_atom];
                    q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[current_atom];
                    printf("%s %d\n", q->compounds[current_out].constituents[current_element_in_compound].name, current_element_in_compound);
                    current_element_in_compound++;
                    current_element_in_mixture++;
                    p->constituents[current_out].valency = 0;
                    p->constituents[current_out].present = 0;
                } else {
                    if (p->constituents[last_atom].valency == 0) {
                        q->compounds[current_out].num_constituents = current_element_in_compound;
                        current_out++;
                        q->num_compounds = current_out;
                        last_atom = current_atom;
                        current_element_in_compound = 0;
                        printf("HELLO WORLD\n");
                        q->num_compounds++;
                    } else {
                        p->constituents[current_atom].valency -= p->constituents[last_atom].valency;
                        p->constituents[last_atom].valency = 0;
                        p->constituents[last_atom].present = 0;
                        q->constituents[current_element_in_mixture] = p->constituents[current_atom];
                        q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[current_atom];
                        last_atom = current_atom;
                    }
                }

            }
            /*if (qs == 1)
                break;*/
        }
    }
    q->compounds[current_out].num_constituents = current_element_in_compound;
    return 1;
}
