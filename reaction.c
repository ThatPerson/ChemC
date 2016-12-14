int number_remaining(struct Mixture * p) {
    int i, count = 0;
    for (i = 0; i < p->num_constituents; i++) {
        if (p->constituents[i].present == 1)
            count++;
    }
    return count;
}

int predict_reaction(struct Mixture * p, struct Mixture * q) {
    int i, s;
    // Get all components in initial reaction mixture.
    p->num_constituents = 0;
    for (i = 0; i < p->num_compounds; i++) {
        find_constituents(&p->compounds[i]);

        for (s = 0; s < p->compounds[i].num_constituents; s++) {
            p->constituents[p->num_constituents] = p->compounds[i].constituents[s]; // Get copying function sorted. (actually pointers so shouldn't matter.
            p->constituents[p->num_constituents].valency = atom_valency(&p->constituents[p->num_constituents]);
            p->constituents[p->num_constituents].present = 1;
            p->num_constituents++;

        }
    }
    // loop from 1 to the number of constituent atoms.
  //   // Same scheme as predict bonding to make it easier to deal with. Essentially the algorithm is identical.
    int current_out = 0;
    int n = 0;
    int last_atom = get_electronegative_m(p, n);
    int c_left = p->num_constituents;
    int current_element_in_compound = 0;
    int current_element_in_mixture = 0;
    for (i = 0; i < c_left; i++) {
        VERBOSE("%s ::::: %s\n", p->constituents[i].small, p->constituents[last_atom].small);
        if (strcmp(p->constituents[i].small, p->constituents[last_atom].small) == 0) {
            p->constituents[i].present = 0;
            q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[i];
            q->constituents[current_element_in_mixture] = p->constituents[i];
            current_element_in_mixture++;
            current_element_in_compound++;
            break;
        }
    }
    VERBOSE("%s %d\n", q->compounds[0].constituents[0].name, current_element_in_compound);
    c_left--;
    int current_atom;
    int reset = 0;
    VERBOSE("====================== INTO LOOP %d ======================\n", c_left);
    while (number_remaining(p) > 0) {

        VERBOSE("====================== %d ======================\n", number_remaining(p));
        n = (n==0)?1:0;
        current_atom = get_electronegative_m(p, n);
        for (i = 0; i < p->num_constituents; i++) {
            if (strcmp(p->constituents[i].small, p->constituents[current_atom].small) == 0 && p->constituents[i].present == 1) {
                VERBOSE("Adding %s to compound. Valency left = %d \n", p->constituents[current_atom].small, p->constituents[last_atom].valency);
                c_left--;
                if (reset == 1) {
                    // Add it regardless.
                    VERBOSE("RESET ALGO\n");
                    last_atom = current_atom;
                    q->compounds[current_out].constituents[0] = p->constituents[i];
                    p->constituents[i].present = 0;
                    current_element_in_compound++;
                    reset = 0;
                } else {

                    if (p->constituents[last_atom].valency >= p->constituents[i].valency) {
                        VERBOSE("ADDITION WITHOUT REMOVAL ALGO\n");
                        p->constituents[last_atom].valency -= p->constituents[i].valency;
                        VERBOSE("NEW VALENCY %d :::: %s %d\n", p->constituents[last_atom].valency, p->constituents[i].small, p->constituents[i].valency);
                        q->constituents[current_element_in_mixture] = p->constituents[i];
                        q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[i];
                        current_element_in_compound++;
                        current_element_in_mixture++;
                        p->constituents[i].valency = 0;
                        p->constituents[i].present = 0;
                    } else {
                        if (p->constituents[last_atom].valency == 0) {
                            VERBOSE("VALENCY 0\n");
                            q->compounds[current_out].num_constituents = current_element_in_compound;
                            p->constituents[last_atom].present = 0;
                            current_out++;
                            q->num_compounds = current_out + 1;
                            last_atom = i;
                            reset = 1;
                            n = 0;
                            current_element_in_compound = 0;
                            VERBOSE("Compound Finished.\n");
                            q->num_compounds++;
                          } else {
                            VERBOSE("REMOVAL ALGO\n");
                            p->constituents[i].valency -= p->constituents[last_atom].valency;
                            p->constituents[last_atom].valency = 0;
                            p->constituents[last_atom].present = 0;
                            q->constituents[current_element_in_mixture] = p->constituents[i];
                            q->compounds[current_out].constituents[current_element_in_compound] = p->constituents[i];
                            last_atom = i;
                          }
                      }
                }

                break;

            }
            /*if (qs == 1)
                break;*/
        }
    }
    q->num_compounds = current_out + 1;
    q->compounds[current_out].num_constituents = current_element_in_compound;
    return 1;
}
