

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
    int present; // Only relevant to prediction algorithms.
};

struct Compound {
    struct Element constituents[50];
    char name[50];
    float molar;
    int num_constituents;
    int present;
};

struct Bond {
    struct Element atoms[2];
    int num_bonds; /* Number of electrons involved in bond? Could then get bond order? */
    int num_electrons;
    int bond_order;
};

struct Mixture {
    struct Compound compounds[50];
    int num_compounds;
    int num_constituents;
    struct Element constituents[50];
};

/* Ignore below */
int hideme ( const char * format, ... ) {
  return 1;
}
#define VERBOSE hideme
