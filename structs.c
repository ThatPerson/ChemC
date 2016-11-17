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

struct Compound {
    struct Element constituents[50];
    char name[50];
    float molar;
    int num_constituents;
};
