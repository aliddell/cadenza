/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * blueharvest.c: Main file for Blue Harvest
 */
#include "blueharvest.h"

int main(int argc, char *argv[]) {
    int i, num_var, num_poly, num_points;

    /* polynomial system */
    polynomial_system F;
    
    /* settings (for aC functions) */
    configurations S;

    /* constant vector v_i */
    rational_complex_vector v_rational;
    complex_vector v_float;
    void *v = NULL;

    /* t_i */
    mpq_t *t_rational = NULL;
    mpf_t *t_float = NULL;
    void *t = NULL;

    /* vector w_i */
    rational_complex_vector *w_rational = NULL;
    complex_vector *w_float = NULL;
    void *w = NULL;

    /* init random seed */
    srand(time(NULL));

    /* get command-line arguments before anything else happens */
    getargs(argc, argv);
    if (verbosity > BH_LACONIC)
        prog_info();

    /* do this before checking filenames, duh */
    if (help_flag) {
        usage();
        exit(BH_EXIT_SUCCESS);
    }

    /* ensure filenames are properly defined */
    if (sysfile == NULL) {
        print_error("You need to define a system file.");
        usage();
        exit(BH_EXIT_BADFILE);
    } else if (pointsfile == NULL) {
        print_error("You need to define a points file.");
        usage();
        exit(BH_EXIT_BADFILE);
    }

    if (verbosity > BH_CHATTY)
        display_config();

    /* set void and function pointers */
    set_function_pointers();
    if (arithmetic_type == BH_USE_FLOAT) {
        v = (void *) &v_float;
    } else {
        v = (void *) &v_rational;
    }

    /* set configurations */
    S.arithmeticType = arithmetic_type;
    S.startingPrecision = default_precision;
    S.algorithm = 1; /* real distinct certify */

    read_system_file(sysfile, &F, v);
    num_var = F.numVariables;
    num_poly = F.numPolynomials;

    num_points = read_points_file(pointsfile, &t, &w, num_var);

    /* test system entered correctly; temporary */
    w_rational = (rational_complex_vector *) w;
    t_rational = (mpq_t *) t;

    print_system_rational(&F);

    /* test vector entered correctly; temporary */
    puts("");
    printf("[");
    for (i = 0; i < num_var - 1; i++) {
        gmp_printf("%Qd + %Qdi, ", v_rational->coord[i]->re, v_rational->coord[i]->im);
    }
    gmp_printf("%Qd + %Qdi]\n", v_rational->coord[i]->re, v_rational->coord[i]->im);

    /* test points entered correctly; temporary */
    puts("");
    for (i = 0; i < num_points; i++) {
        gmp_printf("%Qd", t_rational[i]);
        print_points_rational(&w_rational[i], num_var);
    }

    test_system(&F, &S, v, t, w, num_points);

    /* clean up */
    free_system((void *) &F, v);
    free_vector(w, t, num_points);

    exit(BH_EXIT_SUCCESS);
}

/******************************************************************************
 * free all dynamically-allocated variables in the rational polynomial system *
 ******************************************************************************/
void free_system_rational(void *system, void *v) {
    int i, j;
    polynomial_system *F = (polynomial_system *) system;
    rational_complex_vector *v_rational = (rational_complex_vector *) v;

    for (i = 0; i < F->numPolynomials; i++) {
        polynomial p = F->polynomials[i];
        int num_terms = p.numTerms;
        for (j = 0; j < num_terms; j++)
            free(p.exponents[j]);

        free(p.exponents);

        for (j = 0; j < num_terms; j++)
            clear_rational_number(p.coeff[j]);

        free(p.coeff);
        mpq_clear(p.norm_sqr);
    }

    mpq_clear(F->norm_sqr);
    free(F->polynomials);
    clear_rational_vector(*v_rational);
}

/*********************************************************************
 * free all dynamically-allocated variables in the polynomial system *
 *********************************************************************/
void free_system_float(void *system, void *v) {
    int i, j;
    polynomial_system *F = (polynomial_system *) system;
    complex_vector *v_float = (complex_vector *) v;

    for (i = 0; i < F->numPolynomials; i++) {
        polynomial p = F->polynomials[i];
        int num_terms = p.numTerms;
        for (j = 0; j < num_terms; j++)
            free(p.exponents[j]);

        free(p.exponents);

        for (j = 0; j < num_terms; j++)
            clear_number(p.coeff[j]);

        mpq_clear(p.norm_sqr);
        free(p.coeff);
    }

    mpq_clear(F->norm_sqr);
    free(F->polynomials);
    clear_vector(*v_float);
}

/**************************************************************************
 * free all dynamically-allocated variables in the rational points vector *
 **************************************************************************/
void free_vector_rational(void *w, void *t, int num_points) {
    int i;
    rational_complex_vector *w_rational = (rational_complex_vector *) w;
    mpq_t *t_rational = (mpq_t *) t;

    for (i = 0; i < num_points; i++)
        mpq_clear(t_rational[i]);
    for (i = 0; i < num_points; i++)
        clear_rational_vector(w_rational[i]);

    free(t_rational);
    free(w_rational);
}

/***********************************************************************
 * free all dynamically-allocated variables in the float points vector *
 ***********************************************************************/
void free_vector_float(void *w, void *t, int num_points) {
    int i;

    complex_vector *w_float = (complex_vector *) w;
    mpf_t *t_float = (mpf_t *) t;

    for (i = 0; i < num_points; i++)
        mpf_clear(t_float[i]);
    for (i = 0; i < num_points; i++)
        clear_vector(w_float[i]);

    free(t_float);
    free(w_float);
}

/*******************************************************
 * parse command-line args and set flags and filenames *
 *******************************************************/
void getargs(int argc, char *argv[]) {
    /* define default values */
    help_flag = 0;
    verbosity = BH_LACONIC;
    arithmetic_type = BH_USE_RATIONAL;
    default_precision = MPFR_PREC_MIN;
    sysfile = NULL;
    pointsfile = NULL;

    int c = 0;

    while (c != -1) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"chatty",     no_argument, &verbosity, BH_CHATTY},
            {"verbose",    no_argument, &verbosity, BH_VERBOSE},
            {"loquacious", no_argument, &verbosity, BH_LOQUACIOUS},
            {"help",       no_argument, &help_flag, 1},
            {"float",      no_argument, &arithmetic_type, BH_USE_FLOAT},
            {"rational",   no_argument, &arithmetic_type, BH_USE_RATIONAL},
            /* These options don't set a flag. */
            {"points",     required_argument, 0, 'p'},
            {"system",     required_argument, 0, 's'},
            {"precision",  required_argument, 0, 'm'},
            {0, 0, 0, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvfqp:s:", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0) break;

                printf ("option %s", long_options[option_index].name);

                if (optarg) printf (" with arg %s", optarg);

                puts("\n");
                break;

            case 'f':
                arithmetic_type = BH_USE_FLOAT;
                break;

            case 'h':
                help_flag = 1;
                break;

            case 'm':
                default_precision = atoi(optarg);
                break;

            case 'p':
                pointsfile = optarg;
                break;

            case 'q':
                arithmetic_type = BH_USE_RATIONAL;
                break;

            case 's':
                sysfile = optarg;
                break;

            case 'v':
                verbosity++;
                if (verbosity > BH_LOQUACIOUS)
                    verbosity = BH_LOQUACIOUS;
                break;
        }
    }
}

/**********************************************************************************
 * set function pointers to use the proper arithmetic type with a minimum of fuss *
 **********************************************************************************/
void set_function_pointers() {

    /* floating point functions and data types */
    if (arithmetic_type == BH_USE_FLOAT) {
        read_system_file = &read_system_file_float;
        read_points_file = &read_points_file_float;
        test_system = &test_system_float;
        free_system = &free_system_float;
        free_vector = &free_vector_float;
    }

    /* rational functions */
    else {
        read_system_file = &read_system_file_rational;
        read_points_file = &read_points_file_rational;
        test_system = &test_system_rational;
        free_system = &free_system_rational;
        free_vector = &free_vector_rational;
    }
}
