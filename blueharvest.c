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
    int i, j, k, num_vars, num_poly;
    polynomial_system ps;
    rational_complex_vector *rational_vec = NULL;
    complex_vector *float_vec = NULL;
    void *generic_vec = NULL;

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

    set_arithmetic_type(rational_vec, float_vec, generic_vec);

    read_system_file(sysfile, &ps);
    num_vars = ps.numVariables;
    num_poly = ps.numPolynomials;

    int num_points = read_points_file(pointsfile, &generic_vec, num_vars);
    rational_vec = (rational_complex_vector *) generic_vec;

    printf("system ps has %d variables and %d polynomials\n", num_vars, num_poly);

    /* test system entered correctly; temporary */
    char variables[] = "xyzwuvabcdjkmnpqrs";
    for (i = 0; i < ps.numPolynomials; i++) {
        polynomial p = ps.polynomials[i];
        for (j = 0; j < p.numTerms - 1; j++) {
            gmp_printf("(%Qd + %Qdi)", p.coeff[j]->re, p.coeff[j]->im);
            for (k = 0; k < p.numVariables; k++)
                printf("%c^%d", variables[k], p.exponents[j][k]);
            printf(" + ");
        }
        gmp_printf("(%Qd + %Qdi)", p.coeff[j]->re, p.coeff[j]->im);
        for (k = 0; k < ps.numVariables; k++)
            printf("%c^%d", variables[k], p.exponents[j][k]);
        printf("\n");
    }

    /* test points entered correctly; temporary */
    printf("ps has %d test points\n", num_points);
    for (i = 0; i < num_points; i++) {
        for (j = 0; j < num_vars; j++) {
            gmp_printf("%Qd + %Qdi\n", rational_vec[i]->coord[j]->re, rational_vec[i]->coord[j]->im);
        }
    }

    /* clean up */
    void *system = &ps;
    free_system(system);
    free_vector(generic_vec, num_points);

    exit(0);
}

/*********************************************************************
 * free all dynamically-allocated variables in the polynomial system *
 *********************************************************************/
void free_system_float(void *system) {
    int i, j;
    polynomial_system *ps = (polynomial_system *) system;

    for (i = 0; i < (*ps).numPolynomials; i++) {
        polynomial p = (*ps).polynomials[i];
        int num_terms = p.numTerms;
        for (j = 0; j < num_terms; j++)
            free(p.exponents[j]);

        free(p.exponents);

        for (j = 0; j < num_terms; j++)
            clear_number(p.coeff[j]);

        free(p.coeff);
    }

    free((*ps).polynomials);
}


/******************************************************************************
 * free all dynamically-allocated variables in the rational polynomial system *
 ******************************************************************************/
void free_system_rational(void *system) {
    int i, j;
    polynomial_system *ps = (polynomial_system *) system;

    for (i = 0; i < (*ps).numPolynomials; i++) {
        polynomial p = (*ps).polynomials[i];
        int num_terms = p.numTerms;
        for (j = 0; j < num_terms; j++)
            free(p.exponents[j]);

        free(p.exponents);

        for (j = 0; j < num_terms; j++)
            clear_rational_number(p.coeff[j]);

        free(p.coeff);
    }

    free((*ps).polynomials);
}

/**************************************************************************
 * free all dynamically-allocated variables in the rational points vector *
 **************************************************************************/
void free_vector_rational(void *vec, int num_points) {
    rational_complex_vector *vector = (rational_complex_vector *) vec;
    int i;

    for (i = 0; i < num_points; i++)
        clear_rational_vector(vector[i]);

    free(vector);
}

/***********************************************************************
 * free all dynamically-allocated variables in the float points vector *
 ***********************************************************************/
void free_vector_float(void *vec, int num_points) {
    complex_vector *vector = (complex_vector *) vec;
    int i;

    for (i = 0; i < num_points; i++)
        clear_vector(vector[i]);

    free(vector);
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

                printf ("\n");
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
void set_arithmetic_type(rational_complex_vector *rational_vec, complex_vector *float_vec, void *generic_vec) {

    /* floating point functions and data types */
    if (arithmetic_type == BH_USE_FLOAT) {
        read_system_file = &read_system_file_float;
        read_points_file = &read_points_file_float;
        generic_vec = (void *) float_vec;
        free_system = &free_system_float;
        free_vector = &free_vector_float;
    }

    /* rational functions */
    else {
        read_system_file = &read_system_file_rational;
        read_points_file = &read_points_file_rational;
        generic_vec = (void *) rational_vec;
        free_system = &free_system_rational;
        free_vector = &free_vector_rational;
    }
}
