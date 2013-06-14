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
    int i, j, k;

    /* get command-line arguments before anything else happens */

    getargs(argc, argv);
    if (verbosity > BH_TACITURN)
        prog_info();

    /* do this before checking filenames, duh */
    if (help_flag) {
        usage();
        exit(0);
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

    set_function_pointers();

    polynomial_system ps = read_system_file(sysfile);

    printf("system ps has %d variables and %d polynomials\n", ps.numVariables, ps.numPolynomials);

    /* test system entered correctly */
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

    void *system = &ps;

    free_system(system);

    exit(0);
}

/*********************************************************************
 * free all dynamically-allocated variables in the polynomial system *
 *********************************************************************/
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

/*******************************************************
 * parse command-line args and set flags and filenames *
 *******************************************************/
void getargs(int argc, char *argv[]) {
    /* define default values */
    help_flag = 0;
    verbosity = BH_TACITURN;
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
void set_function_pointers() {

    /* floating point functions */
    if (arithmetic_type == BH_USE_FLOAT) {
        read_system_file = &read_system_file_float;
        free_system = &free_system_float;
    }

    /* rational functions */
    else {
        read_system_file = &read_system_file_rational;
        free_system = &free_system_rational;
    }
}
