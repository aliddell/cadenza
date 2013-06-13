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

    prog_info(stderr);

    /***********************************************************
     * get command-line arguments before anything else happens *
     ***********************************************************/

    getargs(argc, argv);

    /*****************************************
     * ensure filenames are properly defined *
     *****************************************/
    if (sysfile == NULL) {
        print_error("You need to define a system file.");
        usage();
        exit(BH_EXIT_BADFILE);
    } else if (pointsfile == NULL) {
        print_error("You need to define a points file.");
        usage();
        exit(BH_EXIT_BADFILE);
    }

    if (help_flag) {
        usage();
        exit(0);
    }

    if (verbose_flag == 1)
        display_config();

    polynomial_system ps = read_system_file(sysfile);

    printf("system ps has %d variables and %d polynomials\n", ps.numVariables, ps.numPolynomials);

    /* test system entered correctly */
    char variables[] = "xyzwuvabcdjkmnpqrs";
    for (i = 0; i < ps.numPolynomials; i++) {
        polynomial p = ps.polynomials[i];
        for (j = 0; j < p.numTerms - 1; j++) {
            gmp_printf("(%Q + %Qi)", p.coeff[j]->re, p.coeff[j]->im);
            for (k = 0; k < p.numVariables; k++)
                printf("%c^%d", variables[j], p.exponents[j][k]);
            printf(" + ");
        }
        gmp_printf("(%Q + %Qi)", p.coeff[j]->re, p.coeff[j]->im);
        for (k = 0; k < ps.numVariables; k++)
            printf("%c^%d", variables[j], p.exponents[j][k]);
        printf("\n");
    }

    free_system(ps);

    exit(0);
}

/*********************************************************************
 * free all dynamically-allocated variables in the polynomial system *
 *********************************************************************/
void free_system(polynomial_system ps) {
    int i, j;

    for (i = 0; i < ps.numPolynomials; i++) {
        polynomial p = ps.polynomials[i];
        int num_vars = p.numVariables;
        int num_terms = p.numTerms;
        for (j = 0; j < num_terms; j++)
            free(p.exponents[j]);

        free(p.exponents);

        if (arithmetic_type == BH_USE_FLOAT)
            for (j = 0; j < num_terms; j++)
                mpf_clear(p.coeff[j]);
        else
            for (j = 0; j < num_terms; j++)
                mpq_clear(p.coeff[j]);

        free(p.coeff);
    }

    free(ps.polynomials);
}

/*******************************************************
 * parse command-line args and set flags and filenames *
 *******************************************************/
void getargs(int argc, char *argv[]) {
    /* define default values */
    help_flag = 0;
    verbose_flag = 0;
    arithmetic_type = BH_USE_RATIONAL;
    default_precision = MPFR_PREC_MIN;
    sysfile = NULL;
    pointsfile = NULL;

    int c = 0;

    while (c != -1) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"verbose", no_argument,      &verbose_flag, 1},
            {"help",   no_argument,       &help_flag, 1},
            {"float",     no_argument,    &arithmetic_type, BH_USE_FLOAT},
            {"rational",  no_argument,    &arithmetic_type, BH_USE_RATIONAL},
            /* These options don't set a flag. */
            {"points",    required_argument, 0, 'p'},
            {"system",    required_argument, 0, 's'},
            {"precision", required_argument, 0, 'm'},
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
                verbose_flag = 1;
                break;
        }
    }
}
