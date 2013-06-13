/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io_rational.c: Process system and points file using rational arithmetic
 */
#include "blueharvest.h"

/***********************************************
 * read the polynomial/exponential system file *
 ***********************************************/
polynomial_system read_system_file_float(char *filename) {
    int i, res, num_var, num_poly;
    polynomial_system system;

    /* sanity-check the system file */
    errno = 0;
    FILE *sysfile = fopen(filename, "r");
    if (sysfile == NULL) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Couldn't open system file %s: %s.", filename, strerror(errno));

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADFILE);
    }

    /* gather number of variables and polynomials */
    errno = 0;
    res = fscanf(sysfile, "%d %d", &num_var, &num_poly);
    if (res == EOF || res == 0) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    system.numVariables = num_var;
    system.numPolynomials = num_poly;
    system.polynomials = malloc(num_poly * sizeof(polynomial));

    /* read in each polynomial piece by piece */

    for (i = 0; i < num_poly; i++) {
        system.polynomials[i] = parse_polynomial_float(sysfile, filename, num_var);
    }

    /* clean up the file */
    errno = 0;
    res = fclose(sysfile);

    if (res == EOF) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Couldn't close %s: %s.", filename, strerror(errno));

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADREAD);
    }

    return system;
}

/***********************************
* parse polynomial from input file *
************************************/
polynomial parse_polynomial_float(FILE *sysfile, char *filename, int num_var) {
    int num_terms, res, i, j;
    polynomial p;

    p.numVariables = num_var;

    /* first get the number of terms */
    res = fscanf(sysfile, "%d", &num_terms);

    p.exponents = malloc(num_terms * sizeof(int*));

    for (i = 0; i < num_terms; i++)
        p.exponents[i] = malloc(num_var * sizeof(int));
    p.coeff = malloc(num_terms * sizeof(rational_complex_number));

    if (res == EOF || res == 0) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    p.numTerms = num_terms;

    /* get the exponents and coefficients for each term */
    for (i = 0; i < num_terms; i++) {
        /* get exponent for each variable in the term */
        for (j = 0; j < num_var; j++) {
            res = fscanf(sysfile, "%d", &p.exponents[i][j]);
            if (res == EOF || res == 0) {
                char *error_string = malloc(80 * sizeof(char));
                snprintf(error_string, (size_t) 80, "Error reading %s: unexpected EOF", filename);

                print_error(error_string);
                free(error_string);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char *error_string = malloc(80 * sizeof(char));
                snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

                print_error(error_string);
                free(error_string);
                exit(BH_EXIT_BADPARSE);
            }

        }

        /* get (real & imag) coefficients for the term
         * takes extra parsing since these will be fractions */
        char str_coeff_real[50], str_coeff_imag[50];
        res = fscanf(sysfile, "%s %s", str_coeff_real, str_coeff_imag);

        if (res == EOF || res == 0) {
            char *error_string = malloc(80 * sizeof(char));
            snprintf(error_string, (size_t) 80, "Error reading %s: unexpected EOF", filename);

            print_error(error_string);
            free(error_string);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char *error_string = malloc(80 * sizeof(char));
            snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            free(error_string);
            exit(BH_EXIT_BADPARSE);
        }

        /* get mpq_t coefficients for real + imag */
        initialize_number(p.coeff[i]);
        parse_coeff_float(str_coeff_real, str_coeff_imag, p.coeff[i]);

    }

    return p;
}

/*************************************************************
 * parse string into numerator and denominator, return mpf_t *
 *************************************************************/
void parse_coeff_float(char *str_coeff_real, char *str_coeff_imag, complex_number c) {
    int res;

    errno = 0;

    res = mpf_set_str(c->re, str_coeff_real, 0);

    /* error! */
    if (res == -1) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_real);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    res = mpf_set_str(c->im, str_coeff_imag, 0);

    /* error! */
    if (res == -1) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_imag);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    }
}
