/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io.c: Input/output functions for Blue Harvest
 */
#include "blueharvest.h"

/********************************************************
 * print a helpful message about command-line arguments *
 ********************************************************/
void usage() {
    fprintf(stderr, "usage: blueharvest [-h|--help] [-v|--verbose] [-f|--float] [-q|--rational] --precision INT --system FILENAME --points FILENAME\n");
}

/********************************************
 * uniformly print error messages to stderr *
 ********************************************/
void print_error(char *msg) {
    fprintf(stderr, "\nERROR: %s\n\n", msg);
}

/**********************************************
 * print a helpful message about this program *
 **********************************************/
void prog_info(FILE *OUT) {
    char *compile_date = malloc((BH_MAX_DATECHAR+1) * sizeof(char));
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    size_t res = strftime(compile_date, (size_t) BH_MAX_DATECHAR, "%b %d, %Y", timeinfo);
    if (res == 0)
        compile_date = "Jan 1, 1970";

    fprintf(OUT, "\n");
    fprintf(OUT, "\t%s v%s (built %s) (compiled %s)\n", BH_PROGRAM_NAME, BH_VERSION, BH_BUILD_DATE, compile_date);
    fprintf(OUT, "\t%s\n", BH_AUTHORS);
    fprintf(OUT, "\tGMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());

    free(compile_date);
}

void display_config() {
    char str_arithmetic_type[15];
    if (arithmetic_type == BH_USE_FLOAT)
        strcpy(str_arithmetic_type, "floating point");
    else
        strcpy(str_arithmetic_type, "rational");

    printf("\tComputing using %s arithmetic\n", str_arithmetic_type);
    printf("\tPolynomial system file is %s\n", sysfile);
    printf("\tPoint set file is %s\n", pointsfile);
}

/***********************************************
 * read the polynomial/exponential system file *
 ***********************************************/
polynomial_system read_system_file(char *filename) {
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
    if (res == EOF) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: EOF", filename);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_READERR);
    } else if (errno == EILSEQ) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_READERR);
    }

    system.numVariables = num_var;
    system.numPolynomials = num_poly;
    system.polynomials = malloc(num_poly * sizeof(polynomial));

    /* read in each polynomial piece by piece */

    for (i = 0; i < num_poly; i++) {
        system.polynomials[i] = parse_polynomial(sysfile, filename, num_var);
    }

    return system;
}

/***********************************
* parse polynomial from input file *
************************************/
polynomial parse_polynomial(FILE *sysfile, char *filename, int num_var) {
    int num_terms, res, i, j;
    polynomial p;

    p.numVariables = num_var;
    p.exponents = malloc(num_terms * sizeof(int*));
    for (i = 0; i < num_terms; i++)
        p.exponents[i] = malloc(num_var * sizeof(int));
    p.coeff = malloc(num_terms * sizeof(rational_complex_number));

    /* first get the number of terms */
    res = fscanf(sysfile, "%d", &num_terms);

    if (res == EOF) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: EOF", filename);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_READERR);
    } else if (errno == EILSEQ) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_READERR);
    }

    p.numTerms = num_terms;

    /* get the exponents and coefficients for each term */
    for (i = 0; i < num_terms; i++) {
        /* get exponent for each variable in the term */
        for (j = 0; j < num_var; j++) {
            res = fscanf(sysfile, "%d", &p.exponents[i][j]);
            if (res == EOF) {
                char *error_string = malloc(80 * sizeof(char));
                snprintf(error_string, (size_t) 80, "Error reading %s: EOF", filename);

                print_error(error_string);
                free(error_string);
                exit(BH_EXIT_READERR);
            } else if (errno == EILSEQ) {
                char *error_string = malloc(80 * sizeof(char));
                snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

                print_error(error_string);
                free(error_string);
                exit(BH_EXIT_READERR);
            }

        }

        /* get (real & imag) coefficients for the term
         * takes extra parsing since these can be fractions */
        char *str_coeff_real, *str_coeff_imag;
        res = fscanf(sysfile, "%s %s", str_coeff_real, str_coeff_imag);

        if (res == EOF) {
            char *error_string = malloc(80 * sizeof(char));
            snprintf(error_string, (size_t) 80, "Error reading %s: EOF", filename);

            print_error(error_string);
            free(error_string);
            exit(BH_EXIT_READERR);
        } else if (errno == EILSEQ) {
            char *error_string = malloc(80 * sizeof(char));
            snprintf(error_string, (size_t) 80, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            free(error_string);
            exit(BH_EXIT_READERR);
        }

        if (arithmetic_type == BH_USE_FLOAT) {
            /* get mpf_t coefficients for real + imag */
            initialize_number(p.coeff[i]);
            parse_coeff_float(str_coeff_real, str_coeff_imag, &p.coeff[i]);
        } else {
            /* get mpq_t coefficients for real + imag */
            initialize_rational_number(p.coeff[i]);
            parse_coeff_rational(str_coeff_real, str_coeff_imag, &p.coeff[i]);
        }

    }

    return p;
}

/*************************************************************
 * parse string into numerator and denominator, return mpq_t *
 *************************************************************/
void parse_coeff_rational(char *str_coeff_real, char *str_coeff_imag, rational_complex_number *c) {
    long int num;
    unsigned long int denom;
    int res;
    char *tok;

    errno = 0;

    tok = strtok(str_coeff_real, "/");
    num = strtol(tok, NULL, 0);
    if (num == 0) {
        mpq_set_ui(c[0]->re, 0, 1);
    } else if (errno == ERANGE) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_real);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    } else {
        tok = strtok(NULL, "/");
        if (tok == NULL)
            denom = 1;
        else {
            denom = strtoul(tok, NULL, 0);
            if (denom == 0) {
                print_error("Division by zero in system file.");
                exit(BH_EXIT_BADPARSE);
            } else if (errno == ERANGE) {
                char *error_string = malloc(80 * sizeof(char));
                snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_real);

                print_error(error_string);
                free(error_string);
                exit(BH_EXIT_BADPARSE);
            }
        }
        mpq_set_si(c[0]->re, num, denom);
    }

    tok = strtok(str_coeff_imag, "/");
    num = strtol(tok, NULL, 0);
    if (num == 0) {
        mpq_set_ui(c[0]->im, 0, 1);
    } else if (errno == ERANGE) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_imag);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    } else {
        tok = strtok(NULL, "/");
        if (tok == NULL)
            denom = 1;
        else {
            denom = strtoul(tok, NULL, 0);
            if (denom == 0) {
                print_error("Division by zero in system file.");
                exit(BH_EXIT_BADPARSE);
            } else if (errno == ERANGE) {
                char *error_string = malloc(80 * sizeof(char));
                snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_imag);

                print_error(error_string);
                free(error_string);
                exit(BH_EXIT_BADPARSE);
            }
        }
        mpq_set_si(c[0]->im, num, denom);
    }
}

/*************************************************************
 * parse string into numerator and denominator, return mpf_t *
 *************************************************************/
void parse_coeff_float(char *str_coeff_real, char *str_coeff_imag, complex_number *c) {
    int res;

    errno = 0;

    res = mpf_set_str(c[0]->re, str_coeff_real, 0);

    /* error! */
    if (res == -1) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_real);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    res = mpf_set_str(c[0]->im, str_coeff_imag, 0);

    /* error! */
    if (res == -1) {
        char *error_string = malloc(80 * sizeof(char));
        snprintf(error_string, (size_t) 80, "Invalid coefficient: %s", str_coeff_imag);

        print_error(error_string);
        free(error_string);
        exit(BH_EXIT_BADPARSE);
    }
}
