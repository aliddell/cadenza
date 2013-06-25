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
void read_system_file_rational(char *filename, polynomial_system *system, void *v) {
    int i, res, num_var, num_poly;
    rational_complex_vector *rational_v = (rational_complex_vector *) v;

    /* sanity-check the system file */
    errno = 0;
    FILE *sysfile = fopen(filename, "r");
    if (sysfile == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Couldn't open system file %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADFILE);
    }

    /* gather number of variables */
    errno = 0;
    res = fscanf(sysfile, "%d", &num_var);
    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    num_poly = num_var; /* square system */

    system->numVariables = num_var;
    system->numPolynomials = num_poly;
    system->maximumDegree = 0;
    system->isReal = 0;
    mpq_init(system->norm_sqr);
    system->numExponentials = 0;
    system->polynomials = malloc(num_poly * sizeof(polynomial));

    /* read in each polynomial piece by piece */
    for (i = 0; i < num_poly; i++) {
        system->polynomials[i] = parse_polynomial_rational(sysfile, filename, num_var);
        if (system->polynomials[i].degree > system->maximumDegree)
            system->maximumDegree = system->polynomials[i].degree;
    }

    /* read in v */
    initialize_rational_vector(*rational_v, num_var);

    for (i = 0; i < num_var; i++) {
        /* get (real & imag) parts for each v_i */
        errno = 0;
        char str_real[50], str_imag[50];
        res = fscanf(sysfile, "%s %s", str_real, str_imag);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF", filename);

            print_error(error_string);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            exit(BH_EXIT_BADPARSE);
        }

        /* get mpq_t values for real + imag */
        parse_complex_rational(str_real, str_imag, (*rational_v)->coord[i]);
    }

    /* clean up the file */
    errno = 0;
    res = fclose(sysfile);

    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Couldn't close %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    }

    v = (void *) rational_v;
}

/***********************************
* parse polynomial from input file *
************************************/
polynomial parse_polynomial_rational(FILE *sysfile, char *filename, int num_var) {
    int num_terms, res, i, j, max_degree;
    polynomial p;

    p.numVariables = num_var;

    /* first get the number of terms */
    res = fscanf(sysfile, "%d", &num_terms);

    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    p.exponents = malloc(num_terms * sizeof(int*));

    for (i = 0; i < num_terms; i++)
        p.exponents[i] = malloc(num_var * sizeof(int));
    p.coeff = malloc(num_terms * sizeof(rational_complex_number));
    p.numTerms = num_terms;
    p.degree = 0;
    p.isReal = 0;

    /* get the exponents and coefficients for each term */
    for (i = 0; i < num_terms; i++) {
        max_degree = 0;
        /* get exponent for each variable in the term */
        for (j = 0; j < num_var; j++) {
            res = fscanf(sysfile, "%d", &p.exponents[i][j]);
            if (res == EOF || res == 0) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF", filename);

                print_error(error_string);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

                print_error(error_string);
                exit(BH_EXIT_BADPARSE);
            }

            max_degree += p.exponents[i][j];
        }

        /* check if degree is greater than p.degree */
        if (max_degree > p.degree)
            p.degree = max_degree;

        /* get (real & imag) coefficients for the term
         * takes extra parsing since these will be fractions */
        char str_coeff_real[50], str_coeff_imag[50];
        res = fscanf(sysfile, "%s %s", str_coeff_real, str_coeff_imag);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF", filename);

            print_error(error_string);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            exit(BH_EXIT_BADPARSE);
        }

        /* get mpq_t coefficients for real + imag */
        initialize_rational_number(p.coeff[i]);
        parse_complex_rational(str_coeff_real, str_coeff_imag, p.coeff[i]);
    }

    mpq_init(p.norm_sqr);
    norm_sqr_polynomial(p.norm_sqr, &p);
    return p;
}

/*************************************************************
 * parse string into numerator and denominator, return mpq_t *
 *************************************************************/
void parse_complex_rational(char *str_real, char *str_imag, rational_complex_number c) {
    long int num;
    unsigned long int denom;
    char *tok;

    errno = 0;

    tok = strtok(str_real, "/");
    num = strtol(tok, NULL, 0);
    if (num == 0) {
        mpq_set_ui(c->re, 0, 1);
    } else if (errno == ERANGE) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Invalid coefficient: %s", str_real);

        print_error(error_string);
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
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Invalid coefficient: %s", str_real);

                print_error(error_string);
                exit(BH_EXIT_BADPARSE);
            }
        }
        mpq_set_si(c->re, num, denom);
    }

    tok = strtok(str_imag, "/");
    num = strtol(tok, NULL, 0);
    if (num == 0) {
        mpq_set_ui(c->im, 0, 1);
    } else if (errno == ERANGE) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Invalid coefficient: %s", str_imag);

        print_error(error_string);
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
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Invalid coefficient: %s", str_imag);

                print_error(error_string);
                exit(BH_EXIT_BADPARSE);
            }
        }
        mpq_set_si(c->im, num, denom);
    }
}

int read_points_file_rational(char *filename, void **t, void **w, int num_var) {
    if (verbosity > BH_VERBOSE)
        fprintf(stderr, "reading points file %s\n", filename);

    int i, j, num_points, res;
    mpq_t *t_rational;
    rational_complex_vector *w_rational;

    /* sanity-check the points file */
    errno = 0;
    FILE *pfile = fopen(filename, "r");
    if (pfile == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Couldn't open points file %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADFILE);
    }

    /* get number of points */
    errno = 0;

    res = fscanf(pfile, "%d", &num_points);

    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    /* check that we have enought points to test */
    if (num_points < 2) {
        char error_string[] = "You must define 2 or more points to test";

        print_error(error_string);
        exit(BH_EXIT_BADDEF);
    }

    t_rational = malloc(num_points * sizeof(mpq_t));
    w_rational = malloc(num_points * sizeof(rational_complex_vector));

    for (i = 0; i < num_points; i++) {
        /* read in each t */
        mpq_init(t_rational[i]);
        initialize_rational_vector(w_rational[i], num_var);

        errno = 0;
        res = gmp_fscanf(pfile, "%Qd", t_rational[i]);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF.", filename);

            print_error(error_string);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            exit(BH_EXIT_BADPARSE);
        }

        /* check if 0 < t < 1 */
        if (mpq_cmp_ui(t_rational[i], 0, 1) <= 0 || mpq_cmp_ui(t_rational[i], 1, 1) >= 0) {
            /* error here */
        }

        for (j = 0; j < num_var; j++) {
            /* get real and imag points */
            char str_point_real[50], str_point_imag[50];
            res = fscanf(pfile, "%s %s", str_point_real, str_point_imag);

            if (res == EOF || res == 0) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: unexpected EOF.", filename);

                print_error(error_string);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Error reading %s: %s.", filename, strerror(errno));

                print_error(error_string);
                exit(BH_EXIT_BADPARSE);
            }

            parse_complex_rational(str_point_real, str_point_imag, w_rational[i]->coord[j]);
        }
    }

    /* clean up the file */
    errno = 0;
    res = fclose(pfile);

    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH + 1, "Couldn't close %s: %s.", filename, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }


    *w = (void *) w_rational;
    *t = (void *) t_rational;
    return num_points;
}
