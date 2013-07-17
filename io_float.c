/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io_float.c: Process system and points file using floating-point arithmetic
 */
#include "blueharvest.h"

/***********************************************
 * read the polynomial/exponential system file *
 ***********************************************/
void read_system_file_float(char *filename, polynomial_system *system, void *v) {
    int i, res, num_var, num_poly;
    complex_vector *v_rational = (complex_vector *) v;

    /* sanity-check the system file */
    errno = 0;
    FILE *sysfile = fopen(filename, "r");
    if (sysfile == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open system file %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADFILE);
    }

    /* gather number of variables and polynomials */
    errno = 0;

    res = fscanf(sysfile, "%d", &num_var);

    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    num_poly = num_var; /* square system */

    system->numVariables = num_var;
    system->numPolynomials = num_poly;
    system->maximumDegree = 0;
    system->isReal = 1;
    mpq_init(system->norm_sqr);
    system->numExponentials = 0;
    system->polynomials = malloc(num_var * sizeof(polynomial));

    /* read in each polynomial piece by piece */

    for (i = 0; i < num_poly; i++) {
        system->polynomials[i] = parse_polynomial_float(sysfile, filename, num_var);
        if (system->polynomials[i].degree > system->maximumDegree)
            system->maximumDegree = system->polynomials[i].degree;
    }

    /* read in v */
    initialize_vector(*v_rational, num_var);
    for (i = 0; i < num_var; i++) {
        /* get (real & imag) coefficients for the term */
        errno = 0;
        char str_real[50], str_imag[50];
        res = fscanf(sysfile, "%s %s", str_real, str_imag);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", filename);

            print_error(error_string);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            exit(BH_EXIT_BADPARSE);
        }

        /* get mpf_t coefficients for real + imag */
        parse_complex_float(str_real, str_imag, (*v_rational)->coord[i]);
    }

    /* clean up the file */
    errno = 0;
    res = fclose(sysfile);

    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    }

    v = (void *) v_rational;
}

/***********************************
* parse polynomial from input file *
************************************/
polynomial parse_polynomial_float(FILE *sysfile, char *filename, int num_var) {
    int num_terms, res, i, j, max_degree;
    polynomial p;

    p.numVariables = num_var;

    /* first get the number of terms */
    res = fscanf(sysfile, "%d", &num_terms);

    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

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
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", filename);

                print_error(error_string);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

                print_error(error_string);
                exit(BH_EXIT_BADPARSE);
            }

            max_degree += p.exponents[i][j];
        }

        /* check if degree is greater than p.degree */
        if (max_degree > p.degree)
            p.degree = max_degree;

        /* get (real & imag) coefficients for the term */
        char str_coeff_real[50], str_coeff_imag[50];
        res = fscanf(sysfile, "%s %s", str_coeff_real, str_coeff_imag);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", filename);

            print_error(error_string);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            exit(BH_EXIT_BADPARSE);
        }

        /* get mpf_t coefficients for real + imag */
        initialize_number(p.coeff[i]);
        parse_complex_float(str_coeff_real, str_coeff_imag, p.coeff[i]);
    }

    mpq_init(p.norm_sqr);
    norm_sqr_polynomial(p.norm_sqr, &p);
    return p;
}

/*************************************************************
 * parse string into numerator and denominator, return mpf_t *
 *************************************************************/
void parse_complex_float(char *str_real, char *str_imag, complex_number c) {
    int res;

    errno = 0;

    res = mpf_set_str(c->re, str_real, 0);

    /* error! */
    if (res == -1) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Invalid coefficient: %s", str_real);

        print_error(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    res = mpf_set_str(c->im, str_imag, 0);

    /* error! */
    if (res == -1) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Invalid coefficient: %s", str_imag);

        print_error(error_string);
        exit(BH_EXIT_BADPARSE);
    }
}

int read_points_file_float(char *filename, void **t, void **w, int num_var) {
    int i, j, num_points, res;
    mpf_t *t_float;
    complex_vector *w_float;

    /* sanity-check the points file */
    errno = 0;
    FILE *pfile = fopen(filename, "r");
    if (pfile == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open points file %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADFILE);
    }

    /* get number of points */
    errno = 0;

    res = fscanf(pfile, "%d", &num_points);

    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", filename);

        print_error(error_string);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

        print_error(error_string);
        exit(BH_EXIT_BADPARSE);
    }

    /* check that we have enough points to test */
    if (num_points < 2) {
        char error_string[] = "You must define 2 or more points to test";

        print_error(error_string);
        exit(BH_EXIT_BADDEF);
    }

    t_float = malloc(num_points * sizeof(mpf_t));
    w_float = malloc(num_points * sizeof(complex_vector));

    for (i = 0; i < num_points; i++) {
        /* read in each t */
        mpf_init(t_float[i]);
        initialize_vector(w_float[i], num_var);

        errno = 0;
        res = gmp_fscanf(pfile, "%Fd", t_float[i]);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF.", filename);

            print_error(error_string);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

            print_error(error_string);
            exit(BH_EXIT_BADPARSE);
        }

        /* check if 0 < t < 1 */
        if (mpf_cmp_ui(t_float[i], 0) < 0 || mpf_cmp_ui(t_float[i], 1) > 0) {
            char error_string[BH_TERMWIDTH];
            gmp_snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Value for t not between 0 and 1: %Fd", t_float[i]);

            print_error(error_string);
            exit(BH_EXIT_BADDEF);
        }

        for (j = 0; j < num_var; j++) {
            /* get real and imag points */
            char str_point_real[50], str_point_imag[50];
            res = fscanf(pfile, "%s %s", str_point_real, str_point_imag);

            if (res == EOF || res == 0) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF.", filename);

                print_error(error_string);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", filename, strerror(errno));

                print_error(error_string);
                exit(BH_EXIT_BADPARSE);
            }

            parse_complex_float(str_point_real, str_point_imag, w_float[i]->coord[j]);
        }
    }

    /* clean up the file */
    errno = 0;
    res = fclose(pfile);

    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", filename, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }

    *w = (void *) w_float;
    *t = (void *) t_float;
    return num_points;
}

/********************************
 * print a test point to stdout *
 ********************************/
void print_points_float(complex_vector points, FILE *OUT) {
    int i, num_var = points->size;

    for (i = 0; i < num_var; i++) {
        gmp_printf("\t\t[%Fd + %Fdi]\n", points->coord[i]->re, points->coord[i]->im);
    }

}

/***************************************
 * print a polynomial system to stdout *
 ***************************************/
void print_system_float(polynomial_system *system, FILE *OUT) {
    int i, j, k;
    char variables[] = "xyzwuvabcdjkmnpqrs";

    for (i = 0; i < system->numPolynomials; i++) {
        polynomial p = system->polynomials[i];
        for (j = 0; j < p.numTerms - 1; j++) {
            gmp_printf("(%Fd + %Fdi)", p.coeff[j]->re, p.coeff[j]->im);
            for (k = 0; k < p.numVariables; k++)
                printf("%c^%d", variables[k], p.exponents[j][k]);
            printf(" + ");
        }
        gmp_printf("(%Fd + %Fdi)", p.coeff[j]->re, p.coeff[j]->im);
        for (k = 0; k < system->numVariables; k++)
            printf("%c^%d", variables[k], p.exponents[j][k]);
        puts("");
    }
}
