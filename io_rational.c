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

/**********************************************
 * function to compare mpq_t needed for qsort *
 **********************************************/
int compare_mpq(const void *a, const void *b) {
    const mpq_t *qa = (const mpq_t *) a;
    const mpq_t *qb = (const mpq_t *) b;

    return (int) mpq_cmp(*qa, *qb);
}

/********************************************************************
 * sort the t array from lowest to highest and adjust w accordingly *
 ********************************************************************/
void sort_points_rational(mpq_t *t, rational_complex_vector *w, int num_points) {
    int i, j, k, *indices = malloc(num_points * sizeof(int));
    mpq_t *t_cpy = malloc(num_points * sizeof(mpq_t));
    rational_complex_vector *w_cpy = malloc(num_points * sizeof(rational_complex_vector));

    /* make copies of t and w */
    for (i = 0; i < num_points; i++) {
        mpq_init(t_cpy[i]);
        mpq_set(t_cpy[i], t[i]);

        initialize_rational_vector(w_cpy[i], w[i]->size);
        copy_rational_vector(w_cpy[i], w[i]);
    }

    /* sort t, leave t_cpy intact */
    qsort((void *) t, (size_t) num_points, sizeof(mpq_t), &compare_mpq);

    /* collect the indices of the unsorted t-values as they relate to the sorted ones */
    for (i = 0; i < num_points; i++) {
        for (j = 0; j < num_points; j++) {
            if (mpq_cmp(t_cpy[j], t[i]) == 0) {
                int found = 0;
                for (k = 0; k < i; k++)
                    if (j == indices[k])
                        found = 1;
                if (!found) {
                    indices[i] = j;
                    break;
                }
            }
        }
    }

    /* make a `sorted' copy of w */
    for (i = 0; i < num_points; i++)
        copy_rational_vector(w_cpy[i], w[indices[i]]);

    /* copy the sorted copy back into w */
    for (i = 0; i < num_points; i++)
        copy_rational_vector(w[i], w_cpy[i]);

    /* clean up */
    for (i = 0; i < num_points; i++) {
        mpq_clear(t_cpy[i]);
        clear_rational_vector(w_cpy[i]);
    }

    free(indices);
    indices = NULL;
    free(t_cpy);
    t_cpy = NULL;
    free(w_cpy);
    w_cpy = NULL;
}

/***********************************
 * read the polynomial system file *
 ***********************************/
void read_system_file_rational(polynomial_system *system, void *v) {
    int i, res, num_var, num_poly;
    rational_complex_vector *v_rational = (rational_complex_vector *) v;

    /* sanity-check the system file */
    errno = 0;
    FILE *sysfh = fopen(sysfile, "r");
    if (sysfh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open system file %s: %s.", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    /* gather number of variables */
    errno = 0;
    res = fscanf(sysfh, "%d", &num_var);
    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", sysfile, strerror(errno));

        print_error(error_string, stderr);
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
    system->exponentials = NULL;

    /* read in each polynomial piece by piece */
    for (i = 0; i < num_poly; i++) {
        system->polynomials[i] = parse_polynomial_rational(sysfh, num_var);
        if (system->polynomials[i].degree > system->maximumDegree)
            system->maximumDegree = system->polynomials[i].degree;
    }

    /* read in v */
    initialize_rational_vector(*v_rational, num_var);

    for (i = 0; i < num_var; i++) {
        /* get (real & imag) parts for each v_i */
        errno = 0;

        res = gmp_fscanf(sysfh, "%Qd %Qd", (*v_rational)->coord[i]->re, (*v_rational)->coord[i]->im);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }

    }

    /* clean up the file */
    errno = 0;

    res = fclose(sysfh);
    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    }

    v = (void *) v_rational;
}

/***********************************
* parse polynomial from input file *
************************************/
polynomial parse_polynomial_rational(FILE *sysfh, int num_var) {
    int num_terms, res, i, j, max_degree;
    polynomial p;

    p.numVariables = num_var;

    /* first get the number of terms */
    res = fscanf(sysfh, "%d", &num_terms);

    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", sysfile, strerror(errno));

        print_error(error_string, stderr);
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
            res = fscanf(sysfh, "%d", &p.exponents[i][j]);
            if (res == EOF || res == 0) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", sysfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADPARSE);
            }

            max_degree += p.exponents[i][j];
        }

        /* check if degree is greater than p.degree */
        if (max_degree > p.degree)
            p.degree = max_degree;

        initialize_rational_number(p.coeff[i]);

        errno = 0;

        res = gmp_fscanf(sysfh, "%Qd %Qd", p.coeff[i]->re, p.coeff[i]->im);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }
    }

    mpq_init(p.norm_sqr);
    norm_sqr_polynomial(p.norm_sqr, &p);

    return p;
}

/********************
 * read points file *
 ********************/
int read_points_file_rational(void **t, void **w, int num_var) {
    int i, j, num_points, res;
    mpq_t *t_rational;
    rational_complex_vector *w_rational;

    /* sanity-check the points file */
    errno = 0;
    FILE *pointsfh = fopen(pointsfile, "r");
    if (pointsfh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open points file %s: %s.", pointsfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    /* get number of points */
    errno = 0;

    res = fscanf(pointsfh, "%d", &num_points);

    if (res == EOF || res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", pointsfile);

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", pointsfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADPARSE);
    }

    /* check that we have enought points to test */
    if (num_points < 2) {
        char error_string[] = "You must define 2 or more points to test";

        print_error(error_string, stderr);
        exit(BH_EXIT_BADDEF);
    }

    t_rational = malloc(num_points * sizeof(mpq_t));
    w_rational = malloc(num_points * sizeof(rational_complex_vector));

    for (i = 0; i < num_points; i++) {
        /* read in each t */
        mpq_init(t_rational[i]);
        initialize_rational_vector(w_rational[i], num_var);

        errno = 0;
        res = gmp_fscanf(pointsfh, "%Qd", t_rational[i]);

        if (res == EOF || res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF.", pointsfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", pointsfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }

        /* check if 0 < t < 1 */
        if (mpq_cmp_ui(t_rational[i], 0, 1) < 0 || mpq_cmp_ui(t_rational[i], 1, 1) > 0) {
            char error_string[BH_TERMWIDTH];
            gmp_snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Value for t not between 0 and 1: %Qd", t_rational[i]);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADDEF);
        }

        for (j = 0; j < num_var; j++) {
            /* get real and imag points */
            res = gmp_fscanf(pointsfh, "%Qd %Qd", w_rational[i]->coord[j]->re, w_rational[i]->coord[j]->im);

            if (res == EOF || res == 0) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF.", pointsfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s.", pointsfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADPARSE);
            }
        }
    }

    /* clean up the file */
    errno = 0;

    res = fclose(pointsfh);
    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", pointsfile, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }

    sort_points_rational(t_rational, w_rational, num_points);

    *w = (void *) w_rational;
    *t = (void *) t_rational;

    return num_points;
}

/**************************************************************
 * print the file input back to stdout for debugging purposes *
 **************************************************************/
void print_back_input_rational(polynomial_system *system, void *v, void *t, void *w, int num_points) {
    int i, num_var = system->numVariables;

    rational_complex_vector *v_rational = (rational_complex_vector *) v;
    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *w_rational = (rational_complex_vector *) w;

    puts("F:");
    print_system_rational(system, stdout);

    puts("\n");

    puts("v:");
    printf("[");
    for (i = 0; i < num_var - 1; i++) {
        gmp_printf("%Qd + %Qdi, ", (*v_rational)->coord[i]->re, (*v_rational)->coord[i]->im);
    }
    gmp_printf("%Qd + %Qdi]\n", (*v_rational)->coord[i]->re, (*v_rational)->coord[i]->im);

    puts("\n");

    puts("(t_i, w_i)");
    for (i = 0; i < num_points; i++) {
        gmp_printf("%Qd, ", t_rational[i]);
        print_points_rational(w_rational[i], stdout);
    }
}

/*********************************
 * print a test point to outfile *
 *********************************/
void print_points_rational(rational_complex_vector points, FILE *outfile) {
    int i, num_var = points->size;

    fprintf(outfile, "[");
    for (i = 0; i < num_var - 1; i++) {
        if (mpq_cmp_ui(points->coord[i]->im, 0, 1) == 0)
            gmp_fprintf(outfile, "%Qd, ", points->coord[i]->re);
        else
            gmp_fprintf(outfile, "%Qd + %Qdi, ", points->coord[i]->re, points->coord[i]->im);
    }
    if (mpq_cmp_ui(points->coord[i]->im, 0, 1) == 0)
        gmp_fprintf(outfile, "%Qd]\n", points->coord[i]->re);
    else
        gmp_fprintf(outfile, "%Qd + %Qdi]\n", points->coord[i]->re, points->coord[i]->im);
}

/****************************************
 * print a polynomial system to outfile *
 ****************************************/
void print_system_rational(polynomial_system *system, FILE *outfile) {
    int i, j, k;

    for (i = 0; i < system->numPolynomials; i++) {
        polynomial p = system->polynomials[i];
        for (j = 0; j < p.numTerms - 1; j++) {
            /* real part is nonzero */
            if (mpq_cmp_ui(p.coeff[j]->re, 0, 1) != 0) {
                /* imaginary part is zero */
                if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0) {
                    /* real coefficient is not 1 */
                    if (mpq_cmp_ui(p.coeff[j]->re, 1, 1) != 0)
                        gmp_fprintf(outfile, "%Qd", p.coeff[j]->re);
                } 
                /* imaginary part is 1 */
                else if (mpq_cmp_ui(p.coeff[j]->im, 1, 1) == 1) {
                    gmp_fprintf(outfile, "(%Qd + i)", p.coeff[j]->re);
                }
                /* imaginary part is neither 0 nor 1 */
                else {
                    gmp_fprintf(outfile, "(%Qd + %Qdi)", p.coeff[j]->re, p.coeff[j]->im);
                }
            }
            /* real part is zero */
            else {
                /* imaginary part is zero too */
                if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0)
                    continue;
                /* imaginary part is 1 */
                else if (mpq_cmp_ui(p.coeff[j]->im, 1, 1) == 0)
                    fprintf(outfile, "i");
                /* imaginary part is neither 0 nor 1 */
                else
                    gmp_fprintf(outfile, "%Qdi", p.coeff[j]->im);
            }

            for (k = 0; k < p.numVariables; k++) {
                if (p.exponents[j][k] == 1)
                    fprintf(outfile, "(x_%d)", k);
                else if (p.exponents[j][k] != 0)
                    fprintf(outfile, "(x_%d)^%d", k, p.exponents[j][k]);
            }
            fprintf(outfile, " + ");
        }

        /* last term */
        if (mpq_cmp_ui(p.coeff[j]->re, 0, 1) != 0) {
            /* imaginary part is zero */
            if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0) {
                /* real coefficient is not 1 */
                if (mpq_cmp_ui(p.coeff[j]->re, 1, 1) != 0)
                    gmp_fprintf(outfile, "%Qd", p.coeff[j]->re);
            } 
            /* imaginary part is 1 */
            else if (mpq_cmp_ui(p.coeff[j]->im, 1, 1) == 1) {
                gmp_fprintf(outfile, "(%Qd + i)", p.coeff[j]->re);
            }
            /* imaginary part is neither 0 nor 1 */
            else {
                gmp_fprintf(outfile, "(%Qd + %Qdi)", p.coeff[j]->re, p.coeff[j]->im);
            }
        }
        /* real part is zero */
        else {
            /* imaginary part is zero too */
            if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0)
                fprintf(outfile, "0");
            /* imaginary part is 1 */
            else if (mpq_cmp_ui(p.coeff[j]->im, 1, 1) == 0)
                fprintf(outfile, "i");
            /* imaginary part is neither 0 nor 1 */
            else
                gmp_fprintf(outfile, "%Qdi", p.coeff[j]->im);
        }

        for (k = 0; k < p.numVariables; k++) {
            if (p.exponents[j][k] == 1)
                fprintf(outfile, "(x_%d)", k);
            else if (p.exponents[j][k] != 0)
                fprintf(outfile, "(x_%d)^%d", k, p.exponents[j][k]);
        }
        fputs("\n", outfile);
    }
}

/**************************************
 * print continuous intervals to file *
 **************************************/
void fprint_continuous_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector w_left, rational_complex_vector w_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int i, res;

    errno = 0;

    FILE *contfh = fopen(BH_FCONT, "a");
    if (contfh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s.", BH_FCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    for (i = 0; i < BH_TERMWIDTH; i++)
        fprintf(contfh, "=");

    fprintf(contfh, "\n");
    gmp_fprintf(contfh, "t interval: [%Qd, %Qd]\n", t_left, t_right);
    gmp_fprintf(contfh, "x_left: ");
    print_points_rational(w_left, contfh);

    gmp_fprintf(contfh, "alpha^2 (x_left): %Qd\n", alpha_sqr_left);
    gmp_fprintf(contfh, "beta^2 (x_left): %Qd\n", beta_sqr_left);
    gmp_fprintf(contfh, "gamma^2 (x_left): %Qd\n", gamma_sqr_left);

    gmp_fprintf(contfh, "x_right: ");
    print_points_rational(w_right, contfh);

    gmp_fprintf(contfh, "alpha^2 (x_right): %Qd\n", alpha_sqr_right);
    gmp_fprintf(contfh, "beta^2 (x_right): %Qd\n", beta_sqr_right);
    gmp_fprintf(contfh, "gamma^2 (x_right): %Qd\n", gamma_sqr_right);

    errno = 0;

    res = fclose(contfh);
    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", BH_FCONT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}

/*****************************************
 * print discontinuous intervals to file *
 *****************************************/
void fprint_discontinuous_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector w_left, rational_complex_vector w_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int i, res;

    errno = 0;

    FILE *discfh = fopen(BH_FDISCONT, "a");
    if (discfh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s.", BH_FDISCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    for (i = 0; i < BH_TERMWIDTH; i++)
        fprintf(discfh, "=");

    fprintf(discfh, "\n");
    gmp_fprintf(discfh, "t interval: [%Qd, %Qd]\n", t_left, t_right);
    gmp_fprintf(discfh, "x_left: ");
    print_points_rational(w_left, discfh);

    gmp_fprintf(discfh, "alpha^2 (x_left): %Qd\n", alpha_sqr_left);
    gmp_fprintf(discfh, "beta^2 (x_left): %Qd\n", beta_sqr_left);
    gmp_fprintf(discfh, "gamma^2 (x_left): %Qd\n", gamma_sqr_left);

    gmp_fprintf(discfh, "x_right: ");
    print_points_rational(w_right, discfh);

    gmp_fprintf(discfh, "alpha^2 (x_right): %Qd\n", alpha_sqr_right);
    gmp_fprintf(discfh, "beta^2 (x_right): %Qd\n", beta_sqr_right);
    gmp_fprintf(discfh, "gamma^2 (x_right): %Qd\n", gamma_sqr_right);

    errno = 0;

    res = fclose(discfh);
    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", BH_FDISCONT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}

/*****************************************************
 * print solutions in input format to an output file *
 *****************************************************/
void fprint_solutions_rational(void *t, void *w, int num_points) {
    int i, j, res;

    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *w_rational = (rational_complex_vector *) w;

    FILE *outfile = fopen(BH_FPTSOUT, "w");

    errno = 0;

    if (outfile == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s.", BH_FPTSOUT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(outfile, "%d\n", num_points);

    for (i = 0; i < num_points; i++) {
        gmp_fprintf(outfile, "%Qd\n", t_rational[i]);

        for (j = 0; j < w_rational[i]->size; j++)
            gmp_fprintf(outfile, "%Qd\t%Qd\n", w_rational[i]->coord[j]->re, w_rational[i]->coord[j]->im);
    }

    errno = 0;

    res = fclose(outfile);

    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", BH_FPTSOUT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}
