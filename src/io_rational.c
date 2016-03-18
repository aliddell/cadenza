/*
 * Cadenza
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io_rational.c: Process system and points file using rational arithmetic
 */
#include "cadenza.h"

/**********************************************
 * function to compare mpq_t needed for qsort *
 **********************************************/
int compare_mpq(const void *a, const void *b) {
    const mpq_t *qa = (const mpq_t *) a;
    const mpq_t *qb = (const mpq_t *) b;

    if (sort_order == BH_ASCENDING)
        return (int) mpq_cmp(*qa, *qb);
    else
        return (int) mpq_cmp(*qb, *qa);
}

/********************************************************************
 * sort the t array from highest to lowest and adjust w accordingly *
 ********************************************************************/
void sort_points_rational(mpq_t *t, rational_complex_vector *x, int num_points) {
    int i, j, k, *indices = malloc(num_points * sizeof(int));
    mpq_t *t_cpy = malloc(num_points * sizeof(mpq_t));
    rational_complex_vector *x_cpy = malloc(num_points * sizeof(rational_complex_vector));

    /* make copies of t and w */
    for (i = 0; i < num_points; i++) {
        mpq_init(t_cpy[i]);
        mpq_set(t_cpy[i], t[i]);

        initialize_rational_vector(x_cpy[i], x[i]->size);
        copy_rational_vector(x_cpy[i], x[i]);
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
        copy_rational_vector(x_cpy[i], x[indices[i]]);

    /* copy the sorted copy back into w */
    for (i = 0; i < num_points; i++)
        copy_rational_vector(x[i], x_cpy[i]);

    /* clean up */
    for (i = 0; i < num_points; i++) {
        mpq_clear(t_cpy[i]);
        clear_rational_vector(x_cpy[i]);
    }

    free(indices);
    indices = NULL;
    free(t_cpy);
    t_cpy = NULL;
    free(x_cpy);
    x_cpy = NULL;
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
        snprintf(error_string, (size_t) termwidth, "Couldn't open system file `%s': %s", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    /* gather number of variables */
    errno = 0;
    res = fscanf(sysfh, "%d", &num_var);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", sysfile);

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (res == 0) {
        snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", sysfile, strerror(errno));

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
        system->polynomials[i] = parse_polynomial(sysfh, num_var);
        if (system->polynomials[i].degree > system->maximumDegree)
            system->maximumDegree = system->polynomials[i].degree;
    }

    /* read in v */
    initialize_rational_vector(*v_rational, num_var);

    for (i = 0; i < num_var; i++) {
        /* get (real & imag) parts for each v_i */
        errno = 0;

        res = gmp_fscanf(sysfh, "%Qd %Qd", (*v_rational)->coord[i]->re, (*v_rational)->coord[i]->im);

        if (res == EOF) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", sysfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (res == 0) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }
        mpq_canonicalize((*v_rational)->coord[i]->re);
        mpq_canonicalize((*v_rational)->coord[i]->im);
    }

    /* clean up the file */
    errno = 0;

    res = fclose(sysfh);
    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    }

    v = (void *) v_rational;
}

/********************
 * read points file *
 ********************/
int read_points_file_rational(void **t, void **x, int num_var) {
    int i, j, num_points, res;
    mpq_t *t_rational;
    rational_complex_vector *x_rational;

    /* sanity-check the points file */
    errno = 0;
    FILE *pointsfh = fopen(pointsfile, "r");
    if (pointsfh == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open points file `%s': %s", pointsfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    /* get number of points */
    errno = 0;

    res = fscanf(pointsfh, "%d", &num_points);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", pointsfile);

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (res == 0) {
        snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", pointsfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", pointsfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADPARSE);
    }

    /* check that we have enought points to test */
    if (num_points < 2) {
        char err_string[] = "You must define 2 or more points to test";

        print_error(err_string, stderr);
        exit(BH_EXIT_BADDEF);
    }

    t_rational = malloc(num_points * sizeof(mpq_t));
    x_rational = malloc(num_points * sizeof(rational_complex_vector));

    for (i = 0; i < num_points; i++) {
        /* read in each t */
        mpq_init(t_rational[i]);
        initialize_rational_vector(x_rational[i], num_var);

        errno = 0;
        res = gmp_fscanf(pointsfh, "%Qd", t_rational[i]);

        if (res == EOF) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", pointsfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (res == 0) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': are you using floating point input?", pointsfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", pointsfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }

        mpq_canonicalize(t_rational[i]);

        /* check if 0 <= t <= 1 */
        if (mpq_cmp_ui(t_rational[i], 0, 1) < 0 || mpq_cmp_ui(t_rational[i], 1, 1) > 0) {
            gmp_snprintf(error_string, (size_t) termwidth, "Value for t not between 0 and 1: %Qd", t_rational[i]);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADDEF);
        }

        for (j = 0; j < num_var; j++) {
            errno = 0;
            /* get real and imag points */
            res = gmp_fscanf(pointsfh, "%Qd %Qd", x_rational[i]->coord[j]->re, x_rational[i]->coord[j]->im);

            if (res == EOF) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", pointsfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (res == 0) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': are you using floating point input?", pointsfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", pointsfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADPARSE);
            }

            mpq_canonicalize(x_rational[i]->coord[j]->re);
            mpq_canonicalize(x_rational[i]->coord[j]->im);
        }
    }

    /* clean up the file */
    errno = 0;

    res = fclose(pointsfh);
    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", pointsfile, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }

    sort_points_rational(t_rational, x_rational, num_points);

    *x = (void *) x_rational;
    *t = (void *) t_rational;

    return num_points;
}

/**************************************************************
 * print the file input back to stdout for debugging purposes *
 **************************************************************/
void fprint_input_rational(FILE *outfile, polynomial_system *system, void *v, void *t, void *x, int num_points) {
    int i, num_var = system->numVariables;

    rational_complex_vector *v_rational = (rational_complex_vector *) v;
    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *x_rational = (rational_complex_vector *) x;

    fputs("F:\n", outfile);
    print_system(outfile, system);

    fputs("\n", outfile);

    fputs("v:\n", outfile);
    fprintf(outfile, "[");
    for (i = 0; i < num_var - 1; i++) {
        gmp_fprintf(outfile, "%Qd + I * %Qd, ", (*v_rational)->coord[i]->re, (*v_rational)->coord[i]->im);
    }
    gmp_fprintf(outfile, "%Qd + I * %Qd]\n", (*v_rational)->coord[i]->re, (*v_rational)->coord[i]->im);

    fputs("\n", outfile);

    fputs("(t, x)\n", outfile);
    for (i = 0; i < num_points; i++) {
        gmp_fprintf(outfile, "%Qd, ", t_rational[i]);
        print_points_rational(outfile, x_rational[i]);
    }
}

/*********************************
 * print a test point to outfile *
 *********************************/
void print_points_rational(FILE *outfile, rational_complex_vector points) {
    int i, num_var = points->size;

    fprintf(outfile, "[");
    for (i = 0; i < num_var - 1; i++) {
        if (mpq_cmp_ui(points->coord[i]->im, 0, 1) == 0)
            gmp_fprintf(outfile, "%Qd, ", points->coord[i]->re);
        else
            gmp_fprintf(outfile, "%Qd + I * %Qd, ", points->coord[i]->re, points->coord[i]->im);
    }
    if (mpq_cmp_ui(points->coord[i]->im, 0, 1) == 0)
        gmp_fprintf(outfile, "%Qd]\n", points->coord[i]->re);
    else
        gmp_fprintf(outfile, "%Qd + I * %Qd]\n", points->coord[i]->re, points->coord[i]->im);
}

/*****************************************
 * print interval information to outfile *
 *****************************************/
void fprint_interval_rational(FILE *outfile, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int i;

    for (i = 0; i < BH_TERMWIDTH; i++)
        fprintf(outfile, "=");

    fprintf(outfile, "\n");
    gmp_fprintf(outfile, "t interval: [%Qd, %Qd]\n", t_left, t_right);
    gmp_fprintf(outfile, "x_left: ");
    print_points_rational(outfile, x_left);

    gmp_fprintf(outfile, "alpha^2 (x_left): %Qd\n", alpha_sqr_left);
    gmp_fprintf(outfile, "beta^2 (x_left): %Qd\n", beta_sqr_left);
    gmp_fprintf(outfile, "gamma^2 (x_left): %Qd\n", gamma_sqr_left);

    gmp_fprintf(outfile, "x_right: ");
    print_points_rational(outfile, x_right);

    gmp_fprintf(outfile, "alpha^2 (x_right): %Qd\n", alpha_sqr_right);
    gmp_fprintf(outfile, "beta^2 (x_right): %Qd\n", beta_sqr_right);
    gmp_fprintf(outfile, "gamma^2 (x_right): %Qd\n", gamma_sqr_right);
}

/**************************************
 * print continuous intervals to file *
 **************************************/
void fprint_continuous_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int res;

    errno = 0;

    FILE *fh = fopen(BH_FCONT, "a");

    if (fh == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprint_interval_rational(fh, t_left, t_right, x_left, x_right, alpha_sqr_left, alpha_sqr_right, beta_sqr_left, beta_sqr_right, gamma_sqr_left, gamma_sqr_right);

    errno = 0;

    res = fclose(fh);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", BH_FCONT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}

/*****************************************
 * print discontinuous intervals to file *
 *****************************************/
void fprint_discontinuous_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int res;

    errno = 0;

    FILE *fh = fopen(BH_FDISCONT, "a");

    if (fh == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FDISCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprint_interval_rational(fh, t_left, t_right, x_left, x_right, alpha_sqr_left, alpha_sqr_right, beta_sqr_left, beta_sqr_right, gamma_sqr_left, gamma_sqr_right);

    errno = 0;

    res = fclose(fh);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", BH_FDISCONT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}

/*********************************************
 * print uncertain intervals to summary file *
 *********************************************/
void fprint_uncertain_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int res;

    errno = 0;
    FILE *fh = fopen(BH_FUNCERTIFIED, "a");

    if (fh == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FUNCERTIFIED, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprint_interval_rational(fh, t_left, t_right, x_left, x_right, alpha_sqr_left, alpha_sqr_right, beta_sqr_left, beta_sqr_right, gamma_sqr_left, gamma_sqr_right);

    errno = 0;

    res = fclose(fh);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", BH_FUNCERTIFIED, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}

/*****************************************************
 * print solutions in input format to an output file *
 *****************************************************/
void fprint_solutions_rational(void *t, void *x, int num_points) {
    int i, j, res;

    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *x_rational = (rational_complex_vector *) x;

    FILE *outfile = fopen(BH_FPTSOUT, "w");

    errno = 0;

    if (outfile == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FPTSOUT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(outfile, "%d\n", num_points);

    for (i = 0; i < num_points; i++) {
        gmp_fprintf(outfile, "%Qd\n", t_rational[i]);

        for (j = 0; j < x_rational[i]->size; j++)
            gmp_fprintf(outfile, "%Qd\t%Qd\n", x_rational[i]->coord[j]->re, x_rational[i]->coord[j]->im);
    }

    errno = 0;

    res = fclose(outfile);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", BH_FPTSOUT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}
