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

/**********************************************
 * function to compare mpf_t needed for qsort *
 **********************************************/
int compare_mpf(const void *a, const void *b) {
    const mpf_t *fa = (const mpf_t *) a;
    const mpf_t *fb = (const mpf_t *) b;

    return (int) mpf_cmp(*fa, *fb);
}

/********************************************************************
 * sort the t array from lowest to highest and adjust w accordingly *
 ********************************************************************/
void sort_points_float(mpf_t *t, complex_vector *w, int num_points) {
    int i, j, k, *indices = malloc(num_points * sizeof(int));
    mpf_t *t_cpy = malloc(num_points * sizeof(mpf_t));
    complex_vector *w_cpy = malloc(num_points * sizeof(complex_vector));

    /* make copies of t and w */
    for (i = 0; i < num_points; i++) {
        mpf_init(t_cpy[i]);
        mpf_set(t_cpy[i], t[i]);

        initialize_vector(w_cpy[i], w[i]->size);
        copy_vector(w_cpy[i], w[i]);
    }

    /* sort t, leave t_cpy intact */
    qsort((void *) t, (size_t) num_points, sizeof(mpf_t), &compare_mpf);

    /* collect the indices of the unsorted t-values as they relate to the sorted ones */
    for (i = 0; i < num_points; i++) {
        for (j = 0; j < num_points; j++) {
            if (mpf_cmp(t_cpy[j], t[i]) == 0) {
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
        copy_vector(w_cpy[i], w[indices[i]]);

    /* copy the sorted copy back into w */
    for (i = 0; i < num_points; i++)
        copy_vector(w[i], w_cpy[i]);

    /* clean up */
    for (i = 0; i < num_points; i++) {
        mpf_clear(t_cpy[i]);
        clear_vector(w_cpy[i]);
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
void read_system_file_float(polynomial_system *system, void *v) {
    int i, res, num_var = 0, num_poly = 0, base = 10;

    complex_vector *v_float = (complex_vector *) v;

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
        system->polynomials[i] = parse_polynomial(sysfh, num_var);
        if (system->polynomials[i].degree > system->maximumDegree)
            system->maximumDegree = system->polynomials[i].degree;
    }

    initialize_vector(*v_float, num_var);

    for (i = 0; i < num_var; i++) {
        /* get (real & imag) parts for each v_i */
        errno = 0;

        res = mpf_inp_str((*v_float)->coord[i]->re, sysfh, base);
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

        res = mpf_inp_str((*v_float)->coord[i]->im, sysfh, base);
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

    v = (void *) v_float;
}

/********************
 * read points file *
 ********************/
int read_points_file_float(void **t, void **w, int num_var) {
    int i, j, num_points, res, base = 10;;

    mpf_t *t_float;

    complex_vector *w_float;

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

    t_float = malloc(num_points * sizeof(mpf_t));
    w_float = malloc(num_points * sizeof(complex_vector));

    for (i = 0; i < num_points; i++) {
        /* read in each t */
        mpf_init(t_float[i]);
        initialize_vector(w_float[i], num_var);

        errno = 0;

        res = mpf_inp_str(t_float[i], pointsfh, base);
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

        /* check if 0 < t < 1 */
        if (mpf_cmp_ui(t_float[i], 0) < 0 || mpf_cmp_ui(t_float[i], 1) > 0) {
            char error_string[BH_TERMWIDTH];
            mpfr_snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Value for t not between 0 and 1: %.Rf", t_float[i]);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADDEF);
        }

        for (j = 0; j < num_var; j++) {
            /* get real and imag points */
            res = mpf_inp_str(w_float[i]->coord[j]->re, pointsfh, base);
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

            res = mpf_inp_str(w_float[i]->coord[j]->im, pointsfh, base);
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

    sort_points_float(t_float, w_float, num_points);

    *w = (void *) w_float;
    *t = (void *) t_float;

    return num_points;
}

/**************************************************************
 * print the file input back to stdout for debugging purposes *
 **************************************************************/
void print_back_input_float(polynomial_system *system, void *v, void *t, void *w, int num_points) {
    int i, num_var = system->numVariables;

    complex_vector *v_float = (complex_vector *) v;
    mpf_t *t_float = (mpf_t *) t;
    complex_vector *w_float = (complex_vector *) w;

    puts("F:");
    print_system_float(system, stdout);

    puts("\n");

    puts("v:");
    printf("[");
    for (i = 0; i < num_var - 1; i++) {
        mpfr_printf("%.Rf + %.Rfi, ", (*v_float)->coord[i]->re, (*v_float)->coord[i]->im);
    }
    mpfr_printf("%.Rf + %.Rfi]\n", (*v_float)->coord[i]->re, (*v_float)->coord[i]->im);

    puts("\n");

    puts("(t_i, w_i)");
    for (i = 0; i < num_points; i++) {
        mpfr_printf("%.Rf, ", t_float[i]);
        print_points_float(w_float[i], stdout);
    }
}

/*********************************
 * print a test point to outfile *
 *********************************/
void print_points_float(complex_vector points, FILE *outfile) {
    int i, num_var = points->size;

    fprintf(outfile, "[");
    for (i = 0; i < num_var - 1; i++) {
        if (mpf_cmp_ui(points->coord[i]->im, 0) == 0)
            mpfr_fprintf(outfile, "%.Rf, ", points->coord[i]->re);
        else
            mpfr_fprintf(outfile, "%.Rf + %.Rfi, ", points->coord[i]->re, points->coord[i]->im);
    }
    if (mpf_cmp_ui(points->coord[i]->im, 0) == 0)
        mpfr_fprintf(outfile, "%.Rf]\n", points->coord[i]->re);
    else
        mpfr_fprintf(outfile, "%.Rf + %.Rfi]\n", points->coord[i]->re, points->coord[i]->im);
}

/****************************************
 * print a polynomial system to outfile *
 ****************************************/
void print_system_float(polynomial_system *system, FILE *outfile) {
    int i, j, k;

    mpf_t re;
    mpf_t im;
    mpf_init(re);
    mpf_init(im);

    for (i = 0; i < system->numPolynomials; i++) {
        polynomial p = system->polynomials[i];
        for (j = 0; j < p.numTerms - 1; j++) {
            /* real part is nonzero */
            if (mpq_cmp_ui(p.coeff[j]->re, 0, 1) != 0) {
                /* imaginary part is zero */
                if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0) {
                    /* real coefficient is not 1 */
                    if (mpq_cmp_ui(p.coeff[j]->re, 1, 1) != 0) {
                        mpf_set_q(re, p.coeff[j]->re);
                        mpfr_fprintf(outfile, "%.Rf", re);
                    } else {
                        mpfr_fprintf(outfile, "1");
                    }
                } 
                /* imaginary part is 1 */
                else if (mpq_cmp_ui(p.coeff[j]->im, 1, 1) == 1) {
                    mpf_set_q(re, p.coeff[j]->re);
                    mpfr_fprintf(outfile, "(%.Rf + i)", re);
                }
                /* imaginary part is neither 0 nor 1 */
                else {
                    mpf_set_q(re, p.coeff[j]->re);
                    mpf_set_q(im, p.coeff[j]->im);
                    mpfr_fprintf(outfile, "(%.Rf + %.Rfi)", re, im);
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
                else {
                    mpf_set_q(im, p.coeff[j]->im);
                    mpfr_fprintf(outfile, "%.Rfi", im);
                }
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
                if (mpq_cmp_ui(p.coeff[j]->re, 1, 1) != 0) {
                    mpf_set_q(re, p.coeff[j]->re);
                    mpfr_fprintf(outfile, "%.Rf", re);
                } else {
                    mpfr_fprintf(outfile, "1");
                }
            } 
            /* imaginary part is 1 */
            else if (mpq_cmp_ui(p.coeff[j]->im, 1, 1) == 1) {
                mpf_set_q(re, p.coeff[j]->re);
                mpfr_fprintf(outfile, "(%.Rf + i)", re);
            }
            /* imaginary part is neither 0 nor 1 */
            else {
                mpf_set_q(re, p.coeff[j]->re);
                mpf_set_q(im, p.coeff[j]->im);
                mpfr_fprintf(outfile, "(%.Rf + %.Rfi)", re, im);
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
            else {
                mpf_set_q(im, p.coeff[j]->im);
                mpfr_fprintf(outfile, "%.Rfi", im);
            }
        }

        for (k = 0; k < p.numVariables; k++) {
            if (p.exponents[j][k] == 1)
                fprintf(outfile, "(x_%d)", k);
            else if (p.exponents[j][k] != 0)
                fprintf(outfile, "(x_%d)^%d", k, p.exponents[j][k]);
        }
        fputs("\n", outfile);
    }

    mpf_clear(re);
    mpf_clear(im);
}

/**************************************
 * print continuous intervals to file *
 **************************************/
void fprint_continuous_float(mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right) {
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
    mpfr_fprintf(contfh, "t interval: [%.Rf, %.Rf]\n", t_left, t_right);
    mpfr_fprintf(contfh, "x_left: ");
    print_points_float(w_left, contfh);

    mpfr_fprintf(contfh, "alpha (x_left): %.Rf\n", alpha_left);
    mpfr_fprintf(contfh, "beta (x_left): %.Rf\n", beta_left);
    mpfr_fprintf(contfh, "gamma (x_left): %.Rf\n", gamma_left);

    mpfr_fprintf(contfh, "x_right: ");
    print_points_float(w_right, contfh);

    mpfr_fprintf(contfh, "alpha (x_right): %.Rf\n", alpha_right);
    mpfr_fprintf(contfh, "beta (x_right): %.Rf\n", beta_right);
    mpfr_fprintf(contfh, "gamma (x_right): %.Rf\n", gamma_right);

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
void fprint_discontinuous_float(mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right) {
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
    mpfr_fprintf(discfh, "t interval: [%.Rf, %.Rf]\n", t_left, t_right);
    mpfr_fprintf(discfh, "x_left: ");
    print_points_float(w_left, discfh);

    mpfr_fprintf(discfh, "alpha (x_left): %.Rf\n", alpha_left);
    mpfr_fprintf(discfh, "beta (x_left): %.Rf\n", beta_left);
    mpfr_fprintf(discfh, "gamma (x_left): %.Rf\n", gamma_left);

    mpfr_fprintf(discfh, "x_right: ");
    print_points_float(w_right, discfh);

    mpfr_fprintf(discfh, "alpha (x_right): %.Rf\n", alpha_right);
    mpfr_fprintf(discfh, "beta (x_right): %.Rf\n", beta_right);
    mpfr_fprintf(discfh, "gamma (x_right): %.Rf\n", gamma_right);

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
void fprint_solutions_float(void *t, void *w, int num_points) {
    int i, j, res;

    mpf_t *t_float = (mpf_t *) t;
    complex_vector *w_float = (complex_vector *) w;

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
        mpfr_fprintf(outfile, "%.Rf\n", t_float[i]);

        for (j = 0; j < w_float[i]->size; j++)
            mpfr_fprintf(outfile, "%.Rf\t%.Rf\n", w_float[i]->coord[j]->re, w_float[i]->coord[j]->im);
    }

    errno = 0;

    res = fclose(outfile);

    if (res == EOF) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't close %s: %s.", BH_FPTSOUT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}
