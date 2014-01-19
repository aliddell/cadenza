/*
 * Cadenza
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io_float.c: Process system and points file using floating-point arithmetic
 */
#include "cadenza.h"

/**********************************************
 * function to compare mpf_t needed for qsort *
 **********************************************/
int compare_mpf(const void *a, const void *b) {
    const mpf_t *fa = (const mpf_t *) a;
    const mpf_t *fb = (const mpf_t *) b;

    if (sort_order == BH_ASCENDING)
        return (int) mpf_cmp(*fa, *fb);
    else
        return (int) mpf_cmp(*fb, *fa);
}

/********************************************************************
 * sort the t array from highest to lowest and adjust w accordingly *
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

    initialize_vector(*v_float, num_var);

    for (i = 0; i < num_var; i++) {
        /* get (real & imag) parts for each v_i */
        errno = 0;
        res = mpf_inp_str((*v_float)->coord[i]->re, sysfh, base);

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

        errno = 0;
        res = mpf_inp_str((*v_float)->coord[i]->im, sysfh, base);

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
    }

    /* clean up the file */
    errno = 0;
    res = fclose(sysfh);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    }

    v = (void *) v_float;
}

/********************
 * read points file *
 ********************/
int read_points_file_float(void **t, void **w, int num_var) {
    int i, j, num_points, res, base = 10;

    mpf_t *t_float;

    complex_vector *w_float;

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

    t_float = malloc(num_points * sizeof(mpf_t));
    w_float = malloc(num_points * sizeof(complex_vector));

    for (i = 0; i < num_points; i++) {
        /* read in each t */
        mpf_init(t_float[i]);
        initialize_vector(w_float[i], num_var);

        errno = 0;
        res = mpf_inp_str(t_float[i], pointsfh, base);

        if (res == EOF) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", pointsfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (res == 0) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': are you using rational input?", pointsfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", pointsfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }

        /* check if 0 < t < 1 */
        if (mpf_cmp_ui(t_float[i], 0) < 0 || mpf_cmp_ui(t_float[i], 1) > 0) {
            mpfr_snprintf(error_string, (size_t) termwidth, "Value for t not between 0 and 1: %.15Re", t_float[i]);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADDEF);
        }

        /* get real and imag points */
        for (j = 0; j < num_var; j++) {
            errno = 0;
            res = mpf_inp_str(w_float[i]->coord[j]->re, pointsfh, base);

            if (res == EOF) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", pointsfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (res == 0) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': are you using rational input?", pointsfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", pointsfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADPARSE);
            }

            errno = 0;
            res = mpf_inp_str(w_float[i]->coord[j]->im, pointsfh, base);

            if (res == EOF) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': unexpected EOF", pointsfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (res == 0) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': are you using rational input?", pointsfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                snprintf(error_string, (size_t) termwidth, "Error reading `%s': %s", pointsfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADPARSE);
            }
        }
    }

    /* clean up the file */
    errno = 0;
    res = fclose(pointsfh);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", pointsfile, strerror(errno));
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
void fprint_input_float(FILE *outfile, polynomial_system *system, void *v, void *t, void *w, int num_points) {
    int i, num_var = system->numVariables;

    complex_vector *v_float = (complex_vector *) v;
    mpf_t *t_float = (mpf_t *) t;
    complex_vector *w_float = (complex_vector *) w;

    fputs("F:\n", outfile);
    print_system(outfile, system);

    fputs("\n", outfile);

    fputs("v:\n", outfile);
    fprintf(outfile, "[");
    for (i = 0; i < num_var - 1; i++) {
        mpfr_fprintf(outfile, "%.15Re + I * %.15Re, ", (*v_float)->coord[i]->re, (*v_float)->coord[i]->im);
    }
    mpfr_fprintf(outfile, "%.15Re + I * %.15Re]\n", (*v_float)->coord[i]->re, (*v_float)->coord[i]->im);

    fputs("\n", outfile);

    fputs("(t_i, w_i)\n", outfile);
    for (i = 0; i < num_points; i++) {
        mpfr_fprintf(outfile, "%.15Re, ", t_float[i]);
        print_points_float(outfile, w_float[i]);
    }
}

/*********************************
 * print a test point to outfile *
 *********************************/
void print_points_float(FILE *outfile, complex_vector points) {
    int i, num_var = points->size;

    fprintf(outfile, "[");
    for (i = 0; i < num_var - 1; i++) {
        if (mpf_cmp_ui(points->coord[i]->im, 0) == 0)
            mpfr_fprintf(outfile, "%.15Re, ", points->coord[i]->re);
        else
            mpfr_fprintf(outfile, "%.15Re + I * %.15Re, ", points->coord[i]->re, points->coord[i]->im);
    }
    if (mpf_cmp_ui(points->coord[i]->im, 0) == 0)
        mpfr_fprintf(outfile, "%.15Re]\n", points->coord[i]->re);
    else
        mpfr_fprintf(outfile, "%.15Re + I * %.15Re]\n", points->coord[i]->re, points->coord[i]->im);
}

/*****************************************
 * print interval information to outfile *
 *****************************************/
void fprint_interval_float(FILE *outfile, mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right) {
    int i;

    for (i = 0; i < BH_TERMWIDTH; i++)
        fprintf(outfile, "=");

    fprintf(outfile, "\n");
    mpfr_fprintf(outfile, "t interval: [%.15Re, %.15Re]\n", t_left, t_right);
    mpfr_fprintf(outfile, "x_left: ");
    print_points_float(outfile, w_left);

    mpfr_fprintf(outfile, "alpha (x_left): %.15Re\n", alpha_left);
    mpfr_fprintf(outfile, "beta (x_left): %.15Re\n", beta_left);
    mpfr_fprintf(outfile, "gamma (x_left): %.15Re\n", gamma_left);

    mpfr_fprintf(outfile, "x_right: ");
    print_points_float(outfile, w_right);

    mpfr_fprintf(outfile, "alpha (x_right): %.15Re\n", alpha_right);
    mpfr_fprintf(outfile, "beta (x_right): %.15Re\n", beta_right);
    mpfr_fprintf(outfile, "gamma (x_right): %.15Re\n", gamma_right);
}

/**************************************
 * print continuous intervals to file *
 **************************************/
void fprint_continuous_float(mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right) {
    int res;

    errno = 0;
    FILE *fh = fopen(BH_FCONT, "a");

    if (fh == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprint_interval_float(fh, t_left, t_right, w_left, w_right, alpha_left, alpha_right, beta_left, beta_right, gamma_left, gamma_right);

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
void fprint_discontinuous_float(mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right) {
    int res;

    errno = 0;
    FILE *fh = fopen(BH_FDISCONT, "a");

    if (fh == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FDISCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprint_interval_float(fh, t_left, t_right, w_left, w_right, alpha_left, alpha_right, beta_left, beta_right, gamma_left, gamma_right);

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
void fprint_uncertain_float(mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right) {
    int res;

    errno = 0;
    FILE *fh = fopen(BH_FUNCERTIFIED, "a");

    if (fh == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FUNCERTIFIED, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(fh, "Uncertain of interval:\n");
    fprint_interval_float(fh, t_left, t_right, w_left, w_right, alpha_left, alpha_right, beta_left, beta_right, gamma_left, gamma_right);
    fputs("\n", fh);

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
void fprint_solutions_float(void *t, void *w, int num_points) {
    int i, j, res;

    mpf_t *t_float = (mpf_t *) t;
    complex_vector *w_float = (complex_vector *) w;

    errno = 0;
    FILE *outfile = fopen(BH_FPTSOUT, "w");

    if (outfile == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open output file `%s': %s", BH_FPTSOUT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(outfile, "%d\n", num_points);

    for (i = 0; i < num_points; i++) {
        mpfr_fprintf(outfile, "%.15Re\n", t_float[i]);

        for (j = 0; j < w_float[i]->size; j++)
            mpfr_fprintf(outfile, "%.15Re\t%.15Re\n", w_float[i]->coord[j]->re, w_float[i]->coord[j]->im);
    }

    errno = 0;
    res = fclose(outfile);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Couldn't close `%s': %s", BH_FPTSOUT, strerror(errno));
        exit(BH_EXIT_BADREAD);
    }
}
