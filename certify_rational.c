/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * certify_rational.c: Rational-valued functions related to certification
 */
#include "blueharvest.h"

/***************************************
 * print a polynomial system to a file *
 ***************************************/
void print_system_file_rational(polynomial_system *F, int counter) {
    int i, j, k;
    char filename[BH_MAX_FILENAME];
    snprintf(filename, (size_t) BH_MAX_FILENAME + 1, "system%d", counter);

    /* check for errors here later */
    FILE *OUT = fopen(filename, "w");

    /* print number of variables, polynomials */
    fprintf(OUT, "%d\t%d\n\n", F->numVariables, F->numPolynomials);
    
    /* print each monomial degree and coefficient */
    for (i = 0; i < F->numPolynomials; i++) {
        polynomial p = F->polynomials[i];
        fprintf(OUT, "%d\n", p.numTerms);
        for (j = 0; j < p.numTerms; j++) {
            for (k = 0; k < p.numVariables; k++) {
                gmp_fprintf(OUT, "%d\t", p.exponents[j][k]);
            }

            gmp_fprintf(OUT, "%Qd\t%Qd\n", p.coeff[j]->re, p.coeff[j]->im);
        }
    }

    fclose(OUT);
}

/***********************************
 * print a set of points to a file *
 ***********************************/
void print_points_file_rational(rational_complex_vector w, int num_var, int counter) {
    int i;
    char filename[BH_MAX_FILENAME];
    snprintf(filename, (size_t) BH_MAX_FILENAME + 1, "points%d", counter);

    /* check for errors here later */
    FILE *OUT = fopen(filename, "w");

    /* print number of points */
    fprintf(OUT, "1\n");
    
    /* print each monomial degree and coefficient */
    for (i = 0; i < num_var; i++)
        gmp_fprintf(OUT, "%Qd\t%Qd\n", w->coord[i]->re, w->coord[i]->im);

    fclose(OUT);
}
/**************************
 * apply f + tv for all t *
 **************************/
void deform_rational(polynomial_system *system, configurations *config, void *v, void *t, void *w, int num_points) {
    int i, j, k, h, num_var, num_terms;
    rational_complex_vector *v_rational = (rational_complex_vector *) v;
    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *w_rational = (rational_complex_vector *) w;

    num_var = system->numVariables;

    /* for each t */
    for (i = 0; i < num_points; i++) {
        polynomial_system Fn;
        Fn.numVariables = num_var;
        Fn.numPolynomials = num_var;
        Fn.maximumDegree = system->maximumDegree;
        Fn.isReal = 0;
        Fn.numExponentials = 0;
        Fn.polynomials = malloc(num_var * sizeof(polynomial));

        /* for each polynomial in F */
        for (j = 0; j < num_var; j++) {
            /* copy system->p[j] */
            polynomial p;
            num_terms = system->polynomials[j].numTerms + 1;
            p.numTerms = num_terms;
            p.numVariables = num_var;
            p.degree = system->polynomials[j].degree;

            /* copy exponents from system->polynomials[j] */
            p.exponents = malloc(num_terms * sizeof(int *));
            for (k = 0; k < num_terms - 1; k++) {
                p.exponents[k] = malloc(num_var * sizeof(int));
                for (h = 0; h < num_var; h++)
                    p.exponents[k][h] = system->polynomials[j].exponents[k][h];
            }
            
            /* add 0's for exponents in last (constant) term */
            p.exponents[k] = malloc(num_var * sizeof(int));
            for (h = 0; h < num_var; h++) {
                p.exponents[k][h] = 0;
            }

            /* copy the coefficients */
            p.coeff = malloc(num_terms * sizeof(rational_complex_number));

            for (k = 0; k < num_terms - 1; k++) {
                initialize_rational_number(p.coeff[k]);
                set_rational_number(p.coeff[k], system->polynomials[j].coeff[k]);
            }

            /* now add tv */
            initialize_rational_number(p.coeff[k]);
            mpq_mul(p.coeff[k]->re, t_rational[i], (*v_rational)->coord[i]->re);
            mpq_mul(p.coeff[k]->im, t_rational[i], (*v_rational)->coord[i]->im);

            Fn.polynomials[j] = p;
        }

        /*
        print_system_file_rational(&Fn, i + 1);
        print_points_file_rational(w_rational[i], num_var, i + 1);
        */

        classify_points_rational(1, &w_rational[i], &Fn, &config);
    }
}
