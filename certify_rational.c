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

void deform_rational(polynomial_system *system, void *v, void *t, void *w, int num_points) {
    int i, j, k, h, num_var, num_terms;
    rational_complex_vector *v_rational = (rational_complex_vector *) v;
    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_number *w_rational = (rational_complex_number *) w;

    num_var = (*system).numVariables;
    char variables[] = "xyzwuvabcdjkmnpqrs";
    /* for each t */
    for (i = 0; i < num_points; i++) {
        polynomial_system Fn;
        Fn.numVariables = num_var;
        Fn.numPolynomials = num_var;
        Fn.polynomials = malloc(num_var * sizeof(polynomial));

        printf("i: %d\n", i);
        for (j = 0; j < num_var; j++) {
            /* copy system->p[j] */
            polynomial p;
            num_terms = (*system).polynomials[j].numTerms + 1;
            p.numTerms = num_terms;
            p.numVariables = num_var;

            /* copy exponents from system->polynomials[j] */
            p.exponents = malloc(num_terms * sizeof(int *));
            for (k = 0; k < num_terms - 1; k++) {
                p.exponents[k] = malloc(num_var * sizeof(int));
                for (h = 0; h < num_var; h++)
                    p.exponents[k][h] = (*system).polynomials[j].exponents[k][h];
            }
            
            /* add 0's for exponents in last (constant) term */
            for (h = 0; h < num_var; h++) {
                p.exponents[k] = malloc(num_var * sizeof(int));
                p.exponents[k][h] = 0;
            }

            /* copy the coefficients */
            p.coeff = malloc(num_terms * sizeof(rational_complex_number));

            for (k = 0; k < num_terms - 1; k++) {
                initialize_rational_number(p.coeff[k]);
                set_rational_number(p.coeff[k], (*system).polynomials[j].coeff[k]);
            }

            /* now add tv */
            initialize_rational_number(p.coeff[k]);
            mpq_mul(p.coeff[k], t_rational[i], (*v_rational[i]).coord[j]->re);
            mpq_mul(p.coeff[k], t_rational[i], (*v_rational[i]).coord[j]->im);

            for (k = 0; k < p.numTerms - 1; k++) {
                gmp_printf("(%Qd + %Qdi)", p.coeff[j]->re, p.coeff[j]->im);
                for (h = 0; h < p.numVariables - 1; h++)
                    printf("%c^%d", variables[h], p.exponents[k][h]);
                printf(" + ");
            }
            gmp_printf("(%Qd + %Qdi)", p.coeff[k]->re, p.coeff[k]->im);
            for (h = 0; h < (*system).numVariables; h++)
                printf("%c^%d", variables[h], p.exponents[k][h]);
            puts("");
        }
    }
}
