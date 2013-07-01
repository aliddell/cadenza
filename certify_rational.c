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

/**********************************************
 * test for the continuity of a given segment *
 **********************************************/
int segment_is_continuous_rational(mpq_t t_1, mpq_t t_2, rational_complex_vector w_1, rational_complex_vector w_2) {
    /* this will be a real test eventually */
    return rand() % 2;
}

/**************************************************************************************
 * subdivide a segment (t_1, w_1) <=> (t_2, w_2) for ((t_1 + t_2)/2, (w_1 + w_2) / 2) *
 **************************************************************************************/
void subdivide_segment_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector w_left, rational_complex_vector w_right, mpq_t *t_mid, rational_complex_vector *w_mid, int num_var) {
    if (verbosity > BH_VERBOSE)
        printf("subdividing...\n");
    int i;
    mpq_t one_half_rational;
    mpq_init(one_half_rational);
    mpq_set_ui(one_half_rational, 1, 2);
    rational_complex_number one_half_complex;
    initialize_rational_number(one_half_complex);
    mpq_set_ui(one_half_complex->re, 1, 2);
    mpq_set_ui(one_half_complex->im, 0, 1);

    /* t_mid */
    mpq_add(*t_mid, t_left, t_right);
    mpq_mul(*t_mid, *t_mid, one_half_rational);

    /* w_mid */
    for (i = 0; i < num_var; i++) {
        add_rational((*w_mid)->coord[i], w_left->coord[i], w_right->coord[i]);
        multiply_rational((*w_mid)->coord[i], (*w_mid)->coord[i], one_half_complex);
    }

    mpq_clear(one_half_rational);
    clear_rational_number(one_half_complex);
}

long int get_max_num_points(mpq_t *beta_min) {
    if (mpq_cmp_ui(*beta_min, 0, 1) == 0)
        return 50; /* arbitrary */
    else {
        mpfr_t beta_min_float;
        mpfr_t log_beta;
        mpfr_init(log_beta);
        mpfr_init_set_q(beta_min_float, beta_min, MPFR_RNDU);
        mpfr_log10(log_beta, beta_min_float, MPFR_RNDU);
        long int retval = -1 * mpfr_get_si(log_beta, MPFR_RNDU);
        return retval;
    }
}

/*******************************
 * apply f + tv for a single t *
 *******************************/
void apply_tv_rational(polynomial_system *base, polynomial_system *F, mpq_t t, rational_complex_vector v) {
    int i, j, k, num_terms, num_var, num_poly;
    /* get information from base */
    num_var = base->numVariables;
    num_poly = base->numPolynomials;

    F->numVariables = num_var;
    F->numPolynomials = num_poly;
    F->maximumDegree = base->maximumDegree;
    F->isReal = 0;
    F->numExponentials = 0;
    F->polynomials = malloc(num_var * sizeof(polynomial));
    F->exponentials = NULL;
    mpq_init(F->norm_sqr);

    /* for each polynomial in F */
    for (i = 0; i < num_var; i++) {
        num_terms = base->polynomials[i].numTerms + 1;

        /* copy base->p[i] */
        polynomial p;
        p.numTerms = num_terms;
        p.numVariables = num_var;
        p.degree = base->polynomials[i].degree;

        /* copy exponents from base->polynomials[i] */
        p.exponents = malloc(num_terms * sizeof(int *));
        for (j = 0; j < num_terms - 1; j++) {
            p.exponents[j] = malloc(num_var * sizeof(int));
            for (k = 0; k < num_var; k++)
                p.exponents[j][k] = base->polynomials[i].exponents[j][k];
        }
        
        /* add 0's for exponents in last (constant) term */
        p.exponents[j] = malloc(num_var * sizeof(int));
        for (k = 0; k < num_var; k++)
            p.exponents[j][k] = 0;

        /* copy the coefficients */
        p.coeff = malloc(num_terms * sizeof(rational_complex_number));

        for (j = 0; j < num_terms - 1; j++) {
            initialize_rational_number(p.coeff[j]);
            set_rational_number(p.coeff[j], base->polynomials[i].coeff[j]);
        }

        /* now add tv */
        initialize_rational_number(p.coeff[j]);
        mpq_mul(p.coeff[j]->re, t, v->coord[i]->re);
        mpq_mul(p.coeff[j]->im, t, v->coord[i]->im);

        mpq_init(p.norm_sqr);
        norm_sqr_polynomial(p.norm_sqr, &p);

        F->polynomials[i] = p;
    }
}

/***********************************************************************************
 * test that each w is in the quadratic convergence basin and check for continuity *
 ***********************************************************************************/
void test_pairwise_rational(polynomial_system *system, rational_complex_vector *v, mpq_t t_left, mpq_t t_right, rational_complex_vector w_left, rational_complex_vector w_right, int num_var, int iter) {
    /* test if we've been through this 20 times or what */
    if (iter > BH_TOLERANCE) {
        char error_string[] = "Number of subdivisions exceeds tolerance";
        print_error(error_string);
        exit(BH_EXIT_INTOLERANT);
    }

    int w_left_solution, w_right_solution, seg_continuous;
    polynomial_system F1, F2;
    mpq_t alpha_left, beta_left, gamma_left, alpha_right, beta_right, gamma_right, beta_min;
    mpq_init(alpha_left);
    mpq_init(beta_left);
    mpq_init(gamma_left);
    mpq_init(alpha_right);
    mpq_init(beta_right);
    mpq_init(gamma_right);

    apply_tv_rational(system, &F1, t_left, *v);
    apply_tv_rational(system, &F2, t_right, *v);

    w_left_solution = get_alpha_beta_gamma_rational(w_left, &F1, &alpha_left, &beta_left, &gamma_left);
    w_right_solution = get_alpha_beta_gamma_rational(w_right, &F2, &alpha_right, &beta_right, &gamma_right);

    mpq_set_min(beta_min, beta_left, beta_right);
    long int retval = get_max_num_points(&beta_min);

    if (w_left_solution && w_right_solution) {
        seg_continuous = segment_is_continuous_rational(t_left, t_right, w_left, w_right);

        if (!seg_continuous) {
            mpq_t t_mid;
            mpq_init(t_mid);
            rational_complex_vector w_mid;
            initialize_rational_vector(w_mid, num_var);
            subdivide_segment_rational(t_left, t_right, w_left, w_right, &t_mid, &w_mid, num_var);

            /* recurse! */
            test_pairwise_rational(system, v, t_left, t_mid, w_left, w_mid, num_var, iter + 1);
            test_pairwise_rational(system, v, t_mid, t_right, w_mid, w_right, num_var, iter + 1);
            mpq_clear(t_mid);
            clear_rational_vector(w_mid);
        }
    } else {
        char error_string[] = "Points not in convergence basin; aborting";
        print_error(error_string);

        exit(BH_EXIT_NOCONVERGE);
    }

    /* free up stuff */
    clear_polynomial_system(&F1);
    clear_polynomial_system(&F2);
    mpq_clear(alpha_left);
    mpq_clear(beta_left);
    mpq_clear(gamma_left);
    mpq_clear(alpha_right);
    mpq_clear(beta_right);
    mpq_clear(gamma_right);
    mpq_clear(beta_min);
}

/* here's what we really want to do:
 * take each t_i, t_{i+1} and get the alpha, beta, gamma values for each
 * then test that each is in convergence zone
 * if so, great, test for continuity btw w_i and w_{i+1}
 * if not, subdivide and retest
 */
void test_system_rational(polynomial_system *system, configurations *config, void *v, void *t, void *w, int num_points) {
    int i, num_var = system->numVariables;

    rational_complex_vector *v_rational = (rational_complex_vector *) v;
    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *w_rational = (rational_complex_vector *) w;

    /* for each t_i, t_{i+1} */
    /* this is the part that needs to get parallelized */
    for (i = 0; i < num_points - 1; i++) {
        test_pairwise_rational(system, v_rational, t_rational[i], t_rational[i+1], w_rational[i], w_rational[i+1], num_var, 1);
    }
}

/*************************************
 * get alpha, beta, gamma values, &c *
 *************************************/
int get_alpha_beta_gamma_rational(rational_complex_vector points, polynomial_system *F, mpq_t *alpha, mpq_t *beta, mpq_t *gamma) {
    int rV, numApproxSolns, num_var = F->numVariables;

    rational_point_struct P;
    /* setup point struct and determine if it is an approximate solution */
    initialize_rational_point_struct(&P, num_var);

    /* set to active */
    P.isActive = 1;

    /* copy rational point */
    copy_rational_vector(P.origX, points);
    copy_rational_vector(P.x, points);

    /* compute ||x||_2^2 */
    norm_sqr_rational_vector(P.norm_sqr_x, P.x);

    /* compute alpha^2, beta^2, & gamma^2 */
    rV = compute_alpha_beta_gamma_sqr_rational(P.Nx, P.alpha_sqr, P.beta_sqr, P.gamma_sqr, F, P.x);

    mpq_set(*alpha, P.alpha_sqr);
    mpq_set(*beta, P.beta_sqr);
    mpq_set(*gamma, P.gamma_sqr);
    /* holdover from when we defined alpha, beta and gamma as complex numbers
    set_rational_number(*alpha, P.alpha_sqr);
    set_rational_number(*beta, P.beta_sqr);
    set_rational_number(*gamma, P.gamma_sqr);
    */

    /* check to see if we have successfully computed alpha_sqr, beta_sqr, & gamma_sqr */
    if (rV == EXACT_SOLUTION_LU_ERROR) { // exact solution
      //numApproxSolns += P.isApproxSoln = 1;
      return 1;
    } else if (rV) { // error 
      P.isApproxSoln = 0; // unknown 
      return 0;
    } else { // determine if alpha is small enough to be an approximate solution
      numApproxSolns += P.isApproxSoln = determine_approximate_solution_rational(P.alpha_sqr);
      return numApproxSolns;
    }

    clear_rational_point_struct(&P);
}
