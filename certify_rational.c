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


/* here's what we really want to do:
 * take each t_i, t_{i+1} and get the alpha, beta, gamma values for each
 * then test that each is in convergence zone
 * if so, great, test for continuity btw w_i and w_{i+1}
 * if not, subdivide and retest
 */
void test_pairwise_rational(polynomial_system *system, configurations *config, void *v, void *t, void *w, int num_points) {
    int i, num_var = system->numVariables;
    rational_complex_vector *v_rational = (rational_complex_vector *) v;
    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *w_rational = (rational_complex_vector *) w;

    /* for each t_i, t_{i+1} */
    /* this is the part that needs to get parallelized */
    for (i = 0; i < num_points - 1; i++) {
        int w1_is_solution, w2_is_solution;
        polynomial_system F1, F2;
        rational_complex_number alpha1, beta1, gamma1, alpha2, beta2, gamma2;
        initialize_rational_number(alpha1);
        initialize_rational_number(beta1);
        initialize_rational_number(gamma1);
        initialize_rational_number(alpha2);
        initialize_rational_number(beta2);
        initialize_rational_number(gamma2);

        apply_tv_rational(system, &F1, t_rational[i], w_rational[i]);
        apply_tv_rational(system, &F2, t_rational[i+1], w_rational[i+1]);

        w1_is_solution = get_alpha_beta_gamma(w_rational[i], &F1, &alpha1, &beta1, &gamma1);
        w2_is_solution = get_alpha_beta_gamma(w_rational[i+1], &F2, &alpha2, &beta2, &gamma2);

        print_system_rational(&F1);
        print_points_rational(&w_rational[i], num_var);
        gmp_printf("alpha: %Qd\nbeta: %Qd\ngamma: %Qd\n", alpha1->re, beta1->re, gamma1->re);

        if (w1_is_solution)
            printf("P1 is an approximate solution\n\n");
        else 
            printf("P1 is not an approximate solution\n\n");

        print_system_rational(&F2);
        print_points_rational(&w_rational[i+1], num_var);
        gmp_printf("alpha: %Qd\nbeta: %Qd\ngamma: %Qd\n", alpha2->re, beta2->re, gamma2->re);

        if (w2_is_solution)
            printf("P2 is an approximate solution\n\n");
        else 
            printf("P2 is not an approximate solution\n\n");

        /* need to free F1/F2 */
    }
}

/*************************************
 * get alpha, beta, gamma values, &c *
 *************************************/
int get_alpha_beta_gamma(rational_complex_vector points, polynomial_system *F, rational_complex_number *alpha, rational_complex_number *beta, rational_complex_number *gamma) {
    int i, rV, numApproxSolns, num_var = F->numVariables;

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

    /* compute alpha^2, beta^2, & gamma^2 (and save to original values of alpha^2, beta^2, & gamma^2) */
    rV = compute_alpha_beta_gamma_sqr_rational(P.Nx, P.alpha_sqr, P.beta_sqr, P.gamma_sqr, F, P.x);
    set_rational_number(*alpha, P.alpha_sqr);
    set_rational_number(*beta, P.beta_sqr);
    set_rational_number(*gamma, P.gamma_sqr);

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
}


/*
    if (S->algorithm >= 1)
    { // now that we have approximate solutions, isolate them
    printf("Isolating %d approximate solution%s.\n\n", numApproxSolns, numApproxSolns == 1 ? "" : "s");
    numDistinctSolns = isolate_approximate_solutions_rational(numPoints, Points_struct, F);

    // now that we have distinct ones, determine which ones are real
    if (S->algorithm >= 2 && F->isReal)
    { // print message and do the analysis
      printf("Classifying %d distinct approximate solution%s.\n\n", numDistinctSolns, numDistinctSolns == 1 ? "" : "s");
      if (S->realityTest)
      { // use global approach
        numRealSolns = classify_real_points_global_rational(numPoints, Points_struct, F);
      }
      else
      { // use local approach
        numRealSolns = classify_real_points_rational(numPoints, Points_struct, F);
      }
    }
    }

    // refine the solutions
    refine_points_rational(numPoints, Points_struct, F, S->refineDigits);

    // print the data out
    classify_rational_output(numPoints, Points_struct, numApproxSolns, numDistinctSolns, numRealSolns, F->isReal, S, 0, F);

    // clear Points_struct
    for (i = 0; i < numPoints; i++)
    clear_rational_point_struct(&Points_struct[i]);
    free(Points_struct);
    Points_struct = NULL;

    return;
}
*/
