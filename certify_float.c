/*
 * Cadenza
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * certify_float.c: Floating-point functions related to certification
 */
#include "cadenza.h"

/****************************
 * get an mpq_t from mpfr_t *
 ****************************/
void mpq_set_r(mpq_t *q, mpfr_t  r) {
    int i, base = 10;
    mpfr_exp_t expptr;
    char *str_mp = mpfr_get_str(NULL, &expptr, base, 0, r, MPFR_RNDN);

    if (expptr < 0) {
        int exp_int = -1 * (int) expptr;

        exp_int += strlen(str_mp);

        if (str_mp[0] == '-')
            exp_int--;

        mpq_t denominator;
        mpq_t one_tenth;

        mpq_init(one_tenth);
        mpq_init(denominator);

        mpq_set_str(*q, str_mp, base);
        mpq_set_ui(one_tenth, 1, 10);
        mpq_set_ui(denominator, 1, 1);

        for (i = 0; i < exp_int; i++)
            mpq_mul(denominator, denominator, one_tenth);

        mpq_mul(*q, *q, denominator);

        mpq_clear(denominator);
        mpq_clear(one_tenth);
    } else {
        int exp_int = (int) expptr;

        if (str_mp[0] == '-')
            exp_int++;

        exp_int = strlen(str_mp) - exp_int;

        mpq_set_str(*q, str_mp, base);
        mpq_t denominator;
        mpq_t one_tenth;

        mpq_init(one_tenth);
        mpq_init(denominator);

        mpq_set_str(*q, str_mp, base);
        mpq_set_ui(one_tenth, 1, 10);
        mpq_set_ui(denominator, 1, 1);

        for (i = 0; i < exp_int; i++)
            mpq_mul(denominator, denominator, one_tenth);

        mpq_mul(*q, *q, denominator);

        mpq_clear(denominator);
        mpq_clear(one_tenth);
    }

    mpfr_free_str(str_mp);
}

/*****************************
 * compute ||J(f)^-1 * v||^2 *
 *****************************/
void compute_norm_Jv_float(polynomial_system *F, complex_vector v, complex_vector w, mpf_t *norm_Jv) {
    int retval = 0, *rowswaps = NULL;

    mpf_t pivot_tol, pivot_drop_tol;
    mpf_init(pivot_tol);
    mpf_init(pivot_drop_tol);

    determine_pivot_tolerances(pivot_tol, pivot_drop_tol, default_precision);

    complex_vector f;
    /* X = J^-1 * v */
    complex_vector X;

    complex_matrix J;
    complex_matrix LU;

    initialize_vector(f, 0);
    initialize_vector(X, 0);

    initialize_matrix(J, 0, 0);
    initialize_matrix(LU, 0, 0);

    eval_polynomial_system(f, J, F, w, default_precision);

    retval = LUdecomp(LU, &rowswaps, J, pivot_tol, pivot_drop_tol, default_precision);

    /* success */
    if (!retval) {
        LUsolve_vector(X, LU, rowswaps, v, default_precision);
        norm_vector(*norm_Jv, X);
    } else {
        mpf_set_ui(*norm_Jv, 0); /* this may need to be something else */
    }

    mpf_clear(pivot_tol);
    mpf_clear(pivot_drop_tol);
    clear_vector(f);
    clear_vector(X);
    clear_matrix(J);
    clear_matrix(LU);

    free(rowswaps);
    rowswaps = NULL;
}

/*********************************************
 * return 1 if interval certifies continuous *
 *********************************************/
int is_continuous_float(complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, polynomial_system F_left, polynomial_system F_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t gamma_left, mpf_t gamma_right) {
    int retval = 0, left_success = 0, right_success = 0;

    mpf_t alpha_star;
    mpf_t t_dist;

    mpf_init(alpha_star);
    mpf_init(t_dist);

    mpf_set_ui(alpha_star, 4);
    mpf_div_ui(alpha_star, alpha_star, 100);

    mpf_sub(t_dist, t_right, t_left); /* t_2 -t_1 */

    /* we want either of magic_number_left or _right to be less than alpha_star */
    mpf_t magic_number_left;
    mpf_t norm_Jv_left;
    mpf_t magic_number_right;
    mpf_t norm_Jv_right;

    mpf_init(magic_number_left);
    mpf_init(norm_Jv_left);
    mpf_init(magic_number_right);
    mpf_init(norm_Jv_right);

    compute_norm_Jv_float(&F_left, v, w_left, &norm_Jv_left);
    compute_norm_Jv_float(&F_right, v, w_right, &norm_Jv_right);

    /* ||J^-1 * v|| */
    mpf_set(magic_number_left, norm_Jv_left);
    mpf_set(magic_number_right, norm_Jv_right);

    /* ||J^-1 * v|| * g */
    mpf_mul(magic_number_left, magic_number_left, gamma_left);
    mpf_mul(magic_number_right, magic_number_right, gamma_right);

    /* ||J^-1 * v|| * g * |t_1 - t_2|*/
    mpf_mul(magic_number_left, magic_number_left, t_dist);
    mpf_mul(magic_number_right, magic_number_right, t_dist);

    /* a + ||J^-1 * v|| * g * |t_1 - t_2|*/
    mpf_add(magic_number_left, alpha_left, magic_number_left);
    mpf_add(magic_number_right, alpha_right, magic_number_right);

    /* run the tests */
    left_success = (mpf_cmp(magic_number_left, alpha_star) <= 0);
    right_success = (mpf_cmp(magic_number_right, alpha_star) <= 0);

    /* test fails for both left and right */
    if (left_success == 0 && right_success == 0)
        retval = 0;
    else {
        mpf_t gamma_inv;
        mpf_t points_dist;
        mpf_t magic_number;

        mpf_init(gamma_inv);
        mpf_init(points_dist);
        mpf_init(magic_number);

        complex_vector points_dist_vector;
        initialize_vector(points_dist_vector, w_left->size);

        /* test succeeds for both left and right */
        if (left_success && right_success) {
            mpf_set_min(gamma_inv, gamma_left, gamma_right); /* choose the lowest value for gamma */
            mpf_ui_div(gamma_inv, 1, gamma_inv); /* set gamma to 1/gamma */
        } else if (left_success) { /* test succeeds for left only */
            mpf_ui_div(gamma_inv, 1, gamma_left); /* set gamma to 1/gamma */
        } else { /* test succeeds for right only */
            mpf_ui_div(gamma_inv, 1, gamma_right); /* set gamma to 1/gamma */
        }

        /* magic number is now to be 79/(1000*gamma) */
        mpf_set_ui(magic_number, 79);
        mpf_div_ui(magic_number, magic_number, 1000);
        mpf_mul(magic_number, magic_number, gamma_inv);

        subtract_vector(points_dist_vector, w_left, w_right);
        norm_vector(points_dist, points_dist_vector); /* ||x_1 - x_2|| */

        if (mpf_cmp(points_dist, magic_number) <= 0) {
            retval = 1;
        } else {
            retval = 0;
        }

        mpf_clear(points_dist);
        mpf_clear(magic_number);
        mpf_clear(gamma_inv);

        clear_vector(points_dist_vector);
    }

    mpf_clear(alpha_star);
    mpf_clear(t_dist);
    mpf_clear(magic_number_left);
    mpf_clear(magic_number_right);
    mpf_clear(norm_Jv_left);
    mpf_clear(norm_Jv_right);

    return retval;
}

/************************************************
 * return 1 if interval certifies discontinuous *
 ************************************************/
int is_not_continuous_float(complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, polynomial_system F_left, polynomial_system F_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t gamma_left, mpf_t gamma_right) {
    int retval = 0, beta_retval = 0, left_success = 0, right_success = 0;

    mpf_t alpha_star;
    mpf_t t_dist;

    mpf_init(alpha_star);
    mpf_init(t_dist);

    mpf_set_ui(alpha_star, 4);
    mpf_div_ui(alpha_star, alpha_star, 100);

    mpf_sub(t_dist, t_right, t_left); /* t_2 -t_1 */

    /* we want either of magic_number_left or _right to be less than alpha_star */
    mpf_t magic_number_left;
    mpf_t norm_Jv_left;
    mpf_t magic_number_right;
    mpf_t norm_Jv_right;

    mpf_init(magic_number_left);
    mpf_init(norm_Jv_left);
    mpf_init(magic_number_right);
    mpf_init(norm_Jv_right);

    compute_norm_Jv_float(&F_left, v, w_left, &norm_Jv_left);
    compute_norm_Jv_float(&F_right, v, w_right, &norm_Jv_right);

    /* ||J^-1 * v|| */
    mpf_set(magic_number_left, norm_Jv_left);
    mpf_set(magic_number_right, norm_Jv_right);

    /* ||J^-1 * v|| * g */
    mpf_mul(magic_number_left, magic_number_left, gamma_left);
    mpf_mul(magic_number_right, magic_number_right, gamma_right);

    /* ||J^-1 * v|| * g * |t_1 - t_2|*/
    mpf_mul(magic_number_left, magic_number_left, t_dist);
    mpf_mul(magic_number_right, magic_number_right, t_dist);

    /* a + ||J^-1 * v|| * g * |t_1 - t_2|*/
    mpf_add(magic_number_left, alpha_left, magic_number_left);
    mpf_add(magic_number_right, alpha_right, magic_number_right);

    /* run the tests */
    left_success = (mpf_cmp(magic_number_left, alpha_star) <= 0);
    right_success = (mpf_cmp(magic_number_right, alpha_star) <= 0);

    /* test fails for both left and right */
    if (left_success == 0 && right_success == 0)
        retval = 0;
    else {
        complex_number beta_complex1;
        complex_number beta_complex2;

        mpf_t beta1;
        mpf_t beta2;
        mpf_t beta_sum;
        mpf_t points_dist;

        initialize_number(beta_complex1);
        initialize_number(beta_complex2);

        mpf_init(beta1);
        mpf_init(beta2);
        mpf_init(beta_sum);
        mpf_init(points_dist);

        /* test succeeds for at least left (possibly right, but use left anyway) */
        if (left_success) {
            beta_retval = compute_beta(beta_complex1, &F_left, w_left, default_precision);

            if (!beta_retval)
                mpf_set(beta1, beta_complex1->re);

            beta_retval = compute_beta(beta_complex2, &F_left, w_right, default_precision);

            if (!beta_retval)
                mpf_set(beta2, beta_complex2->re);
        } else { /* test succeeds for right only */
            beta_retval = compute_beta(beta_complex1, &F_right, w_left, default_precision);

            if (!beta_retval)
                mpf_set(beta1, beta_complex1->re);

            beta_retval = compute_beta(beta_complex2, &F_right, w_right, default_precision);

            if (!beta_retval)
                mpf_set(beta2, beta_complex2->re);
        }
        
        /* 2 * (beta1 + beta2) */
        mpf_add(beta_sum, beta1, beta2);
        mpf_mul_ui(beta_sum, beta_sum, 2);

        complex_vector points_dist_vector;
        initialize_vector(points_dist_vector, w_left->size);

        subtract_vector(points_dist_vector, w_left, w_right);
        norm_vector(points_dist, points_dist_vector); /* ||x_1 - x_2|| */

        if (mpf_cmp(beta_sum, points_dist) <= 0)
            retval = 1;
        else
            retval = 0;

        clear_number(beta_complex1);
        clear_number(beta_complex2);

        clear_vector(points_dist_vector);

        mpf_clear(beta1);
        mpf_clear(beta2);
        mpf_clear(beta_sum);
        mpf_clear(points_dist);
    }

    mpf_clear(alpha_star);
    mpf_clear(t_dist);
    mpf_clear(magic_number_left);
    mpf_clear(magic_number_right);
    mpf_clear(norm_Jv_left);
    mpf_clear(norm_Jv_right);

    return retval;
}

/**********************************************
 * test for the continuity of a given segment *
 **********************************************/
int test_continuity_float(complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, polynomial_system F_left, polynomial_system F_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t gamma_left, mpf_t gamma_right) {
    if (is_continuous_float(v, t_left, t_right, w_left, w_right, F_left, F_right, alpha_left, alpha_right, gamma_left, gamma_right))
        return 1;
    else if (is_not_continuous_float(v, t_left, t_right, w_left, w_right, F_left, F_right, alpha_left, alpha_right, gamma_left, gamma_right))
        return -1;

    return 0;
}

/**************************************************************************************
 * subdivide a segment (t_1, w_1) <=> (t_2, w_2) for ((t_1 + t_2)/2, (w_1 + w_2) / 2) *
 **************************************************************************************/
void subdivide_segment_float(polynomial_system *base, complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, mpf_t *t_mid, complex_vector *w_mid, int num_var) {
    if (verbosity > BH_VERBOSE)
        printf("subdividing... ");

    int i, retval = 0, *rowswaps = NULL;

    mpf_t beta;
    mpf_init(beta);

    complex_number one_half_complex;
    initialize_number(one_half_complex);

    complex_number beta_complex;
    initialize_number(beta_complex);

    mpf_set_d(one_half_complex->re, 0.5);
    mpf_set_ui(one_half_complex->im, 0);

    mpf_set(beta_complex->re, beta);
    mpf_set_ui(beta_complex->im, 0);

    complex_matrix LU;
    initialize_matrix(LU, 0, 0);

    polynomial_system F_mid;

    complex_vector new_point;

    initialize_vector(new_point, num_var);

    mpf_t pivot_tol, pivot_drop_tol;
    mpf_init(pivot_tol);
    mpf_init(pivot_drop_tol);

    setPrec_number(beta_complex, default_precision);
    determine_pivot_tolerances(pivot_tol, pivot_drop_tol, default_precision);

    /* t_mid */
    mpf_add(*t_mid, t_left, t_right);
    mpf_div_ui(*t_mid, *t_mid, 2);

    /* w_mid */
    for (i = 0; i < num_var; i++) {
        add((*w_mid)->coord[i], w_left->coord[i], w_right->coord[i]);
        multiply((*w_mid)->coord[i], (*w_mid)->coord[i], one_half_complex);
    }

    /* apply Newton to w_mid */
    apply_tv_float(base, &F_mid, *t_mid, v);
    retval = newton_iteration(new_point, beta_complex, LU, &rowswaps, &F_mid, *w_mid, pivot_tol, pivot_drop_tol, default_precision);
    if (retval == ERROR_LU_DECOMP) {
        mpfr_fprintf(stderr, "singularity suspected at t = %.15Re near point ", *t_mid);
        print_points_float(stderr, new_point);
        fprintf(stderr, ERROR_MESSAGE);
        exit(BH_EXIT_OTHER);
    }

    copy_vector(*w_mid, new_point);

    if (verbosity > BH_VERBOSE)
        mpfr_printf("new intervals are [%.15Re, %.15Re] and [%.15Re, %.15Re]\n", t_left, *t_mid, *t_mid, t_right);

    mpf_clear(beta);
    mpf_clear(pivot_tol);
    mpf_clear(pivot_drop_tol);
    clear_number(one_half_complex);
    clear_number(beta_complex);
    clear_matrix(LU);
    clear_polynomial_system(&F_mid);
    clear_vector(new_point);

    free(rowswaps);
    rowswaps = NULL;
}

/*******************************
 * apply f + tv for a single t *
 *******************************/
void apply_tv_float(polynomial_system *base, polynomial_system *F, mpf_t t, complex_vector v) {
    int i, j, k, num_terms, num_var, num_poly;

    /* convert t to t_rational */
    mpq_t t_rational;
    mpq_init(t_rational);
    mpq_set_r(&t_rational, t);

    /* convert v to v_rational */
    rational_complex_vector v_rational;
    initialize_rational_vector(v_rational, v->size);
    for (i = 0; i < v->size; i++) {
        mpq_set_r(&v_rational->coord[i]->re, v->coord[i]->re);
        mpq_set_r(&v_rational->coord[i]->im, v->coord[i]->im);
    }

    /* get information from base */
    num_var = base->numVariables;
    num_poly = base->numPolynomials;

    errno = 0;

    F->numVariables = num_var;
    F->numPolynomials = num_poly;
    F->maximumDegree = base->maximumDegree;
    F->isReal = 0;
    F->numExponentials = 0;
    F->polynomials = malloc(num_var * sizeof(polynomial));
    if (F->polynomials == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
        print_error(error_string, stderr);

        exit(BH_EXIT_MEMORY);
    }

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
        errno = 0;

        p.exponents = malloc(num_terms * sizeof(int *));
        if (p.exponents == NULL) {
            snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
            print_error(error_string, stderr);

            exit(BH_EXIT_MEMORY);
        }

        for (j = 0; j < num_terms - 1; j++) {
            p.exponents[j] = malloc(num_var * sizeof(int));
            if (p.exponents[j] == NULL) {
                snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
                print_error(error_string, stderr);

                exit(BH_EXIT_MEMORY);
            }

            for (k = 0; k < num_var; k++)
                p.exponents[j][k] = base->polynomials[i].exponents[j][k];
        }
        
        /* add 0's for exponents in last (constant) term */
        p.exponents[j] = malloc(num_var * sizeof(int));
        if (p.exponents[j] == NULL) {
            snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
            print_error(error_string, stderr);

            exit(BH_EXIT_MEMORY);
        }

        for (k = 0; k < num_var; k++)
            p.exponents[j][k] = 0;

        /* copy the coefficients */
        p.coeff = malloc(num_terms * sizeof(rational_complex_number));
        if (p.coeff == NULL) {
            snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
            print_error(error_string, stderr);

            exit(BH_EXIT_MEMORY);
        }

        for (j = 0; j < num_terms - 1; j++) {
            initialize_rational_number(p.coeff[j]);
            set_rational_number(p.coeff[j], base->polynomials[i].coeff[j]);
        }

        /* now add tv */
        initialize_rational_number(p.coeff[j]);
        mpq_mul(p.coeff[j]->re, t_rational, v_rational->coord[i]->re);
        mpq_mul(p.coeff[j]->im, t_rational, v_rational->coord[i]->im);

        mpq_init(p.norm_sqr);
        norm_sqr_polynomial(p.norm_sqr, &p);

        F->polynomials[i] = p;
    }

    /* clean up */
    mpq_clear(t_rational);
    clear_rational_vector(v_rational);
}

/*********************************
 * get alpha, beta, gamma values *
 *********************************/
int compute_abg_float(complex_vector points, polynomial_system *F, mpf_t *alpha, mpf_t *beta, mpf_t *gamma) {
    int res = 0, retval = 0, num_approx_solns = 0, num_var = F->numVariables;

    point_struct P;
    /* setup point struct and determine if it is an approximate solution */
    initialize_point_struct(&P, num_var);

    /* set to active */
    P.isActive = 1;

    /* copy point */
    copy_vector(P.origX, points);
    copy_vector(P.x, points);

    /* compute ||x||_2 */
    norm_vector(P.norm_x, P.x);

    /* compute alpha, beta, & gamma */
    res = compute_alpha_beta_gamma(P.Nx, P.alpha, P.beta, P.gamma, F, P.x, default_precision);

    mpf_set(*alpha, P.alpha->re);
    mpf_set(*beta, P.beta->re);
    mpf_set(*gamma, P.gamma->re);

    /* check to see if we have successfully computed alpha, beta, & gamma */
    if (res == EXACT_SOLUTION_LU_ERROR) {
        /* exact solution */
        retval = 1;
    } else if (res) {
        /* unknown */
        retval = 0;
    } else {
        /* determine if alpha is small enough to be an approximate solution */
        num_approx_solns += P.isApproxSoln = determine_approximate_solution(P.alpha);
        retval = num_approx_solns;
    }

    clear_point_struct(&P);

    return retval;
}

/***********************************************************************************
 * test that each w is in the quadratic convergence basin and check for continuity *
 ***********************************************************************************/
void test_pairwise_float(polynomial_system *system, complex_vector *v, mpf_t t_left, mpf_t t_right, complex_vector w_left, complex_vector w_right, int num_var, int iter, mpf_t **t_final, complex_vector **w_final, complex_vector **sing, int *tested, int *succeeded, int *failed, int *num_sing, int check_left) {
    int w_left_solution, w_right_solution, seg_continuous;
    polynomial_system F_left, F_right;
    mpf_t alpha_left, beta_left, gamma_left, alpha_right, beta_right, gamma_right, beta_min;
    mpf_init(alpha_left);
    mpf_init(beta_left);
    mpf_init(gamma_left);
    mpf_init(alpha_right);
    mpf_init(beta_right);
    mpf_init(gamma_right);
    mpf_init(beta_min);

    apply_tv_float(system, &F_left, t_left, *v);
    apply_tv_float(system, &F_right, t_right, *v);

    w_left_solution = compute_abg_float(w_left, &F_left, &alpha_left, &beta_left, &gamma_left);
    w_right_solution = compute_abg_float(w_right, &F_right, &alpha_right, &beta_right, &gamma_right);

    /* test if we've been through this too many times */
    if (iter > subd_tolerance) {
        *tested += 1;

        if (verbosity > BH_LACONIC)
            mpfr_printf("unsure whether segment [%.15Re, %.15Re] is continuous\n", t_left, t_right);
        fprint_uncertain_float(t_left, t_right, w_left, w_right, alpha_left, alpha_right, beta_left, beta_right, gamma_left, gamma_right);

        /* alert user to singularities */
        if (check_left) {
            if (mpfr_inf_p(gamma_left)) {
                mpfr_fprintf(stderr, "singularity found for t = %.15Re at point ", t_left);
                print_points_float(stderr, w_left);
                *num_sing += 1;
            }
        }
        if (mpfr_inf_p(gamma_right)) {
            mpfr_fprintf(stderr, "singularity found for t = %.15Re at point ", t_right);
            print_points_float(stderr, w_right);
            *num_sing += 1;
        }

        return;
    }

    mpf_set_min(beta_min, beta_left, beta_right);

    /* perform Newton iterations on non-compliant points */
    int newton_counter = 1;
    while (!w_left_solution || !w_right_solution) {
        if (newton_counter == newton_tolerance) {
            char err_string[] = "Points not in convergence basin; aborting";
            print_error(err_string, stderr);

            exit(BH_EXIT_OTHER);
        }

        int retval, *rowswaps = NULL;
        complex_number beta_complex;
        complex_vector new_point;
        complex_matrix LU;

        /* perform Newton iterations if w_left is not a solution */
        if (!w_left_solution) {
            if (verbosity > BH_VERBOSE) {
                if (newton_counter % 10 == 1 && newton_counter % 100 != 11)
                    fprintf(stderr, "performing %dst Newton iteration on w_left", newton_counter);
                else if (newton_counter % 10 == 2 && newton_counter % 100 != 12)
                    fprintf(stderr, "performing %dnd Newton iteration on w_left", newton_counter);
                else if (newton_counter % 10 == 3 && newton_counter % 100 != 13)
                    fprintf(stderr, "performing %drd Newton iteration on w_left", newton_counter);
                else
                    fprintf(stderr, "performing %dth Newton iteration on w_left", newton_counter);
                
                mpfr_fprintf(stderr, " (interval: [%.15Re, %.15Re])\n", t_left, t_right);
            }

            mpf_t pivot_tol, pivot_drop_tol;
            mpf_init(pivot_tol);
            mpf_init(pivot_drop_tol);

            initialize_number(beta_complex);
            initialize_vector(new_point, num_var);
            initialize_matrix(LU, 0, 0);

            setPrec_number(beta_complex, default_precision);
            determine_pivot_tolerances(pivot_tol, pivot_drop_tol, default_precision);

            retval = newton_iteration(new_point, beta_complex, LU, &rowswaps, &F_left, w_left, pivot_tol, pivot_drop_tol, default_precision);
            if (retval == ERROR_LU_DECOMP) {
                mpfr_fprintf(stderr, "singularity suspected at t = %.15Re near point ", t_left);
                print_points_float(stderr, w_left);
                fprintf(stderr, ERROR_MESSAGE);
                exit(BH_EXIT_OTHER);
            }

            copy_vector(w_left, new_point);

            w_left_solution = compute_abg_float(w_left, &F_left, &alpha_left, &beta_left, &gamma_left);

            mpf_clear(pivot_tol);
            mpf_clear(pivot_drop_tol);

            clear_number(beta_complex);
            clear_vector(new_point);
            clear_matrix(LU);

            free(rowswaps);
            rowswaps = NULL;
        }

        /* perform Newton iterations if w_right is not a solution */
        if (!w_right_solution) {
            if (verbosity > BH_VERBOSE) {
                if (newton_counter % 10 == 1 && newton_counter % 100 != 11)
                    fprintf(stderr, "performing %dst Newton iteration on w_right", newton_counter);
                else if (newton_counter % 10 == 2 && newton_counter % 100 != 12)
                    fprintf(stderr, "performing %dnd Newton iteration on w_right", newton_counter);
                else if (newton_counter % 10 == 3 && newton_counter % 100 != 13)
                    fprintf(stderr, "performing %drd Newton iteration on w_right", newton_counter);
                else
                    fprintf(stderr, "performing %dth Newton iteration on w_right", newton_counter);

                mpfr_fprintf(stderr, " (interval: [%.15Re, %.15Re])\n", t_left, t_right);
            }

            mpf_t pivot_tol, pivot_drop_tol;
            mpf_init(pivot_tol);
            mpf_init(pivot_drop_tol);

            initialize_number(beta_complex);
            initialize_vector(new_point, num_var);
            initialize_matrix(LU, 0, 0);

            setPrec_number(beta_complex, default_precision);
            determine_pivot_tolerances(pivot_tol, pivot_drop_tol, default_precision);

            retval = newton_iteration(new_point, beta_complex, LU, &rowswaps, &F_right, w_right, pivot_tol, pivot_drop_tol, default_precision);
            if (retval == ERROR_LU_DECOMP) {
                mpfr_fprintf(stderr, "singularity suspected at t = %.15Re near point ", t_right);
                print_points_float(stderr, w_right);
                fprintf(stderr, ERROR_MESSAGE);
                exit(BH_EXIT_OTHER);
            }

            copy_vector(w_right, new_point);

            w_right_solution = compute_abg_float(w_right, &F_right, &alpha_right, &beta_right, &gamma_right);

            mpf_clear(pivot_tol);
            mpf_clear(pivot_drop_tol);

            clear_number(beta_complex);
            clear_vector(new_point);
            clear_matrix(LU);

            free(rowswaps);
            rowswaps = NULL;
        }

        newton_counter++;
    }

    if (verbosity > BH_CHATTY)
        mpfr_printf("testing segment [%.15Re, %.15Re]\n", t_left, t_right);

    seg_continuous = test_continuity_float(*v, t_left, t_right, w_left, w_right, F_left, F_right, alpha_left, alpha_right, gamma_left, gamma_right);

    if (seg_continuous == 0) {
        if (verbosity > BH_CHATTY)
            mpfr_printf("unsure whether segment [%.15Re, %.15Re] is continuous\n", t_left, t_right);

        mpf_t t_mid;
        mpf_init(t_mid);
        complex_vector w_mid;
        initialize_vector(w_mid, num_var);

        if (mpfr_cmp(t_left, t_right) == 0) {
            snprintf(error_string, (size_t) termwidth, "Segment is a point; cannot continue\n");
            print_error(error_string, stderr);
            exit(BH_EXIT_MEMORY);
        }

        subdivide_segment_float(system, *v, t_left, t_right, w_left, w_right, &t_mid, &w_mid, num_var);

        /* recurse! */
        test_pairwise_float(system, v, t_left, t_mid, w_left, w_mid, num_var, iter + 1, t_final, w_final, sing, tested, succeeded, failed, num_sing, 0);
        test_pairwise_float(system, v, t_mid, t_right, w_mid, w_right, num_var, iter + 1, t_final, w_final, sing, tested, succeeded, failed, num_sing, 0);

        mpf_clear(t_mid);
        clear_vector(w_mid);
    } else if (seg_continuous == 1) {
        if (verbosity > BH_LACONIC)
            mpfr_printf("segment [%.15Re, %.15Re] is continuous\n", t_left, t_right);

        fprint_continuous_float(t_left, t_right, w_left, w_right, alpha_left, alpha_right, beta_left, beta_right, gamma_left, gamma_right);

        *tested += 1;
        *succeeded += 1;

        /* alert user to singularities */
        if (check_left) {
            if (mpfr_inf_p(gamma_left)) {
                mpfr_fprintf(stderr, "singularity found for t = %.15Re at point ", t_left);
                print_points_float(stderr, w_left);
            }
        }
        if (mpfr_inf_p(gamma_right)) {
            mpfr_fprintf(stderr, "singularity found for t = %.15Re at point ", t_right);
            print_points_float(stderr, w_right);
        }

        /* copy t and w values into t_final and w_final respectively */
        errno = 0;

        mpf_t *t_final_tmp = realloc(*t_final, (size_t) (*succeeded + 1) * sizeof(mpf_t));
        if (t_final_tmp == NULL) {
            snprintf(error_string, (size_t) termwidth, "Couldn't realloc: %s\n", strerror(errno));
            print_error(error_string, stderr);
            exit(BH_EXIT_MEMORY);
        } else {
            *t_final = t_final_tmp;
            mpf_init((*t_final)[*succeeded]);
            mpf_set((*t_final)[*succeeded], t_right);
        }

        errno = 0;

        complex_vector *w_final_tmp = realloc(*w_final, (size_t) (*succeeded + 1) * sizeof(complex_vector));
        if (w_final_tmp == NULL) {
            snprintf(error_string, (size_t) termwidth, "Couldn't realloc: %s\n", strerror(errno));
            print_error(error_string, stderr);
            exit(BH_EXIT_MEMORY);
        } else {
            *w_final = w_final_tmp;
            initialize_vector((*w_final)[*succeeded], w_right->size);
            copy_vector((*w_final)[*succeeded], w_right);
        }
    } else {
        if (verbosity > BH_LACONIC)
            mpfr_printf("segment [%.15Re, %.15Re] is not continuous\n", t_left, t_right);

        fprint_discontinuous_float(t_left, t_right, w_left, w_right, alpha_left, alpha_right, beta_left, beta_right, gamma_left, gamma_right);

        /* alert user to singularities */
        if (check_left) {
            if (mpfr_inf_p(gamma_left)) {
                mpfr_fprintf(stderr, "singularity found for t = %.15Re at point ", t_left);
                print_points_float(stderr, w_left);
            }
        }
        if (mpfr_inf_p(gamma_right)) {
            mpfr_fprintf(stderr, "singularity found for t = %.15Re at point ", t_right);
            print_points_float(stderr, w_right);
        }

        *tested += 1;
        *failed += 1;
    }

    /* free up stuff */
    clear_polynomial_system(&F_left);
    clear_polynomial_system(&F_right);

    mpf_clear(alpha_left);
    mpf_clear(beta_left);
    mpf_clear(gamma_left);
    mpf_clear(alpha_right);
    mpf_clear(beta_right);
    mpf_clear(gamma_right);
    mpf_clear(beta_min);
}

/*******************
 * certify H(x, t) *
 *******************/
void test_system_float(polynomial_system *system, void *v, void *t, void *w, int num_points, void **t_final, void **w_final, void **sing, int *tested, int *succeeded, int *failed, int *num_sing) {
    int i, num_var = system->numVariables;

    complex_vector *v_float = (complex_vector *) v;
    mpf_t *t_float = (mpf_t *) t;
    complex_vector *w_float = (complex_vector *) w;

    mpf_t **t_final_float = (mpf_t **) t_final;
    complex_vector **w_final_float = (complex_vector **) w_final;
    complex_vector **sing_float = (complex_vector **) sing;

    /* set up t_final and w_final */
    *t_final_float = malloc(sizeof(mpf_t));
    if (*t_final_float == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
        print_error(error_string, stderr);

        exit(BH_EXIT_MEMORY);
    }

    mpf_init((*t_final_float)[0]);
    mpf_set((*t_final_float)[0], t_float[0]);

    *w_final_float = malloc(sizeof(complex_vector));
    if (*w_final_float == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
        print_error(error_string, stderr);

        exit(BH_EXIT_MEMORY);
    }

    initialize_vector((*w_final_float)[0], w_float[0]->size);
    copy_vector((*w_final_float)[0], w_float[0]);

    /* for each t_i, t_{i+1} */
    /* this is the part that needs to get parallelized */
    for (i = 0; i < num_points - 1; i++) {
        test_pairwise_float(system, v_float, t_float[i], t_float[i+1], w_float[i], w_float[i+1], num_var, 1, t_final_float, w_final_float, sing_float, tested, succeeded, failed, num_sing, (i == 0));
    }
}
