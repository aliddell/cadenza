/*
 * Cadenza
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * certify_rational.c: Rational-valued functions related to certification
 */
#include "cadenza.h"

/*********************************************************************
 * takes X^2, Y^2, Z^2, returns 1 if X + Y <= Z (where X, Y, Z >= 0) *
 *********************************************************************/
int lte_sqr_rational(mpq_t X_sqr, mpq_t Y_sqr, mpq_t Z_sqr) {
    /* first inequality: X^2 + Y^2 - Z^2 <= 0 */
    mpq_t first;
    mpq_init(first);
    mpq_add(first, X_sqr, Y_sqr);
    mpq_sub(first, first, Z_sqr);

    /* fails */
    if (mpq_cmp_ui(first, 0, 1) >= 0) {
        mpq_clear(first);
        return 0;
    }
    mpq_clear(first);

    /* second inequality: (Z^2 - (X^2 + Y^2))^2 >= 4X^2Y^2 */
    mpq_t second_left;
    mpq_t second_right;
    mpq_t four;

    mpq_init(second_left);
    mpq_init(second_right);
    mpq_init(four);

    mpq_set_ui(four, 4, 1);

    mpq_add(second_left, X_sqr, Y_sqr);
    mpq_sub(second_left, Z_sqr, second_left); /* Z^2 - (X^2 + Y^2) */
    mpq_mul(second_left, second_left, second_left); /* (Z^2 - (X^2 + Y^2))^2 */

    mpq_mul(second_right, X_sqr, Y_sqr); /* X^2Y^2 */
    mpq_mul(second_right, four, second_right); /* 4X^2Y^2 */

    /* succeeds */
    if (mpq_cmp(second_left, second_right) >= 0) {
        mpq_clear(second_left);
        mpq_clear(second_right);
        mpq_clear(four);

        return 1;
    }

    mpq_clear(second_left);
    mpq_clear(second_right);
    mpq_clear(four);

    return 0;
}

/*****************************
 * compute ||J(f)^-1 * v||^2 *
 *****************************/
void compute_norm_sqr_Jv_rational(polynomial_system *F, rational_complex_vector v, rational_complex_vector x, mpq_t *norm_sqr_Jv) {
    int retval = 0, *rowswaps = NULL;

    rational_complex_vector f;
    /* X = J^-1 * v */
    rational_complex_vector X;

    rational_complex_matrix J;
    rational_complex_matrix LU;

    initialize_rational_vector(f, 0);
    initialize_rational_vector(X, 0);

    initialize_rational_matrix(J, 0, 0);
    initialize_rational_matrix(LU, 0, 0);

    eval_polynomial_system_rational(f, J, F, x);

    retval = LUdecomp_rational(LU, &rowswaps, J);

    /* success */
    if (!retval) {
        LUsolve_rational_vector(X, LU, rowswaps, v);
        norm_sqr_rational_vector(*norm_sqr_Jv, X);
    } else {
        mpq_set_ui(*norm_sqr_Jv, 0, 1); /* this may need to be something else */
    }

    clear_rational_vector(f);
    clear_rational_vector(X);
    clear_rational_matrix(J);
    clear_rational_matrix(LU);

    free(rowswaps);
    rowswaps = NULL;
}

/*********************************************
 * return 1 if interval certifies continuous *
 *********************************************/
int is_continuous_rational(rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int retval = 0, left_success = 0, right_success = 0;

    mpq_t alpha_star_sqr;
    mpq_t t_dist_sqr;

    mpq_init(alpha_star_sqr);
    mpq_init(t_dist_sqr);

    mpq_set_ui(alpha_star_sqr, 4, 100);
    mpq_mul(alpha_star_sqr, alpha_star_sqr, alpha_star_sqr);

    mpq_sub(t_dist_sqr, t_right, t_left); /* t_2 -t_1 */
    mpq_mul(t_dist_sqr, t_dist_sqr, t_dist_sqr); /* |t_2 - t_1|^2 */

    /* we want either of magic_number_left or _right to be less than alpha_star */
    mpq_t magic_number_left;
    mpq_t norm_sqr_Jv_left;
    mpq_t magic_number_right;
    mpq_t norm_sqr_Jv_right;

    mpq_init(magic_number_left);
    mpq_init(norm_sqr_Jv_left);
    mpq_init(magic_number_right);
    mpq_init(norm_sqr_Jv_right);

    compute_norm_sqr_Jv_rational(&F_left, v, x_left, &norm_sqr_Jv_left);
    compute_norm_sqr_Jv_rational(&F_right, v, x_right, &norm_sqr_Jv_right);

    /* ||J^-1 * v||^2 */
    mpq_set(magic_number_left, norm_sqr_Jv_left);
    mpq_set(magic_number_right, norm_sqr_Jv_right);

    /* ||J^-1 * v||^2 * g^2 */
    mpq_mul(magic_number_left, magic_number_left, gamma_sqr_left);
    mpq_mul(magic_number_right, magic_number_right, gamma_sqr_right);

    /* ||J^-1 * v||^2 * g^2 * |t_1 - t_2|^2*/
    mpq_mul(magic_number_left, magic_number_left, t_dist_sqr);
    mpq_mul(magic_number_right, magic_number_right, t_dist_sqr);

    /* run the tests */
    left_success = lte_sqr_rational(alpha_sqr_left, magic_number_left, alpha_star_sqr);
    right_success = lte_sqr_rational(alpha_sqr_right, magic_number_right, alpha_star_sqr);

    /* test fails for both left and right */
    if (left_success == 0 && right_success == 0)
        retval = 0;
    else {
        mpq_t gamma_sqr_inv;
        mpq_t points_dist_sqr;
        mpq_t magic_number_sqr;

        mpq_init(gamma_sqr_inv);
        mpq_init(points_dist_sqr);
        mpq_init(magic_number_sqr);

        rational_complex_vector points_dist;
        initialize_rational_vector(points_dist, x_left->size);

        /* test succeeds for both left and right */
        if (left_success && right_success) {
            /* singularities segfault on the Mac */
            mpz_t denom_left;
            mpz_t denom_right;
            mpz_init(denom_left);
            mpz_init(denom_right);

            mpq_get_den(denom_left, gamma_sqr_left);
            mpq_get_den(denom_right, gamma_sqr_right);
            if (mpz_cmp_ui(denom_left, 0) == 0) {
                if (mpz_cmp_ui(denom_right, 0) == 0)
                    mpq_set(gamma_sqr_inv, gamma_sqr_left); /* might as well pick left since both are singularities */
                else
                    mpq_set(gamma_sqr_inv, gamma_sqr_left); /* gamma_sqr_right wins without a fight */
            } else if (mpz_cmp_ui(denom_right, 0) == 0) {
                mpq_set(gamma_sqr_inv, gamma_sqr_left); /* gamma_sqr_left wins without a fight */
            } else {
                mpq_set_min(gamma_sqr_inv, gamma_sqr_left, gamma_sqr_right); /* choose the lowest value for (both finite) gamma */
            }
            mpq_inv(gamma_sqr_inv, gamma_sqr_inv); /* set gamma to 1/gamma */
        } else if (left_success) { /* test succeeds for left only */
            mpq_inv(gamma_sqr_inv, gamma_sqr_left);
        } else { /* test succeeds for right only */
            mpq_inv(gamma_sqr_inv, gamma_sqr_right); 
        }

        /* magic number is now to be 79/(1000*gamma) */
        mpq_set_ui(magic_number_sqr, 79, 1000);
        mpq_mul(magic_number_sqr, magic_number_sqr, magic_number_sqr); /* (79 / 1000)^2 */
        mpq_mul(magic_number_sqr, magic_number_sqr, gamma_sqr_inv); /* (79 / (1000 * gamma))^2 */

        subtract_rational_vector(points_dist, x_left, x_right);
        norm_sqr_rational_vector(points_dist_sqr, points_dist); /* ||x_1 - x_2||^2 */

        if (mpq_cmp(points_dist_sqr, magic_number_sqr) <= 0) {
            retval = 1;
        } else {
            retval = 0;
        }

        mpq_clear(points_dist_sqr);
        mpq_clear(magic_number_sqr);
        mpq_clear(gamma_sqr_inv);

        clear_rational_vector(points_dist);
    }

    mpq_clear(alpha_star_sqr);
    mpq_clear(t_dist_sqr);
    mpq_clear(magic_number_left);
    mpq_clear(magic_number_right);
    mpq_clear(norm_sqr_Jv_left);
    mpq_clear(norm_sqr_Jv_right);

    return retval;
}

/************************************************
 * return 1 if interval certifies discontinuous *
 ************************************************/
int is_not_continuous_rational(rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    int retval = 0, beta_retval = 0, left_success = 0, right_success = 0;

    mpq_t alpha_star_sqr;
    mpq_t t_dist_sqr;

    mpq_init(alpha_star_sqr);
    mpq_init(t_dist_sqr);

    mpq_set_ui(alpha_star_sqr, 4, 100);
    mpq_mul(alpha_star_sqr, alpha_star_sqr, alpha_star_sqr);

    mpq_sub(t_dist_sqr, t_right, t_left); /* t_2 -t_1 */
    mpq_mul(t_dist_sqr, t_dist_sqr, t_dist_sqr); /* |t_2 - t_1|^2 */

    /* we want either of magic_number_left or _right to be less than alpha_star */
    mpq_t magic_number_left;
    mpq_t norm_sqr_Jv_left;
    mpq_t magic_number_right;
    mpq_t norm_sqr_Jv_right;

    mpq_init(magic_number_left);
    mpq_init(norm_sqr_Jv_left);
    mpq_init(magic_number_right);
    mpq_init(norm_sqr_Jv_right);

    compute_norm_sqr_Jv_rational(&F_left, v, x_left, &norm_sqr_Jv_left);
    compute_norm_sqr_Jv_rational(&F_right, v, x_right, &norm_sqr_Jv_right);

    /* ||J^-1 * v||^2 */
    mpq_set(magic_number_left, norm_sqr_Jv_left);
    mpq_set(magic_number_right, norm_sqr_Jv_right);

    /* ||J^-1 * v||^2 * g^2 */
    mpq_mul(magic_number_left, magic_number_left, gamma_sqr_left);
    mpq_mul(magic_number_right, magic_number_right, gamma_sqr_right);

    /* ||J^-1 * v||^2 * g^2 * |t_1 - t_2|^2*/
    mpq_mul(magic_number_left, magic_number_left, t_dist_sqr);
    mpq_mul(magic_number_right, magic_number_right, t_dist_sqr);

    /* run the tests */
    left_success = lte_sqr_rational(alpha_sqr_left, magic_number_left, alpha_star_sqr);
    right_success = lte_sqr_rational(alpha_sqr_right, magic_number_right, alpha_star_sqr);

    /* test fails for both left and right */
    if (left_success == 0 && right_success == 0)
        retval = 0;
    else {
        rational_complex_number beta_sqr_complex1;
        rational_complex_number beta_sqr_complex2;

        mpq_t beta_sqr1;
        mpq_t beta_sqr2;
        mpq_t four;
        mpq_t points_dist_sqr;

        initialize_rational_number(beta_sqr_complex1);
        initialize_rational_number(beta_sqr_complex2);

        mpq_init(beta_sqr1);
        mpq_init(beta_sqr2);
        mpq_init(four);
        mpq_init(points_dist_sqr);

        mpq_set_ui(four, 4, 1);

        /* test succeeds for at least left (possibly right, but use left anyway) */
        if (left_success) {
            beta_retval = compute_beta_sqr_rational(beta_sqr_complex1, &F_left, x_left);

            if (!beta_retval)
                mpq_set(beta_sqr1, beta_sqr_complex1->re);

            beta_retval = compute_beta_sqr_rational(beta_sqr_complex2, &F_left, x_right);

            if (!beta_retval)
                mpq_set(beta_sqr2, beta_sqr_complex2->re);
        } else { /* test succeeds for right only */
            beta_retval = compute_beta_sqr_rational(beta_sqr_complex1, &F_right, x_left);

            if (!beta_retval)
                mpq_set(beta_sqr1, beta_sqr_complex1->re);

            beta_retval = compute_beta_sqr_rational(beta_sqr_complex2, &F_right, x_right);

            if (!beta_retval)
                mpq_set(beta_sqr2, beta_sqr_complex2->re);
        }
        
        /* 4* beta^2 */
        mpq_mul(beta_sqr1, four, beta_sqr1);
        mpq_mul(beta_sqr2, four, beta_sqr2);

        rational_complex_vector points_dist;
        initialize_rational_vector(points_dist, x_left->size);

        subtract_rational_vector(points_dist, x_left, x_right);
        norm_sqr_rational_vector(points_dist_sqr, points_dist); /* ||x_1 - x_2||^2 */

        if (lte_sqr_rational(beta_sqr1, beta_sqr2, points_dist_sqr))
            retval = 1;
        else
            retval = 0;

        clear_rational_number(beta_sqr_complex1);
        clear_rational_number(beta_sqr_complex2);

        clear_rational_vector(points_dist);

        mpq_clear(four);
        mpq_clear(beta_sqr1);
        mpq_clear(beta_sqr2);
        mpq_clear(points_dist_sqr);
    }

    mpq_clear(alpha_star_sqr);
    mpq_clear(t_dist_sqr);
    mpq_clear(magic_number_left);
    mpq_clear(magic_number_right);
    mpq_clear(norm_sqr_Jv_left);
    mpq_clear(norm_sqr_Jv_right);

    return retval;
}

/**********************************************
 * test for the continuity of a given segment *
 **********************************************/
int test_continuity_rational(rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right) {
    if (is_continuous_rational(v, t_left, t_right, x_left, x_right, F_left, F_right, alpha_sqr_left, alpha_sqr_right, gamma_sqr_left, gamma_sqr_right))
        return 1;
    else if (is_not_continuous_rational(v, t_left, t_right, x_left, x_right, F_left, F_right, alpha_sqr_left, alpha_sqr_right, gamma_sqr_left, gamma_sqr_right))
        return -1;

    return 0;
}

/**************************************************************************************
 * subdivide a segment (t_1, w_1) <=> (t_2, w_2) for ((t_1 + t_2)/2, (w_1 + w_2) / 2) *
 **************************************************************************************/
void subdivide_segment_rational(polynomial_system *base, rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t *t_mid, rational_complex_vector *x_mid, int num_var) {
    if (verbosity > BH_VERBOSE)
        printf("subdividing... ");

    int i, x_left_solution, x_mid_solution, x_right_solution, retval = 0, *rowswaps = NULL;

    mpq_t one_half_rational;
    mpq_t alpha_sqr;
    mpq_t beta_sqr;
    mpq_t gamma_sqr;
    mpq_t alpha0_sqr;

    mpq_init(one_half_rational);
    mpq_init(alpha_sqr);
    mpq_init(beta_sqr);
    mpq_init(gamma_sqr);
    mpq_init(alpha0_sqr);

    mpq_set_ui(one_half_rational, 1, 2);
    mpq_set_ui(alpha0_sqr, 12321, 495616); /* 12321/495616 < ((13 - 3 * sqrt(17))/4)^2 */

    rational_complex_number one_half_complex;
    initialize_rational_number(one_half_complex);
    mpq_set_ui(one_half_complex->re, 1, 2);
    mpq_set_ui(one_half_complex->im, 0, 1);

    rational_complex_number beta_sqr_complex;
    initialize_rational_number(beta_sqr_complex);

    mpq_set(beta_sqr_complex->re, beta_sqr);
    mpq_set_ui(beta_sqr_complex->im, 0, 1);

    rational_complex_matrix LU;
    initialize_rational_matrix(LU, 0, 0);

    polynomial_system F_left;
    polynomial_system F_mid;
    polynomial_system F_right;

    /* t_left and t_right are fixed, this will not change */
    apply_tv_rational(base, &F_left, *t_mid, v);
    apply_tv_rational(base, &F_right, *t_mid, v);

    rational_complex_vector new_point;

    initialize_rational_vector(new_point, num_var);

    /* perform a Newton iteration on x_left */
    retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_left, x_left);
    if (retval == ERROR_LU_DECOMP) {
        fprintf(stderr, "singularity suspected near point ");
        print_points_rational(stderr, new_point);
        fprintf(stderr, ERROR_MESSAGE);
        exit(BH_EXIT_OTHER);
    }
    copy_rational_vector(x_left, new_point);

    x_left_solution = compute_abg_sqr_rational(x_left, &F_left, &alpha_sqr, &beta_sqr, &gamma_sqr);

    i = 1;
    /* perform Newton iterations on x_left until within tolerance */
    while (i < newton_tolerance && mpq_cmp(alpha_sqr, alpha0_sqr) >= 0) {
        retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_left, x_left);
        if (retval == ERROR_LU_DECOMP) {
            fprintf(stderr, "singularity suspected near point ");
            print_points_rational(stderr, new_point);
            fprintf(stderr, ERROR_MESSAGE);
            exit(BH_EXIT_OTHER);
        }
        copy_rational_vector(x_left, new_point);

        x_left_solution = compute_abg_sqr_rational(x_left, &F_left, &alpha_sqr, &beta_sqr, &gamma_sqr);

        i++;
    }

    /* perform a Newton iteration on x_right */
    retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_right, x_right);
    if (retval == ERROR_LU_DECOMP) {
        fprintf(stderr, "singularity suspected near point ");
        print_points_rational(stderr, new_point);
        fprintf(stderr, ERROR_MESSAGE);
        exit(BH_EXIT_OTHER);
    }
    copy_rational_vector(x_right, new_point);

    x_right_solution = compute_abg_sqr_rational(x_right, &F_right, &alpha_sqr, &beta_sqr, &gamma_sqr);

    i = 1;
    /* perform Newton iterations on x_right until within tolerance */
    while (i < newton_tolerance && mpq_cmp(alpha_sqr, alpha0_sqr) >= 0) {
        retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_right, x_right);
        if (retval == ERROR_LU_DECOMP) {
            fprintf(stderr, "singularity suspected near point ");
            print_points_rational(stderr, new_point);
            fprintf(stderr, ERROR_MESSAGE);
            exit(BH_EXIT_OTHER);
        }
        copy_rational_vector(x_right, new_point);

        x_right_solution = compute_abg_sqr_rational(x_right, &F_right, &alpha_sqr, &beta_sqr, &gamma_sqr);

        i++;
    }

    /* t_mid */
    mpq_add(*t_mid, t_left, t_right);
    mpq_mul(*t_mid, *t_mid, one_half_rational);

    /* x_mid */
    for (i = 0; i < num_var; i++) {
        add_rational((*x_mid)->coord[i], x_left->coord[i], x_right->coord[i]);
        multiply_rational((*x_mid)->coord[i], (*x_mid)->coord[i], one_half_complex);
    }

    /* apply Newton to x_mid */
    apply_tv_rational(base, &F_mid, *t_mid, v);
    retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_mid, *x_mid);
    if (retval == ERROR_LU_DECOMP) {
        fprintf(stderr, "singularity suspected near point ");
        print_points_rational(stderr, new_point);
        fprintf(stderr, ERROR_MESSAGE);
        exit(BH_EXIT_OTHER);
    }

    copy_rational_vector(*x_mid, new_point);

    x_mid_solution = compute_abg_sqr_rational(*x_mid, &F_right, &alpha_sqr, &beta_sqr, &gamma_sqr);

    i = 1;
    /* perform Newton iterations on x_mid until within tolerance */
    while (i < newton_tolerance && mpq_cmp(alpha_sqr, alpha0_sqr) >= 0) {
        retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_right, *x_mid);
        if (retval == ERROR_LU_DECOMP) {
            fprintf(stderr, "singularity suspected near point ");
            print_points_rational(stderr, new_point);
            fprintf(stderr, ERROR_MESSAGE);
            exit(BH_EXIT_OTHER);
        }
        copy_rational_vector(*x_mid, new_point);

        x_mid_solution = compute_abg_sqr_rational(*x_mid, &F_right, &alpha_sqr, &beta_sqr, &gamma_sqr);

        i++;
    }

    if (verbosity > BH_VERBOSE)
        gmp_printf("new intervals are [%Qd, %Qd] and [%Qd, %Qd]\n", t_left, *t_mid, *t_mid, t_right);

    mpq_clear(one_half_rational);
    mpq_clear(alpha_sqr);
    mpq_clear(beta_sqr);
    mpq_clear(gamma_sqr);
    mpq_clear(alpha0_sqr);
    clear_rational_number(one_half_complex);
    clear_rational_number(beta_sqr_complex);
    clear_rational_matrix(LU);
    clear_polynomial_system(&F_mid);
    clear_rational_vector(new_point);

    free(rowswaps);
    rowswaps = NULL;
}

/*******************************
 * apply f + tv for a single t *
 *******************************/
void apply_tv_rational(polynomial_system *base, polynomial_system *F, mpq_t t, rational_complex_vector v) {
    int i, j, k, num_terms, num_var, num_poly;
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
        mpq_mul(p.coeff[j]->re, t, v->coord[i]->re);
        mpq_mul(p.coeff[j]->im, t, v->coord[i]->im);

        mpq_init(p.norm_sqr);
        norm_sqr_polynomial(p.norm_sqr, &p);

        F->polynomials[i] = p;
    }
}

/***************************************
 * get alpha^2, beta^2, gamma^2 values *
 ***************************************/
int compute_abg_sqr_rational(rational_complex_vector points, polynomial_system *F, mpq_t *alpha_sqr, mpq_t *beta_sqr, mpq_t *gamma_sqr) {
    int res = 0, retval = 0, num_approx_solns = 0, num_var = F->numVariables;

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
    res = compute_alpha_beta_gamma_sqr_rational(P.Nx, P.alpha_sqr, P.beta_sqr, P.gamma_sqr, F, P.x);

    mpq_set(*alpha_sqr, P.alpha_sqr->re);
    mpq_set(*beta_sqr, P.beta_sqr->re);
    mpq_set(*gamma_sqr, P.gamma_sqr->re);

    //copy_rational_vector(points, P.Nx);

    /* check to see if we have successfully computed alpha_sqr, beta_sqr, & gamma_sqr */
    if (res == EXACT_SOLUTION_LU_ERROR) {
        /* exact solution */
        retval = 1;
    } else if (res) {
        /* unknown */
        retval = 0;
    } else {
        /* determine if alpha_sqr is small enough to be an approximate solution */
        num_approx_solns += P.isApproxSoln = determine_approximate_solution_rational(P.alpha_sqr);
        retval = num_approx_solns;
    }

    clear_rational_point_struct(&P);

    return retval;
}

/***********************************************************************************
 * test that each w is in the quadratic convergence basin and check for continuity *
 ***********************************************************************************/
void test_pairwise_rational(polynomial_system *system, rational_complex_vector *v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, int num_var, int iter, mpq_t **t_final, rational_complex_vector **x_final, rational_complex_vector **sing, int *tested, int *succeeded, int *failed, int *num_sing, int check_left) {
    int x_left_solution, x_right_solution, seg_continuous;
    polynomial_system F_left, F_right;
    mpq_t alpha_sqr_left, beta_sqr_left, gamma_sqr_left, alpha_sqr_right, beta_sqr_right, gamma_sqr_right, beta_sqr_min;
    mpq_init(alpha_sqr_left);
    mpq_init(beta_sqr_left);
    mpq_init(gamma_sqr_left);
    mpq_init(alpha_sqr_right);
    mpq_init(beta_sqr_right);
    mpq_init(gamma_sqr_right);
    mpq_init(beta_sqr_min);

    apply_tv_rational(system, &F_left, t_left, *v);
    apply_tv_rational(system, &F_right, t_right, *v);

    x_left_solution = compute_abg_sqr_rational(x_left, &F_left, &alpha_sqr_left, &beta_sqr_left, &gamma_sqr_left);
    x_right_solution = compute_abg_sqr_rational(x_right, &F_right, &alpha_sqr_right, &beta_sqr_right, &gamma_sqr_right);

    /* test if we've been through too many times */
    if (iter > subd_tolerance) {
        *tested += 1;

        if (verbosity > BH_LACONIC)
            gmp_printf("unsure whether segment [%Qd, %Qd] is continuous\n", t_left, t_right);
        fprint_uncertain_rational(t_left, t_right, x_left, x_right, alpha_sqr_left, alpha_sqr_right, beta_sqr_left, beta_sqr_right, gamma_sqr_left, gamma_sqr_right);

        /* alert user to singularities */
        mpz_t denom;
        mpz_init(denom);

        if (check_left) {
            mpq_get_den(denom, gamma_sqr_left);
            if (mpz_cmp_ui(denom, 0) == 0) {
                gmp_fprintf(stderr, "singularity found for t = %Qd at point ", t_left);
                print_points_rational(stderr, x_left);
            }
        }
        mpq_get_den(denom, gamma_sqr_right);
        if (mpz_cmp_ui(denom, 0) == 0) {
            gmp_fprintf(stderr, "singularity found for t = %Qd at point ", t_right);
            print_points_rational(stderr, x_right);
        }

        mpz_clear(denom);

        return;
    }

    mpq_set_min(beta_sqr_min, beta_sqr_left, beta_sqr_right);

    /* perform Newton iterations on non-compliant points */
    int newton_counter = 1;
    while (!x_left_solution || !x_right_solution) {
        if (newton_counter == newton_tolerance) {
            char err_string[] = "Points not in convergence basin; aborting";
            print_error(err_string, stderr);

            exit(BH_EXIT_OTHER);
        }

        int retval, *rowswaps = NULL;
        rational_complex_number beta_sqr_complex;
        rational_complex_vector new_point;
        rational_complex_matrix LU;

        /* perform Newton iterations if x_left is not a solution */
        if (!x_left_solution) {
            if (verbosity > BH_VERBOSE) {
                if (newton_counter % 10 == 1 && newton_counter % 100 != 11)
                    fprintf(stderr, "performing %dst Newton iteration on x_left", newton_counter);
                else if (newton_counter % 10 == 2 && newton_counter % 100 != 12)
                    fprintf(stderr, "performing %dnd Newton iteration on x_left", newton_counter);
                else if (newton_counter % 10 == 3 && newton_counter % 100 != 13)
                    fprintf(stderr, "performing %drd Newton iteration on x_left", newton_counter);
                else
                    fprintf(stderr, "performing %dth Newton iteration on x_left", newton_counter);
                
                gmp_fprintf(stderr, " (interval: [%Qd, %Qd])\n", t_left, t_right);
            }

            initialize_rational_number(beta_sqr_complex);
            initialize_rational_vector(new_point, num_var);
            initialize_rational_matrix(LU, 0, 0);

            retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_left, x_left);
            if (retval == ERROR_LU_DECOMP) {
                fprintf(stderr, "singularity suspected near point ");
                print_points_rational(stderr, new_point);
                fprintf(stderr, ERROR_MESSAGE);
                exit(BH_EXIT_OTHER);
            }

            copy_rational_vector(x_left, new_point);

            x_left_solution = compute_abg_sqr_rational(x_left, &F_left, &alpha_sqr_left, &beta_sqr_left, &gamma_sqr_left);

            clear_rational_number(beta_sqr_complex);
            clear_rational_vector(new_point);
            clear_rational_matrix(LU);

            free(rowswaps);
            rowswaps = NULL;
        }

        /* perform Newton iterations if x_right is not a solution */
        if (!x_right_solution) {
            if (verbosity > BH_VERBOSE) {
                if (newton_counter % 10 == 1 && newton_counter % 100 != 11)
                    fprintf(stderr, "performing %dst Newton iteration on x_right", newton_counter);
                else if (newton_counter % 10 == 2 && newton_counter % 100 != 12)
                    fprintf(stderr, "performing %dnd Newton iteration on x_right", newton_counter);
                else if (newton_counter % 10 == 3 && newton_counter % 100 != 13)
                    fprintf(stderr, "performing %drd Newton iteration on x_right", newton_counter);
                else
                    fprintf(stderr, "performing %dth Newton iteration on x_right", newton_counter);

                gmp_fprintf(stderr, " (interval: [%Qd, %Qd])\n", t_left, t_right);
            }

            initialize_rational_number(beta_sqr_complex);
            initialize_rational_vector(new_point, num_var);
            initialize_rational_matrix(LU, 0, 0);

            retval = newton_iteration_rational(new_point, beta_sqr_complex, LU, &rowswaps, &F_right, x_right);
            if (retval == ERROR_LU_DECOMP) {
                fprintf(stderr, "singularity suspected near point ");
                print_points_rational(stderr, new_point);
                fprintf(stderr, ERROR_MESSAGE);
                exit(BH_EXIT_OTHER);
            }

            copy_rational_vector(x_right, new_point);

            x_right_solution = compute_abg_sqr_rational(x_right, &F_right, &alpha_sqr_right, &beta_sqr_right, &gamma_sqr_right);

            clear_rational_number(beta_sqr_complex);
            clear_rational_vector(new_point);
            clear_rational_matrix(LU);

            free(rowswaps);
            rowswaps = NULL;
        }

        newton_counter++;
    }

    if (verbosity > BH_CHATTY)
        gmp_printf("testing segment [%Qd, %Qd]\n", t_left, t_right);

    seg_continuous = test_continuity_rational(*v, t_left, t_right, x_left, x_right, F_left, F_right, alpha_sqr_left, alpha_sqr_right, gamma_sqr_left, gamma_sqr_right);

    if (seg_continuous == 0) {
        mpq_t t_mid;
        mpq_init(t_mid);
        rational_complex_vector x_mid;
        initialize_rational_vector(x_mid, num_var);

        subdivide_segment_rational(system, *v, t_left, t_right, x_left, x_right, &t_mid, &x_mid, num_var);

        /* recurse! */
        test_pairwise_rational(system, v, t_left, t_mid, x_left, x_mid, num_var, iter + 1, t_final, x_final, sing, tested, succeeded, failed, num_sing, 0);
        test_pairwise_rational(system, v, t_mid, t_right, x_mid, x_right, num_var, iter + 1, t_final, x_final, sing, tested, succeeded, failed, num_sing, 0);
        mpq_clear(t_mid);
        clear_rational_vector(x_mid);
    } else if (seg_continuous == 1) {
        if (verbosity > BH_LACONIC)
            gmp_printf("segment [%Qd, %Qd] is continuous\n", t_left, t_right);
        fprint_continuous_rational(t_left, t_right, x_left, x_right, alpha_sqr_left, alpha_sqr_right, beta_sqr_left, beta_sqr_right, gamma_sqr_left, gamma_sqr_right);
        *tested += 1;
        *succeeded += 1;

        /* alert user to singularities */
        mpz_t denom;
        mpz_init(denom);

        if (check_left) {
            mpq_get_den(denom, gamma_sqr_left);
            if (mpz_cmp_ui(denom, 0) == 0) {
                gmp_fprintf(stderr, "singularity found for t = %Qd at point ", t_left);
                print_points_rational(stderr, x_left);
            }
        }
        mpq_get_den(denom, gamma_sqr_right);
        if (mpz_cmp_ui(denom, 0) == 0) {
            gmp_fprintf(stderr, "singularity found for t = %Qd at point ", t_right);
            print_points_rational(stderr, x_right);
        }

        mpz_clear(denom);

        /* copy t and w values into t_final and x_final respectively */
        errno = 0;

        mpq_t *t_final_tmp = realloc(*t_final, (size_t) (*succeeded + 1) * sizeof(mpq_t));
        if (t_final_tmp == NULL) {
            snprintf(error_string, (size_t) termwidth, "Couldn't realloc: %s\n", strerror(errno));
            print_error(error_string, stderr);
            exit(BH_EXIT_MEMORY);
        } else {
            *t_final = t_final_tmp;
            mpq_init((*t_final)[*succeeded]);
            mpq_set((*t_final)[*succeeded], t_right);
        }

        errno = 0;

        rational_complex_vector *x_final_tmp = realloc(*x_final, (size_t) (*succeeded + 1) * sizeof(rational_complex_vector));
        if (x_final_tmp == NULL) {
            snprintf(error_string, (size_t) termwidth, "Couldn't realloc: %s\n", strerror(errno));
            print_error(error_string, stderr);
            exit(BH_EXIT_MEMORY);
        } else {
            *x_final = x_final_tmp;
            initialize_rational_vector((*x_final)[*succeeded], x_right->size);
            copy_rational_vector((*x_final)[*succeeded], x_right);
        }
    } else {
        if (verbosity > BH_LACONIC)
            gmp_printf("segment [%Qd, %Qd] is not continuous\n", t_left, t_right);
        fprint_discontinuous_rational(t_left, t_right, x_left, x_right, alpha_sqr_left, alpha_sqr_right, beta_sqr_left, beta_sqr_right, gamma_sqr_left, gamma_sqr_right);

        /* alert user to singularities */
        mpz_t denom;
        mpz_init(denom);

        if (check_left) {
            mpq_get_den(denom, gamma_sqr_left);
            if (mpz_cmp_ui(denom, 0) == 0) {
                gmp_fprintf(stderr, "singularity found for t = %Qd at point ", t_left);
                print_points_rational(stderr, x_left);
            }
        }
        mpq_get_den(denom, gamma_sqr_right);
        if (mpz_cmp_ui(denom, 0) == 0) {
            gmp_fprintf(stderr, "singularity found for t = %Qd at point ", t_right);
            print_points_rational(stderr, x_right);
        }

        mpz_clear(denom);

        *tested += 1;
        *failed += 1;
    }

    /* free up stuff */
    clear_polynomial_system(&F_left);
    clear_polynomial_system(&F_right);

    mpq_clear(alpha_sqr_left);
    mpq_clear(beta_sqr_left);
    mpq_clear(gamma_sqr_left);
    mpq_clear(alpha_sqr_right);
    mpq_clear(beta_sqr_right);
    mpq_clear(gamma_sqr_right);
    mpq_clear(beta_sqr_min);
}

/*******************
 * certify H(x, t) *
 *******************/
void test_system_rational(polynomial_system *system, void *v, void *t, void *x, int num_points, void **t_final, void **x_final, void **sing, int *tested, int *succeeded, int *failed, int *num_sing) {
    int i, num_var = system->numVariables;

    rational_complex_vector *v_rational = (rational_complex_vector *) v;
    mpq_t *t_rational = (mpq_t *) t;
    rational_complex_vector *x_rational = (rational_complex_vector *) x;

    mpq_t **t_final_rational = (mpq_t **) t_final;
    rational_complex_vector **x_final_rational = (rational_complex_vector **) x_final;
    rational_complex_vector **sing_rational = (rational_complex_vector **) sing;

    /* set up t_final and x_final */
    *t_final_rational = malloc(sizeof(mpq_t));
    if (*t_final_rational == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
        print_error(error_string, stderr);

        exit(BH_EXIT_MEMORY);
    }

    mpq_init((*t_final_rational)[0]);
    mpq_set((*t_final_rational)[0], t_rational[0]);

        *x_final_rational = malloc(sizeof(rational_complex_vector));
    if (*x_final_rational == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't alloc: %s\n", strerror(errno));
        print_error(error_string, stderr);

        exit(BH_EXIT_MEMORY);
    }

    initialize_rational_vector((*x_final_rational)[0], x_rational[0]->size);
    copy_rational_vector((*x_final_rational)[0], x_rational[0]);

    /* for each t_i, t_{i+1} */
    /* this is the part that needs to get parallelized */
    for (i = 0; i < num_points - 1; i++) {
        test_pairwise_rational(system, v_rational, t_rational[i], t_rational[i+1], x_rational[i], x_rational[i+1], num_var, 1, t_final_rational, x_final_rational, sing_rational, tested, succeeded, failed, num_sing, (i == 0));
    }
}
