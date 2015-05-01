#include "cadenza.h"

int send_complex_number(complex_number c, int to) {
    int sizes[2], SENDSIZE=0, SENDRE=1, SENDIM=2;
    char rebuf[BH_MAX_STRING];
    char imbuf[BH_MAX_STRING];
    
    sizes[0] = mpfr_sprintf(rebuf, "%Re", c->re);
    sizes[1] = mpfr_sprintf(imbuf, "%Re", c->im);

    MPI_Send(sizes, 2, MPI_INT, to, SENDSIZE, MPI_COMM_WORLD);
    MPI_Send(rebuf, sizes[0], MPI_CHAR, to, SENDRE, MPI_COMM_WORLD);
    MPI_Send(imbuf, sizes[1], MPI_CHAR, to, SENDIM, MPI_COMM_WORLD);

    return 0;
}

int recv_complex_number(complex_number c, int from) {
    int sizes[2], SENDSIZE=0, SENDRE=1, SENDIM=2;
    char rebuf[BH_MAX_STRING];
    char imbuf[BH_MAX_STRING];
    MPI_Status status;

    MPI_Recv(sizes, 2, MPI_INT, from, SENDSIZE, MPI_COMM_WORLD, &status);
    MPI_Recv(rebuf, sizes[0], MPI_CHAR, from, SENDRE, MPI_COMM_WORLD, &status);
    MPI_Recv(imbuf, sizes[1], MPI_CHAR, from, SENDIM, MPI_COMM_WORLD, &status);

    mpf_set_str(c->re, rebuf, 0);
    mpf_set_str(c->im, imbuf, 0);

    return 0;
}

int send_complex_vector(complex_vector v, int to) {
    int i, vals[3], SENDVALS=0, SENDCOORD=1;

    vals[0] = v->alloc_size;
    vals[1] = v->curr_prec;
    vals[2] = v->size;
    MPI_Send(vals, 3, MPI_INT, to, SENDVALS, MPI_COMM_WORLD);

    for (i=0; i<vals[0]; i++) {
        MPI_Send(&i, 1, MPI_INT, to, SENDCOORD, MPI_COMM_WORLD);
        send_complex_number(v->coord[i], to);
    }
}

int recv_complex_vector(complex_vector v, int from) {
    int i, index, vals[3], SENDVALS=0, SENDCOORD=1;
    MPI_Status status;

    complex_number tmp;
    initialize_number(tmp);

    MPI_Recv(vals, 3, MPI_INT, from, SENDVALS, MPI_COMM_WORLD, &status);
    initialize_vector(v, vals[0]);
    for (i=0; i<vals[0]; i++) {
        MPI_Recv(&index, 1, MPI_INT, from, SENDCOORD, MPI_COMM_WORLD, &status);
        recv_complex_number(tmp, from);
        set_number(v->coord[index], tmp);
    }

    clear_number(tmp);
}

int send_rational_complex_number(rational_complex_number c, int to) {
    int sizes[2], SENDSIZE=0, SENDRE=1, SENDIM=2;
    char rebuf[BH_MAX_STRING];
    char imbuf[BH_MAX_STRING];

    sizes[0] = mpfr_sprintf(rebuf, "%Qd", c->re);
    sizes[1] = mpfr_sprintf(imbuf, "%Qd", c->im);

    MPI_Send(sizes, 2, MPI_INT, to, SENDSIZE, MPI_COMM_WORLD);
    MPI_Send(rebuf, sizes[0], MPI_CHAR, to, SENDRE, MPI_COMM_WORLD);
    MPI_Send(imbuf, sizes[1], MPI_CHAR, to, SENDIM, MPI_COMM_WORLD);

    return 0;
}

int recv_rational_complex_number(rational_complex_number c, int from) {
    int sizes[2], SENDSIZE=0, SENDRE=1, SENDIM=2;
    char rebuf[BH_MAX_STRING];
    char imbuf[BH_MAX_STRING];
    MPI_Status status;

    mpq_t re, im;
    mpq_init(re);
    mpq_init(im);

    MPI_Recv(sizes, 2, MPI_INT, from, SENDSIZE, MPI_COMM_WORLD, &status);
    MPI_Recv(rebuf, sizes[0], MPI_CHAR, from, SENDRE, MPI_COMM_WORLD, &status);
    MPI_Recv(imbuf, sizes[1], MPI_CHAR, from, SENDIM, MPI_COMM_WORLD, &status);

    // I don't know why this is necessary but it is
    rebuf[sizes[0]] = '\0';
    imbuf[sizes[1]] = '\0';

    mpq_set_str(re, rebuf, 0);
    mpq_set_str(im, imbuf, 0);

    mpq_set(c->re, re);
    mpq_set(c->im, im);

    return 0;
}

int send_polynomial(polynomial *p, int to) {
    int i, j, c=0, vals[5], SENDVALS=0, SENDNORM2=1, SENDINDEX=2, SENDEXP=4, tv=p->numTerms*p->numVariables;
    int *exponents = malloc(tv*sizeof(int));
    char buf[BH_MAX_STRING];

    vals[0] = p->numVariables;
    vals[1] = p->numTerms;
    vals[2] = p->degree;
    vals[3] = p->isReal;
    vals[4] = mpfr_sprintf(buf, "%Qd", p->norm_sqr);

    for (i=0; i<p->numTerms; i++) {
        for (j=0; j<p->numVariables; j++) {
            exponents[c++] = p->exponents[i][j];
        }
    }

    MPI_Send(vals, 5, MPI_INT, to, SENDVALS, MPI_COMM_WORLD);
    MPI_Send(buf, vals[4], MPI_CHAR, to, SENDNORM2, MPI_COMM_WORLD);
    MPI_Send(exponents, tv, MPI_INT, to, SENDEXP, MPI_COMM_WORLD);
    for (i=0; i<p->numTerms; i++) {
        MPI_Send(&i, 1, MPI_INT, to, SENDINDEX, MPI_COMM_WORLD);
        send_rational_complex_number(p->coeff[i], to);
    }

    free(exponents);
}

int recv_polynomial(polynomial *p, int from) {
    int i, j, c=0, index, vals[5], SENDVALS=0, SENDNORM2=1, SENDINDEX=2, SENDEXP=4, tv;
    int *exponents = NULL;
    char buf[BH_MAX_STRING];
    MPI_Status status;

    mpq_init(p->norm_sqr);

    MPI_Recv(vals, 5, MPI_INT, from, SENDVALS, MPI_COMM_WORLD, &status);

    p->numVariables = vals[0];
    p->numTerms = vals[1];
    p->degree = vals[2];
    p->isReal = vals[3];

    MPI_Recv(buf, vals[4], MPI_CHAR, from, SENDNORM2, MPI_COMM_WORLD, &status);
    mpq_set_str(p->norm_sqr, buf, 0);

    tv = vals[0]*vals[1];
    exponents = malloc(tv*sizeof(int));
    MPI_Recv(exponents, tv, MPI_INT, from, SENDEXP, MPI_COMM_WORLD, &status);

    p->exponents = (int**)malloc(p->numTerms*sizeof(int*));
    for (i=0; i<p->numTerms; i++) {
        p->exponents[i] = (int *)malloc(p->numVariables*sizeof(int));
        for (j=0; j<p->numVariables; j++) {
            p->exponents[i][j] = exponents[c++];
        }
    }

    p->coeff = malloc(p->numTerms*sizeof(rational_complex_number));
    for (i=0; i<p->numTerms; i++) {
        MPI_Recv(&index, 1, MPI_INT, from, SENDINDEX, MPI_COMM_WORLD, &status);
        initialize_rational_number(p->coeff[index]);
        recv_rational_complex_number(p->coeff[index], from);
    }

    free(exponents);
}

int send_polynomial_system(polynomial_system *F, int to) {
    int i, vals[6], SENDVALS=0, SENDNORM2=1, SENDINDEX=2;
    char buf[BH_MAX_STRING];

    vals[0] = F->numVariables;
    vals[1] = F->numPolynomials;
    vals[2] = F->maximumDegree;
    vals[3] = F->isReal;
    vals[4] = F->numExponentials;
    vals[5] = mpfr_sprintf(buf, "%Qd", F->norm_sqr);

    if (vals[4] > 0) {
        fprintf(stderr, "we don't deal in exponentials in this here program\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Send(vals, 6, MPI_INT, to, SENDVALS, MPI_COMM_WORLD);
    MPI_Send(buf, vals[5], MPI_CHAR, to, SENDNORM2, MPI_COMM_WORLD);

    for (i=0; i<F->numPolynomials; i++) {
        MPI_Send(&i, 1, MPI_INT, to, SENDINDEX, MPI_COMM_WORLD);
        send_polynomial(&(F->polynomials[i]), to);
    }
}

int recv_polynomial_system(polynomial_system *F, int from) {
    int i, index, vals[6], SENDVALS=0, SENDNORM2=1, SENDINDEX=2;
    char buf[BH_MAX_STRING];
    MPI_Status status;

    mpq_init(F->norm_sqr);

    MPI_Recv(vals, 6, MPI_INT, from, SENDVALS, MPI_COMM_WORLD, &status);

    F->numVariables = vals[0];
    F->numPolynomials = vals[1];
    F->maximumDegree = vals[2];
    F->isReal = vals[3];
    F->numExponentials = vals[4];

    MPI_Recv(buf, vals[5], MPI_CHAR, from, SENDNORM2, MPI_COMM_WORLD, &status);
    mpq_set_str(F->norm_sqr, buf, 0);

    F->polynomials = malloc(F->numPolynomials*sizeof(polynomial));
    for (i=0; i<F->numPolynomials; i++) {
        MPI_Recv(&index, 1, MPI_INT, from, SENDINDEX, MPI_COMM_WORLD, &status);
        recv_polynomial(&(F->polynomials[i]), from);
    }

    F->exponentials = NULL;
}
