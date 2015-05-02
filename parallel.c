#include "cadenza.h"

int send_complex_number(complex_number c, int to) {
    int sizes[2], SENDSIZE=0, SENDRE=1, SENDIM=2;
    long e1, e2;
    char *re = NULL;
    char *im = NULL;
    char rebuf[BH_MAX_STRING], imbuf[BH_MAX_STRING];
    
    re = mpf_get_str(NULL, &e1, 10, 0, c->re);
    im = mpf_get_str(NULL, &e2, 10, 0, c->im);
    
    if (re[0] == '-') {
        sprintf(rebuf, "-0.%se%ld", &re[1], e1);
    } else {
        sprintf(rebuf, "0.%se%ld", re, e1);
    }
    
    if (im[0] == '-') {
        sprintf(imbuf, "-0.%se%ld", &im[1], e2);
    } else {
        sprintf(imbuf, "0.%se%ld", im, e2);
    }
    
//     printf("sending restring is %s\n", rebuf);
//     printf("sending imstring is %s\n", imbuf);
    
    sizes[0] = strlen(rebuf) + 1;
    sizes[1] = strlen(imbuf) + 1;

    MPI_Send(sizes, 2, MPI_INT, to, SENDSIZE, MPI_COMM_WORLD);
    MPI_Send(rebuf, sizes[0], MPI_CHAR, to, SENDRE, MPI_COMM_WORLD);
    MPI_Send(imbuf, sizes[1], MPI_CHAR, to, SENDIM, MPI_COMM_WORLD);

    free(re);
    free(im);
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

    mpf_set_str(c->re, rebuf, 10);
    mpf_set_str(c->im, imbuf, 10);

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
    char *rebuf = NULL;
    char *imbuf = NULL;

    rebuf = mpq_get_str(NULL, 10, c->re);
    imbuf = mpq_get_str(NULL, 10, c->im);
    
    sizes[0] = strlen(rebuf) + 1;
    sizes[1] = strlen(imbuf) + 1;

    MPI_Send(sizes, 2, MPI_INT, to, SENDSIZE, MPI_COMM_WORLD);
    MPI_Send(rebuf, sizes[0], MPI_CHAR, to, SENDRE, MPI_COMM_WORLD);
    MPI_Send(imbuf, sizes[1], MPI_CHAR, to, SENDIM, MPI_COMM_WORLD);

    free(rebuf); free(imbuf);
    return 0;
}

int recv_rational_complex_number(rational_complex_number c, int from) {
    int sizes[2], SENDSIZE=0, SENDRE=1, SENDIM=2;
    char rebuf[BH_MAX_STRING];
    char imbuf[BH_MAX_STRING];
    MPI_Status status;

    MPI_Recv(sizes, 2, MPI_INT, from, SENDSIZE, MPI_COMM_WORLD, &status);
    MPI_Recv(rebuf, sizes[0], MPI_CHAR, from, SENDRE, MPI_COMM_WORLD, &status);
    MPI_Recv(imbuf, sizes[1], MPI_CHAR, from, SENDIM, MPI_COMM_WORLD, &status);
    
    mpq_set_str(c->re, rebuf, 10);
    mpq_set_str(c->im, imbuf, 10);
    
    return 0;
}

int send_polynomial(polynomial *p, int to) {
    int i, j, c=0, vals[5], SENDVALS=0, SENDNORM2=1, SENDINDEX=2, SENDEXP=4, tv=p->numTerms*p->numVariables;
    int *exponents = malloc(tv*sizeof(int));
    char *normbuf = NULL;

    vals[0] = p->numVariables;
    vals[1] = p->numTerms;
    vals[2] = p->degree;
    vals[3] = p->isReal;
    
    normbuf = mpq_get_str(NULL, 10, p->norm_sqr);
    vals[4] = strlen(normbuf);

    for (i=0; i<p->numTerms; i++) {
        for (j=0; j<p->numVariables; j++) {
            exponents[c++] = p->exponents[i][j];
        }
    }

    MPI_Send(vals, 5, MPI_INT, to, SENDVALS, MPI_COMM_WORLD);
    MPI_Send(normbuf, vals[4], MPI_CHAR, to, SENDNORM2, MPI_COMM_WORLD);
    MPI_Send(exponents, tv, MPI_INT, to, SENDEXP, MPI_COMM_WORLD);
    for (i=0; i<p->numTerms; i++) {
        MPI_Send(&i, 1, MPI_INT, to, SENDINDEX, MPI_COMM_WORLD);
        send_rational_complex_number(p->coeff[i], to);
    }

    free(exponents);
    free(normbuf);
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
    mpq_set_str(p->norm_sqr, buf, 10);

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
    char *normbuf;

    vals[0] = F->numVariables;
    vals[1] = F->numPolynomials;
    vals[2] = F->maximumDegree;
    vals[3] = F->isReal;
    vals[4] = F->numExponentials;
    
    normbuf = mpq_get_str(NULL, 10, F->norm_sqr);
    vals[5] = strlen(normbuf);

    if (vals[4] > 0) {
        fprintf(stderr, "we don't deal in exponentials in this here program\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Send(vals, 6, MPI_INT, to, SENDVALS, MPI_COMM_WORLD);
    MPI_Send(normbuf, vals[5], MPI_CHAR, to, SENDNORM2, MPI_COMM_WORLD);

    for (i=0; i<F->numPolynomials; i++) {
        MPI_Send(&i, 1, MPI_INT, to, SENDINDEX, MPI_COMM_WORLD);
        send_polynomial(&(F->polynomials[i]), to);
    }
    
    free(normbuf);
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
