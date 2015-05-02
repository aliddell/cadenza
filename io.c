/*
 * Cadenza
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io.c: Generic input/output functions for Cadenza
 */
#include "cadenza.h"

/********************************************************
 * print a helpful message about command-line arguments *
 ********************************************************/
void usage() {
    fputs("usage: cadenza ARGS\n", stderr);
    fputs("\nif you do not specify a config file with `-c FILENAME', the mandatory arguments are:\n", stderr);
    fputs("\t-s || --system FILENAME (designates FILENAME as system file)\n", stderr);
    fputs("\t-p || --points FILENAME (designates FILENAME as points file)\n", stderr);
    fputs("\noptional arguments are:\n", stderr);
    fputs("\t-h || --help (displays this help dialog)\n", stderr);
    fputs("\t-f || --float (use floating point arithmetic (default is rational))\n", stderr);
    fputs("\t-a || --ascending (certify intervals from t=0 to t=1 (default is descending))\n", stderr);
    fputs("\t-c || --config FILENAME (use options in FILENAME (command-line arguments override))\n", stderr);
    fputs("\t-m || --precision INT (sets floating-point precision to INT (default is 32))\n", stderr);
    fputs("\t-n || --newtons INT (sets maximum number of Newton iterations to INT (default is 20))\n", stderr);
    fputs("\t-d || --subdivisions INT (sets maximum number of segment subdivisions to INT (default is 100))\n", stderr);
}

/*****************************************
 * uniformly print error messages to outfile *
 *****************************************/
void print_error(char *msg, FILE *outfile) {
    fprintf(outfile, "\nERROR: %s\n\n", msg);
    fflush(outfile);
}

/*****************************************************
 * print a helpful message about this program to outfile *
 *****************************************************/
void prog_info(FILE *outfile) {
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    fprintf(outfile, "%s v%s (built %s) (compiled %s)\n", BH_PROGRAM_NAME, BH_VERSION, BH_BUILD_DATE, __DATE__);
    fprintf(outfile, "%s\n", BH_AUTHORS);
    fprintf(outfile, "GMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
}

/******************************************
 * print the details of the configuration *
 ******************************************/
void display_config(FILE *outfile) {
    char precision_str[termwidth];

//     if (arithmetic_type == BH_USE_FLOAT)
//         //snprintf(precision_str, (size_t) termwidth, "%d-bit precision floating point", default_precision);
//         sprintf(precision_str, "%d-bit precision floating point", default_precision);
//     else
//         strcpy(precision_str, "exact rational");
// 
//     fprintf(outfile, "Computing using %s arithmetic\n", precision_str);
    fprintf(outfile, "Polynomial system file is %s\n", sysfile);
    fprintf(outfile, "Point set file is %s\n", pointsfile);
    if (strcmp(configfile, "") != 0)
        fprintf(outfile, "Configuration file is %s\n", configfile);
    fprintf(outfile, "Maximum Newton iterations: %d\n", newton_tolerance);
    fprintf(outfile, "Maximum segment subdivisions: %d\n\n", subd_tolerance);

    //free(precision_str);
}

/*************************************************************
 * read the configuration file, if given, and apply settings *
 *************************************************************/
void read_config_file() {
    int verb = BH_LACONIC, arith_type = BH_USE_RATIONAL, def_prec = BH_PREC_MIN,  newtol = BH_NEWT_TOLERANCE, subtol = BH_SUB_TOLERANCE;
    char sysf[BH_MAX_FILENAME] = "";
    char pointsf[BH_MAX_FILENAME] = "";
    char line[BH_MAX_FILENAME];

    FILE *config = fopen(configfile, "r");
    if (config == NULL) {
        snprintf(error_string, (size_t) termwidth, "Couldn't open config file `%s': %s", configfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    /* read information from the file and copy it to local values */
    while (fgets(line, BH_TERMWIDTH, config) != NULL) {
        char *key = strtok(line, " \t=\n");
        char *val = strtok(NULL, " \t=\n");

        if (strcmp(key, "verbosity") == 0) {
            errno = 0;

            if (strcmp(val, "chatty") == 0)
                verb = BH_CHATTY;
            else if (strcmp(val, "verbose") == 0)
                verb = BH_VERBOSE;
            else if (strcmp(val, "loquacious") == 0)
                verb = BH_LOQUACIOUS;
            else {
                verb = (int) strtol(val, NULL, 10); /* also accepts numerical values */
                if (verb > BH_LOQUACIOUS)
                    verb = BH_LOQUACIOUS;
                if (errno == ERANGE)
                    verb = BH_LACONIC;
            }
        } else if (strcmp(key, "arithmetic") == 0) {
            if (strcmp(val, "float") == 0)
                arith_type = BH_USE_FLOAT;
            else if (strcmp(val, "rational") == 0)
                arith_type = BH_USE_RATIONAL;
        } else if (strcmp(key, "precision") == 0) {
            errno = 0;
            
            def_prec = (int) strtol(val, NULL, 10);

            if (errno == ERANGE) {
                def_prec = BH_PREC_MIN;
            } else if (def_prec < MPFR_PREC_MIN) {
                def_prec = MPFR_PREC_MIN;
            }
        } else if (strcmp(key, "newtons") == 0) {
            errno = 0;
            
            newtol = (int) strtol(val, NULL, 10);

            if (errno == ERANGE || newtol < 1)
                newtol = BH_NEWT_TOLERANCE;
        } else if (strcmp(key, "subdivisions") == 0) {
            errno = 0;

            subtol = (int) strtol(val, NULL, 10);

            if (errno == ERANGE || subtol < 1)
                subtol = BH_SUB_TOLERANCE;
        } else if (strcmp(key, "sysfile") == 0) {
            strcpy(sysf, val);
        } else if (strcmp(key, "pointsfile") == 0) {
            strcpy(pointsf, val);
        } else if (strcmp(key, "sort") == 0) {
            if (strcmp(val, "ascending") == 0) {
                sort_order = BH_ASCENDING;
            }
        }
    }

    fclose(config);

    /* command-line flags override config file options */
    if (verbosity == BH_UNSET)
        verbosity = verb;
    if (arithmetic_type == BH_UNSET)
        arithmetic_type = arith_type;
    if (default_precision == BH_UNSET)
        default_precision = def_prec;
    if (newton_tolerance == BH_UNSET)
        newton_tolerance = newtol;
    if (subd_tolerance == BH_UNSET)
        subd_tolerance = subtol;
    if (strcmp(sysfile, "") == 0)
        strcpy(sysfile, sysf);
    if (strcmp(pointsfile, "") == 0)
        strcpy(pointsfile, pointsf);
}

/***********************************
* parse polynomial from input file *
************************************/
polynomial parse_polynomial(FILE *sysfh, int num_var) {
    int num_terms, res, i, j, max_degree;
    polynomial p;

    p.numVariables = num_var;

    /* first get the number of terms */
    errno = 0;
    res = fscanf(sysfh, "%d", &num_terms);

    if (res == EOF) {
        snprintf(error_string, (size_t) termwidth, "Error reading %s: unexpected EOF", sysfile);

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (res == 0) {
        snprintf(error_string, (size_t) termwidth, "Error reading %s: %s", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        snprintf(error_string, (size_t) termwidth, "Error reading %s: %s", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADPARSE);
    }

    p.exponents = malloc(num_terms * sizeof(int*));

    for (i = 0; i < num_terms; i++)
        p.exponents[i] = malloc(num_var * sizeof(int));
    p.coeff = malloc(num_terms * sizeof(rational_complex_number));
    p.numTerms = num_terms;
    p.degree = 0;
    p.isReal = 0;

    /* get the exponents and coefficients for each term */
    for (i = 0; i < num_terms; i++) {
        max_degree = 0;
        /* get exponent for each variable in the term */
        for (j = 0; j < num_var; j++) {
            errno = 0;
            res = fscanf(sysfh, "%d", &p.exponents[i][j]);

            if (res == EOF) {
                snprintf(error_string, (size_t) termwidth, "Error reading %s: unexpected EOF", sysfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (res == 0) {
                snprintf(error_string, (size_t) termwidth, "Error reading %s: %s", sysfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                snprintf(error_string, (size_t) termwidth, "Error reading %s: %s", sysfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADPARSE);
            }

            max_degree += p.exponents[i][j];
        }

        /* check if degree is greater than p.degree */
        if (max_degree > p.degree)
            p.degree = max_degree;

        initialize_rational_number(p.coeff[i]);

        errno = 0;
        res = gmp_fscanf(sysfh, "%Qd %Qd", p.coeff[i]->re, p.coeff[i]->im);

        if (res == EOF) {
            snprintf(error_string, (size_t) termwidth, "Error reading %s: unexpected EOF", sysfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (res == 0) {
            snprintf(error_string, (size_t) termwidth, "Error reading %s: %s", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            snprintf(error_string, (size_t) termwidth, "Error reading %s: %s", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }
    }

    mpq_init(p.norm_sqr);
    norm_sqr_polynomial(p.norm_sqr, &p);

    return p;
}

/**********************************************************
 * sum the exponents, return 1 if sum is 0 or 0 otherwise *
 **********************************************************/
int is_constant(int *exponents, int num_var) {
    int i, sum = 0;
    for (i = 0; i < num_var; i++) {
        sum += exponents[i];
    }

    return (sum == 0);
}

void fprint_monomial(FILE *outfile, int *exponents, int num_var) {
    int i, first_term = 1;

    for (i = 0; i < num_var; i++) {
        if (exponents[i] > 1) {
            if (first_term) {
                fprintf(outfile, "x%d^%d", i, exponents[i]);
                first_term = 0;
            } else
                fprintf(outfile, " * x%d^%d", i, exponents[i]);
        } else if (exponents[i] > 0) {
            if (first_term) {
                fprintf(outfile, "x%d", i);
                first_term = 0;
            } else
                fprintf(outfile, " * x%d", i);
        }
    }
}

/****************************************
 * print a polynomial system to outfile *
 ****************************************/
void print_system(FILE *outfile, polynomial_system *system) {
    int i, j;
    mpz_t r_num, r_denom, i_num, i_denom;
    mpz_init(r_num);
    mpz_init(r_denom);
    mpz_init(i_num);
    mpz_init(i_denom);

    for (i = 0; i < system->numPolynomials; i++) {
        polynomial p = system->polynomials[i];
        for (j = 0; j < p.numTerms; j++) {
            mpq_get_num(r_num, p.coeff[j]->re);
            mpq_get_den(r_denom, p.coeff[j]->re);
            mpq_get_num(i_num, p.coeff[j]->im);
            mpq_get_den(i_denom, p.coeff[j]->im);

            /* term is zero; this shouldn't normally happen */
            if (mpq_cmp_ui(p.coeff[j]->re, 0, 1) == 0 && mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0) {
                if (j == p.numTerms - 1)
                    gmp_fprintf(outfile, "0");
                else
                    continue;
            }
            /* constant term */
            else if (is_constant(p.exponents[j], p.numVariables)) {
                if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0)
                    gmp_fprintf(outfile, "%Qd", p.coeff[j]->re);
                else if (mpq_cmp_ui(p.coeff[j]->re, 0, 1) == 0)
                    gmp_fprintf(outfile, "I * %Qd", p.coeff[j]->im);
                else 
                    gmp_fprintf(outfile, "%Qd + I * %Qd", p.coeff[j]->re, p.coeff[j]->im);
            }
            /* coefficient is one */
            else if (mpq_cmp_ui(p.coeff[j]->re, 1, 1) == 0 && mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0) {
                fprint_monomial(outfile, p.exponents[j], p.numVariables);
            }
            /* coefficient is an integer */
            else if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0 && mpz_cmp_ui(r_denom, 1) == 0) {
                gmp_fprintf(outfile, "%Zd * ", r_num);
                fprint_monomial(outfile, p.exponents[j], p.numVariables);
            }
            /* coefficient is rational */
            else if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0 && mpz_cmp_ui(r_denom, 1) != 0) {
                gmp_fprintf(outfile, "(%Zd * ", r_num);
                fprint_monomial(outfile, p.exponents[j], p.numVariables);
                gmp_fprintf(outfile, ")/%Zd", r_denom);
            }
            /* coefficient is a Gaussian integer with no real part */
            else if (mpq_cmp_ui(p.coeff[j]->re, 0, 1) == 0 && mpz_cmp_ui(i_denom, 1) == 0) {
                gmp_fprintf(outfile, "I * %Zd * ", r_num);
                fprint_monomial(outfile, p.exponents[j], p.numVariables);
            }
            /* coefficient is a complex rational with no real part */
            else if (mpq_cmp_ui(p.coeff[j]->im, 0, 1) == 0 && mpz_cmp_ui(r_denom, 1) != 0) {
                gmp_fprintf(outfile, "I * (%Zd * ", r_num);
                fprint_monomial(outfile, p.exponents[j], p.numVariables);
                gmp_fprintf(outfile, ")/%Zd", r_denom);
            }
            /* anything else */
            else {
                gmp_fprintf(outfile, "(%Qd + I * %Qd) * ", p.coeff[j]->re, p.coeff[j]->im);
                fprint_monomial(outfile, p.exponents[j], p.numVariables);
            }

            if (j < p.numTerms - 1)
                fprintf(outfile, " + ");
            else
                fprintf(outfile, "\n");
        }
    }

    mpz_clear(r_num);
    mpz_clear(r_denom);
    mpz_clear(i_num);
    mpz_clear(i_denom);
}

/*********************************
 * create empty files for output *
 *********************************/
void initialize_output_files(polynomial_system *system, void *v, void *t, void *w, int num_points) {
    FILE *fh;

    errno = 0;

    fh = fopen(BH_FSUMMARY, "w");
    if (fh == NULL) {
        sprintf(error_string, "Couldn't open output file %s: %s", BH_FSUMMARY, strerror(errno));

        print_error(error_string, stderr);
        MPI_Abort(MPI_COMM_WORLD, BH_EXIT_BADFILE);
    }

    fprintf(fh, "Summary:\n\n");
    prog_info(fh);
    display_config(fh);
    fprint_input(fh, system, v, t, w, num_points);
    fputs("\n", fh);

    fclose(fh);

    errno = 0;

    fh = fopen(BH_FCONT, "w");
    if (fh == NULL) {
        sprintf(error_string, (size_t) termwidth, "Couldn't open output file %s: %s", BH_FCONT, strerror(errno));

        print_error(error_string, stderr);
        MPI_Abort(MPI_COMM_WORLD, BH_EXIT_BADFILE);
    }

    fprintf(fh, "The following intervals have been certified continuous:\n");

    fclose(fh);

    errno = 0;

    fh = fopen(BH_FDISCONT, "w");
    if (fh == NULL) {
        sprintf(error_string, "Couldn't open output file %s: %s", BH_FDISCONT, strerror(errno));

        print_error(error_string, stderr);
        MPI_Abort(MPI_COMM_WORLD, BH_EXIT_BADFILE);
    }

    fprintf(fh, "The following intervals have been certified discontinuous:\n");

    fclose(fh);

    errno = 0;

    fh = fopen(BH_FUNCERTIFIED, "w");
    if (fh == NULL) {
        sprintf(error_string, "Couldn't open output file %s: %s", BH_FUNCERTIFIED, strerror(errno));

        print_error(error_string, stderr);
        MPI_Abort(MPI_COMM_WORLD, BH_EXIT_BADFILE);
    }

    fprintf(fh, "The following intervals could not be certified:\n");

    fclose(fh);
}

/*****************************
 * print out a short summary *
 *****************************/
void summarize(int tested, int continuous, int discontinuous, int singularities) {
    puts("\nSUMMARY\n");
    printf("Number of intervals tested: %d\n", tested);
    printf("Number of intervals continuous: %d\n", continuous);
    printf("Number of intervals discontinuous: %d\n", discontinuous);
    if (singularities > 0)
        printf("Number of singularities found: %d\n", singularities);

    puts("\nThe following files have been created:");
    printf("%s:\t\tA general summary\n", BH_FSUMMARY);
    printf("%s:\t\tA summary of continuous intervals\n", BH_FCONT);
    printf("%s:\tA summary of discontinuous intervals\n", BH_FDISCONT);
    printf("%s:\tA summary of uncertified intervals\n", BH_FUNCERTIFIED);

    if (tested == continuous)
        printf("%s:\t\tAn updated points file\n", BH_FPTSOUT);

    /* no intervals were continuous; this should never happen */
    if (continuous == 0) {
        FILE *fh = fopen(BH_FCONT, "a");
        fputs("\nNo intervals were certified continuous\n", fh);
        fclose(fh);

        fh = NULL;
    }

    /* no intervals were discontinous */
    if (discontinuous == 0) {
        FILE *fh = fopen(BH_FDISCONT, "a");
        fputs("\nNo intervals were certified discontinuous\n", fh);
        fclose(fh);

        fh = NULL;
    }

    /* all intervals were certified */
    if (continuous + discontinuous == tested) {
        FILE *fh = fopen(BH_FUNCERTIFIED, "a");
        fputs("\nAll intervals were certified\n", fh);
        fclose(fh);

        fh = NULL;
    }
}
