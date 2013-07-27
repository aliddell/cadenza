/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io.c: Generic input/output functions for Blue Harvest
 */
#include "blueharvest.h"

/********************************************************
 * print a helpful message about command-line arguments *
 ********************************************************/
void usage() {
    fputs("usage: blueharvest [-h|--help] [-v|--verbose] [-f|--float] [-q|--rational] [--precision INT] [--config FILENAME] --system FILENAME --points FILENAME\n", stderr);
}

/*****************************************
 * uniformly print error messages to outfile *
 *****************************************/
void print_error(char *msg, FILE *outfile) {
    fprintf(outfile, "\nERROR: %s\n\n", msg);
}

/*****************************************************
 * print a helpful message about this program to outfile *
 *****************************************************/
void prog_info(FILE *outfile) {
    char compile_date[BH_MAX_DATECHAR];
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    size_t res = strftime(compile_date, (size_t) BH_MAX_DATECHAR + 1, "%b %d, %Y", timeinfo);
    if (res == 0)
        strcpy(compile_date, "Jan 1, 1970");

    fprintf(outfile, "%s v%s (built %s) (compiled %s)\n", BH_PROGRAM_NAME, BH_VERSION, BH_BUILD_DATE, compile_date);
    fprintf(outfile, "%s\n", BH_AUTHORS);
    fprintf(outfile, "GMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
}

/******************************************
 * print the details of the configuration *
 ******************************************/
void display_config(FILE *outfile) {
    char precision_str[BH_TERMWIDTH];

    if (arithmetic_type == BH_USE_FLOAT)
        snprintf(precision_str, (size_t) BH_TERMWIDTH+1, "%d-bit precision floating point", default_precision);
    else
        strcpy(precision_str, "exact rational");

    fprintf(outfile, "Computing using %s arithmetic\n", precision_str);
    fprintf(outfile, "Polynomial system file is %s\n", sysfile);
    fprintf(outfile, "Point set file is %s\n", pointsfile);
    if (configfile != NULL)
        fprintf(outfile, "Configuration file is %s\n", configfile);
    fprintf(outfile, "Maximum Newton iterations: %d\n", newton_tolerance);
    fprintf(outfile, "Maximum segment subdivisions: %d\n\n", subd_tolerance);
}

/*************************************************************
 * read the configuration file, if given, and apply settings *
 *************************************************************/
void read_config_file() {
    int verb = BH_LACONIC, arith_type = BH_USE_RATIONAL, def_prec = MPFR_PREC_MIN,  newtol = BH_NEWT_TOLERANCE, subtol = BH_SUB_TOLERANCE, errno;
    char sysf[BH_MAX_FILENAME] = "";
    char pointsf[BH_MAX_FILENAME] = "";
    char line[BH_TERMWIDTH];

    FILE *config = fopen(configfile, "r");
    if (config == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open config file `%s': %s", configfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    /* read information from the file and copy it to local values */
    while (fgets(line, BH_TERMWIDTH, config) != NULL) {
        char *key = strtok(line, "=\n");
        char *val = strtok(NULL, "=\n");

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

            if (errno == ERANGE || def_prec < MPFR_PREC_MIN)
                def_prec = MPFR_PREC_MIN;
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
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (res == 0) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s", sysfile, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADREAD);
    } else if (errno == EILSEQ) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s", sysfile, strerror(errno));

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
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (res == 0) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s", sysfile, strerror(errno));

                print_error(error_string, stderr);
                exit(BH_EXIT_BADREAD);
            } else if (errno == EILSEQ) {
                char error_string[BH_TERMWIDTH];
                snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s", sysfile, strerror(errno));

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
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: unexpected EOF", sysfile);

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (res == 0) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADREAD);
        } else if (errno == EILSEQ) {
            char error_string[BH_TERMWIDTH];
            snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Error reading %s: %s", sysfile, strerror(errno));

            print_error(error_string, stderr);
            exit(BH_EXIT_BADPARSE);
        }
    }

    mpq_init(p.norm_sqr);
    norm_sqr_polynomial(p.norm_sqr, &p);

    return p;
}

/*********************************
 * create empty files for output *
 *********************************/
void initialize_output_files(polynomial_system *system, void *v, void *t, void *w, int num_points) {
    FILE *fh;

    errno = 0;

    fh = fopen(BH_FSUMMARY, "w");
    if (fh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s", BH_FSUMMARY, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
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
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s", BH_FCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(fh, "The following segments are continuous:\n");

    fclose(fh);

    errno = 0;

    fh = fopen(BH_FDISCONT, "w");
    if (fh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s", BH_FDISCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(fh, "The following segments are discontinuous:\n");

    fclose(fh);
}

/*****************************
 * print out a short summary *
 *****************************/
void summarize(int tested, int succeeded, int failed) {
    puts("\nSUMMARY\n");
    printf("Number of intervals tested: %d\n", tested);
    printf("Number of intervals continuous: %d\n", succeeded);
    printf("Number of intervals discontinuous: %d\n", failed);

    puts("\nThe following files have been created:");
    printf("%s:\t\tA general summary\n", BH_FSUMMARY);
    printf("%s:\t\tA summary of continuous intervals\n", BH_FCONT);
    printf("%s:\tA summary of discontinuous intervals\n", BH_FDISCONT);

    if (tested == succeeded)
        printf("%s:\t\tAn updated points file\n", BH_FPTSOUT);
}
