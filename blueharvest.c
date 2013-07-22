/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * blueharvest.c: Main file for Blue Harvest
 */
#include "blueharvest.h"

int main(int argc, char *argv[]) {
    int i, num_var = 0, num_points = 0, tested = 0, succeeded = 0, failed = 0;

    /* polynomial system */
    polynomial_system F;
    
    /* vector v */
    rational_complex_vector v_rational;
    complex_vector v_float;
    void *v = NULL;

    /* initial vector t */
    mpq_t *t_rational = NULL;
    mpf_t *t_float = NULL;
    void *t = NULL;

    /* final vector t */
    void *t_final = NULL;

    /* initial vector w */
    rational_complex_vector *w_rational = NULL;
    complex_vector *w_float = NULL;
    void *w = NULL;

    /* final vector w */
    void *w_final = NULL;

    /* init random seed */
    srand(time(NULL));

    /* get command-line arguments before anything else happens */
    getargs(argc, argv);

    if (verbosity > BH_LACONIC)
        prog_info(stderr);

    /* do this before checking filenames, duh */
    if (help_flag) {
        usage();
        exit(BH_EXIT_SUCCESS);
    }

    /* ensure filenames are properly defined */
    if (strcmp(sysfile, "unset") == 0) {
        print_error("You need to define a system file.", stderr);
        usage();
        exit(BH_EXIT_BADFILE);
    } else if (strcmp(pointsfile, "unset") == 0) {
        print_error("You need to define a points file.", stderr);
        usage();
        exit(BH_EXIT_BADFILE);
    }

    if (verbosity > BH_CHATTY)
        display_config(stderr);

    /* set void and function pointers */
    set_function_pointers();
    if (arithmetic_type == BH_USE_FLOAT) {
        v = (void *) &v_float;
    } else {
        v = (void *) &v_rational;
    }

    read_system_file(&F, v);
    num_var = F.numVariables;

    num_points = read_points_file(&t, &w, num_var);

    w_rational = (rational_complex_vector *) w;
    t_rational = (mpq_t *) t;
    /* print out the system, vector and points */
    if (verbosity > BH_CHATTY) {

        fputs("F:\n", stderr);
        print_system_rational(&F, stderr);

        fputs("\n", stderr);

        fputs("v:\n", stderr);
        fprintf(stderr, "[");
        for (i = 0; i < num_var - 1; i++) {
            gmp_fprintf(stderr, "%Qd + %Qdi, ", v_rational->coord[i]->re, v_rational->coord[i]->im);
        }
        gmp_fprintf(stderr, "%Qd + %Qdi]\n", v_rational->coord[i]->re, v_rational->coord[i]->im);

        fputs("\n", stderr);

        fputs("(t_i, w_i)\n", stderr);
        for (i = 0; i < num_points; i++) {
            gmp_fprintf(stderr, "%Qd, ", t_rational[i]);
            print_points_rational(w_rational[i], stderr);
        }
    }

    initialize_output_files();
    test_system(&F, v, t, w, num_points, &t_final, &w_final, &tested, &succeeded, &failed);

    /* print an output file only if all intervals certified continuous */
    if (tested == succeeded)
        fprint_solutions(t_final, w_final, succeeded + 1);
    
    summarize(tested, succeeded, failed);

    /* clean up */
    clear_polynomial_system(&F);
    free_v(v);
    free_t(t, num_points);
    free_t(t_final, succeeded + 1);
    free_w(w, num_points);
    free_w(w_final, succeeded + 1);

    exit(BH_EXIT_SUCCESS);
}

/*******************************************************
 * parse command-line args and set flags and filenames *
 *******************************************************/
void getargs(int argc, char *argv[]) {
    /* define default values */
    help_flag = 0;
    verbosity = BH_UNSET;
    arithmetic_type = BH_UNSET;
    default_precision = BH_UNSET;
    newton_tolerance = BH_UNSET;
    strcpy(sysfile, "unset");
    strcpy(pointsfile, "unset");
    strcpy(configfile, "unset");

    int c = 0;

    while (c != -1) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"chatty",     no_argument, &verbosity, BH_CHATTY},
            {"verbose",    no_argument, &verbosity, BH_VERBOSE},
            {"loquacious", no_argument, &verbosity, BH_LOQUACIOUS},
            {"help",       no_argument, &help_flag, 1},
            {"float",      no_argument, &arithmetic_type, BH_USE_FLOAT},
            {"rational",   no_argument, &arithmetic_type, BH_USE_RATIONAL},
            /* These options don't set a flag. */
            {"points",     required_argument, 0, 'p'},
            {"system",     required_argument, 0, 's'},
            {"precision",  required_argument, 0, 'm'},
            {"config",     required_argument, 0, 'c'},
            {"tolerance",  required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvfqp:s:m:c:t:", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0) break;

                printf ("option %s", long_options[option_index].name);

                if (optarg) printf (" with arg %s", optarg);

                puts("\n");
                break;

            case 'c':
                strcpy(configfile, optarg);
                break;

            case 'f':
                arithmetic_type = BH_USE_FLOAT;
                break;

            case 'h':
                help_flag = 1;
                break;

            case 'm':
                default_precision = atoi(optarg);
                break;

            case 'p':
                strcpy(pointsfile, optarg);
                break;

            case 'q':
                arithmetic_type = BH_USE_RATIONAL;
                break;

            case 's':
                strcpy(sysfile, optarg);
                break;

            case 't':
                newton_tolerance = atoi(optarg);
                break;

            case 'v':
                if (verbosity == BH_UNSET)
                    verbosity = BH_LACONIC;

                verbosity++;
                if (verbosity > BH_LOQUACIOUS)
                    verbosity = BH_LOQUACIOUS;
                break;
        }
    }

    if (strcmp(configfile, "unset") != 0)
        read_config_file();

    /* if after reading flags and config file, options are still unset */
    if (verbosity == BH_UNSET)
        verbosity = BH_LACONIC;
    if (arithmetic_type == BH_UNSET)
        arithmetic_type = BH_USE_RATIONAL;
    if (default_precision == BH_UNSET)
        default_precision = MPFR_PREC_MIN;
    if (newton_tolerance == BH_UNSET)
        newton_tolerance = BH_NEWT_TOLERANCE;
}

/**********************************************************************************
 * set function pointers to use the proper arithmetic type with a minimum of fuss *
 **********************************************************************************/
void set_function_pointers() {

    /* floating point functions and data types */
    if (arithmetic_type == BH_USE_FLOAT) {
        read_system_file = &read_system_file_float;
        read_points_file = &read_points_file_float;
        test_system = &test_system_float;
        fprint_solutions = &fprint_solutions_float;
    }

    /* rational functions */
    else {
        read_system_file = &read_system_file_rational;
        read_points_file = &read_points_file_rational;
        test_system = &test_system_rational;
        fprint_solutions = &fprint_solutions_rational;
    }
}

/************************************
 * free [rational_]complex_vector v *
 ************************************/
void free_v(void *v) {
    if (arithmetic_type == BH_USE_FLOAT) {
        complex_vector *v_float = (complex_vector *) v;
        clear_vector(*v_float);
        v_float = NULL;
        v = NULL;
    } else {
        rational_complex_vector *v_rational = (rational_complex_vector *) v;
        clear_rational_vector(*v_rational);
        v_rational = NULL;
        v = NULL;
    }
}

/*******************
 * free mp[qf]_t t *
 *******************/
void free_t(void *t, int num_points) {
    int i;

    if (arithmetic_type == BH_USE_FLOAT) {
        mpf_t *t_float = (mpf_t *) t;

        for (i = 0; i < num_points; i++)
            mpf_clear(t_float[i]);

        free(t_float);
        t_float = NULL;
        t = NULL;
    } else {
        mpq_t *t_rational = (mpq_t *) t;

        for (i = 0; i < num_points; i++)
            mpq_clear(t_rational[i]);

        free(t_rational);
        t_rational = NULL;
        t = NULL;
    }
}

/*************************************
 * free [rational_]complex_vector *w *
 *************************************/
void free_w(void *w, int num_points) {
    int i;

    if (arithmetic_type == BH_USE_FLOAT) {
        complex_vector *w_float = (complex_vector *) w;

        for (i = 0; i < num_points; i++)
            clear_vector(w_float[i]);

        free(w_float);
        w_float = NULL;
        w = NULL;
    } else {
        rational_complex_vector *w_rational = (rational_complex_vector *) w;

        for (i = 0; i < num_points; i++)
            clear_rational_vector(w_rational[i]);

        free(w_rational);
        w_rational = NULL;
        w = NULL;
    }
}
