/*
 * Cadenza
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * cadenza.c: Main file for Cadenza
 */
#include "cadenza.h"

int main(int argc, char *argv[]) {
    int num_var = 0, num_points = 0, num_sing = 0, tested = 0, succeeded = 0, failed = 0;

    /* set termwidth for maximum string length */
    termwidth = set_termwidth();

    error_string = malloc(termwidth * sizeof(char));

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

    /* singularities array */
    void *sing = NULL;

    /* init random seed */
    srand(time(NULL));

    /* get command-line arguments before anything else happens */
    getargs(argc, argv);

    /* tell em who we are */
    fputs("\n", stderr);
    prog_info(stderr);

    /* do this before checking filenames */
    if (help_flag) {
        usage();
        exit(BH_EXIT_SUCCESS);
    }

    /* ensure filenames are properly defined */
    if (strcmp(sysfile, "") == 0) {
        print_error("You need to define a system file", stderr);
        usage();
        exit(BH_EXIT_BADFILE);
    } else if (strcmp(pointsfile, "") == 0) {
        print_error("You need to define a points file", stderr);
        usage();
        exit(BH_EXIT_BADFILE);
    }

    /* display the options the user has selected */
    if (verbosity > BH_LACONIC) {
        fputs("\n", stderr);
        display_config(stderr);
    }

    /* set void and function pointers */
    set_function_pointers();
    if (arithmetic_type == BH_USE_FLOAT) {
        mpf_set_default_prec(default_precision);
        v = (void *) &v_float;
    } else {
        v = (void *) &v_rational;
    }

    read_system_file(&F, v);
    num_var = F.numVariables;

    num_points = read_points_file(&t, &w, num_var);

    /* print out the system, vector and points */
    if (verbosity > BH_CHATTY) {
        fprint_input(stdout, &F, v, t, w, num_points);
    }

    initialize_output_files(&F, v, t, w, num_points);
    test_system(&F, v, t, w, num_points, &t_final, &w_final, &sing, &tested, &succeeded, &failed, &num_sing);

    /* print an output file only if all intervals certified continuous */
    if (tested == succeeded)
        fprint_solutions(t_final, w_final, succeeded + 1);
    
    summarize(tested, succeeded, failed, num_sing);

    /* clean up */
    clear_polynomial_system(&F);
    free_v(v);
    free_t(t, num_points);
    free_t(t_final, succeeded + 1);
    free_w(w, num_points);
    free_w(w_final, succeeded + 1);
    free_w(sing, num_sing);
    free(error_string);

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
    subd_tolerance = BH_UNSET;
    strcpy(sysfile, "");
    strcpy(pointsfile, "");
    strcpy(configfile, "");

    int c = 0;

    while (c != -1) {
        static struct option long_options[] = {
            /* These options set a flag */
            {"laconic",    no_argument, &verbosity, BH_LACONIC},
            {"chatty",     no_argument, &verbosity, BH_CHATTY},
            {"verbose",    no_argument, &verbosity, BH_VERBOSE},
            {"loquacious", no_argument, &verbosity, BH_LOQUACIOUS},
            {"help",       no_argument, &help_flag, 1},
            {"float",      no_argument, &arithmetic_type, BH_USE_FLOAT},
            {"rational",   no_argument, &arithmetic_type, BH_USE_RATIONAL},
            /* These options set a global */
            {"points",       required_argument, 0, 'p'},
            {"system",       required_argument, 0, 's'},
            {"precision",    required_argument, 0, 'm'},
            {"newtons",       required_argument, 0, 'n'},
            {"subdivisions", required_argument, 0, 'd'},
            {"config",       required_argument, 0, 'c'},
            {0, 0, 0, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvfqp:s:m:n:d:c:", long_options, &option_index);
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

            case 'n':
                newton_tolerance = atoi(optarg);
                break;

            case 'd':
                subd_tolerance = atoi(optarg);
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

    if (strcmp(configfile, "") != 0)
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
    if (subd_tolerance == BH_UNSET)
        subd_tolerance = BH_SUB_TOLERANCE;
}

/**************************
 * set the terminal width *
 **************************/
int set_termwidth() {
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);

    return w.ws_col;
}

/**********************************************************************************
 * set function pointers to use the proper arithmetic type with a minimum of fuss *
 **********************************************************************************/
void set_function_pointers() {

    /* floating point functions and data types */
    if (arithmetic_type == BH_USE_FLOAT) {
        read_system_file = &read_system_file_float;
        read_points_file = &read_points_file_float;
        fprint_input = &fprint_input_float;
        test_system = &test_system_float;
        fprint_solutions = &fprint_solutions_float;
    }

    /* rational functions */
    else {
        read_system_file = &read_system_file_rational;
        read_points_file = &read_points_file_rational;
        fprint_input = &fprint_input_rational;
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

/********************
 * free mp[qf]_t *t *
 ********************/
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
