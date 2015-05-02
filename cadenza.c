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
    int i, j, num_var = 0, num_points = 0, num_sing = 0, tested = 0, succeeded = 0, failed = 0;
    int rank, size, tmpsize, intbuf[2];
    char tmpbuf[BH_MAX_STRING];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    /* initial vector x */
    rational_complex_vector *x_rational = NULL;
    complex_vector *x_float = NULL;
    void *x = NULL;

    /* final vector x */
    void *x_final = NULL;

    /* singularities array */
    void *sing = NULL;

    /* init random seed */
    srand(time(NULL));

    /* get command-line arguments before anything else happens */
    getargs(argc, argv);

    /* tell em who we are */
//     if (rank == 0) { 
//         fputs("\n", stderr);
//         prog_info(stderr);
//     }

    /* do this before checking filenames */
    if (help_flag) {
        if (rank == 0)
            usage();
        MPI_Finalize();
        exit(BH_EXIT_SUCCESS);
    }

    /* ensure filenames are properly defined */
    if (rank == 0) {
        if (strcmp(sysfile, "") == 0) {
            print_error("You need to define a system file", stderr);
            usage();
            MPI_Abort(MPI_COMM_WORLD, BH_EXIT_BADFILE);
        } else if (strcmp(pointsfile, "") == 0) {
            print_error("You need to define a points file", stderr);
            usage();
            MPI_Abort(MPI_COMM_WORLD, BH_EXIT_BADFILE);
        }
    }

    /* display the options the user has selected */
    if (rank == 0 && verbosity > BH_LACONIC) {
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

    if (rank == 0) {
        read_system_file(&F, v);
        num_var = F.numVariables;
        num_points = read_points_file(&t_float, &x_float, num_var);
        initialize_output_files(&F, v, t_float, x_float, num_points);
        intbuf[0] = num_var;
        intbuf[1] = num_points;
    } 
    
    MPI_Bcast(intbuf, 2, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        num_var = intbuf[0];
        num_points = intbuf[1];
        x_float = malloc(num_points * sizeof(complex_vector));
        t_float = malloc(num_points * sizeof(mpf_t));
    }

    /* print out the system, vector and points */
    if (rank == 0 && verbosity > BH_CHATTY) {
        fprint_input(stdout, &F, v, t, x, num_points);
    }
    
    /* pass F to all processes */
    if (rank == 0) {
        for (i=1; i<size; i++) {
            send_polynomial_system(&F, i);
            send_complex_vector(v, i);
             for (j=0; j<num_points; j++) {
                 send_complex_vector(x_float[j], i);
                 tmpsize = mpfr_sprintf(tmpbuf, "%Rf", t_float[j]);
                 MPI_Send(&tmpsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                 MPI_Send(tmpbuf, tmpsize, MPI_CHAR, i, 1, MPI_COMM_WORLD);
             }
        }
    } else {
        recv_polynomial_system(&F, 0);
        recv_complex_vector(v, 0);
        for (j=0; j<num_points; j++) {
             recv_complex_vector(x_float[j], 0);
             MPI_Recv(&tmpsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
             MPI_Recv(tmpbuf, tmpsize, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status);
             tmpbuf[tmpsize] = '\0';
             mpf_init(t_float[j]);
             mpf_set_str(t_float[j], tmpbuf, 0);
         }
    }
    
    x = (void *) x_float;
    t = (void *) t_float;
//     if (rank == 1)
//         mpfr_printf("rank %d; t_float[0] = %Rf\n", rank, t_float[0]);

//     if (rank == 0) {
//         printf("rank 0\n");
//         print_points_float(stdout, x_float[0]);
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//     if (rank == 1) {
//         printf("rank 1\n");
//         print_points_float(stdout, x_float[0]);
//     }
    test_system(&F, v, t, x, num_points, &t_final, &x_final, &sing, &tested, &succeeded, &failed, &num_sing);

    /* print an output file only if all intervals certified continuous */
//     if (rank == 0 && tested == succeeded)
//         fprint_solutions(t_final, x_final, succeeded + 1);
//     
//     if (rank == 0)
//         summarize(tested, succeeded, failed, num_sing);

    /* clean up */
    clear_polynomial_system(&F);
    free_v(v);
    free_t(t, num_points);
    //free_t(t_final, succeeded + 1);
    free_x(x, num_points);
    //free_x(x_final, succeeded + 1);
    //free_x(sing, num_sing);
    free(error_string);

    MPI_Finalize();

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
    sort_order = BH_DESCENDING;
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
            {"ascending",  no_argument, &sort_order, BH_ASCENDING},
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

        c = getopt_long (argc, argv, "ahvfqp:s:m:n:d:c:", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0) break;

                printf ("option %s", long_options[option_index].name);

                if (optarg) printf (" with arg %s", optarg);

                puts("\n");
                break;

            case 'a':
                sort_order = BH_ASCENDING;
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
        default_precision = BH_PREC_MIN;
    if (newton_tolerance == BH_UNSET)
        newton_tolerance = BH_NEWT_TOLERANCE;
    if (subd_tolerance == BH_UNSET)
        subd_tolerance = BH_SUB_TOLERANCE;

    /* set significant digits */
    if (default_precision < 32)
        sigdig = 15;
    else
        sigdig = default_precision / 3;
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
 * free [rational_]complex_vector *x *
 *************************************/
void free_x(void *x, int num_points) {
    int i;

    if (arithmetic_type == BH_USE_FLOAT) {
        complex_vector *x_float = (complex_vector *) x;

        for (i = 0; i < num_points; i++)
            clear_vector(x_float[i]);

        free(x_float);
        x_float = NULL;
        x = NULL;
    } else {
        rational_complex_vector *x_rational = (rational_complex_vector *) x;

        for (i = 0; i < num_points; i++)
            clear_rational_vector(x_rational[i]);

        free(x_rational);
        x_rational = NULL;
        x = NULL;
    }
}
