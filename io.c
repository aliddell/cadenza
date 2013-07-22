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

    fputs("\n", outfile);
    fprintf(outfile, "\t%s v%s (built %s) (compiled %s)\n", BH_PROGRAM_NAME, BH_VERSION, BH_BUILD_DATE, compile_date);
    fprintf(outfile, "\t%s\n", BH_AUTHORS);
    fprintf(outfile, "\tGMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
}

/*************************************************************
 * read the configuration file, if given, and apply settings *
 *************************************************************/
void read_config_file() {
    int verb = BH_LACONIC;
    int arith_type = BH_USE_RATIONAL;
    int def_prec = MPFR_PREC_MIN;
    int newtol = BH_NEWT_TOLERANCE;
    char sysf[BH_MAX_FILENAME] = "unset";
    char pointsf[BH_MAX_FILENAME] = "unset";
    char line[BH_TERMWIDTH];

    FILE *config = fopen(configfile, "r");

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
            errno = 0;

            if (strcmp(val, "rational") == 0)
                arith_type = BH_USE_RATIONAL;
            else if (strcmp(val, "float") == 0)
                arith_type = BH_USE_FLOAT;
        } else if (strcmp(key, "precision") == 0) {
            errno = 0;
            
            def_prec = (int) strtol(val, NULL, 10);

            if (errno == ERANGE || def_prec < MPFR_PREC_MIN)
                def_prec = MPFR_PREC_MIN;
        } else if (strcmp(key, "tolerance") == 0) {
            errno = 0;
            
            newtol = (int) strtol(val, NULL, 10);

            if (errno == ERANGE || newtol < 1)
                newtol = BH_NEWT_TOLERANCE;

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
    if (strcmp(sysfile, "unset") == 0)
        strcpy(sysfile, sysf);
    if (strcmp(pointsfile, "unset") == 0)
        strcpy(pointsfile, pointsf);
}

/******************************************
 * print the details of the configuration *
 ******************************************/
void display_config(FILE *outfile) {
    char str_arithmetic_type[15];
    if (arithmetic_type == BH_USE_FLOAT)
        strcpy(str_arithmetic_type, "floating point");
    else
        strcpy(str_arithmetic_type, "rational");

    fprintf(outfile, "\tComputing using %d-bit precision %s arithmetic\n", default_precision, str_arithmetic_type);
    fprintf(outfile, "\tPolynomial system file is %s\n", sysfile);
    fprintf(outfile, "\tPoint set file is %s\n", pointsfile);
    if (configfile != NULL)
        fprintf(outfile, "\tConfiguration file is %s\n\n", configfile);
    else
        puts("");
}

/*********************************
 * create empty files for output *
 *********************************/
void initialize_output_files() {
    errno = 0;
    FILE *fh = fopen(BH_FDISCONT, "w");
    if (fh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s.", BH_FDISCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(fh, "The following segments are discontinuous:\n");

    fclose(fh);

    errno = 0;

    fh = fopen(BH_FCONT, "w");
    if (fh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s.", BH_FCONT, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(fh, "The following segments are continuous:\n");

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
    printf("%s:\t\tA summary of continuous intervals\n", BH_FCONT);
    printf("%s:\tA summary of discontinuous intervals\n", BH_FDISCONT);

    if (tested == succeeded)
        printf("%s:\t\tAn updated points file\n", BH_FPTSOUT);
}
