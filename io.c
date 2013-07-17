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
    fputs("usage: blueharvest [-h|--help] [-v|--verbose] [-f|--float] [-q|--rational] --precision INT --system FILENAME --points FILENAME\n", stderr);
}

/*****************************************
 * uniformly print error messages to OUT *
 *****************************************/
void print_error(char *msg, FILE *OUT) {
    fprintf(OUT, "\nERROR: %s\n\n", msg);
}

/*****************************************************
 * print a helpful message about this program to OUT *
 *****************************************************/
void prog_info(FILE *OUT) {
    char compile_date[BH_MAX_DATECHAR];
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    size_t res = strftime(compile_date, (size_t) BH_MAX_DATECHAR + 1, "%b %d, %Y", timeinfo);
    if (res == 0)
        strcpy(compile_date, "Jan 1, 1970");

    fputs("\n", OUT);
    fprintf(OUT, "\t%s v%s (built %s) (compiled %s)\n", BH_PROGRAM_NAME, BH_VERSION, BH_BUILD_DATE, compile_date);
    fprintf(OUT, "\t%s\n", BH_AUTHORS);
    fprintf(OUT, "\tGMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
}

/******************************************
 * print the details of the configuration *
 ******************************************/
void display_config(FILE *OUT) {
    char str_arithmetic_type[15];
    if (arithmetic_type == BH_USE_FLOAT)
        strcpy(str_arithmetic_type, "floating point");
    else
        strcpy(str_arithmetic_type, "rational");

    fprintf(OUT, "\tComputing using %d-bit precision %s arithmetic\n", default_precision, str_arithmetic_type);
    fprintf(OUT, "\tPolynomial system file is %s\n", sysfile);
    fprintf(OUT, "\tPoint set file is %s\n", pointsfile);
    if (configfile != NULL)
        fprintf(OUT, "\tConfiguration file is %s\n\n", configfile);
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

    errno = 0;

    fh = fopen(BH_FSUMMARY, "w");
    if (fh == NULL) {
        char error_string[BH_TERMWIDTH];
        snprintf(error_string, (size_t) BH_TERMWIDTH+1, "Couldn't open output file %s: %s.", BH_FSUMMARY, strerror(errno));

        print_error(error_string, stderr);
        exit(BH_EXIT_BADFILE);
    }

    fprintf(fh, "Summary:\n");

    fclose(fh);
}

