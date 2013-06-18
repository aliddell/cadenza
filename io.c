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

/********************************************
 * uniformly print error messages to stderr *
 ********************************************/
void print_error(char *msg) {
    fprintf(stderr, "\nERROR: %s\n\n", msg);
}

/********************************************************
 * print a helpful message about this program to stderr *
 ********************************************************/
void prog_info() {
    char compile_date[BH_MAX_DATECHAR];
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    size_t res = strftime(compile_date, (size_t) BH_MAX_DATECHAR + 1, "%b %d, %Y", timeinfo);
    if (res == 0)
        strcpy(compile_date, "Jan 1, 1970");

    fputs("\n", stderr);
    fprintf(stderr, "\t%s v%s (built %s) (compiled %s)\n", BH_PROGRAM_NAME, BH_VERSION, BH_BUILD_DATE, compile_date);
    fprintf(stderr, "\t%s\n", BH_AUTHORS);
    fprintf(stderr, "\tGMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
}

void display_config() {
    char str_arithmetic_type[15];
    if (arithmetic_type == BH_USE_FLOAT)
        strcpy(str_arithmetic_type, "floating point");
    else
        strcpy(str_arithmetic_type, "rational");

    fprintf(stderr, "\tComputing using %d-bit precision %s arithmetic\n", default_precision, str_arithmetic_type);
    fprintf(stderr, "\tPolynomial system file is %s\n", sysfile);
    fprintf(stderr, "\tPoint set file is %s\n", pointsfile);
}
