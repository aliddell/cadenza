/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io.c: Input/output functions for Blue Harvest
 */
#include "blueharvest.h"

/********************************************************
 * print a helpful message about command-line arguments *
 ********************************************************/
void usage() {
    fprintf(stderr, "usage: blueharvest [-h|--help] [-v|--verbose] [-f|--float] [-q|--rational] --system FILENAME --points FILENAME\n");
}

/********************************************
 * uniformly print error messages to stderr *
 ********************************************/
void print_error(char *msg) {
    fprintf(stderr, "\nERROR: %s\n\n", msg);
}

/**********************************************
 * print a helpful message about this program *
 **********************************************/
void prog_info(FILE *OUT) {
    char *compile_date = malloc((BH_MAX_DATECHAR+1) * sizeof(char));
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    size_t res = strftime(compile_date, (size_t) BH_MAX_DATECHAR, "%b %d, %Y", timeinfo);
    if (res == 0)
        compile_date = "Jan 1, 1970";

    fprintf(OUT, "\n");
    fprintf(OUT, "\t%s v%s (built %s) (compiled %s)\n", BH_PROGRAM_NAME, BH_VERSION, BH_BUILD_DATE, compile_date);
    fprintf(OUT, "\t%s\n", BH_AUTHORS);
    fprintf(OUT, "\tGMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());

    free(compile_date);
}

void display_config() {
    char str_arithmetic_type[15];
    if (arithmetic_type == BH_USE_FLOAT)
        strcpy(str_arithmetic_type, "floating point");
    else
        strcpy(str_arithmetic_type, "rational");

    printf("\tComputing using %s arithmetic\n", str_arithmetic_type);
    printf("\tPolynomial system file is %s\n", sysfile);
    printf("\tPoint set file is %s\n", pointsfile);
}

/***********************************
 * read the polynomial system file *
 ***********************************/
void read_poly_file(FILE *fp) {

}

/******************************
* parse polynomial from input *
******************************/
polynomial parse_polynomial(char *input) {
    polynomial p;

    return p;
}
