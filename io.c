/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * io.c: Input/output functions for Blue Harvest
 */
#include "alphaCertified.h"

void getargs(int argc, char* argv[], short* verbose_flag, short* help_flag, short* arithmetic_type) {
    int c = 0;

    while (c != -1) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"verbose", no_argument,      verbose_flag, 1},
            {"brief",   no_argument,      verbose_flag, 0},
            {"help",   no_argument,       help_flag, 1},
            {"float",     no_argument,    0, 'f'},
            {"rational",  no_argument,    0, 'q'},
            /* These options don't set a flag.
            *  We distinguish them by their indices. */
            {"points",    required_argument, 0, 'p'},
            {"system",    required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hf:m:p:i:n:s:", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;

                printf ("option %s", long_options[option_index].name);

                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 'f':
                filename = optarg;
                break;

            case 'n':
                noise = atof(optarg);
                break;

            case 'm':
                mut_rate = atof(optarg);
                break;
        }
    }
}


/******************************
* parse polynomial from input *
******************************/
polynomial parse_polynomial(char* input) {
    input
}
