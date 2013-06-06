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

int
main(int argc, char* argv[])
{
    short verbose_flag, help_flag, arithmetic_type;

    getargs(argc, argv, &verbose_flag, &help_flag);

    if (argc > 1) {
        return 0;
    }

    return 0;
}

