/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   loadSettings.c: Load the configurations from the settings file
*/

#include "alphaCertified.h"

void load_default_settings(configurations *S)
/***************************************************************\
* USAGE: load default settings into S                           *
\***************************************************************/
{
  S->arithmeticType = 0;
  S->startingPrecision = 96;
  S->algorithm = 2;
  S->refineDigits = 0;
  S->numRandomSystems = 2;
  S->randomDigits = 10;
  S->randomSeed = (unsigned int) time(NULL);
  S->newtonOnly = 0;
  S->newtonIts = 2;
  S->realityCheck = 1;
  S->realityTest = 0;

  return;
}

void move_to_eol(FILE *IN)
/***************************************************************\
* USAGE: Move to EOL (or EOF)                                   *
\***************************************************************/
{ 
  char ch;

  do
  {
    ch = fgetc(IN);
  } while (ch != '\n' && ch != EOF);
  
  return;
}

void update_settings_item(configurations *S, char *str, FILE *IN)
/***************************************************************\
* USAGE: determine which item str describes and update it       *
\***************************************************************/
{
  int tempInt;

  if (!strcmp(str, "ALGORITHM:"))
  { // read in the algorithm number
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (0 <= tempInt && tempInt <= 2)
    { // update algorithm
      S->algorithm = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'ALGORITHM' must be between 0 and 2.  The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "ARITHMETICTYPE:"))
  { // read in the arithmetic type
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (0 <= tempInt && tempInt <= 1)
    { // update arithmetic type
      S->arithmeticType = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'ARITHMETICTYPE' must be either 0 or 1.  The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "PRECISION:"))
  { // read in the precision
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (tempInt >= 64)
    { // update precision
      S->startingPrecision = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'PRECISION' must be >= 64.  The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "REFINEDIGITS:"))
  { // read in the number of digits
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (tempInt >= 0)
    { // update digits
      S->refineDigits = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'REFINEDIGITS' must be a nonnegative integer.  The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "NUMRANDOMSYSTEMS:"))
  { // read in the number of random systems
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (tempInt >= 2)
    { // update number of random systems
      S->numRandomSystems = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'NUMRANDOMSYSTEMS' must be >= 2.  The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "RANDOMDIGITS:"))
  { // read in the number of digits for random systems
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (tempInt > 0)
    { // update number of digits
      S->randomDigits = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'RANDOMDIGITS' must be positive.  The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "RANDOMSEED:"))
  { // read in the random seed
    unsigned int randSeed = 0;
    tempInt = fscanf(IN, "%u", &randSeed);

    // error checking
    if (tempInt > 0)
    { // update random seed
      S->randomSeed = randSeed;
    }
    else
    { // print error message
      printf("NOTE: 'RANDOMSEED' was not read in properly.\n");
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "NEWTONONLY:"))
  { // read in if performing only Newton iterations
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (0 <= tempInt && tempInt <= 1)
    { // update newton only
      S->newtonOnly = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'NEWTONONLY' must be either 0 or 1.  The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "NUMITERATIONS:"))
  { // read in number of newton iterations
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (tempInt > 0)
    { // update newtonIts
      S->newtonIts = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'NUMITERATIONS' must be a positive integer. The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "REALITYCHECK:"))
  { // read in the reality check configuration
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (tempInt == -1 || tempInt == 0 || tempInt == 1)
    { // update realityCheck
      S->realityCheck = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'REALITYCHECK' must be either -1, 0, or 1. The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (!strcmp(str, "REALITYTEST:"))
  { // read in the reality test configuration
    fscanf(IN, "%d", &tempInt);

    // error checking
    if (tempInt == 0 || tempInt == 1)
    { // update realityTest
      S->realityTest = tempInt;
    }
    else
    { // print error message
      printf("NOTE: 'REALITYTEST' must be either 0 or 1. The value %d has been ignored.\n", tempInt);
    }

    // read in the rest of the line
    move_to_eol(IN);
  }
  else if (strcmp(str, "\n") && strcmp(str, "\r\n")) // ignore blank lines
  { // display error message
    printf("NOTE: The following line has been ignored: \"%s", str);

    char ch;
    do
    { // read in next ch
      ch = fgetc(IN);

      if (ch != '\r' && ch != '\n' && ch != EOF)
      { // print ch
        printf("%c", ch);
      }
      else 
      { // move to eol (if not there) and break
        if (ch == '\r')
          move_to_eol(IN);
        break;
      }
    }
    while (1);
    printf("\"\n");
  }

  return;
}

int read_settings_item(configurations *S, FILE *IN)
/***************************************************************\
* USAGE: read a settings item from IN - return 0 if EOF         *
\***************************************************************/
{
  int strSize = 0;
  char ch, *str = NULL;

  // read in the next item
  ch = fgetc(IN);
  if (ch != EOF)
  { // reallocate str
    strSize = 1;
    str = errRealloc(str, strSize * sizeof(char));
    str[0] = ch;

    if (ch == ':' || ch == '\n')
    { // null terminate str
      str = errRealloc(str, (strSize + 1) * sizeof(char));
      str[strSize] = '\0';
      strSize++;
    }
    else
    { // add on to str
      while (1)
      { // read in next ch
        ch = fgetc(IN);

        if (ch != EOF)
        { // add to str
          str = errRealloc(str, (strSize + 1) * sizeof(char));
          str[strSize] = ch;
          strSize++;
        }
  
        // see if to break
        if (ch == ':' || ch == '\n' || ch == EOF)
        { // null terminate str
          str = errRealloc(str, (strSize + 1) * sizeof(char));
          str[strSize] = '\0';
          strSize++;
 
          // exit loop
          break;
        }
      }
    }

    // determine if this item is valid and update S
    update_settings_item(S, str, IN);
  }

  return (ch != EOF);
}

void load_settings_file(configurations *S, FILE *IN)
/***************************************************************\
* USAGE: load S from IN                                         *
\***************************************************************/
{
  // load default settings
  load_default_settings(S);

  // loop until complete
  while (read_settings_item(S,IN));
  
  return;
}

void load_settings(configurations *S, char *fileName)
/***************************************************************\
* USAGE: load S from fileName, if exists                        *
\***************************************************************/
{
  FILE *IN = fopen(fileName, "r");

  // see if file exists
  if (IN == NULL)
  { // load default settings
    load_default_settings(S);
  }
  else
  { // load settings from IN
    load_settings_file(S, IN);

    // close file
    fclose(IN);
  }
  printf("\n");

  return;
}


