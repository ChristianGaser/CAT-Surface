/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "CAT_ArgvCheck.h"

/**
 * \brief Report leftover command-line arguments that look like options.
 *
 * ParseArgv() leaves an argument it does not recognize in argv instead of
 * failing, so a mistyped or removed flag is consumed as a file name later and
 * the tool fails with a confusing message about that "file". Call this after
 * ParseArgv() and before the positional arguments are read: it writes one
 * message per offending argument to stderr and returns how many it found, so
 * the caller can print its usage and exit.
 *
 * A negative number is a value, not an option -- several tools take one
 * positionally, e.g. the -1 index of CAT_SurfSeparatePolygons or the -0.5
 * extent of CAT_SurfCentral2Pial -- and so is a lone "-".
 *
 * \param argc       (in) argument count left by ParseArgv()
 * \param argv       (in) argument vector left by ParseArgv()
 * \param executable (in) name to print the messages under, usually argv[0]
 * \return number of arguments that look like options, 0 when there are none
 */
int
cat_check_unknown_options(int argc, char *argv[], const char *executable)
{
    const char *name = (executable && executable[0]) ? executable : "CAT";
    int i, n_unknown = 0;

    for (i = 1; i < argc; i++)
    {
        const char *arg = argv[i];
        char *end;

        if (!arg || arg[0] != '-' || arg[1] == '\0')
            continue; /* a file name, or "-" on its own */

        /* a number is a value even though it starts with a minus sign */
        (void)strtod(arg, &end);
        if (*end == '\0')
            continue;

        fprintf(stderr, "%s: unrecognized option \"%s\"\n", name, arg);
        n_unknown++;
    }

    return n_unknown;
}
