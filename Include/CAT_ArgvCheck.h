/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_ARGVCHECK_H_
#define _CAT_ARGVCHECK_H_

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
int cat_check_unknown_options(int argc, char *argv[], const char *executable);

#endif
