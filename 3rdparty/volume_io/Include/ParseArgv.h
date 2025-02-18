/*
 * ParseArgv.h --
 *
 *	Declarations for Tk-related things that are visible
 *	outside of the Tk module itself.
 *
 * Copyright 1989-1992 Regents of the University of California.
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any purpose and without
 * fee is hereby granted, provided that the above copyright
 * notice appear in all copies.  The University of California
 * makes no representations about the suitability of this
 * software for any purpose.  It is provided "as is" without
 * express or implied warranty.
 *
 * This file has been modified to be used only for argv parsing without
 * reference to tk, tcl or X11. Base on tk.h from tk2.3
 * $Header: /private-cvsroot/minc/libsrc/ParseArgv.h,v 6.5.2.1 2004/09/28 20:23:39 bert Exp $ SPRITE (Berkeley)
 */

#include "minc.h"
/*
 * Definitions that allow this header file to be used either with or
 * without ANSI C features like function prototypes.
 */

#undef _ANSI_ARGS_
#if ((defined(__STDC__) || defined(SABER)) && !defined(NO_PROTOTYPE)) || defined(__cplusplus)
#   define _ANSI_ARGS_(x)	x
#else
#   define _ANSI_ARGS_(x)	()
#endif

/*
 * Structure used to specify how to handle argv options.
 */

typedef struct {
    char *key;		/* The key string that flags the option in the
			 * argv array. */
    int type;		/* Indicates option type;  see below. */
    char *src;		/* Value to be used in setting dst;  usage
			 * depends on type. */
    char *dst;		/* Address of value to be modified;  usage
			 * depends on type. */
    char *help;		/* Documentation message describing this option. */
} ArgvInfo;

/*
 * Legal values for the type field of a ArgvInfo: see the user
 * documentation for details.
 */

#define ARGV_CONSTANT		15
#define ARGV_INT			16
#define ARGV_STRING			17
#define ARGV_LONG           100
#define ARGV_REST			19
#define ARGV_FLOAT			20
#define ARGV_FUNC			21
#define ARGV_GENFUNC			22
#define ARGV_HELP			23
#define ARGV_VERSION                    24
#define ARGV_VERINFO                    25
#define ARGV_END			27

/*
 * Flag bits for passing to ParseArgv:
 */

#define ARGV_NO_DEFAULTS		0x1
#define ARGV_NO_LEFTOVERS		0x2
#define ARGV_NO_ABBREV		0x4
#define ARGV_DONT_SKIP_FIRST_ARG	0x8
#define ARGV_NO_PRINT 0x10	/* bert- value was 0x16, which was bogus */

/*
 *--------------------------------------------------------------
 *
 * Exported procedures and variables.
 *
 *--------------------------------------------------------------
 */
#if defined(__cplusplus)
extern "C" {
#endif

MNCAPI int ParseArgv _ANSI_ARGS_((int *argcPtr, char **argv,
                                  ArgvInfo *argTable, int flags));

#if defined(__cplusplus)
}
#endif
