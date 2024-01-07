#ifndef BICAPI
#if defined(_MSC_VER)
/* If we are building on the Microsoft C compiler, we want to
 * explicitly export all public functions from the DLL
 */
#define BICAPI __declspec(dllexport)
#else
#define BICAPI
#endif /* _MSC_VER not defined */
#endif /* BICAPI not defined */

#include "bicpl.h"

