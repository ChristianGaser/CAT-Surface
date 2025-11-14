/* Safe allocation helpers for CAT-Surface
 * Centralizes NULL checking and error reporting for malloc/calloc/fopen.
 * Use SAFE_MALLOC(T,n), SAFE_CALLOC(T,n), SAFE_FOPEN(path,mode) to get
 * checked pointers that abort with a diagnostic on failure.
 */
#ifndef CAT_SAFE_ALLOC_H
#define CAT_SAFE_ALLOC_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void *cat_safe_malloc(size_t bytes, const char *type, size_t count,
                       const char *file, int line);
void *cat_safe_calloc(size_t count, size_t size, const char *type,
                       size_t total, const char *file, int line);
FILE *cat_safe_fopen(const char *path, const char *mode,
                      const char *file, int line);

#define SAFE_MALLOC(T, N) ( (T*) cat_safe_malloc(sizeof(T) * (size_t)(N), #T, (size_t)(N), __FILE__, __LINE__) )
#define SAFE_CALLOC(T, N) ( (T*) cat_safe_calloc((size_t)(N), sizeof(T), #T, (size_t)(N), __FILE__, __LINE__) )
#define SAFE_FOPEN(PATH, MODE) cat_safe_fopen((PATH), (MODE), __FILE__, __LINE__)

#ifdef __cplusplus
}
#endif

#endif /* CAT_SAFE_ALLOC_H */
