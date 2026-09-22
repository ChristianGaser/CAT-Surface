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

/**
 * \brief Safely allocate memory and exit if allocation fails.
 *
 * \param bytes (in)  number of bytes to allocate
 * \param type  (in)  string describing data type (e.g., "float"); used in error messages
 * \param count (in)  number of elements (e.g., array length); used in error messages
 * \param file  (in)  source filename for error reporting (use __FILE__)
 * \param line  (in)  source line number for error reporting (use __LINE__)
 * \return             pointer to allocated memory (never NULL; exits on failure)
 */
void *cat_safe_malloc(size_t bytes, const char *type, size_t count,
                      const char *file, int line);
/**
 * \brief Safely allocate and zero-initialize memory, exit if allocation fails.
 *
 * \param count (in)  number of elements to allocate
 * \param size  (in)  size of each element in bytes
 * \param type  (in)  string describing data type (e.g., "float"); used in error messages
 * \param total (in)  total count (redundant with count; used in error messages for clarity)
 * \param file  (in)  source filename for error reporting (use __FILE__)
 * \param line  (in)  source line number for error reporting (use __LINE__)
 * \return             pointer to allocated and zero-initialized memory (never NULL; exits on failure)
 */
void *cat_safe_calloc(size_t count, size_t size, const char *type, size_t total,
                      const char *file, int line);
/**
 * \brief Safely open a file and exit if open fails.
 *
 * \param path  (in)  file path to open (relative or absolute)
 * \param mode  (in)  fopen()-style mode string ("r", "w", "rb", etc.)
 * \param file  (in)  source filename for error reporting (use __FILE__)
 * \param line  (in)  source line number for error reporting (use __LINE__)
 * \return             FILE pointer (never NULL; exits on failure)
 */
FILE *cat_safe_fopen(const char *path, const char *mode, const char *file,
                     int line);

#define SAFE_MALLOC(T, N) ( (T*) cat_safe_malloc(sizeof(T) * (size_t)(N), #T, (size_t)(N), __FILE__, __LINE__) )
#define SAFE_CALLOC(T, N) ( (T*) cat_safe_calloc((size_t)(N), sizeof(T), #T, (size_t)(N), __FILE__, __LINE__) )
#define SAFE_FOPEN(PATH, MODE) cat_safe_fopen((PATH), (MODE), __FILE__, __LINE__)

#ifdef __cplusplus
}
#endif

#endif /* CAT_SAFE_ALLOC_H */
