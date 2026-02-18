#include "CAT_SafeAlloc.h"

/**
 * \brief Safely allocate memory and exit if allocation fails.
 *
 * Attempts to allocate memory using malloc(). If allocation fails, prints a detailed
 * error message with file/line information and terminates the program.
 *
 * \param bytes (in)  number of bytes to allocate
 * \param type  (in)  string describing data type (e.g., "float"); used in error messages
 * \param count (in)  number of elements (e.g., array length); used in error messages
 * \param file  (in)  source filename for error reporting (use __FILE__)
 * \param line  (in)  source line number for error reporting (use __LINE__)
 * \return             pointer to allocated memory (never NULL; exits on failure)
 */
void *cat_safe_malloc(size_t bytes, const char *type, size_t count,
                       const char *file, int line) {
    void *ptr = malloc(bytes);
    if (!ptr) {
        fprintf(stderr, "OOM allocating %zu bytes (%s[%zu]) at %s:%d\n", bytes, type, count, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/**
 * \brief Safely allocate and zero-initialize memory, exit if allocation fails.
 *
 * Attempts to allocate and zero-initialize memory using calloc(). If allocation fails,
 * prints a detailed error message with file/line information and terminates the program.
 * The allocated memory is guaranteed to be initialized to zero.
 *
 * \param count (in)  number of elements to allocate
 * \param size  (in)  size of each element in bytes
 * \param type  (in)  string describing data type (e.g., "float"); used in error messages
 * \param total (in)  total count (redundant with count; used in error messages for clarity)
 * \param file  (in)  source filename for error reporting (use __FILE__)
 * \param line  (in)  source line number for error reporting (use __LINE__)
 * \return             pointer to allocated and zero-initialized memory (never NULL; exits on failure)
 */
void *cat_safe_calloc(size_t count, size_t size, const char *type,
                       size_t total, const char *file, int line) {
    void *ptr = calloc(count, size);
    if (!ptr) {
        fprintf(stderr, "OOM allocating %zu bytes (%s[%zu]) at %s:%d\n", count*size, type, total, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/**
 * \brief Safely open a file and exit if open fails.
 *
 * Attempts to open a file using fopen(). If the open operation fails,
 * prints a detailed error message with file/line information and terminates the program.
 *
 * \param path  (in)  file path to open (relative or absolute)
 * \param mode  (in)  fopen()-style mode string ("r", "w", "rb", etc.)
 * \param file  (in)  source filename for error reporting (use __FILE__)
 * \param line  (in)  source line number for error reporting (use __LINE__)
 * \return             FILE pointer (never NULL; exits on failure)
 */
FILE *cat_safe_fopen(const char *path, const char *mode,
                      const char *file, int line) {
    FILE *fp = fopen(path, mode);
    if (!fp) {
        fprintf(stderr, "Failed to open file '%s' (mode %s) at %s:%d\n", path, mode, file, line);
        exit(EXIT_FAILURE);
    }
    return fp;
}
