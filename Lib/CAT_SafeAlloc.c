#include "CAT_SafeAlloc.h"

void *cat_safe_malloc(size_t bytes, const char *type, size_t count,
                       const char *file, int line) {
    void *ptr = malloc(bytes);
    if (!ptr) {
        fprintf(stderr, "OOM allocating %zu bytes (%s[%zu]) at %s:%d\n", bytes, type, count, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void *cat_safe_calloc(size_t count, size_t size, const char *type,
                       size_t total, const char *file, int line) {
    void *ptr = calloc(count, size);
    if (!ptr) {
        fprintf(stderr, "OOM allocating %zu bytes (%s[%zu]) at %s:%d\n", count*size, type, total, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

FILE *cat_safe_fopen(const char *path, const char *mode,
                      const char *file, int line) {
    FILE *fp = fopen(path, mode);
    if (!fp) {
        fprintf(stderr, "Failed to open file '%s' (mode %s) at %s:%d\n", path, mode, file, line);
        exit(EXIT_FAILURE);
    }
    return fp;
}
