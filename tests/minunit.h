#ifndef MINUNIT_H
#define MINUNIT_H
#include <stdio.h>
static int tests_run = 0;
static int tests_failed = 0;
#define MU_ASSERT(message, test) do { if (!(test)) { printf("%s\n", message); tests_failed++; } } while (0)
#define MU_RUN_TEST(test) do { test(); tests_run++; } while (0)
#endif
