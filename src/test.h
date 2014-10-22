/*
 * @file
 * Test module
 */

#ifndef MBTRACK_TEST_H
#define MBTRACK_TEST_H

#include <stdbool.h>
#include <math.h>

typedef struct test_results
{
unsigned n;
unsigned errors;
unsigned warnings;
} test_results_t;

extern test_results_t test_results;

/**
 * Initialize test module (set assert count and error count to 0)
 */
void test_init(bool enable_fprint);

/**
 * Show test summary, return 0 if OKAY, return -1 if error(s)
 */
int test_backend();

/* Assertions */

bool
_assert(bool expression, const char * format, ...);

/**
 * General assertion, record an error if expression is not verified,
 * and display a custom error message provided in format.
 * Note: format is a printf-style format string, and can be followed by
 * additional parameters
 */
#define ASSERT(expression, format, ...) \
  _assert(expression, format, __VA_ARGS__)

/**
 * Assert expression is true, record an error if false
 */
#define ASSERT_TRUE(expression) \
  _assert(expression, "%s is not true", #expression)

/**
 * Assert expression is false, record an error if true
 */
#define ASSERT_FALSE(expression) \
  _assert(!expression, "%s is not false", #expression)

/**
 * Assert a equal b, record an error if not
 */
#define ASSERT_EQUAL(a, b) \
  _assert((a) == (b), "%s is not equal to %s", #a, #b)

/**
 * Assert a not equal b, record an error if equal
 */
#define ASSERT_NOT_EQUAL(a, b) \
  _assert((a) != (b), "%s is equal to %s", #a, #b)

/**
 * Assert a equal b, with an error margin of e
 */
#define ASSERT_EQUAL_FLOAT(a, b, e) \
  _assert(fabs(((double) a) - ((double) b)) <= e, "%s is not equal to %s (using a tolerance of %s)", #a, #b, #e)

/**
 * Assert string a equal string b, using strcmp, record an error if not
 */
#define ASSERT_EQUAL_STR(a, b) \
  _assert(strcmp(a, b) == 0, "%s is not equal to %s", a, b)

#endif /* MBTRACK_TEST_H */
