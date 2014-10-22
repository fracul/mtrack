#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
#include "test.h"

/* Global variable for test results */
test_results_t test_results;
bool test_enable_fprint;

void test_init(bool enable_fprint)
{
  test_results.n = 0;
  test_results.errors = 0;
  test_results.warnings = 0;
  test_enable_fprint = enable_fprint;
}

int test_backend()
{
  if(test_results.errors == 0)
  {
    printf("Test PASSED\n  Asserts: %d\n  Errors: %d\n  Warnings: %d\n",
           test_results.n, test_results.errors, test_results.warnings);
    return 0;
  }
  else
  {
    fprintf(stderr, "Test FAILED\n  Asserts: %d\n  Errors: %d\n  Warnings: %d\n",
           test_results.n, test_results.errors, test_results.warnings);
    return -1;
  }
}

bool
_assert(bool expression, const char * format, ...)
{
  va_list args;
  va_start(args, format);
  
  (test_results.n)++;
  if(expression)
  {
    return true;
  }
  else
  {
    (test_results.errors)++;
    if(test_enable_fprint)
    {
      fprintf(stderr, "Test %d failed ! ", test_results.n);
      vfprintf(stderr, format, args);
      fprintf(stderr, "\n");
    }
    return false;
  }
  
  va_end(args);
}

