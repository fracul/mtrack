#include <stdio.h>
#include "input.h"
#include "test.h"

int main(int argc, char ** argv)
{
  test_init(true);
  
  FILE * fp = fopen("test_input.input", "r");
  ASSERT_TRUE(fp != NULL);
  if(fp != NULL)
  {
    double cx, a, b, c, d, e, f, g, h, fortytwo, clight;
    double tmp;
    
    int n,m;

    double bogus;
    int bogusi;
    
    char str[128];

    /* Bogus lines shouldn't be accepted */
    ASSERT_FALSE(config_get_double(fp, "", "bogus1", &bogus));
    ASSERT_FALSE(config_get_double(fp, "bogus_land", "bogus1", &bogus));
    ASSERT_FALSE(config_get_double(fp, "bogus_land", "bogus2", &bogus));
    ASSERT_FALSE(config_get_double(fp, "bogus_land", "bogus3", &bogus));
    ASSERT_FALSE(config_get_int(fp, "bogus_land", "bogus4", &bogusi));
    
    /* Commented lines should'nt be read */
    
    ASSERT_FALSE(config_get_double(fp, "general", "cx", &cx));
    
    /*
     * Section testing
     */
     
    /* Looking in the bad section */
    ASSERT_FALSE(config_get_double(fp, "physic", "a", &a));
    ASSERT_FALSE(config_get_double(fp, "general", "clight", &clight));
    
    /* Looking for section existance */
    ASSERT_TRUE(config_have_section(fp, "general"));
    ASSERT_TRUE(config_have_section(fp, "physic"));
    ASSERT_FALSE(config_have_section(fp, "ghost"));
    
    /* Looking if original values preserved */
    tmp = -1;
    config_get_double(fp, "general", "DoNotExist", &tmp);
    ASSERT_EQUAL(tmp, -1);
    
    /* Ensure section uniqueness is ensured */
    ASSERT_FALSE(config_have_section(fp, "twin"));
    
    /*
     * Doubles that should be read
     */
    
    ASSERT_TRUE(config_get_double(fp, "general", "a", &a));
    ASSERT_EQUAL_FLOAT(a, 1.1, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "b", &b));
    ASSERT_EQUAL_FLOAT(b, 2.2, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "c", &c));
    ASSERT_EQUAL_FLOAT(c, 3.3, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "d", &d));
    ASSERT_EQUAL_FLOAT(d, 4.4, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "e", &e));
    ASSERT_EQUAL_FLOAT(e, 5.5, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "f", &f));
    ASSERT_EQUAL_FLOAT(f, 6.6, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "g", &g));
    ASSERT_EQUAL_FLOAT(g, 7.7, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "h", &h));
    ASSERT_EQUAL_FLOAT(h, 8.8, 1e-16);
    
    ASSERT_TRUE(config_get_double(fp, "general", "fortytwo", &fortytwo));
    ASSERT_EQUAL_FLOAT(fortytwo, 42.0, 1e-16); 
    
    ASSERT_TRUE(config_get_double(fp, "physic", "clight", &clight));
    ASSERT_EQUAL_FLOAT(clight, 2.997925e+08, 1e-16); 
    
    /* String */
    
    if(config_get_str(fp, "general", "message", str))
    ASSERT_EQUAL_STR(str, "hello World!");
    
    if(config_get_str(fp, "general", "messagewident", str))
    ASSERT_EQUAL_STR(str, "hello World!");
    
    
    /* Integers that should be read */
    
    ASSERT_TRUE(config_get_int(fp, "general", "n", &n));
    ASSERT_EQUAL(n, -2); 
    
    ASSERT_TRUE(config_get_int(fp, "general", "m", &m));
    ASSERT_EQUAL(m, 5); 
  }

  if(argc == 2)
  {
    FILE * fp2 = fopen(argv[1], "r");
    char identifier[256];
    char line[256];
    char trash[256];
    double n_double;
    int n_int;
    while(fp2!= NULL && !feof(fp2)) if(fgets(line, 256, fp2) != NULL)
    {
      if(sscanf(line, " [%[^\\]]] ", identifier) == 1)
        printf("Section %s\n", identifier);
      if(sscanf(line, " %[^ \t] %[=:] %d ", identifier, trash, &n_int) == 3)
        printf("Variable %s = %d\n", identifier, n_int);
      if(sscanf(line, " %[^ \t] %[=:] %lf ", identifier, trash, &n_double) == 3)
        printf("Variable %s = %g\n", identifier, n_double);
    }
    if(fp2 != NULL) fclose(fp2);
  }

  return test_backend();
}
