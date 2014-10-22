/*
 * This tool compare integer/float values from two files, look for differences
 */

#include <math.h>
#include <stdio.h>

int main(int argc, char ** argv)
{
  if(argc != 3)
  {
    fprintf(stderr, "ERROR: wrong arguments\n");
    if(argc > 0)
      fprintf(stderr, "Usage: %s [data_file_1] [data_file_2]\n", argv[0]);
    else
      fprintf(stderr, "Usage: ./compare [data_file_1] [data_file_2]\n");
    return -1;
  }
  else
  {
    if(strcmp(argv[1], argv[2]) == 0)
      fprintf(stderr, "WARNING: Comparing %s with itself.\n", argv[1]);
    FILE * fp1 = fopen(argv[1], "r");
    if(fp1 == NULL)
    {
      fprintf(stderr, "ERROR: cannot read file %s\n", argv[1]);
      return -1;
    }
    FILE * fp2 = fopen(argv[2], "r");
    if(fp2 == NULL)
    {
      fprintf(stderr, "ERROR: cannot read file %s\n", argv[2]);
      return -1;
    }
    double maxe = 0.0; /* Maximum absolute error */
    double maxre = 0.0; /* Maximum relative error */
    unsigned count = 0; /* Number counter */
    unsigned diffcount = 0; /* Difference counter */
    while(!feof(fp1) && !feof(fp2))
    {
      double n1;
      int c1 = fscanf(fp1, "%lf", &n1);
      fgetc(fp1);
      double n2;
      int c2 = fscanf(fp2, "%lf", &n2);
      fgetc(fp2);
      if(c1 != c2) /* No float here, only for one file */
      {
        fprintf(stderr, "WARNING: Could not read a value in one of the files\n");
        count++;
        diffcount++;
      }
      else if(c1 != 1) /* No float here in any of the file */
      {
      }
      else /* Two float to compare */
      {
        count++;
        if(n1 != n2)
        {
          diffcount++;
          if(fabs(n1 - n2) > maxe)
          {
            maxe = fabs(n1 - n2);
          }
          double re = 0.0;
          if(n1 == 0)
          {
            re = fabs(n1 - n2)/n2;
          }
          else
          {
            re = fabs(n1 - n2)/n1;
          }
          if(re > maxre)
          {
            maxre = re;
          }
        }
      }
    }
    fclose(fp1);
    fclose(fp2);
    printf("Values read: %u\n", count);
    printf("Differences of value: %u\n", diffcount);
    printf("Maximum absolute error: %g\n", maxe);
    printf("Maximum relative error: %g\n", maxre);
    return 0;
  }
}
