////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
// Base include for an R extension in C
#include <R.h>
// Many useful R functions (pchisq, ppois, etc.)
#include <Rmath.h>
// Other useful R functions (sort, etc.)
#include <R_ext/Utils.h>
// Standard C libraries (file I/O, etc.)
#include <stdlib.h>
#include <stdio.h>
////////////////////////////////////////////////////////////////////////////////
#ifndef LOCALOPT_FISHER
#define LOCALOPT_FISHER
// -----------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
// This function performs the "domainogram" computation, that is
// multi-resolution statistics based on Fisher's combination of pre-computed
// measurement p-values, and identifies an exhaustive set of localy optimal
// domains.
// The result is a multi-resolution segmentation of the given sequence of
// measurements p-values associated to a genomic profile.
//
// Required computation time:
// time = O(n_data x (n_data-1)/2)
//
// Specifically required computation memory (excluding storage of all inputs):
// memory = n_data x sizeof(double) + w_max x sizeof(double)
//
// NOTE: C array indexes start at 0 while R array indexes start at 1.
// Thus, when comparing the C versus R implementation of the same algoritm,
// differences with i and w index handling may be a little bit confusing.
////////////////////////////////////////////////////////////////////////////////
void localopt_fisher(
    double  *Yi,         // INPUT: sequence of measurements
    int	    *n_data,     // INPUT: number of measurements
    char	  **name, 		 // INPUT: name used for the output file
    int	    *w_max,	     // INPUT: maximal scale to explore as number of measurements
    int	    *w_min,	     // INPUT: minimal domain size as number of measurements
    double  *gamma,		   // INPUT: p-value improvement factor
    double	*threshold,	 // INPUT: p-value threshold for each scale w
    double	*Pw_min,	   // OUTPUT: best p-value for each scale w
    int	    *n_optimums  // OUTPUT: number of optimal domains
) {
  int n     = *n_data;
  int wmax  = *w_max;
  int wmin  = (*w_min);
  double lg = log(*gamma);

  Rprintf("[ localopt_fisher ] name = %s\tn.probes = %i\twmax = %i\twmin = %i\tgamma = %lf\n", *name, n, *w_max, (*w_min)+1, *gamma);

  // -------------------------------------------------------------------------
  // Output file initialization
  char filename[1024];
  *filename = '\0';
  sprintf(filename, "%s_RawSegments.txt", *name);
  FILE * optimums = fopen(filename, "w");
  fprintf(optimums, "[ localopt_fisher ] name = %s\tn.probes = %i\twmax = %i\twmin = %i\tgamma = %lf\n", *name, n, *w_max, (*w_min)+1, *gamma);
  fprintf(optimums, "i\tw\tPiw\n");

  // -------------------------------------------------------------------------
  // Optimization initialization
  (*n_optimums) = 0;
  int progress=0, n_iterations=0;
  int n_compute = wmax * (wmax - 1) / 2 + wmax * (n - wmax);
  Rprintf("|-----------------------------------------------------------------------------|\r");

  double S1 = Yi[0];
  double Fi[n];
  for(int i = 0; i < n; i++)
    Fi[i] = lg + pchisq(-2.0 * Yi[i], 1 * 2, 0, 1); // Fi[i] = lg + Yi[i]

  // -------------------------------------------------------------------------
  // Optimization walkthrough
  for(int w=1; w < wmax; w++) {

    int wn = n - (w+1) + 1;

    // Compute windowed p-value Piw for the first window
    S1 += Yi[w];
    double Siw = S1;
    double Piw = pchisq(-2.0 * Siw, w * 2, 0, 1);

    // Update minimal p-value Piw for the whole scale w
    if(Fi[1] < Fi[0]) Fi[0] = Fi[1];

    if(Piw < Fi[0]) {
      Fi[0] = Piw;
      // We found a local optimum, does it match segmentation constraints?
      if((w > wmin) && (Piw < threshold[w])) {
        fprintf(optimums, "%i\t%i\t%lf\n", w, w+1, Piw);
        (*n_optimums)++;
      }
    }

    // We initialize minimal Piw for whole scale w
    Pw_min[w] = Piw;

    for(int i=1; i <= wn; i++) {

      // Rprintf("i = %i\tw = %i\n", i, w);

      // Compute windowed p-value Piw
      Siw += Yi[i + w] - Yi[i - 1];
      Piw = pchisq(-2.0 * Siw, w * 2, 0, 1);

      // Update locally optimal p-value Piw
      if(Fi[i + 1] < Fi[i]) Fi[i] = Fi[i + 1];

      if(Piw < Fi[i]) {
        Fi[i] = Piw;
        // We found a local optimum, does it match segmentation constraints?
        if((w >= wmin) && (Piw < threshold[w])) {
          fprintf(optimums, "%i\t%i\t%lf\n", i+w+1, w+1, Piw);
          (*n_optimums)++;
        }
      }

      // Update minimal p-value Piw for the whole scale w
      if(Piw < Pw_min[w]) Pw_min[w] = Piw;

      n_iterations++;
      // Print computation progress
      if(floor(n_iterations * 78.0 / n_compute)==progress) {
        Rprintf("|");
        progress++;
      }
    }
  }
  fclose(optimums);

  Rprintf("\n");
}
/////////////////////////////////////////////////////////////////////////////////
#endif //  LOCALOPT_FISHER
