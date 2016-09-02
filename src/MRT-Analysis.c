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
#ifndef MULTI_RESOLUTION_ANALYSIS
#define MULTI_RESOLUTION_ANALYSIS
// -----------------------------------------------------------------------------
// Optimization algorithms
#include "localopt.fisher.h"
// -----------------------------------------------------------------------------
#define MAX_LINE_SIZE 4096

////////////////////////////////////////////////////////////////////////////////
// Data structure supporting the hierarchical analysis performed by function
// MRT_Analysis that follows generation of a multi-resolution segmentation by
// function localopt_fisher
////////////////////////////////////////////////////////////////////////////////
typedef struct MultiResolutionTree {

  struct MultiResolutionTree * up;	  // Container
  struct MultiResolutionTree * down;  // Leftmost contained domain
  struct MultiResolutionTree * left;  // Next domain left
  struct MultiResolutionTree * right; // Next domain right
  int id;
  int start;
  int end;
  int wmin;
  int wmax;
  double P;
  int i;
  int w;

} MRT;

// -----------------------------------------------------------------------------
void MRT_free(MRT * domain) {

  // Depth recursion
  if(domain->down != NULL) MRT_free(domain->down);

  // Breadth recursion
  if(domain->right != NULL) MRT_free(domain->right);

  free(domain);
}

// -----------------------------------------------------------------------------
void MRT_print_domain(MRT * d) {
  if(d!=NULL)
    Rprintf(
      "id=%i\ti=%i\tstart=%i\tend=%i\twmin=%i\twmax=%i\n",
      d->id, d->i, d->start, d->end, d->wmin, d->wmax
    );
  else
    Rprintf("NULL\n");
}

// -----------------------------------------------------------------------------
void MRT_create_domain(
    MRT * d,
    int id, int start, int end, int wmin, int wmax, double P, int i, int w
) {

  d->up = NULL; d->down = NULL; d->left = NULL; d->right = NULL;
  d->id      = id;
  d->start   = start;
  d->end     = end;
  d->wmin    = wmin;
  d->wmax    = wmax;
  d->P       = P;
  d->i       = i;
  d->w       = w;

}

// -----------------------------------------------------------------------------
void MRT_update_domain(MRT * d, int start, int end, double P, int i, int w) {

  // Rprintf("Update domain\t\tfrom\t");
  // MRT_print_domain(d);

  d->start   = start;
  d->end     = end;
  d->wmax    = w;
  if(P < d->P) {
    d->P = P;
    d->i = i;
    d->w = w;
  }

  // Rprintf("\t\t\tto\t");
  // MRT_print_domain(d);
}

// -----------------------------------------------------------------------------
int MRT_contained_domains(MRT ** first,  int start, int end, MRT ** last) {

  int n = 0;

  // Rprintf("Match hierarchy\t\tn=%i\tstart=%i\tend=%i\n", n, start, end);
  // Rprintf("\t\t\t<-F\t");
  // MRT_print_domain(*first);
  // Rprintf("\t\t\tL->\t");
  // MRT_print_domain(*last);

  if(*first != NULL) {

    // We skip domains that don't overlap with [start:end]
    while(((*first)->right!=NULL) && ((*first)->end < start)) {
      *first = (*first)->right;
    }

    *last = *first;

    // We scan only domains overlaping with [start:end]
    while(((*last)->right!=NULL) && ((*last)->start <= end)) {
      *last = (*last)->right;
      n++;
    }

    // We adjust the match
    if((*last)->start > end) *last = (*last)->left;
    else if((*last)->end > start) n++;

  }

  // Rprintf("Match result\t\tn=%i\tstart=%i\tend=%i\n", n, start, end);
  // Rprintf("\t\t\t<-F\t");
  // MRT_print_domain(*first);
  // Rprintf("\t\t\tL->\t");
  // MRT_print_domain(*last);

  return n;
}


// -----------------------------------------------------------------------------
void MRT_insert_domain(MRT * ground, MRT * domain, MRT * first) {

  // Rprintf("Insert\t\t\t\t");
  // MRT_print_domain(first);

  domain->up   = ground;

  // Is this our very first domain insertion?
  if(ground->down==NULL) {
    // Rprintf("       very first\t\t");
    // MRT_print_domain(domain);

    ground->down = domain;
    return;
  }
  if(first->left == NULL && first->start > domain->end) {
    // Rprintf("       leftmost\t\t\t");
    // MRT_print_domain(domain);

    domain->right = first;
    first->left = domain;

    ground->down = domain;
    return;
  }
  if(first->right == NULL && first->end < domain->start) {
    // Rprintf("       rightmost\t\t");
    // MRT_print_domain(domain);

    domain->left  = first;
    first->right = domain;
    return;
  }
  if(first->start > domain->end) {
    // Rprintf("       inter-left\t\t");
    // MRT_print_domain(domain);

    domain->right = first;
    domain->left  = first->left;

    domain->left->right = domain;
    domain->right->left = domain;
    return;
  }
  if(first->end < domain->start) {
    Rprintf("       inter-right\t\t");
    MRT_print_domain(domain);

    domain->left  = first;
    domain->right = first->right;

    domain->left->right = domain;
    domain->right->left = domain;
    return;
  }
  // Rprintf("       DEAD END !!!!!!!!!!!!!!!!");
  // MRT_print_domain(domain);
}

// -----------------------------------------------------------------------------
void MRT_insert_container(MRT * ground, MRT * domain, MRT * first, MRT * last, int n) {

  // Rprintf("Update hierarchy\t\tn=%i\n", n);
  // Rprintf("\t\t\t+E+\t");
  // MRT_print_domain(domain);
  // Rprintf("\t\t\t<-F\t");
  // MRT_print_domain(first);
  // Rprintf("\t\t\tL->\t");
  // MRT_print_domain(last);

  // Rewire left side  links
  domain->left = first->left;
  if(domain->left != NULL) domain->left->right = domain;
  first->left = NULL;

  // Rewire right side  links
  domain->right = last->right;
  if(domain->right != NULL) domain->right->left = domain;
  last->right = NULL;

  // Rewire up and down links
  domain->up   = ground;
  domain->down = first;
  for(int i=0; i<n; i++) {
    first->up = domain;
    first = first->right;
  }
  if(domain->left == NULL)
    ground->down = domain;
}

// -----------------------------------------------------------------------------
void MRT_fprint_preorder(
    MRT  * domain,
    FILE * all,
    FILE * smallest,
    FILE * largest,
    int  * n_all,		   // total number of domains
    int  * n_smallest, // number of Maximum Resolution Domains
    int  * n_largest	 // number of Maximum Scale Domains
) {
  // We write in three output files:
  //      all      => all domains
  //      smallest => Maximum Resolution Domains
  //      largest  => Maximum Scale Domains
  // Output file columns = container, start, end, wmin, wmax, P, i, w

  fprintf(all, "%i\t%i\t%i\t%i\t%i\t%i\t%lf\t%i\t%i\n",
          domain->id, domain->up->id,
          domain->start, domain->end, domain->wmin, domain->wmax,
          domain->P, domain->i + domain->w - 1, domain->w
  );
  (*n_all)++;

  if(domain->down == NULL) {
    fprintf(smallest, "%i\t%i\t%i\t%i\t%i\t%i\t%lf\t%i\t%i\n",
            domain->id, domain->up->id,
            domain->start, domain->end, domain->wmin, domain->wmax,
            domain->P, domain->i + domain->w - 1, domain->w
    );
    (*n_smallest)++;
  }

  if(domain->up->id == 0) {
    fprintf(largest, "%i\t%i\t%i\t%i\t%i\t%i\t%lf\t%i\t%i\n",
            domain->id, domain->up->id,
            domain->start, domain->end, domain->wmin, domain->wmax,
            domain->P, domain->i + domain->w - 1, domain->w
    );
    (*n_largest)++;
  }

  // Depth recursion
  if(domain->down != NULL)
    MRT_fprint_preorder(domain->down, all, smallest, largest, n_all, n_smallest, n_largest);

  // Breadth recursion
  if(domain->right != NULL)
    MRT_fprint_preorder(domain->right, all, smallest, largest, n_all, n_smallest, n_largest);
}

////////////////////////////////////////////////////////////////////////////////
// This function performs a hierachical analysis of the multi-resolution
// segmentation resulting from function localopt_fisher, extending segmented
// genomic intervals according to rules of internal consistency, and finally
// extracts the set of maximal resolution and maximal scale segmented domains.
////////////////////////////////////////////////////////////////////////////////
void MRT_Analysis (
    int	   *ndata,		  // INPUT: number of signal data
    char   **name,	    // INPUT: name
    int	   *n_all,		  // OUTPUT: total number of domains
    int	   *n_smallest, // OUTPUT: Maximum Resolution Domains
    int	   *n_largest   // OUTPUT: number of Maximum Scale Domains
) {
  int n  = *ndata;

  Rprintf("[ MRT_Analysis ] name = %s\tn = %i\n", *name, n);

  // -------------------------------------------------------------------------
  // Input file initialization
  char filename[MAX_LINE_SIZE];
  *filename = '\0';
  sprintf(filename,"%s_RawSegments.txt",*name);
  FILE * optimums = fopen(filename, "r");

  if((n < 1) || (optimums == NULL)) {
    Rprintf("%s_RawSegments.txt is empty or does not exist\n", filename);
    return;
  }

  char * trash;
  char line[MAX_LINE_SIZE];

  // Skip file header => optimization parameters
  trash = fgets(line, MAX_LINE_SIZE, optimums);
  // Skip file header => column names
  trash = fgets(line, MAX_LINE_SIZE, optimums);

  // -------------------------------------------------------------------------
  // Initialize the segmentation structure

  // We define ground as the virtual domain that includes all possible domains
  MRT ground;
  MRT_create_domain(& ground, 0, 0, 0, 0, 0, 0, 0, 0);

  MRT * first  = NULL;
  MRT * last   = NULL;
  MRT * domain = NULL;


  double P;
  int i, w;

  trash = fgets(line, MAX_LINE_SIZE, optimums); // Read the very first optimum
  sscanf(line, "%i\t%i\t%lf\n", & i, & w, & P);
  i = i - w + 1;

  double	best_P = P;
  int	best_i = i;
  int	last_w = w;

  int start = i;
  int end   = i + w - 1;

  // We initialize the segmentation mask
  int mask[n];
  for(int x = 1; x <= n; x++) mask[x - 1] = 0;
  for(int x = start; x <= end; x++) mask[x - 1] = 1;

  int overlap = FALSE, scale_up = FALSE;
  int id = 0;

  int progress = 1;

  // -------------------------------------------------------------------------
  // Analyze the segmentation structure

  while(!feof(optimums)) {

    // Read and parse next optimum
    trash = fgets(line, MAX_LINE_SIZE, optimums);
    sscanf(line, "%i\t%i\t%lf\n", & i, & w, & P);
    i = i - w + 1;

    // Rprintf("P=%lf\ti=%i\tw=%i\t| last_w=%i\tbest_i=%i\tbest_P=%lf\n", P, i, w, last_w, best_i, best_P);

    Rprintf("\r\t\traw segments = %i\t domains = %i",progress, id);
    progress++;

    // Have we just changed scale?
    scale_up = FALSE;
    if(w != last_w) scale_up = TRUE;

    // Is there an overlap with current domain?
    overlap = FALSE;
    if(i <= end + 1) overlap = TRUE;

    // So can we record current domain?
    if(scale_up || ! overlap) {

      // How current domain matches existing hierarchy?
      int n_contained = MRT_contained_domains(& first, start, end, & last);

      // Do we instanciate?
      if(n_contained!=1) {
        id++;
        domain = malloc(sizeof(MRT));
        MRT_create_domain(domain, id, start, end, last_w, last_w, best_P, best_i, last_w);
      }

      switch(n_contained) {
      case 1:	// We update an existing domain
        MRT_update_domain(first, start, end, best_P, best_i, last_w);
        domain = first;
        break;
      case 0: // We insert a new smallest domain
        MRT_insert_domain(& ground, domain, first);
        break;
      default: // We insert a new container
        MRT_insert_container(& ground, domain, first, last, n_contained);
      break;
      }

      // We update next search for contained domains
      first = domain;

      // We move to next optimum/possible domain
      best_P = P;
      best_i = i;
      start = i;
      end   = i + w - 1;
      // We update segmentation mask							/* THINK ABOUT PERFORMANCE IMPROVEMENT  */
      for(int x = start; x <= end; x++) mask[x - 1] = 1;		/* => We could ommit part of the update */
      // We extend current domain coordinates					/* => We could even avoid to use a mask */
      while((start > 1) && (mask[start - 2] > 0)) start--;
      while((end   < n) && (mask[end]       > 0)) end++;
    }
    else {
      // There is an overlap at constant scale
      // We update centroid and p-value of the current domain
      if(P < best_P) {
        best_P = P;
        best_i = i;
      }
      // We update segmentation mask							/* THINK ABOUT PERFORMANCE IMPROVEMENT */
      for(int x = end; x <= i + w - 1; x++) mask[x - 1] = 1;
      // We update and extend current domain coordinates
      end = i + w - 1;
      while((end   < n) && (mask[end] > 0)) end++;
    }

    if(scale_up)
      first = ground.down;

    last_w = w;
  }
  fclose(optimums);

  // -------------------------------------------------------------------------
  // Save results

  *filename = '\0';
  sprintf(filename,"%s_Domains.txt",*name);
  FILE * all = fopen(filename, "w");
  fprintf(all, "id\tcontainer\tstart\tend\twmin\twmax\tP\ti\tw\n");

  *filename = '\0';
  sprintf(filename,"%s_MaxResolutionDomains.txt",*name);
  FILE * smallest = fopen(filename, "w");
  fprintf(smallest, "id\tcontainer\tstart\tend\twmin\twmax\tP\ti\tw\n");

  *filename = '\0';
  sprintf(filename,"%s_MaxScaleDomains.txt",*name);
  FILE * largest = fopen(filename, "w");
  fprintf(largest, "id\tcontainer\tstart\tend\twmin\twmax\tP\ti\tw\n");

  if(ground.down != NULL) {
    MRT_fprint_preorder(ground.down, all, smallest, largest, n_all, n_smallest, n_largest);
  }

  Rprintf("\tMax. Resolution = %i\t Max. Scale = %i\n", *n_smallest, *n_largest);

  fclose(all);
  fclose(largest);
  fclose(smallest);

  if(ground.down != NULL) {
    MRT_free(ground.down);
  }
}
/////////////////////////////////////////////////////////////////////////////////
#endif //  MULTI_RESOLUTION_ANALYSIS
