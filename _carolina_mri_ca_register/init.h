#ifndef __INIT_H__
#define __INIT_H__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "romp_support.h"

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"
#include "gcamorph.h"
#include "transform.h"
#include "mrisegment.h"
#include "version.h"
// #include "mri_ca_register.help.xml.h" // JR: a che cavolo serve??
#include "mri2.h"
#include "fsinit.h"
#include "ctrpoints.h"

#ifdef FS_CUDA
#include "devicemanagement.h"
#endif
#include "gcamorphtestutils.h"

extern int nozero;
extern int gcam_write_grad ; // defined in gcamorph.c for diags
extern int remove_cerebellum;
extern int remove_lh;
extern int remove_rh;

extern int remove_bright;
extern int map_to_flash;
extern double TRs[MAX_GCA_INPUTS] ;
extern double fas[MAX_GCA_INPUTS] ;
extern double TEs[MAX_GCA_INPUTS] ;

char         *Progname ;
extern GCA_MORPH_PARMS  parms ;

extern float gsmooth_sigma;
extern int ninsertions;
extern int insert_labels[MAX_INSERTIONS] ;
extern int insert_intensities[MAX_INSERTIONS] ;
extern int insert_coords[MAX_INSERTIONS][3] ;
extern int insert_whalf[MAX_INSERTIONS] ;

extern int avgs;  /* for smoothing conditional densities */
extern int read_lta;
extern char *T2_mask_fname;
extern double T2_thresh;
extern char *aparc_aseg_fname;
extern char *mask_fname;
extern char *norm_fname;
extern int renormalize;
extern int renormalize_new;
extern int renormalize_align;
extern int renormalize_align_after;

extern int  renorm_with_histos;

extern char *long_reg_fname;
//extern int inverted_xform = 0 ;

extern char *write_gca_fname;
extern float regularize;
extern float regularize_mean;
extern char *example_T1;
extern char *example_segmentation;
extern int register_wm_flag;

extern double TR;
extern double alpha;
extern double TE;
extern char *tl_fname;

#define MAX_READS 100
extern int nreads;
extern char *read_intensity_fname[MAX_READS] ;
extern char *sample_fname;
extern char *transformed_sample_fname;
extern char *normalized_transformed_sample_fname;
extern char *ctl_point_fname;
extern int novar;
extern int reinit;

extern int use_contrast;
extern float min_prior;
extern int reset;

extern FILE *diag_fp;

extern int translation_only;


extern char *twm_fname;  // file with manually specified temporal lobe white matter points
extern char *renormalization_fname;
extern char *tissue_parms_fname;
extern int center;
extern int nreductions;
extern char *xform_name;
extern int noscale;
extern int transform_loaded;
extern char *gca_mean_fname;
extern TRANSFORM  *transform;
extern char *vf_fname;

extern double blur_sigma;

extern int handle_expanded_ventricles;

extern int do_secondpass_renorm;

#define MM_FROM_EXTERIOR  5  // distance into brain mask to go when erasing super bright CSF voxels
/*
 command line consists of three inputs:
 
 argv[1]  - directory containing 'canonical' brain
 argv[2]  - directory containing brain to be registered
 argv[3]  - directory in which to write out registered brain.
 */

#define NPARMS           12
#define DEFAULT_CTL_POINT_PCT   .25
extern double ctl_point_pct;

extern char *rusage_file; // =NULL

//// INIT
extern char         *gca_fname, *in_fname, *out_fname, fname[STRLEN], **av ;
extern MRI          *mri_inputs, *mri_tmp ;
extern GCA          *gca /*, *gca_tmp, *gca_reduced*/ ;
extern int          ac, nargs, ninputs, input, extra;
extern int          msec, hours, minutes, seconds /*, iter*/ ;
extern int          n_omp_threads;
extern struct timeb start ;
extern GCA_MORPH    *gcam ;
// for GCA Renormalization with Alignment (if called sequentially)
extern float        label_scales[MAX_CMA_LABELS], label_offsets[MAX_CMA_LABELS];
extern float        label_peaks[MAX_CMA_LABELS];
extern int          label_computed[MAX_CMA_LABELS];
extern int          got_scales;

//////////////////////////////
int init_jr(int argc, char *argv[]);
int get_option(int argc, char *argv[]);
void print_help(void);
int write_vector_field(MRI *mri, GCA_MORPH *gcam, char *vf_fname);
int remove_bright_stuff(MRI *mri, GCA *gca, TRANSFORM *transform);


#endif
