#include "init.h"


int nozero = 0 ;
int gcam_write_grad ; // defined in gcamorph.c for diags
int remove_cerebellum = 0 ;
int remove_lh = 0 ;
int remove_rh = 0 ;

int remove_bright =0 ;
int map_to_flash = 0 ;
double TRs[MAX_GCA_INPUTS] ;
double fas[MAX_GCA_INPUTS] ;
double TEs[MAX_GCA_INPUTS] ;

char         *Progname ;
GCA_MORPH_PARMS  parms ;

float gsmooth_sigma = -1 ;
int ninsertions = 0 ;
int insert_labels[MAX_INSERTIONS] ;
int insert_intensities[MAX_INSERTIONS] ;
int insert_coords[MAX_INSERTIONS][3] ;
int insert_whalf[MAX_INSERTIONS] ;

int avgs = 0 ;  /* for smoothing conditional densities */
int read_lta = 0 ;
char *T2_mask_fname = NULL ;
double T2_thresh = 0 ;
char *aparc_aseg_fname = NULL ;
char *mask_fname = NULL ;
char *norm_fname = NULL ;
int renormalize = 0 ;
int renormalize_new = 0 ;
int renormalize_align = 0 ;
int renormalize_align_after = 0 ;

int  renorm_with_histos = 0 ;

char *long_reg_fname = NULL ;
//int inverted_xform = 0 ;

char *write_gca_fname = NULL ;
float regularize = 0 ;
float regularize_mean = 0 ;
char *example_T1 = NULL ;
char *example_segmentation = NULL ;
int register_wm_flag = 0 ;

double TR = -1 ;
double alpha = -1 ;
double TE = -1 ;
char *tl_fname = NULL ;

int nreads = 0 ;
char *read_intensity_fname[MAX_READS] ;
char *sample_fname = NULL ;
char *transformed_sample_fname = NULL ;
char *normalized_transformed_sample_fname = NULL ;
char *ctl_point_fname = NULL ;
int novar = 1 ;
int reinit = 0 ;

int use_contrast = 0 ;
float min_prior = MIN_PRIOR ;
int reset = 0 ;

FILE *diag_fp = NULL ;

int translation_only = 0 ;
int get_option(int argc, char *argv[]) ;
int write_vector_field(MRI *mri, GCA_MORPH *gcam, char *vf_fname) ;
int remove_bright_stuff(MRI *mri, GCA *gca, TRANSFORM *transform) ;
void print_help(void);

char *twm_fname = NULL ;  // file with manually specified temporal lobe white matter points
char *renormalization_fname = NULL ;
char *tissue_parms_fname = NULL ;
int center = 1 ;
int nreductions = 1 ;
char *xform_name = NULL ;
int noscale = 0 ;
int transform_loaded = 0 ;
char *gca_mean_fname = NULL ;
TRANSFORM  *transform = NULL ;
char *vf_fname = NULL ;

double blur_sigma = 0.0f ;

int handle_expanded_ventricles = 0;

int do_secondpass_renorm = 0;

double ctl_point_pct = DEFAULT_CTL_POINT_PCT ;

char *rusage_file; // =NULL

//// INIT

char         *gca_fname, *in_fname, *out_fname, fname[STRLEN], **av ;
MRI          *mri_inputs, *mri_tmp ;
GCA          *gca /*, *gca_tmp, *gca_reduced*/ ;
int          ac, nargs, ninputs, input, extra = 0 ;
int          msec, hours, minutes, seconds /*, iter*/ ;
int          n_omp_threads;
struct timeb start ;
GCA_MORPH    *gcam ;
float        label_scales[MAX_CMA_LABELS], label_offsets[MAX_CMA_LABELS];
float        label_peaks[MAX_CMA_LABELS];
int          label_computed[MAX_CMA_LABELS];
int          got_scales =0;



