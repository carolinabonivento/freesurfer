#include "init.h"

int main(int argc, char *argv[]){
    
    int err;
    
    err = init_jr(argc,argv);
    
    if (err){
        printf("error during init_jr");
        return err;
    }
    /////////////////////////////////////////////////
    // -T transform option
    // transform is loaded at get_opt() with -T using TransformRead()
    // assumed to be vox-to-vox
    if (!transform_loaded)   /* wasn't preloaded */
    {
        transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    }
    else
        // calculate inverse and cache it
    {
        TransformInvert(transform, mri_inputs) ;
    }
    
    
    /////////////////////////////////////////////////
    // -novar option  (default novar = 1)
    if (novar)
    {
        GCAunifyVariance(gca) ;
    }
    
    
    
    /////////////////////////////////////////////////
    // XGRAD or YGRAD or ZGRAD set
    // store (x,y,z)gradient info into mri_inputs
    if (gca->flags & GCA_GRAD)
    {
        int i, start = ninputs ;
        MRI *mri_kernel, *mri_smooth, *mri_grad, *mri_tmp ;
        
        mri_kernel = MRIgaussian1d(1.0, 30) ;
        mri_smooth = MRIconvolveGaussian(mri_inputs, NULL, mri_kernel) ;
        
        if (mri_inputs->type != MRI_FLOAT)
        {
            // change data to float
            mri_tmp = MRISeqchangeType(mri_inputs, MRI_FLOAT, 0, 0, 1) ;
            MRIfree(&mri_inputs) ;
            mri_inputs = mri_tmp ;
        }
        start = ninputs ;
        if (gca->flags & GCA_XGRAD)
        {
            for (i = 0 ; i < ninputs ; i++)
            {
                mri_grad = MRIxSobel(mri_smooth, NULL, i) ;
                MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
                MRIfree(&mri_grad) ;
            }
            start += ninputs ;
        }
        if (gca->flags & GCA_YGRAD)
        {
            for (i = 0 ; i < ninputs ; i++)
            {
                mri_grad = MRIySobel(mri_smooth, NULL, i) ;
                MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
                MRIfree(&mri_grad) ;
            }
            start += ninputs ;
        }
        if (gca->flags & GCA_ZGRAD)
        {
            for (i = 0 ; i < ninputs ; i++)
            {
                mri_grad = MRIzSobel(mri_smooth, NULL, i) ;
                MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
                MRIfree(&mri_grad) ;
            }
            start += ninputs ;
        }
        
        MRIfree(&mri_kernel) ;
        MRIfree(&mri_smooth) ;
    }
    
    
    //////////////////////////////////////////////////////////
    // -regularize val option (default = 0)
    if (regularize > 0)
    {
        GCAregularizeCovariance(gca, regularize) ;
    }
    
    
    
    //////////////////////////////////////////////////////////
    // -X prev.m3d option
    if (xform_name)
    {
        gcam = GCAMread(xform_name) ;
        if (!gcam)
            ErrorExit(ERROR_NOFILE,
                      "%s: could not read transform from %s", Progname, xform_name) ;
        if (long_reg_fname && strcmp(long_reg_fname,"identity.nofile") != 0)
        {
            TRANSFORM *transform_long ;
            
            transform_long = TransformRead(long_reg_fname) ;
            if (transform_long == NULL)
                ErrorExit(ERROR_NOFILE,
                          "%s: could not read longitudinal registration file %s",
                          Progname, long_reg_fname) ;
            
            //      if (inverted_xform)
            // {
            //   TransformInvert(transform_long, mri_inputs) ;
            //   TransformSwapInverse(transform_long) ;
            // }
            TransformInvert(transform_long, mri_inputs);
            GCAMapplyInverseTransform(gcam, transform_long) ;
            TransformFree(&transform_long) ;
            //      GCAMwrite(gcam, "combined_gcam.m3z");
        }
        {
            char fname[STRLEN] ;
            MRI  *mri ;
            sprintf(fname, "%s.invalid.mgz", parms.base_name) ;
            mri = GCAMwriteMRI(gcam, NULL, GCAM_INVALID) ;
            printf("writing %s\n", fname) ;
            MRIwrite(mri, fname) ;
            MRIfree(&mri) ;
            sprintf(fname, "%s.status.mgz", parms.base_name) ;
            mri = GCAMwriteMRI(gcam, NULL, GCAM_STATUS) ;
            printf("writing %s\n", fname) ;
            MRIwrite(mri, fname) ;
            MRIfree(&mri) ;
        }
        printf("NON DOVREBBE SUCCEDERE");
        
    }
    else   // default is to create one
    {
        gcam = GCAMalloc(gca->prior_width, gca->prior_height, gca->prior_depth) ;
    }
    
    
    //////////////////////////////////////////////////////////
    // -debug_voxel Gvx Gvy Gvz option
    //////////////////////////////////////////////////////////
    // -debug_voxel Gvx Gvy Gvz option
    if (Gvx > 0)
    {
        float xf, yf, zf ;
        
        if (xform_name)
        {
            GCAMinvert(gcam, mri_inputs) ;
            GCAMsampleInverseMorph(gcam, Gvx, Gvy, Gvz, &xf, &yf, &zf) ;
        }
        else
        {
            TransformInvert(transform, mri_inputs);
            TransformSample(transform, Gvx, Gvy, Gvz, &xf, &yf, &zf) ;
        }
        
        Gsx = nint(xf) ;
        Gsy = nint(yf) ;
        Gsz = nint(zf) ;
        printf("mapping by transform (%d, %d, %d) --> "
               "(%d, %d, %d) for rgb writing\n",
               Gvx, Gvy, Gvz, Gsx, Gsy, Gsz) ;
    }
    
    
    
    //////////////////////////////////////////////////////////////////
    // GCM initialization
    if (!xform_name)  /* only if the transform wasn't previously created */
        GCAMinit(gcam, mri_inputs, gca, transform,
                 parms.relabel_avgs >= parms.navgs) ;
    else
    {
        // added by xhan
        int x, y, z, n, label, max_n, max_label;
        float max_p;
        GC1D *gc;
        GCA_MORPH_NODE  *gcamn ;
        GCA_PRIOR *gcap;
        
        gcam->ninputs = mri_inputs->nframes ;
        getVolGeom(mri_inputs, &gcam->image);
        GCAsetVolGeom(gca, &gcam->atlas);
        gcam->gca = gca ;
        gcam->spacing = gca->prior_spacing;
        
        // use gca information
        for (x = 0 ; x < gcam->width ; x++)
        {
            for (y = 0 ; y < gcam->height ; y++)
            {
                for (z = 0 ; z < gcam->depth ; z++)
                {
                    gcamn = &gcam->nodes[x][y][z] ;
                    gcap = &gca->priors[x][y][z] ;
                    max_p = 0 ;
                    max_n = -1 ;
                    max_label = 0 ;
                    
                    // find the label which has the max p
                    for (n = 0 ; n < gcap->nlabels ; n++)
                    {
                        label = gcap->labels[n] ;   // get prior label
                        if (label == Gdiag_no)
                        {
                            DiagBreak() ;
                        }
                        if (label >= MAX_CMA_LABEL)
                        {
                            printf("invalid label %d at (%d, %d, %d) in prior volume\n",
                                   label, x, y, z);
                        }
                        if (gcap->priors[n] >= max_p) // update the max_p and max_label
                        {
                            max_n = n ;
                            max_p = gcap->priors[n] ;
                            max_label = gcap->labels[n] ;
                        }
                    }
                    
                    gcamn->label = max_label ;
                    gcamn->n = max_n ;
                    gcamn->prior = max_p ;
                    gc = GCAfindPriorGC(gca, x, y, z, max_label) ;
                    // gc can be NULL
                    gcamn->gc = gc ;
                    gcamn->log_p = 0 ;
                    
                }
            }
        }
        
        GCAMcomputeOriginalProperties(gcam) ;
        if (parms.relabel_avgs >= parms.navgs)
        {
            GCAMcomputeLabels(mri_inputs, gcam) ;
        }
        else
        {
            GCAMcomputeMaxPriorLabels(gcam) ;
        }
    }
    
    
    if (regularize_mean > 0)
    {
        GCAregularizeConditionalDensities(gca, regularize_mean) ;
    }
    
    
    if (parms.write_iterations != 0)
    {
        char fname[STRLEN] ;
        
        if (parms.diag_morph_from_atlas )
        {
            sprintf(fname, "%s_target", parms.base_name) ;
            MRIwriteImageViews(mri_inputs, fname, IMAGE_SIZE) ;
            sprintf(fname, "%s_target.mgz", parms.base_name) ;
            printf("writing target volume to %s...\n", fname) ;
            MRIwrite(mri_inputs, fname) ;
        }
        else
        {
            MRI  *mri_gca ;
            
            mri_gca = MRIalloc(gcam->atlas.width, gcam->atlas.height, gcam->atlas.depth, MRI_FLOAT) ;
            MRIcopyHeader(mri_inputs, mri_gca) ;
            GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
#if 0
            if (mri_gca->nframes > 1)
            {
                MRI *mri_tmp ;
                printf("careg: extracting %dth frame\n", mri_gca->nframes-1) ;
                mri_tmp = MRIcopyFrame(mri_gca, NULL, mri_gca->nframes-1, 0) ;
                MRIfree(&mri_gca) ;
                mri_gca = mri_tmp ;
            }
#endif
            sprintf(fname, "%s_target", parms.base_name) ;
            MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
            sprintf(fname, "%s_target.mgz", parms.base_name) ;
            printf("writing target volume to %s...\n", fname) ;
            MRIwrite(mri_gca, fname) ;
            MRIfree(&mri_gca) ;
        }
    }
    
    if (renormalize_align_after) // 1st morph should be smooth
    {
        parms.tol *= 5 ;
        parms.l_smoothness *= 20 ;
    }
    
    gcamComputeMetricProperties(gcam) ;
    //  GCAMremoveNegativeNodes(gcam, mri_inputs, &parms) ;
    
    GCAMregister(gcam, mri_inputs, &parms) ;
    //  printf("registration complete, removing remaining folds if any exist\n") ;
    //  GCAMremoveNegativeNodes(gcam, mri_inputs, &parms) ;
    
    if (renormalize_align_after)
    {
        int old_diag ;
        TRANSFORM  _transform, *transform = &_transform ;
        
        transform->type = MORPH_3D_TYPE ;
        transform->xform = (void *)gcam ;
        old_diag = Gdiag ;
        if (parms.write_iterations == 0)
        {
            Gdiag &= ~DIAG_WRITE ;
        }
        TransformInvert(transform, mri_inputs) ;
        
        if (Gdiag & DIAG_WRITE)
        {
            char fname[STRLEN] ;
            sprintf(fname, "%s.log", parms.base_name) ;
            parms.log_fp = fopen(fname, "a") ;
        }
        
        // GCA Renormalization with Alignment:
        // check whether or not this is a sequential call
        if (!got_scales)
            // this is the first (and also last) call
            // do not bother passing or receiving scales info
        {
            GCAmapRenormalizeWithAlignment
            (gcam->gca,
             mri_inputs,
             transform,
             parms.log_fp,
             parms.base_name,
             NULL,
             0) ;
            
            if (parms.write_iterations != 0)
            {
                char fname[STRLEN] ;
                MRI  *mri_gca, *mri_tmp ;
                if (parms.diag_morph_from_atlas )
                {
                    sprintf(fname, "%s_target", parms.base_name) ;
                    MRIwriteImageViews(mri_inputs, fname, IMAGE_SIZE) ;
                    sprintf(fname, "%s_target.mgz", parms.base_name) ;
                    printf("writing target volume to %s...\n", fname) ;
                    MRIwrite(mri_inputs, fname) ;
                }
                else
                {
                    //          mri_gca = MRIclone(mri_inputs, NULL) ;
                    mri_gca = MRIalloc(gcam->atlas.width, gcam->atlas.height, gcam->atlas.depth, MRI_FLOAT) ;
                    MRIcopyHeader(mri_inputs, mri_gca) ;
                    GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
                    if (mri_gca->nframes > 1)
                    {
                        printf("careg: extracting %dth frame\n", mri_gca->nframes-1) ;
                        mri_tmp = MRIcopyFrame(mri_gca, NULL, mri_gca->nframes-1, 0) ;
                        MRIfree(&mri_gca) ;
                        mri_gca = mri_tmp ;
                    }
                    sprintf(fname, "%s_target", parms.base_name) ;
                    MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
                    sprintf(fname, "%s_target1.mgz", parms.base_name) ;
                    printf("writing target volume to %s...\n", fname) ;
                    MRIwrite(mri_gca, fname) ;
                    MRIfree(&mri_gca) ;
                }
            }
        }
        else // this is a sequential call, pass scales..
            GCAseqRenormalizeWithAlignment
            (gcam->gca,
             mri_inputs,
             transform,
             parms.log_fp,
             parms.base_name,
             NULL,
             0,
             label_scales,label_offsets,label_peaks,label_computed) ;
        
        got_scales = 1;
        
        
        Gdiag = old_diag ;
        if (write_gca_fname)
        {
            printf("writing normalized gca to %s...\n", write_gca_fname) ;
            GCAwrite(gcam->gca, write_gca_fname) ;
        }
        if (parms.noneg < 2)
        {
            parms.tol /= 5 ;  // reset parameters to previous level
            parms.l_smoothness /= 20 ;
            
            GCAMregister(gcam, mri_inputs, &parms) ;
            
            printf("********************* ALLOWING NEGATIVE NODES IN DEFORMATION"
                   "********************************\n") ;
            parms.noneg = 0 ;
            parms.tol = 0.25 ;
            parms.orig_dt = 1e-6 ;
            parms.navgs = 256 ;
            
            GCAMregister(gcam, mri_inputs, &parms) ;
        }
    }
    
    
    if (parms.l_label > 0)
    {
        GCAMcomputeMaxPriorLabels(gcam) ;  /* start out with max
                                            prior labels again */
        if (reset)
        {
            GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, ORIGINAL_POSITIONS) ;
            GCAMstoreMetricProperties(gcam) ;
        }
        parms.l_label = 0 ;
        printf("***************** morphing with label term set to 0 "
               "*******************************\n") ;
        GCAMregister(gcam, mri_inputs, &parms) ;
    }
    
    
    //record GCA filename to gcam
    strcpy(gcam->atlas.fname, gca_fname);
    printf("writing output transformation to %s...\n", out_fname) ;
    // ...
    // GCAMwrite is used not MORPH3D
    if (GCAMwrite(gcam, out_fname) != NO_ERROR)
    {
        ErrorExit(Gerror, "%s: GCAMwrite(%s) failed", Progname, out_fname) ;
    }
    
    GCAMfree(&gcam) ;
    if (mri_inputs)
    {
        MRIfree(&mri_inputs) ;
    }
    if (diag_fp)
    {
        fclose(diag_fp) ;
    }
    msec = TimerStop(&start) ;
    seconds = nint((float)msec/1000.0f) ;
    minutes = seconds / 60 ;
    hours = minutes / (60) ;
    minutes = minutes % 60 ;
    seconds = seconds % 60 ;
    printf("mri_ca_register took %d hours, %d minutes and %d seconds.\n",
           hours, minutes, seconds) ;
    
    // Print usage stats to the terminal (and a file is specified)
    PrintRUsage(RUSAGE_SELF, "mri_ca_register ", stdout);
    if(rusage_file) WriteRUsage(RUSAGE_SELF, "", rusage_file);
    
    // Output formatted so it can be easily grepped
    printf("FSRUNTIME@ mri_ca_register %7.4f hours %d threads\n",msec/(1000.0*60.0*60.0),n_omp_threads);
    
#ifdef FS_CUDA
    PrintGPUtimers();
#endif

    
    return 0;
}
