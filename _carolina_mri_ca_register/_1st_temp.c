#include "init.h"
#include "time.h"

double cclock()
{
    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}

int main(int argc, char *argv[]){
    
    int err;
    
    double tot=0., start, delta;
    
    err = init_jr(argc,argv);
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
    /////////////////////////////////////////////////
    // -T transform option
    // transform is loaded at get_opt() with -T using TransformRead()
    // assumed to be vox-to-vox
    if (!transform_loaded)   // wasn't preloaded
    {
        transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    }
    else
        // calculate inverse and cache it
    {
        TransformInvert(transform, mri_inputs) ;
    }
    
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
    /////////////////////////////////////////////////
    // -novar option  (default novar = 1)
    if (novar)
    {
        GCAunifyVariance(gca) ;
    }
    
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
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
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
    //////////////////////////////////////////////////////////
    // -regularize val option (default = 0)
    if (regularize > 0)
    {
        GCAregularizeCovariance(gca, regularize) ;
    }
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
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
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
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
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
    //////////////////////////////////////////////////////////////////
    // GCM initialization
    if (!xform_name)  // only if the transform wasn't previously created
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
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
    if (regularize_mean > 0)
    {
        GCAregularizeConditionalDensities(gca, regularize_mean) ;
    }
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
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
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    
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
    start=clock(); delta = clock()-start; tot += delta; printf("init: %f \n",delta);
    printf("TOT: %f \n",tot);
    
    return 0;
}
