/*
 * This is the Measure Objects code
 *
 * Upon arrival in this module, an OBJC is expected to have at least one
 * non-NULL OBJECT1. All coordinates in the OBJECT1s are taken to be in
 * their band's native coordinate system
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "dervish.h"
#include "phSignals.h"
#include "phExtract.h"
#include "phPeaks.h"
#include "phMathUtils.h"
#include "phMeasureObj.h"
#include "phCellFitobj.h"
#include "phOffset.h"
#include "phCosmic.h"
#include "phUtils.h"

#define SAVE_PETRO 1			/* save petro ratio to TEST_INFO? */
#if SAVE_PETRO
#  include "phTestInfo.h"
#endif

/*****************************************************************************/
/*
 * Statics in this file
 */
static void setup_geometry(void);

static void update_fieldstat(OBJC *objc,
			     const FIELDPARAMS *fiparams, FIELDSTAT *fstat);

static void prepare_object(OBJC *objc, const FIELDPARAMS *fiparams);
static void do_measure_objc(OBJC *objc, const FIELDPARAMS *fparams,
			    int fit_models);

static double calc_psf_counts(const REGION *reg, int c, float bkgd);

static void create_object1s(OBJC *objc, const FIELDPARAMS *fparams);

static void setup_aperture_3(const FIELDPARAMS *fiparams);

static void set_region(OBJECT1 *obj1, const REGION *reg,
		       OBJMASK *master_mask, float dcol, float drow);
static void set_sky(OBJECT1 *obj1, int color, const FIELDPARAMS *fiparams);

static void calc_aperture_3(OBJC *objc,
			    int color,
			    const CELL_STATS *prof,
			    const FIELDPARAMS *fiparams);
static void calc_iso_ellipse(OBJC *objc, int c, const CELL_STATS *prof,
			     const FIELDPARAMS *fiparams);
static float calc_photom_sigma(float cts, float neff, const OBJECT1 *obj1,
			       const FRAMEPARAMS *fparams);
static void calc_model_counts(OBJC *objc, int color, const CELL_STATS *cstats,
			      const FIELDPARAMS *fiparams);
static void galaxy_ap_corrections(OBJC *objc, int color,
				  const FIELDPARAMS *fiparams);
static SPLINE *calc_petrosian(OBJC *objc,
			      int color,
			      int ngood,
			      const FIELDPARAMS *fiparams);
static void calc_texture(OBJC *objc, int color, const CELL_STATS *cstats,
			 const FIELDPARAMS *fiparams);
#define DO_PHIL 1
#if DO_PHIL
static void calc_adaptive_moments(OBJC *objc,
			          int color,
			          const FIELDPARAMS *fparams);
#endif
static int set_profiles(OBJECT1 *obj1,
			int color,
			const CELL_STATS *prof,
			const FIELDPARAMS *fiparams);
static void trim_profiles(OBJC *objc, int color, int ngood, float rmax);
static void calc_shape(OBJC *objc,
		       int color,
		       const CELL_STATS *prof,
		       const FIELDPARAMS *fiparams,
		       const SPLINE *cumul_sp);
static void classify_obj1(OBJC *objc, const FIELDPARAMS *fiparams, int color);
static void classify_objc(OBJC *objc, const FIELDPARAMS *fiparams);

static void find_fracPSF(OBJC *objc, int color);
#define FIT_PSF_EXPANSION 0
#if FIT_PSF_EXPANSION
static void fit_psf_expansion(OBJC *objc, int c, const CELL_STATS *prof,
			      const FIELDPARAMS *fiparams);
#endif
static void get_galaxy_ap_correction(FIELDPARAMS *fiparams);

/*****************************************************************************/
/*
 * Note that r[], as_r[], and area[] are set by p_phSetParamsFrame()
 */
const float Gamma = 2.5;		/* control phFindSplineTaut */
static float r[NANN + 1], area[NANN];	/* inner radii and areas of annuli */
static float as_r[NANN + 1];		/* asinh(r[]) */
static float petro_f5_pix = -100;	/* petro_f5 measured in pixels */
static RANDOM *randoms;			/* random numbers */
/*
 * see FRAMEPARAMS for definition of nann_ap_frame
 *
 * ap_corr_run[] is the aperture correction that needs to be
 * applied to the flux within nann_ap_frame to correct
 * it to an aperture containing nann_ap_run annuli, large enough that
 * it's unaffected by seeing variations during the entire run
 * (or survey, for that matter)
 */
static int nann_ap_frame[NCOLOR] = {-1, -1, -1, -1, -1};
static int nann_ap_run[NCOLOR] = {-1, -1, -1, -1, -1};
static float ap_corr_run[NCOLOR] = {-1.0, -1.0, -1.0, -1.0, -1.0};
/*
 * use median not mean profile?
 */
static int median_profs = 0;
/*
 * coefficients for PSF magnitudes
 */
#define NCOEFF_PSF 12
static float coeffs_psf_ss[NCOLOR][NCOEFF_PSF][NCOEFF_PSF];
static float *coeffs_psf_s[NCOLOR][NCOEFF_PSF];
static float **coeffs_psf[NCOLOR] = { NULL, NULL, NULL, NULL, NULL };
					/* values of PSF in a sector */
static float neff_psf[NCOLOR];	/* noise-effective number of pixels */
/*
 * The reconstructed PSF at the current point, and the same PSF sinc centred
 * (as calculated by phProfileExtract)
 */
static const PSF_REG *KLPsf[NCOLOR] = { NULL, NULL, NULL, NULL, NULL };
static const REGION *KLPsfReg[NCOLOR] = { NULL, NULL, NULL, NULL, NULL };
/*
 * Spatially variable aperture corrections within a frame
 */
ACOEFF *deV_ap_correction[NCOLOR] = { NULL, NULL, NULL, NULL, NULL };
ACOEFF *deV_model_ap_correction[NCOLOR] = { NULL, NULL, NULL, NULL, NULL };
ACOEFF *exp_ap_correction[NCOLOR] = { NULL, NULL, NULL, NULL, NULL };
ACOEFF *exp_model_ap_correction[NCOLOR] = { NULL, NULL, NULL, NULL, NULL };

/*****************************************************************************/
/*
 * init and fini functions for this file. The medians structure is used
 * to calculate median colours/Stokes parameters for each field, primarily
 * for QA
 */
static struct {
   int ncolor;				/* number of colours being measured */
   int size;				/* size of the colors[] arrays */
   int nobj;				/* number of objects in colors[] */
   float *fiberColors[NCOLOR];		/* fibre colours in band i, wrt i+1, so
					   if filters are u g r i z the colours
					   will be u-g, g-r, r-i, i-z, 0 */
   float *psfColors[NCOLOR];		/* PSF colours in band i, wrt i+1 */
   float *Q[NCOLOR], *U[NCOLOR];	/* Q and U in band i */
} medians = { 0, 0, 0,
		{ NULL, NULL, NULL, NULL, NULL},
		{ NULL, NULL, NULL, NULL, NULL},
		{ NULL, NULL, NULL, NULL, NULL},
		{ NULL, NULL, NULL, NULL, NULL},
	   };

void
p_phInitMeasureObj(int ncolor,		/* number of colours in use */
		   RANDOM *rand)	/* random numbers */
{
   int i;

   randoms = rand;			/* transfer to global */

   if(medians.size == 0) {
      medians.ncolor = ncolor;
      medians.size = 1000;
      medians.nobj = 0;
      for(i = 0;i < ncolor;i++) {
	 medians.fiberColors[i] = shMalloc(medians.size*sizeof(float));
	 medians.psfColors[i] = shMalloc(medians.size*sizeof(float));
	 medians.Q[i] = shMalloc(medians.size*sizeof(float));
	 medians.U[i] = shMalloc(medians.size*sizeof(float));
      }
   }

   p_phInitEllipseFit();
}

void
p_phFiniMeasureObj(void)
{
   int i;

   randoms = NULL;
   
   if(medians.size > 0) {
      for(i = 0;i < medians.ncolor;i++) {
	 shFree(medians.fiberColors[i]); medians.fiberColors[i] = NULL;
	 shFree(medians.psfColors[i]); medians.psfColors[i] = NULL;
	 shFree(medians.Q[i]); medians.Q[i] = NULL;
	 shFree(medians.U[i]); medians.U[i] = NULL;
      }
      medians.size = medians.nobj = 0;
   }

   p_phFiniEllipseFit();
}

/*
 * Set measure objects up for a complete run
 * 
 * Also check that the FIELDPARAMS are as we expect; return 1 if there's
 * a problem
 */
int
p_phSetParamsRun(const FIELDPARAMS *fiparams)
{
   int color;
   shAssert(fiparams != NULL);
   
   if(fabs(fiparams->fiber_rad - 1.5/fiparams->pixscale) > 1e-5) {
      shError("Expected fiparams->fiber_rad == %g, saw %g",
				  1.5/fiparams->pixscale, fiparams->fiber_rad);
      return(1);
   }
   setup_aperture_3(fiparams);
   petro_f5_pix = fiparams->petro_f5/fiparams->pixscale;
   median_profs = fiparams->median_profs;
/*
 * We really set this in p_phSetParamsFrame(), but for the PSP's convenience
 * we'll set it here too.
 */
   for(color = 0; color < fiparams->ncolor; color++) {
      ap_corr_run[color] = fiparams->frame[color].ap_corr_run;
      nann_ap_run[color] = fiparams->frame[color].nann_ap_run;
      nann_ap_frame[color] = fiparams->frame[color].nann_ap_frame;
   }

   shAssert(medians.size > 0);
   medians.nobj = 0;

   return(0);
}

/*****************************************************************************/
/*
 * First a couple of hyperbolic sine functions. asinh() isn't in the 
 * standard, and sinh(100) is a problem on an alpha...
 */
static float
asinh_ph(float x)
{
   const float amax = 35.0;		/* ~ largest number that we can sinh */

   if(fabs(x) > sinh(amax)) {
      return(x > 0 ? amax : -amax);
   }
/*
 * this is equivalent to a simple
 *   return(log(x + sqrt(1 + x*x)));
 * for all values of x, but is far more numerically stable
 */
   if(x < 0) {
      return(-log(-x + sqrt(1 + x*x)));
   } else {
      return(log(x + sqrt(1 + x*x)));
   }
}

static float
sinh_ph(double x)
{
   const float amax = 35.0;		/* ~ largest number that we can sinh */

   if(fabs(x) > amax) {
      return(sinh(amax));
   } else {
      return(sinh(x));
   }
}

/*****************************************************************************/
/*
 * Here's the setup for a given frame
 */
float
p_phSetParamsFrame(int color,
		   FIELDPARAMS *fiparams,
		   FIELDSTAT *fieldstat)
{
   float neff;				/* effective number of pixels */

   shAssert(color >= 0 && color < NCOLOR);
   shAssert(fiparams != NULL && fiparams->frame[color].psf != NULL);

   setup_geometry();
   
   neff = phPsfCountsSetupFromSincRegion(color, fiparams->frame[color].psf);

   fieldstat->neff_psf[color] = neff_psf[color];
   
   medians.nobj = 0;
/*
 * Unpack the correction needed to take the flux within nann_ap_frame
 * to the flux within nann_ap_run; this is calculated for us by the
 * PSP once per frame
 */
   ap_corr_run[color] = fiparams->frame[color].ap_corr_run;
   nann_ap_run[color] = fiparams->frame[color].nann_ap_run;
   nann_ap_frame[color] = fiparams->frame[color].nann_ap_frame;

   shAssert(nann_ap_run[color] < 100);	/* shush compiler; value isn't used */
/*
 * calculate the model magnitude aperture correction for a star
 * at the centre of the frame if so requested.
 *
 * Note that the canonical band _must_ be set last, as this routine
 * uses the PSF_BASIS for all the bands
 *
 * The fieldstat values are -ve if the aperture correction isn't actually
 * applied to model magnitudes
 */
   if(color == fiparams->ref_band_index) {
      int i;
      
      get_galaxy_ap_correction(fiparams);
      
      for(i = 0; i < fiparams->ncolor; i++) {
	 fieldstat->deV_ap_correction[i] =
				     fiparams->frame[i].deV_ap_correction;
	 fieldstat->deV_ap_correctionErr[i] =
				     fiparams->frame[i].deV_ap_correctionErr;
	 fieldstat->exp_ap_correction[i] =
				     fiparams->frame[i].exp_ap_correction;
	 fieldstat->exp_ap_correctionErr[i] =
				     fiparams->frame[i].exp_ap_correctionErr;
	 
	 if(!fiparams->use_galaxy_ap_correction) {
	    fieldstat->deV_ap_correction[i] *= -1;
	    fieldstat->deV_ap_correctionErr[i] *= -1;
	    fieldstat->exp_ap_correction[i] *= -1;
	    fieldstat->exp_ap_correctionErr[i] *= -1;
	 }
      }
   }

   return(neff);
}

/*
 * And here's the unsetup for a given frame
 */
void
p_phUnsetParamsFrame(int c,
		     FIELDPARAMS *fiparams) /* NOTUSED */
{
   phAcoeffDel(deV_ap_correction[c]); deV_ap_correction[c] = NULL;
   phAcoeffDel(deV_model_ap_correction[c]); deV_model_ap_correction[c] = NULL;
   phAcoeffDel(exp_ap_correction[c]); exp_ap_correction[c] = NULL;
   phAcoeffDel(exp_model_ap_correction[c]); exp_model_ap_correction[c] = NULL;

   phPsfRegDel((PSF_REG *)KLPsf[c]); KLPsf[c] = NULL;
   shRegDel((REGION *)KLPsfReg[c]); KLPsfReg[c] = NULL;
}

/*****************************************************************************/
/*
 * Here's the setup for cells' geometry
 */
static void
setup_geometry(void)
{
   static int first = 1;		/* is this the first call? */

/*
 * Note that area[i] is the area of the i'th annulus and r[i] is the
 * value of the INNER radius (so area[0] = pi*r[1]^2).
 *
 * The cumulative profile cumul[i] is the cumululative flux to r[i], so
 * cumul[0] = 0, and cumul[1] == I0*area[0]
 */
   if(first) {
      const CELL_STATS *cstats = phProfileGeometry();
      int i;
      
      for(i = 0;i <= NANN;i++) {
	 r[i] = cstats->radii[i];
	 as_r[i] = asinh_ph(r[i]);
	 if(i < NANN) {
	    area[i] = cstats->area[i];
	 }
      }
      
      first = 0;
   }

   return;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * phMeasureObjc() takes an OBJC and calculates lots of parameters for
 * the embedded OBJECT1s.
 *
 * We increment the frame's FIELDSTAT structure's fields as appropriate 
 * for the given OBJC; it's less efficient than going through a complete CHAIN
 * of all measured and classified objects at once, but we don't have that
 * luxury with this function.
 */
RET_CODE
phMeasureObjc(OBJC *objc,		/* the OBJC to measure */
	      const FIELDPARAMS *fparams, /* and how to measure it */
	      FIELDSTAT *fieldstat,	/* and a summary of what we found */
	      int remove_atlas_images,	/* replace atlas images by noise? */
	      int bright,		/* are we measuring bright objects? */
	      int fit_models)		/* should we fit models to objects? */
{
   int c, i;
   OBJC *child;
   OBJECT1 *obj1;
   int set_child_atlas_images = 0;	/* did we insert childrens' atlas
					   images into data? */
   int union_flag;			/* Rule Brittania etc. But really the
					   union of all the colour's flags */

   CATCH("phMeasureObjc", SIGABRT, SH_GENERIC_ERROR);

   shAssert(fparams != NULL);
   shAssert(objc != NULL && objc->aimage->master_mask != NULL);
   
   union_flag = 0;				/* OR of flags in all bands */
   for(c = 0;c < objc->ncolor;c++) {
      shAssert(fparams->frame[c].data != NULL);
      shAssert(fparams->frame[c].toGCC != NULL);
      if(objc->color[c] != NULL) {
	 union_flag |= objc->color[c]->flags;
      }
   }
/*
 * If the OBJC was detected as a bright object it'll have acquired
 * atlas images and regions with associated masks. These should all
 * have been deleted, but let's make sure.
 * 
 * It will only have been labelled OBJECT3_MEASURED if it was
 * measured, and _that_ will only have happened if OBJECT3_MEASURE_BRIGHT
 * was set. In fact, if the object wasn't measured it shouldn't be
 * treated as bright (it was probably only detected for the sake of its
 * position, to peak up the astrometry), so turn off OBJECT1_BRIGHT too
 *
 * Furthermore, the deblender should _not_ have been run, so check that too
 */
   shAssert(objc->nchild == 0);
   if(!bright && (union_flag & OBJECT1_BRIGHT)) {
      if(objc->flags3 & OBJECT3_MEASURE_BRIGHT) {
	 shAssert(objc->flags3 & OBJECT3_MEASURED); /* already measured */
      } else {
 	 objc->flags &= ~OBJECT1_BRIGHT;
	 objc->flags3 &= ~OBJECT3_MEASURE_BRIGHT;
	 for(c = 0;c < objc->ncolor;c++) {
	    if((obj1 = objc->color[c]) != NULL) {
	       obj1->flags &= ~OBJECT1_BRIGHT;
	       obj1->flags3 &= ~OBJECT3_MEASURE_BRIGHT;
	    }
	 }
      }
   }
   for(c = 0;c < objc->ncolor;c++) {
      if((obj1 = objc->color[c]) != NULL) {
	 shAssert(obj1->region == NULL);
      }
   }
/*
 * If the object's already been measured we want to remeasure it now that
 * we've subtracted the sky and found faint objects. Keep the original
 * (bright) measurements, and make a new OBJC that's the sibling of the
 * bright one. Note that the bright one is also its parent, which is
 * OK in photo's world, even if frowned upon elsewhere
 */
   if(objc->flags & OBJECT1_BRIGHT) {
      shAssert(!bright);		/* BRIGHT was set in a previous call */

      shAssert(objc->sibbs == NULL);
      objc->sibbs = phObjcNewFromObjc(objc, 1);
      objc->nchild++;
      objc = objc->sibbs;

      objc->flags &= ~(OBJECT1_BRIGHT | OBJECT1_NODEBLEND);
      objc->type = OBJ_UNK;

      for(c = 0; c < objc->ncolor; c++) { /* clear old profile */
	 obj1 = objc->color[c];
	 if(obj1 != NULL) {
	    obj1->flags &= ~(OBJECT1_BRIGHT | OBJECT1_NODEBLEND);
	    obj1->type = OBJ_UNK;
	    obj1->nprof = 0;	
	    for(i = 0; i < NANN; i++) {
	       obj1->profMean[i] = -1000.0;
	       obj1->profMed[i] = -1000.0;
	       obj1->profErr[i] = -1000.0;
	    }
	 }
      }
   } else {
      shAssert(!(objc->flags3 & OBJECT3_MEASURED));
   }
/*
 * restore atlas images, or cut them if they don't exist. We don't want to
 * restore bright atlas images, of course, as they date from before we
 * subtracted the sky
 */
   for(c = 0;c < objc->ncolor;c++) {
      if(objc->aimage->pix[c] == NULL) {
	 phAtlasImageCut(objc, -1, fparams, -1, -1, NULL);
      } else {
	 if(!(objc->flags & OBJECT1_BRIGHT)) {
	    phRegionSetFromAtlasImage(objc->aimage, c,
				      (REGION *)fparams->frame[c].data, 0, 0);
	 }
      }
   }
/*
 * Set the PSF for this part of the frame
 * The error term fparams->frame[c].psf_app_correctionErr comes from the PSP
 * as CALIB1.psf_ap_correctionErr
 */
   for(c = 0;c < objc->ncolor;c++) {
      fparams->frame[c].psf_app_correction =
	phPsfSetAtPoint(c, fparams->frame[c].psfBasis, objc->rowc, objc->colc,
			NULL);
   }   
/*
 * See if we need to deblend; note that we don't deblend while measuring
 * bright objects as we don't yet have a complete peaks list (of course,
 * we'll see those objects again in the faint pass, and deblend them then)
 */
   prepare_object(objc, fparams);

   if(bright && (objc->flags & OBJECT1_BLENDED)) {
     objc->flags |= OBJECT1_NODEBLEND;
   }

   if(!bright && objc->peaks != NULL && objc->peaks->npeak > 1) {
      shAssert(objc->flags & OBJECT1_BLENDED);
/*
 * create children
 *
 * give up if there are too many children (or rather, potential
 * children as some may be discarded during template construction)
 */
      if(phObjcMakeChildren(objc, fparams) == 0) {
	 objc->flags |= OBJECT1_NODEBLEND;
      } else {
	 set_child_atlas_images = 1;
	 
	 if(phObjcDeblend(objc, fparams) < 0) {
	    phObjcChildrenDel(objc);
	 } else {
/*
 * We need to find the brightest galaxy child in each band. Why? Because
 * target selection isn't allowed to look at families. PR 1395
 */
	    OBJECT1 *brightest_galaxy_child[NCOLOR];
	    
	    for(c = 0; c < objc->ncolor; c++) {
	       brightest_galaxy_child[c] = NULL;
	    }
/*
 * measure those children
 */
	    (void)phObjcDescendentNext(objc); /* returns objc */
	    while((child = phObjcDescendentNext(NULL)) != NULL) {
	       phInsertAtlasImage(child, fparams);

	       prepare_object(child, fparams);
	       do_measure_objc(child, fparams, fit_models);
/*
 * If it's still labelled blended, and it isn't moving, this implies
 * that we gave up on further deblending
 */
	       if(child->flags & OBJECT1_BLENDED) {
		  if(child->flags2 & OBJECT2_DEBLENDED_AS_MOVING) {
		     child->flags &= ~OBJECT1_BLENDED;
		  } else {
		     child->flags |= OBJECT1_NODEBLEND;
		  }
	       }
/*
 * Is this the brightest galaxy child? Determine this independently in
 * each band, and don't set it for the OBJC (what flux would I use?)
 */
	       for(c = 0;c < objc->ncolor;c++) {
		  obj1 = child->color[c];
		  if(obj1->type == OBJ_GALAXY) {
		     if(brightest_galaxy_child[c] == NULL ||
			obj1->petroCounts >
				      brightest_galaxy_child[c]->petroCounts) {
			brightest_galaxy_child[c] = obj1;
		     }
		  }
	       }
	    }
/*
 * set flag indicating which child galaxy was the brightest
 */
	    for(c = 0;c < objc->ncolor;c++) {
	       if(brightest_galaxy_child[c] != NULL) {
		  brightest_galaxy_child[c]->flags2 |=
						OBJECT2_BRIGHTEST_GALAXY_CHILD;
	       }
	    }
	 }
      }
   }
/*
 * Measure the properties of the primary object.
 *
 * If we ran the deblender we must copy the parent's pixels into the region
 */
   if(set_child_atlas_images) {
      for(c = 0;c < objc->ncolor;c++) {
	 if(!(objc->flags & OBJECT1_BRIGHT)) {
	    phRegionSetFromAtlasImage(objc->aimage, c,
				      (REGION *)fparams->frame[c].data, 0, 0);
	 }
      }
   }
   do_measure_objc(objc, fparams, fit_models);

   if(bright) {
      objc->flags |= OBJECT1_BRIGHT;
   }
/*
 * Measure all those children's fiber magnitudes, as we want the flux
 * in the _parent_ (PR 1031)
 */
   (void)phObjcDescendentNext(objc); /* returns objc */
   while((child = phObjcDescendentNext(NULL)) != NULL) {
      for(c = 0;c < objc->ncolor;c++) {
	 calc_aperture_3(child, c, NULL, fparams);
      }
   }
/*
 * remove the atlas image from the data frames if we are so requested
 */
   if(remove_atlas_images) {
      phRemoveAtlasImage(objc, fparams);
   }
   
   update_fieldstat(objc, fparams, fieldstat);

   END_CATCH(SIGABRT, SH_GENERIC_ERROR);

   return (SH_SUCCESS);
}

/*****************************************************************************/
/*
 * prepare an object to be measured; find the canonical center,
 * create missing OBJECT1s, and set the REGIONs
 */
static void 
prepare_object(OBJC *objc,
	       const FIELDPARAMS *fiparams)
{
   int c;
   float drow, dcol;			/* offsets from reference colour */
   OBJECT1 *obj1;
   
   phDeblendedObjcRecenter(objc, fiparams);
   create_object1s(objc, fiparams);

   for(c = 0; c < objc->ncolor; c++) {
      obj1 = objc->color[c];
      
      phOffsetDo(fiparams, objc->rowc, objc->colc, 
		    fiparams->ref_band_index, c,
		 NULL, NULL, &drow, NULL, &dcol, NULL);

      if(obj1->flags & OBJECT1_CHILD) {
	 shAssert(obj1->region != NULL);
      } else {
	 shAssert(obj1->region == NULL);

	 set_region(obj1, fiparams->frame[c].data,
					objc->aimage->master_mask, dcol, drow);
      }
   }
}

/*****************************************************************************/
/*
 * Here's the real guts of phMeasureObjc, split out so at be callable for
 * both parents and children
 */
static void
do_measure_objc(OBJC *objc,		/* the object to measure */
		const FIELDPARAMS *fparams, /* gain etc. */
		int fit_models)		/* should we fit models to objects? */
{
   int i,c;
   int cindex[NCOLOR];			/* color indices, in order */
   SPLINE *cumul_sp;			/* spline of asinh(cumulative profile)
					   wrt asinh(r) */
   float drow, dcol;			/* offsets from reference colour */
   const REGION *data;			/* == fparams->frame[c].data */
   static int nann_L = 0;		/* number of annuli for model fitting*/
   int nann_old;			/* number of annuli measured by
					   previous call to phProfileExtract */
   int ngood = 0;			/* number of "good" radial points */
   OBJECT1 *obj1;
   CELL_STATS *prof;			/* extracted profile */
   int rad;				/* radius for profile extraction */
   const SPANMASK *sm;			/* == (SPANMASK *)obj1->region->mask */
/*
 * Set the array cindex[] to process the canonical colour first, then the
 * others in their usual order.
 *
 * If the object's detected in the fiparams->ref_band_index band that'll be
 * canonical; but this doesn't always happen
 * 
 */
   i = 0; cindex[i] = -1;
   for(c = 0;c < objc->ncolor;c++) {
      if(objc->color[c]->flags2 & OBJECT2_CANONICAL_BAND) {
	 shAssert(i == 0);		/* check for only one canonical band */
	 cindex[i++] = fparams->ref_band_index; /* do this band first */
      }
   }
   shAssert(cindex[0] >= 0);		/* there _is_ a canonical band */

   for(c = 0;c < objc->ncolor;c++) {
      if(c != cindex[0]) {
	 cindex[i++] = c;
      }
   }
/*
 * Start measuring each band
 */
   for(i = 0;i < objc->ncolor;i++) {
      static char *dump_filename = NULL; /* write data to this file?
					    For use from gdb */
      const float size = 3;		/* size of area to search for interp*/

      c = cindex[i];
      obj1 = objc->color[c];
      data = fparams->frame[c].data;

      if(dump_filename != NULL) {
	 shRegWriteAsFits((REGION *)data,
			  dump_filename, STANDARD, 2, DEF_NONE, NULL, 0);
	 dump_filename = NULL;
      }
      
      phFitCellColorSet(c, NULL);

      shAssert(obj1 != NULL && data != NULL);
      obj1->flags3 |= OBJECT3_MEASURED;
            
      phOffsetDo(fparams, objc->rowc, objc->colc, 
		 fparams->ref_band_index, c,
		 NULL, NULL, &drow, NULL, &dcol, NULL);

      set_sky(obj1, c, fparams);
/*
 * If the object's detected, look to see if its centre's been interpolated;
 * if it has also check if it's saturated (satur => interp).
 *
 * Otherwise set CR, SATUR, and INTERP flags if necessary; for objects that
 * were detected this has already been done, but if it wasn't detected
 * in its own right, we only now know the appropriate area to check.
 *
 * Note that we lie about the master mask's origin in order to achieve this
 */
      sm = (SPANMASK *)obj1->region->mask;
      shAssert(sm != NULL && sm->cookie == SPAN_COOKIE);
      
      if(obj1->flags & OBJECT1_DETECTED) {
	 if(obj1->flags & OBJECT1_INTERP) { /* see if centre's interpolated */
	    const int c0 = obj1->colc - size/2 + 0.5;
	    const int r0 = obj1->rowc - size/2 + 0.5;
	    
	    if(phRectIntersectMask(sm->masks[S_MASK_INTERP],
				   c0, r0, c0 + size + 0.5, r0 + size + 0.5)) {
	       obj1->flags2 |= OBJECT2_INTERP_CENTER;

	       if(phRectIntersectMask(sm->masks[S_MASK_SATUR],
				   c0, r0, c0 + size + 0.5, r0 + size + 0.5)) {
		  obj1->flags2 |= OBJECT2_SATUR_CENTER;
	       }
	    }
	 }
      } else {
	 const int idrow = (drow < 0) ? -(-drow + 0.5) : drow + 0.5;
	 const int idcol = (dcol < 0) ? -(-dcol + 0.5) : dcol + 0.5;
	 
	 objc->aimage->master_mask->row0 += idrow;
	 objc->aimage->master_mask->row0 += idcol;
	 
	 if(!(obj1->flags & OBJECT1_CR) &&
	    phObjmaskIntersectMask(sm->masks[S_MASK_CR],
						  objc->aimage->master_mask)) {
	    obj1->flags |= OBJECT1_CR;
	 }
	 if(!(obj1->flags & OBJECT1_INTERP) &&
	    phObjmaskIntersectMask(sm->masks[S_MASK_INTERP],
						  objc->aimage->master_mask)) {
	    obj1->flags |= OBJECT1_INTERP;
	 }
	 if(!(obj1->flags & OBJECT1_SATUR) &&
	    phObjmaskIntersectMask(sm->masks[S_MASK_SATUR],
						  objc->aimage->master_mask)) {
	    obj1->flags |= OBJECT1_SATUR;
	 }
	 
	 objc->aimage->master_mask->row0 -= idrow;
	 objc->aimage->master_mask->row0 -= idcol;
      }
/*
 * Does center lie in NOTCHECKED region?  Interesting even for non detections
 * if e.g. half of a CCD isn't working
 */
      {
	 const int c0 = obj1->colc - size/2 + 0.5;
	 const int r0 = obj1->rowc - size/2 + 0.5;

	 if(phRectIntersectMask(sm->masks[S_MASK_NOTCHECKED],
				c0, r0, c0 + size + 0.5, r0 + size + 0.5)) {
	    obj1->flags2 |= OBJECT2_NOTCHECKED_CENTER;
	 }
      }
/*
 * Extract the radial profile
 * |rad| is an index into the array extract.c:anndex[]
 */
      if(c == cindex[0]) {
	 rad = (obj1->flags & OBJECT1_BRIGHT) ? 8 : 6;
      } else {
	 rad = objc->color[cindex[0]]->nprof;
      }
      nann_old = 0;
      for(;;) {
	 prof = phProfileExtract(obj1->id, data,
				 objc->rowc + drow, objc->colc + dcol, -rad,
				 fparams->frame[c].bkgd + SOFT_BIAS,
				 obj1->skyErr, 0);
	 if(prof == NULL || prof->syncreg == NULL) { /* too close to edge */
	    prof = NULL;
	    obj1->flags |= OBJECT1_EDGE | OBJECT1_NOPETRO | OBJECT1_NOPROFILE;
	    obj1->petroRad = petro_f5_pix;

	    break;
	 }

	 if(prof->nannuli == nann_old) {	/* we didn't get more points */
	    break;
	 }
	 nann_old = prof->nannuli;	 

	 ngood = set_profiles(obj1, c, prof, fparams);

	 if(prof->nannuli > 0 && obj1->profMed[0] == MAX_U16) {
	    obj1->flags |= OBJECT1_SATUR; /* this overflow could have happened
					     while sinc shifting */
	 }
	 
	 if(prof->nannuli == NANN) {	/* object is too large */
	    obj1->flags |= OBJECT1_TOO_LARGE;
	    break;
	 } else if(ngood < prof->nannuli) {
	    break;
	 }
	 rad++;				/* go out further */
      }

      if(prof == NULL) {
	 continue;			/* try next colour */
      }
      if(ngood < 6) {
	 ngood = 6;			/* minimum size of profile */
	 shAssert(ngood <= prof->nannuli);
      }
/*
 * If the central intensity is very negative we have got the sky level
 * badly wrong; this'll lead to NaNs and such like; so give up now.
 */
      {
	 float cen_val, noise;
	 cen_val = phProfileMedian(prof, 0, 0, 0, NULL);
	 noise = calc_photom_sigma(cen_val,1,obj1,&fparams->frame[c]);

	 if(noise != noise ||		/* i.e. NaN */
	    cen_val < -100*noise) {
	    obj1->flags |= OBJECT1_BADSKY;

	    obj1->petroRad = petro_f5_pix;
	    obj1->flags |= OBJECT1_NOPETRO;

	    continue;
	 }
      }
/*
 * Set parameters in OBJECT1, measured from the canonical centre
 */
      calc_aperture_3(objc, c, prof, fparams);
      cumul_sp = calc_petrosian(objc, c, ngood, fparams);
/*
 * does some part of the object within the r' Petrosian radius
 * lie off the frame?
 */
      if(r[prof->nannuli_c] < objc->color[cindex[0]]->petroRad) {
	 obj1->flags |= OBJECT1_INCOMPLETE_PROFILE;
      }

      if(obj1->flags & OBJECT1_NOPROFILE) {
	 obj1->flags |= OBJECT1_ELLIPFAINT | OBJECT1_NOSTOKES;
      }

      if(!(obj1->flags & OBJECT1_ELLIPFAINT)) {
	 calc_iso_ellipse(objc, c, prof, fparams);
      }
      if(!(obj1->flags & OBJECT1_NOSTOKES)) {
	 calc_shape(objc, c, prof, fparams, cumul_sp);
      }
      phSplineDel(cumul_sp);
/*
 * If this is the canonical band trim the profile to the central region,
 * both for consistency with the other bands, and because it seems to give
 * better values for the likelihoods
 */
      if(c == cindex[0]) {		/* trim to sinc region */
#define USE_SINC_FOR_MODELS 0
#if USE_SINC_FOR_MODELS
	 prof = phProfileExtract(obj1->id, data,
				 obj1->rowc, obj1->colc, SYNC_REG_SIZE/2,
				 fparams->frame[c].bkgd + SOFT_BIAS,
				 obj1->skyErr, 0);
	 if(prof == NULL || prof->syncreg == NULL) { /* too close to edge */
	    prof = NULL;
	    obj1->flags |= OBJECT1_EDGE;
	    obj1->flags2 |= OBJECT2_LOCAL_EDGE;

	    continue;			/* proceed to next colour */
	 }
#endif

	 if(fit_models) {
	    phFitCellAsPsf(objc, c, prof, fparams, nann_L, NULL, NULL);
	    phFitCellAsDeV(objc, c, prof, fparams, nann_L);
	    phFitCellAsExp(objc, c, prof, fparams, nann_L);
	    if(fparams->use_galaxy_ap_correction) {
	       galaxy_ap_corrections(objc, c, fparams);
	    }
	 }
      }
#define MODEL_CANONICAL_CENTER 0	/* counts_model are about canon. cen.*/
#if MODEL_CANONICAL_CENTER
      if(fit_models) {
	 calc_model_counts(objc, c, prof, fparams); /* use canonical centre */
      }
#endif
/*
 * now quantities measured about the _local_ centre; we only use the central
 * sinc shifted region for this for reasons of efficiency. Of course, in the
 * canonical band we don't have to re-extract the profile
 */
      if(c != cindex[0]) {		/* re-extract at this band's centre */
#if USE_SINC_FOR_MODELS
	 const int keep_profile = 0;
#else
	 const int keep_profile = 1;
#endif
	 prof = phProfileExtract(obj1->id, data,
				 obj1->rowc,obj1->colc, SYNC_REG_SIZE/2,
				 fparams->frame[c].bkgd + SOFT_BIAS,
				 obj1->skyErr, keep_profile);
	 if(prof == NULL || prof->syncreg == NULL) { /* too close to edge */
	    prof = NULL;
	    obj1->flags |= OBJECT1_EDGE;
	    obj1->flags2 |= OBJECT2_LOCAL_EDGE;

	    continue;			/* proceed to next colour */
	 }

	 if(fit_models) {
	    phFitCellAsPsf(objc, c, prof, fparams, nann_L, NULL, NULL);
	    phFitCellAsDeV(objc, c, prof, fparams, nann_L);
	    phFitCellAsExp(objc, c, prof, fparams, nann_L);
	    if(fparams->use_galaxy_ap_correction) {
	       galaxy_ap_corrections(objc, c, fparams);
	    }
	 }
      }

      phPsfCountsFromSincRegion(objc, c, prof, fparams);
      calc_texture(objc, c, prof, fparams);

#if !MODEL_CANONICAL_CENTER
      if(fit_models) {
	 calc_model_counts(objc, c, prof, fparams); /* use canonical centre */
      }
#endif
/*
 * Phil's AM measures
 */
#if DO_PHIL
      calc_adaptive_moments(objc, c, fparams);
#endif
/*
 * fit the core of the object using the PSF's derivatives
 */
#if FIT_PSF_EXPANSION
      fit_psf_expansion(objc, c, prof, fparams);
#endif
/*
 * finally classification etc.
 */
      classify_obj1(objc, fparams, c);
      find_fracPSF(objc, c);

      trim_profiles(objc, c, ngood,
		    fparams->petro_f2*objc->color[cindex[0]]->petroRad);
   }
/*
 * estimate the velocity
 */
   if(!(objc->flags2 & OBJECT2_DEBLENDED_AS_MOVING)) {
      float chisq, chisq_max = 3.0;	/* XXX */
      
      chisq = phVelocityFind(objc, fparams, NULL, NULL, NULL, NULL);
      
      if(chisq > chisq_max) {
	 objc->flags2 |= OBJECT2_BAD_MOVING_FIT;
      }
   }
/*
 * propagate some flag information to the OBJC, and some back down to the
 * OBJECT1s
 */
   classify_objc(objc, fparams);

   for(c = 0;c < objc->ncolor;c++) {
      objc->color[c]->flags |= objc->flags & OBJECT1_BLENDED;

      objc->flags |=
	objc->color[c]->flags & (OBJECT1_EDGE |
				 OBJECT1_BLENDED |
				 OBJECT1_CHILD |
				 OBJECT1_NOPETRO |
				 OBJECT1_MANYPETRO |
				 OBJECT1_INTERP |
				 OBJECT1_CR |
				 OBJECT1_SATUR |
				 OBJECT1_NOTCHECKED |
				 OBJECT1_BINNED1 |
				 OBJECT1_BINNED2 |
				 OBJECT1_BINNED4);

      objc->flags2 |=
	objc->color[c]->flags2 & (OBJECT2_DEBLENDED_AS_MOVING |
				  OBJECT2_NODEBLEND_MOVING |
				  OBJECT2_TOO_FEW_DETECTIONS |
				  OBJECT2_BAD_MOVING_FIT |
				  OBJECT2_STATIONARY |
				  OBJECT2_PEAKS_TOO_CLOSE |
				  OBJECT2_BAD_MOVING_FIT_CHILD |
				  OBJECT2_DEBLEND_UNASSIGNED_FLUX |
				  OBJECT2_INTERP_CENTER |
				  OBJECT2_SATUR_CENTER |
				  OBJECT2_DEBLENDED_AT_EDGE |
				  OBJECT2_DEBLEND_NOPEAK |
				  OBJECT2_PSF_FLUX_INTERP |
				  OBJECT2_TOO_FEW_GOOD_DETECTIONS |
				  OBJECT2_CENTER_OFF_AIMAGE |
				  OBJECT2_DEBLEND_DEGENERATE |
				  OBJECT2_MAYBE_CR |
				  OBJECT2_MAYBE_EGHOST);
      /* these flags are needed for bookkeeping: */
      objc->flags3 |=
	objc->color[c]->flags3 & (OBJECT3_MEASURED |
				  OBJECT3_GROWN_MERGED | 
				  OBJECT3_HAS_CENTER |
				  OBJECT3_MEASURE_BRIGHT);
   }
}

/*****************************************************************************/
/*
 * Given the position in a frame, set the PSF to the local value, including
 * calculating and returning the proper aperture correction
 *
 * N.b. This routine is not responsible for setting the PSF coefficients
 * used for the entire frame, it merely adjusts things for the local PSF
 */
float
phPsfSetAtPoint(int c,			/* which colour are we processing? */
		const PSF_BASIS *basis,	/* describe spatial variation of PSF */
		float rowc, float colc,	/* desired position in frame */
		float *ap_correction_err) /* error in ap. correction, or NULL*/
{
   float ap_correction;			/* the desired aperture correction */
   double apCounts;			/* aperture counts for PSF */
   CELL_STATS *prof;			/* extracted profile */
   const DGPSF *psf;			/* estimate of the PSF */
   float psfCounts;			/* measured psfCounts for PSF */

   shAssert(c >= 0 && c < NCOLOR);
   
   if(basis == NULL) {
      if(ap_correction_err != NULL) {
	 *ap_correction_err = 0.0;
      }
      return(1.0);
   }
/*
 * We have a PSF_BASIS, so use it
 */
   phPsfRegDel((PSF_REG *)KLPsf[c]);
   KLPsf[c] = phPsfKLReconstruct(basis, rowc, colc, TYPE_U16);
   
   psf = phFitPsfFromReg(KLPsf[c]->reg, 0); /* n.b. we don't own this! */
   phFitCellColorSet(c, &psf->coeffs);
/*
 * Now calculate the aperture correction proper to this PSF.
 * phFitPsfFromReg() has already extracted a profile so we don't need
 * to do so again.  We save the sinc-shifted region in KLPsfReg
 */
   prof = phProfileGetLast();		/* set in phFitPsfFromReg */
   shAssert(prof->syncreg != NULL);
   if(KLPsfReg[c] == NULL) {
      char buff[20]; sprintf(buff, "KLPSFReg[%d]", c);
      KLPsfReg[c] = shRegNew(buff,
			   prof->syncreg->nrow, prof->syncreg->ncol, TYPE_U16);
   }
   shRegIntCopy((REGION *)KLPsfReg[c], prof->syncreg);

   if(coeffs_psf[c] == NULL) {
      ap_correction = 1.0;
   } else {
      psfCounts = calc_psf_counts(prof->syncreg, c, 0.0);
      apCounts = phApertureCounts(prof, nann_ap_frame[c]);
      
      ap_correction = apCounts/psfCounts*ap_corr_run[c];
   }

   if(ap_correction_err != NULL) {
      *ap_correction_err = 0.0;
   }

   return(ap_correction);
}

/*****************************************************************************/
/*
 * Two procedures that are only needed as one of them needs access to the
 * random numbers in randoms
 */
void
phRemoveAtlasImage(const OBJC *objc,
		   const FIELDPARAMS *fiparams)
{
   int c;
   float sigma;				/* background s.d. */

   shAssert(objc != NULL && objc->aimage != NULL && fiparams != NULL);
   
   for(c = 0;c < objc->ncolor;c++) {
      shAssert(fiparams->frame[c].data != NULL);

      if(objc->aimage->pix[c] != NULL) {
	 sigma = sqrt(objc->color[c]->sky/fiparams->frame[c].gain +
					     fiparams->frame[c].dark_variance);
	 phRegionSetValFromAtlasImage(objc->aimage, c,
				      (REGION *)fiparams->frame[c].data,
				      SOFT_BIAS + fiparams->frame[c].bkgd,
				      sigma, randoms, 0, 0);
      }
   }
}

void
phInsertAtlasImage(const OBJC *objc,
		   const FIELDPARAMS *fiparams)
{
   int c;

   shAssert(objc != NULL && objc->aimage != NULL && fiparams != NULL);
   
   for(c = 0; c < objc->ncolor; c++) {
      shAssert(fiparams->frame[c].data != NULL);

      phRegionSetFromAtlasImage(objc->aimage, c,
				(REGION *)fiparams->frame[c].data, 0, 0);
   }
}

/*****************************************************************************/
/*
 * Reconstruct the PSF and calculate PSF and model fluxes.
 *
 * Note Well that this routine calls phPsfSetAtPoint(), which is expensive.
 * It should therefore NOT be called for every object in the frame. If ever
 * we need to do that, some repackaging will be needed.
 *
 * This code is not part of phPsfSetAtPoint() as that routine must be called
 * once for each _object_, processing the canonical band first
 */
static void
psf_get_and_fit(OBJC *objc,		/* the object to measure */
		int c,			/* in this band */
		FIELDPARAMS *fiparams,
		float *counts_model_deV,
		float *counts_model_exp)
{
   const PSF_BASIS *basis;		/* specify PSF in desired band */
   int cc = fiparams->ref_band_index;	/* canonical colour */
   FRAMEPARAMS *const fparams = &fiparams->frame[c];
   OBJECT1 *obj1 = objc->color[c];	/* object in this band */
   CELL_STATS *prof;			/* extracted profile */
   float psf_app_correction;		/* the local PSF aperture correction */
/*
 * Have we already measured the canonical band?
 */
   if(c == cc) {
      obj1->flags2 |= OBJECT2_CANONICAL_BAND;
   } else {
      shAssert(objc->color[cc]->flags2 & OBJECT2_CANONICAL_BAND);
   }
/*
 * reconstruct the PSF for this point in the frame and set the
 * PSF aperture correction correctly
 */
   basis = fparams->psfBasis;
   shAssert(basis != NULL);
/*
 * reconstruct the PSF at each of those places, and fit the PSF and models
 */
   psf_app_correction = phPsfSetAtPoint(c, basis, obj1->rowc,obj1->colc, NULL);
/*
 * the profile was extracted for us by phPsfSetAtPoint(), but we only
 * want to use the part in the sinc region
 */
   prof = phProfileGetLast();
   prof->nannuli = prof->nannuli_c = NSYNCANN;
/*
 * Set such fields as are needed by the model fit code, either as values
 * or initial guesses for model parameters
 */
   obj1->profMean[0] = 0;
   obj1->rowcErr = obj1->colcErr = 0; 
   obj1->sky = phBinregionInterpolate(fparams->sky, objc->rowc, objc->colc);
   obj1->sky += fparams->bkgd;
   obj1->U = obj1->Q = 0.0;
   obj1->petroR50 = 1;		/* ~ 1" seeing */
/*
 * Calculate PSF counts and find best-fit galaxy models, then set the
 * exp and deV aperture corrections
 */
   obj1->psfCounts = psf_app_correction*calc_psf_counts(prof->syncreg, c, 0.0);
   
   (void)phFitCellAsDeV(objc, c, prof, fiparams, 0);
   (void)phFitCellAsExp(objc, c, prof, fiparams, 0);
/*
 * Evaluate model fluxes. Note that we have to have already fit the
 * canonical band (as is asserted above); this is arranged by ensuring
 * that this routine is called with the canonical band first
 */
   objc->color[cc]->deV_L = 0; objc->color[cc]->exp_L = 1;
   calc_model_counts(objc, c, prof, fiparams);
   *counts_model_exp = obj1->counts_model;
      
   objc->color[cc]->deV_L = 1; objc->color[cc]->exp_L = 0;
   calc_model_counts(objc, c, prof, fiparams);
   *counts_model_deV = obj1->counts_model;
}   

/*****************************************************************************/
/*
 * Calculate the aperture correction for a star from the flux in the the
 * best fit deV/exp models to the PSF flux.
 */
#define N_GALAP_POINT 9

static void
get_galaxy_ap_correction(FIELDPARAMS *fiparams) /* describe field */
{
   int c;
   int cindex[NCOLOR];			/* colours , in order of processing */
   float counts_deV[NCOLOR][N_GALAP_POINT]; /* counts of deV */
   float counts_exp[NCOLOR][N_GALAP_POINT]; /* counts of exp */
   float counts_model_deV[NCOLOR][N_GALAP_POINT]; /* counts when deV is best*/
   float counts_model_exp[NCOLOR][N_GALAP_POINT]; /* counts when exp is best*/
   float counts_psf[NCOLOR][N_GALAP_POINT]; /* counts of PSF fits */
   const int cc = fiparams->ref_band_index; /* canonical colour */
   FRAMEPARAMS *fparams;		/* == &fiparams->frame[] */
   int i, j;
   int n;				/* number of places to eval PSF */
   int nrow, ncol;			/* size of fiparams->frame[].data */
   int nterm_row = 3;			/* number of terms; in row and  */
   int nterm_col = 3;			/*     column direction. linear == 2 */
   OBJC *objc = phObjcNew(fiparams->ncolor);
   OBJECT1 *obj1;			/* == objc->color[] */
   float rowc[N_GALAP_POINT];		/* where to */
   float colc[N_GALAP_POINT];		/*          evaluate corrections */
/*
 * Do we know the PSF_BASIS?
 */
   if(fiparams->frame[0].psfBasis == NULL) { /* No we don't */
      for(c = 0; c < objc->ncolor; c++) {
	 fparams = &fiparams->frame[c];

	 shAssert(fparams->psfBasis == NULL);
	 
	 deV_ap_correction[c] = NULL;
	 fparams->deV_ap_correction = 1.0;
	 fparams->deV_ap_correctionErr = 0.0; /* XXX */
	 
	 exp_ap_correction[c] = NULL;
	 fparams->exp_ap_correction = 1.0;
	 fparams->exp_ap_correctionErr = 0.0; /* XXX */
	 
	 exp_model_ap_correction[c] = NULL;
	 fparams->exp_model_ap_correction = 1.0;
	 fparams->exp_model_ap_correctionErr = 0.0; /* XXX */
	 
	 deV_model_ap_correction[c] = NULL;
	 fparams->deV_model_ap_correction = 1.0;
	 fparams->deV_model_ap_correctionErr = 0.0; /* XXX */
      }

      return;
   }
/*
 * Choose positions in the frame to evaluate the corrections
 */
   nrow = fiparams->frame[cc].data->nrow;
   ncol = fiparams->frame[cc].data->ncol;

   n = 0;
   rowc[n] = 0.1*nrow; colc[n] = 0.1*ncol; n++;
   rowc[n] = 0.1*nrow; colc[n] = 0.5*ncol; n++;
   rowc[n] = 0.1*nrow; colc[n] = 0.9*ncol; n++;
   rowc[n] = 0.5*nrow; colc[n] = 0.1*ncol; n++;
   rowc[n] = 0.5*nrow; colc[n] = 0.5*ncol; n++;
   rowc[n] = 0.5*nrow; colc[n] = 0.9*ncol; n++;
   rowc[n] = 0.9*nrow; colc[n] = 0.1*ncol; n++;
   rowc[n] = 0.9*nrow; colc[n] = 0.5*ncol; n++;
   rowc[n] = 0.9*nrow; colc[n] = 0.9*ncol; n++;
   shAssert(n <= N_GALAP_POINT);
/*
 * Set the array cindex[] to process the canonical colour first, then the
 * others in their usual order.
 */
   i = 0; cindex[i++] = cc;
   for(c = 0;c < objc->ncolor;c++) {
      shAssert(deV_ap_correction[c] == NULL); /* no corrections are current */

      objc->color[c] = phObject1New();

      if(c != cc) {
	 cindex[i++] = c;
      }
   }
/*
 * evaluate corrections at all those points.  We have to do this position
 * by position as we need to use the model fitted at that position for
 * each band
 */
   for(i = 0; i < n; i++) {
      objc->color[cc]->flags2 &= ~OBJECT2_CANONICAL_BAND;
      objc->rowc = rowc[i]; objc->colc = colc[i];
      
      for(j = 0; j < objc->ncolor; j++) {
	 c = cindex[j];

	 obj1 = objc->color[c];
	 obj1->rowc = objc->rowc; obj1->colc = objc->colc;
	 
	 psf_get_and_fit(objc, c, fiparams,
			 &counts_model_deV[c][i], &counts_model_exp[c][i]);
	 counts_psf[c][i] = obj1->psfCounts;
	 counts_deV[c][i] = obj1->counts_deV;
	 counts_exp[c][i] = obj1->counts_exp;
      }
   }
/*
 * Find aperture corrections
 */
   for(c = 0; c < objc->ncolor; c++) {
      for(i = 0; i < n; i++) {
	 counts_deV[c][i] = counts_psf[c][i]/counts_deV[c][i];
	 counts_exp[c][i] = counts_psf[c][i]/counts_exp[c][i];
	 counts_model_deV[c][i] = counts_psf[c][i]/counts_model_deV[c][i];
	 counts_model_exp[c][i] = counts_psf[c][i]/counts_model_exp[c][i];
      }
/*
 * Fit for the aperture corrections' spatial structure
 */
      fparams = &fiparams->frame[c];

      deV_ap_correction[c] =
	phPolynomialFit(counts_deV[c], NULL, rowc, colc, n,
			nterm_row, nterm_col,
			&fparams->deV_ap_correction,
			&fparams->deV_ap_correctionErr);
      exp_ap_correction[c] =
	phPolynomialFit(counts_exp[c], NULL, rowc, colc, n,
			nterm_row, nterm_col,
			&fparams->exp_ap_correction,
			&fparams->exp_ap_correctionErr);
      deV_model_ap_correction[c] =
	phPolynomialFit(counts_model_deV[c], NULL, rowc, colc, n,
			nterm_row, nterm_col,
			&fparams->deV_model_ap_correction,
			&fparams->deV_model_ap_correctionErr);
      exp_model_ap_correction[c] =
	phPolynomialFit(counts_model_exp[c], NULL, rowc, colc, n,
			nterm_row, nterm_col,
			&fparams->exp_model_ap_correction,
			&fparams->exp_model_ap_correctionErr);
   }
/*
 * clean up
 */
   phObjcDel(objc, 1);
}

/*****************************************************************************/
/*
 * A short utility function to calculate the standard deviation of some
 * photometric quantity
 */
static float
calc_photom_sigma(float cts,		/* number of counts detected */
		  float neff,		/* effective area of object */
		  const OBJECT1 *obj1,	/* the object in question (for sky) */
		  const FRAMEPARAMS *fparams) /* for gain, dark_variance */
{
   float var = (cts + neff*obj1->sky)/fparams->gain + 
			neff*(fparams->dark_variance + pow(obj1->skyErr,2));
   return(var >= 0 ? sqrt(var) : 1e10);
}

/*****************************************************************************/
/*
 * go through an OBJC and create a new OBJECT1 for each
 * color which has a NULL entry.
 */
static void
create_object1s(OBJC *objc,		/* the OBJC in question */
		const FIELDPARAMS *fparams) /* offsets etc. */
{
   int c;
   int blend = 0;			/* is this object a blend? */
   float drow, dcol;			/* offsets from reference colour */
   float drowErr, dcolErr;		/* errors in drow, dcol */
   OBJECT1 *obj1;			/* == objc->color[c] */
   const int ncolor = objc->ncolor;
/*
 * loop through and create the missing OBJECT1s
 */
   shAssert(objc->flags3 & OBJECT3_HAS_CENTER);

   for(c = 0; c < ncolor; c++) {
      if((obj1 = objc->color[c]) != NULL) {
	 if(obj1->flags & OBJECT1_BLENDED) {
	    blend = 1;
	 }
	 continue;
      }
      
      obj1 = objc->color[c] = phObject1New();
      obj1->peaks = phPeaksNew(0);

      phOffsetDo(fparams, objc->rowc, objc->colc, 
		 fparams->ref_band_index, c,
		 NULL, NULL, &drow, &drowErr, &dcol, &dcolErr);

      obj1->colc = objc->colc + dcol;
      obj1->colcErr = sqrt(pow(objc->colcErr,2) + pow(dcolErr,2));
      obj1->rowc = objc->rowc + drow;
      obj1->rowcErr = sqrt(pow(objc->rowcErr,2) + pow(drowErr,2));

      obj1->flags |= OBJECT1_CANONICAL_CENTER;
   }

   if(blend) {
      objc->flags |= OBJECT1_BLENDED;
   }
}

/*****************************************************************************/
/*
 * Set an OBJECT1's region
 */
static void
set_region(OBJECT1 *obj1,		/* an object */
	   const REGION *reg,		/* and the region wherein it dwells */
	   OBJMASK *master_mask,	/* master mask for OBJC */
	   float dcol, float drow	/* offsets from reference colour */
	   )
{
   char buff[20];
   int i;
   int cmin, cmax, rmin, rmax;		/* unpacked from master_mask */
   SPANMASK *regmask;			/* (SPANMASK *)reg->mask */
   SPANMASK *sm;

   regmask = (SPANMASK *)reg->mask;
   shAssert(regmask->cookie == SPAN_COOKIE);

   cmin = master_mask->cmin + dcol + 0.5;
   cmax = master_mask->cmax + dcol + 0.5;
   rmin = master_mask->rmin + drow + 0.5;
   rmax = master_mask->rmax + drow + 0.5;
/*
 * non-detections and colour-to-colour offsets can move the master_mask
 * off the edge of the frame; fix this
 */
   cmin = (cmin < 0) ? 0 : ((cmin > reg->ncol - 1) ? reg->ncol - 1 : cmin);
   cmax = (cmax < 0) ? 0 : ((cmax > reg->ncol - 1) ? reg->ncol - 1 : cmax);
   rmin = (rmin < 0) ? 0 : ((rmin > reg->nrow - 1) ? reg->nrow - 1 : rmin);
   rmax = (rmax < 0) ? 0 : ((rmax > reg->nrow - 1) ? reg->nrow - 1 : rmax);
/*
 * Reset origin of OBJECT1->mask if it exists, now that it's to
 * be an atlas image
 */
   if(obj1->mask == NULL) {
      obj1->mask = phObjmaskNew(1);
   }
#if 0
   phObjmaskResetCorner(obj1->mask, rmin, cmin);
#endif
/*
 * Create OBJECT1->region->mask by trimming reg->mask, and set its origin
 * appropriately
 */
   sprintf(buff,"OBJECT1(%d)",obj1->id);
   obj1->region = shSubRegNew(buff,reg,rmax - rmin + 1, cmax - cmin + 1,
			      rmin,cmin,NO_FLAGS);
   shAssert(obj1->region != NULL);
   
   sm = phSpanmaskNew(obj1->region->nrow, obj1->region->ncol);
   for(i = 0; i < NMASK_TYPES; i++) {
      shChainDel(sm->masks[i]);
      sm->masks[i] = phTrimMaskToRect(regmask->masks[i], cmin,rmin, cmax,rmax);
   }
   phSpanmaskResetCorner(sm, rmin, cmin);
   obj1->region->mask = (MASK *)sm;
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Set the canonical centre for an OBJC. If it's detected in the
 * canonical colour, that's easy; if it isn't choose the band with
 * the largest peak counts. A (fibre?) magnitude would be better,
 * but one isn't available at the colour-dependent centres
 */
void
phObjcCenterCalc(OBJC *objc,		/* the object in question */
		 const FIELDPARAMS *fparams) /* properties of field */
{
   int best;				/* index of band with best centre */
   int canonical_color;			/* canonical colour */
   float drow, dcol;			/* convert to reference colour */
   float drowErr, dcolErr;		/* errors in drow, dcol */
   OBJECT1 *obj1;
   int i;

   shAssert(objc != NULL && fparams != NULL);
   
   canonical_color = fparams->ref_band_index;
   
   if(objc->flags3 & OBJECT3_HAS_CENTER) { /* the objc->position is correct */
      return;
   }

   obj1 = objc->color[canonical_color];
   if(obj1 != NULL && (obj1->flags & OBJECT1_DETECTED) &&
		      !(obj1->flags & OBJECT1_SATUR)) { /* easy! */
      best = canonical_color;
   } else {				/* we have to find the most significant
					   of the other bands */
      int peak = -100000;		/* largest peak seen */

      best = -1;
      for(i = 0;i < objc->ncolor;i++) {
	 obj1 = objc->color[i];

	 if(obj1 == NULL ||
	    (obj1->peaks == NULL || obj1->peaks->npeak == 0) ||
	    (obj1->flags & OBJECT1_SATUR)) {
	    continue;
	 }

	 if(obj1->peaks->peaks[0]->peak > peak) {
	    best = i;
	    peak = obj1->peaks->peaks[0]->peak;
	 }
      }
      if(best == -1) {			/* all peaks must be saturated;
					   so relax that condition */
	 for(i = 0;i < objc->ncolor;i++) {
	    obj1 = objc->color[i];
	    
	    if(obj1 == NULL || obj1->peaks == NULL || obj1->peaks->npeak == 0){
	       continue;
	    }
	    if(obj1->peaks->peaks[0]->peak > peak) {
	       best = i;
	       peak = obj1->peaks->peaks[0]->peak;
	    }
	 }
      }
      shAssert(best >= 0);		/* must have been seen in some band */
      
      obj1 = objc->color[best];
   }
/*
 * set the OBJC centre from the chosen band
 */
   objc->colc = obj1->colc;
   objc->colcErr = obj1->colcErr;
   objc->rowc = obj1->rowc;
   objc->rowcErr = obj1->rowcErr;
/*
 * convert coordinate system to the canonical band
 */
   phOffsetDo(fparams, objc->rowc, objc->colc, 
	      best, fparams->ref_band_index,
	      NULL, NULL, &drow, &drowErr, &dcol, &dcolErr);
   objc->colc += dcol;
   objc->colcErr = sqrt(pow(objc->colcErr,2) + pow(dcolErr,2));
   objc->rowc += drow;
   objc->rowcErr = sqrt(pow(objc->rowcErr,2) + pow(drowErr,2));

   objc->flags3 |= OBJECT3_HAS_CENTER;
/*
 * set the bit which indicates which band is canonical.  This may have
 * already been set for some other band (e.g. for children which were
 * not originally detected in the canonical band, but for which the deblender
 * determines that they were detectable)
 *
 * We cannot simply use obj1, as a saturated detection is acceptable.
 */
   for(i = 0;i < objc->ncolor;i++) {
      if(objc->color[i] != NULL) {
	 objc->color[i]->flags2 &= ~OBJECT2_CANONICAL_BAND;
      }
   }

   if(objc->color[canonical_color] != NULL &&
      (objc->color[canonical_color]->flags & OBJECT1_DETECTED)) {
      obj1 = objc->color[canonical_color];
   }
   obj1->flags2 |= OBJECT2_CANONICAL_BAND;
}

/*****************************************************************************/
/*
 * Set an OBJECT1's sky field
 */
static void
set_sky(OBJECT1 *obj1,
	int color,			/* which is this colour */
	const FIELDPARAMS *fiparams	/* describe FIELD */
	)
{
   float sky, skyErr;			/* old values of sky, skyErr */
   
   if(obj1->flags & OBJECT1_CHILD) {
      shAssert(obj1->skyErr > -999);	/* i.e. sky was set */
      sky = obj1->sky; skyErr = obj1->skyErr;
   } else {
      sky = skyErr = 0;
   }
   
   obj1->sky = phBinregionInterpolate(fiparams->frame[color].sky,
					    obj1->rowc,obj1->colc);
   obj1->skyErr = phBinregionInterpolate(fiparams->frame[color].skyErr,
				      obj1->rowc,obj1->colc);
   if(skyErr >= 0) {
      obj1->sky += sky;
      obj1->skyErr = sqrt(skyErr*skyErr + obj1->skyErr*obj1->skyErr);
   }
}

/*****************************************************************************/
/*
 * Given a region, a chain of OBJMASKs, a position, and a set of weights,
 * return the sum(weights*region) for an object centred at that position,
 * but only include pixels _not_ in the OBJMASKs
 */
static float
get_real_counts(const REGION *data,	/* the data */
		const CHAIN *objmasks,	/* ignore these pixels */
		int rpeak, int cpeak,	/* centre of object */
		const float **coeffs,	/* set of weights */
		int ncoeff,		/* weights is [ncoeff][ncoeff] */
		float *neff_p)		/* effective number of good pixels,
					   or NULL */
{
   float coeff;				/* == coeffs[][] */
   REGION *good = shRegNew("", 2*ncoeff - 1, 2*ncoeff - 1, TYPE_U16);
   int i, j;
   float neff;				/* effective number of good pixels */
   int rowc, colc;			/* centre of object in REGION good */
   U16 **rptr;				/* == good->rows */
   REGION *sub = shSubRegNew("", data, 2*ncoeff - 1, 2*ncoeff - 1,
			     rpeak - ncoeff + 1, cpeak - ncoeff + 1, NO_FLAGS);
   float sum;				/* sum of pixel*weight */
   float val;				/* background-subtracted pixel value */
/*
 * copy region around peak to region good
 */
   shAssert(sub != NULL);		/* i.e. ncoeff <= rpeak etc. */
   shRegIntCopy(good, sub);
   good->row0 = sub->row0; good->col0 = sub->col0;
   shRegDel(sub);
/*
 * Set pixels in objmasks to 0
 */
   phRegionSetValFromObjmaskChain(good, objmasks, 0);
/*
 * and see how much flux there is in the good pixels
 */
   rowc = colc = ncoeff - 1;
   rptr = good->rows;

   neff = sum = 0;
   if((val = rptr[rowc][colc]) != 0) {
      coeff = coeffs[0][0];
      sum += val*coeff;
      neff += coeff;
   }
   for(i = 1; i < ncoeff;i++) {
      coeff = coeffs[0][i];
      if((val = rptr[rowc][colc - i]) != 0) {
	 sum += val*coeff;
	 neff += coeff;
      }
      if((val = rptr[rowc][colc + i]) != 0) {
	 sum += val*coeffs[0][i];
	 neff += coeff;
      }
   }
   for(i = 1; i < ncoeff;i++) {
      coeff = coeffs[i][0];
      if((val = rptr[rowc - i][colc]) != 0) {
	 sum += val*coeff;
	 neff += coeff;;
      }
      if((val = rptr[rowc + i][colc]) != 0) {
	 sum += val*coeff;
	 neff += coeff;;
      }
   }
   for(i = 1; i < ncoeff;i++) {
      for(j = 1; j < ncoeff;j++) {
	 coeff = coeffs[j][i];
	 if((val = rptr[rowc - i][colc - i]) != 0) {
	    sum += val*coeff;
	    neff += coeff;;
	 }
	 if((val = rptr[rowc - i][colc + i]) != 0) {
	    sum += val*coeff;
	    neff += coeff;;
	 }
	 if((val = rptr[rowc + i][colc - i]) != 0) {
	    sum += val*coeff;
	    neff += coeff;;
	 }
	 if((val = rptr[rowc + i][colc + i]) != 0) {
	    sum += val*coeff;
	    neff += coeff;;
	 }
      }
   }
   
   shRegDel(good);

   if(neff_p != NULL) {
      *neff_p = neff;
   }
   return(sum);
}      

/*****************************************************************************/
/*
 * Calculate the 3" aperture magnitude for an object. First here are the
 * appropriate coefficients; note that they are NOT corrected to some
 * canonical seeing (how could they be? They're constants)
 */
#define NCOEFF3 9

static float *coeffs3[NCOEFF3] = { NULL };
static float coeffs3_s[NCOEFF3][NCOEFF3] = {
#if 1					/* no pixellation correction */
   {  0.93,  1.02,  1.01,  1.05,  0.25, -0.06,  0.03, -0.02,  0.01 },
   {  1.02,  1.03,  0.93,  1.02,  0.14, -0.03,  0.02, -0.01,  0.01 },
   {  1.01,  0.93,  1.16,  0.63, -0.05,  0.03, -0.02,  0.01, -0.01 },
   {  1.05,  1.02,  0.63,  0.02,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.25,  0.14, -0.05,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   { -0.06, -0.03,  0.03,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.03,  0.02, -0.02,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   { -0.02, -0.01,  0.01,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.01,  0.01, -0.01,  0.00,  0.00,  0.00,  0.00,  0.00,  0.0},
#else					/* allowing for integration over
					   the pixels */
   {  0.92,  1.03,  0.99,  1.10,  0.22, -0.07,  0.04, -0.02,  0.02 },
   {  0.01,  1.04,  0.90,  1.09,  0.09, -0.03,  0.02, -0.01,  0.01 },
   {  0.02,  0.90,  1.21,  0.64, -0.10,  0.04, -0.02,  0.02, -0.01 },
   {  0.05,  1.04,  0.64, -0.02,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.26,  0.14, -0.07,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.06, -0.04,  0.03,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.03,  0.02, -0.02,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.02, -0.01,  0.01,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
   {  0.01,  0.01, -0.01,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 },
#endif
};
static float neff_3 = 0;		/* noise-equivalent number of pixels */

/*
 * Check that the coefficients are symmetric and properly normalised
 */
static void
setup_aperture_3(const FIELDPARAMS *fiparams)
{
   int i,j;
   float fac;				/* factor to correct sum of weights */
   float sum;

   if(coeffs3[0] == NULL) {		/* set up pointers */
      for(i = 0; i < NCOEFF3; i++) {
	 coeffs3[i] = coeffs3_s[i];
      }
   }
      
   sum = 0.0;
   for(i = 1; i < NCOEFF3;i++) {
      for(j = 0; j < NCOEFF3;j++) {
	 shAssert(coeffs3[j][i] == coeffs3[i][j]);
	 sum += coeffs3[i][j];
      }
   }
   sum = coeffs3[0][0] + 4*sum;	/* should be 1/4*M_PI*(3/pixscale)^2 pix^2*/

   neff_3 = 0.25*M_PI*pow(3/fiparams->pixscale,2.0);
   fac = neff_3/sum;
   
   for(i = 0; i < NCOEFF3;i++) {
      for(j = 0; j < NCOEFF3;j++) {
	 coeffs3[i][j] *= fac;
      }
   }
}

/*****************************************************************************/
/*
 * Actually calculate the 3" aperture
 */
static void
calc_aperture_3(OBJC *objc,		/* the object in question */
		int color,		/* the desired colour */
		const CELL_STATS *prof, /* contains sinced region; or NULL */
		const FIELDPARAMS *fiparams) /* describe FIELD */
{
   int colc, rowc;			/* centre of prof->syncreg */
   const FRAMEPARAMS *fparams = &fiparams->frame[color];
   int i,j;
   OBJECT1 *const obj1 = objc->color[color];
   U16 **rptr;				/* == prof->syncreg->rows */
   float sum;
   int val;				/* value of sum of some pixels */

   if(prof == NULL) {
      float drow, dcol;			/* offsets from reference colour */

      phOffsetDo(fiparams, objc->rowc, objc->colc, fiparams->ref_band_index,
		 color, NULL, NULL, &drow, NULL, &dcol, NULL);

      prof = phProfileExtract(obj1->id, fiparams->frame[color].data,
			      objc->rowc + drow, objc->colc + dcol, 1,
			      fiparams->frame[color].bkgd + SOFT_BIAS,
			      obj1->skyErr, 0);
      if(prof == NULL || prof->syncreg == NULL) {
	 return;
      }
   }

   shAssert(prof->syncreg != NULL);
   shAssert(prof->syncreg->ncol >= 2*NCOEFF3 - 1 &&
	    prof->syncreg->nrow >= 2*NCOEFF3 - 1);
   shAssert(prof->syncreg->type == TYPE_U16 &&
	    prof->syncreg->rows != NULL);

   rowc = prof->syncreg->ncol/2;
   colc = prof->syncreg->nrow/2;
   rptr = prof->syncreg->rows;

   sum = 0;
   sum += rptr[rowc][colc]*coeffs3[0][0];
   for(i = 1; i < NCOEFF3;i++) {
      val = rptr[rowc][colc - i] + rptr[rowc][colc + i] +
	    rptr[rowc - i][colc] + rptr[rowc + i][colc];
      sum += val*coeffs3[0][i];
   }
   for(i = 1; i < NCOEFF3;i++) {
      for(j = 1; j < NCOEFF3;j++) {
	 val = rptr[rowc - j][colc - i] + rptr[rowc - j][colc + i] +
	   rptr[rowc + j][colc - i] + rptr[rowc + j][colc + i];
	 sum += val*coeffs3[j][i];
      }
   }
   sum -= SOFT_BIAS*neff_3;		/* subtract away SOFT_BIAS */

   obj1->fiberCounts = sum - fparams->bkgd*neff_3;
   obj1->fiberCountsErr = calc_photom_sigma(sum, neff_3, obj1, fparams);
/*
 * If some of the pixels have been interpolated over, the number of DN
 * actually detected will be smaller than sum, and we have to allow for
 * this in the noise calculation. So find out how many DN were really
 * detected (slightly sloppily, as we don't sinc shift when doing this)
 */
   if((obj1->flags & OBJECT1_INTERP) && fparams->data != NULL) {
      float flux_correction;		/* correction for interpolated flux */
      float neff;			/* n_effective for measurement */
      const SPANMASK *sm = (SPANMASK *)fparams->data->mask;
      const CHAIN *objmasks = sm->masks[S_MASK_INTERP];

      if(objmasks != NULL) {
	 float drow, dcol;
	 phOffsetDo(fiparams, objc->rowc, objc->colc, 
		    fiparams->ref_band_index, color,
		    NULL, NULL, &drow, NULL, &dcol, NULL);

	 sum = get_real_counts(fparams->data, objmasks,
			       objc->rowc + drow, objc->colc + dcol,
			       (const float **)coeffs3, NCOEFF3, &neff);
	 sum -= SOFT_BIAS*neff;		/* subtract away SOFT_BIAS */
	 obj1->fiberCountsErr = calc_photom_sigma(sum, neff, obj1, fparams);

	 sum -= fparams->bkgd*neff;	/* subtract any non-subtracted sky */
	 if(sum < 1 || obj1->fiberCounts < 1) {
	    flux_correction = 0;	/* can't find a reliable correction */
	 } else {
	    flux_correction = obj1->fiberCounts/sum;
	 }

/*
 * If the flux correction is too small there's a problem; we cannot simply
 * force the correction to be >= 1 as the sloppy non-shifted calculation
 * of get_real_counts() can be wrong by 100%; this is OK as it only affects
 * the errors, not the counts. The choice of 0.5 is XXX
 */
	 if(flux_correction < 0.5 ||
				   (obj1->flags2 & OBJECT2_BAD_COUNTS_ERROR)) {
	    if(obj1->fiberCountsErr < fabs(obj1->fiberCounts)) {
	       obj1->flags2 |= OBJECT2_BAD_COUNTS_ERROR;
	       obj1->fiberCountsErr = fabs(obj1->fiberCounts);
	    }
	 } else {
	    obj1->fiberCountsErr *= flux_correction;
	 }
      }
   }
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Calculate an aperture magnitude from a radial profile. Used in finding
 * PSF aperture corrections.
 */
float
phApertureCounts(const CELL_STATS *prof, /* the extracted profile */
		 int nann)		/* number of annuli to use */
{
   double apCounts;			/* aperture counts for PSF */
   int i;

   shAssert(nann > 0);			/* i.e. it was set */

   setup_geometry();

   /* calculate aperture counts */ 
   apCounts = 0.0;
   for(i = 0; i < nann; i++) {
      apCounts += area[i]*phProfileMean(prof, i, 0, 1, NULL);
   }

   return(apCounts);
}

/*****************************************************************************/
/*
 * PSF magnitudes next. For each frame, we need to calculate the proper weights
 * based on the PSFs in each colour; thereafter the actual calculation of the
 * counts is a simple sum
 *
 * There are two sets of routines to do this; those using the sinc-shifted
 * region, and those using the extracted cell profiles. These two routines
 * are the sinc-shifted region routines, and they are followed by a number
 * that use the profile
 */
float
phPsfCountsSetupFromSincRegion(int color, /* the colour in question */
			       const DGPSF *psf) /* PSF in this band */
{
   float fac;                           /* factor to correct sum of weights */
   float sum, sum2;                     /* sums of P and P^2 */
   float val;                           /* value of PSF at a point */
   int x, y;
   float x2;                            /* == x^2 */
   float xfac_1, xfac_2;                /* == 1/(2*sigx_[12]^2) */
   float y2_1, y2_2;                    /* == y^2/(2*sigy_[12]^2) */
   float b, sigx_1, sigy_1, sigx_2, sigy_2; /* unpacked from psf */

   shAssert(color >= 0 && color < NCOLOR);

   sigx_1 = psf->sigma1_2G; sigy_1 = psf->sigma1_2G;
   sigx_2 = psf->sigma2_2G; sigy_2 = psf->sigma2_2G;
   b = psf->b_2G;
   
   xfac_1 = 1/(2*sigx_1*sigx_1);
   xfac_2 = 1/(2*sigx_2*sigx_2);

   if(coeffs_psf[color] == NULL) {	/* we have to set the pointers */
      int i;
      
      coeffs_psf[color] = coeffs_psf_s[color];
      for(i = 0; i < NCOEFF_PSF; i++) {
	 coeffs_psf_s[color][i] = coeffs_psf_ss[color][i];
      }
   }
   
   sum = sum2 = 0.0;
   for(y = 0;y < NCOEFF_PSF;y++) {
      y2_1 = 0.5*pow(y/sigy_1,2);
      y2_2 = 0.5*pow(y/sigy_2,2);
      for(x = 0; x < NCOEFF_PSF;x++) {
         x2 = x*x;
         val = exp(-(y2_1 + xfac_1*x2)) + b*exp(-(y2_2 + xfac_2*x2));
         coeffs_psf[color][y][x] = val;
         
         if(x > 0) {                    /* avoid double counting */
            if(y > 0) {
               sum += 4*val;
               sum2 += 4*val*val;
            } else {
               sum += 2*val;
               sum2 += 2*val*val;
            }
         } else {
            if(y > 0) {
               sum += 2*val;
               sum2 += 2*val*val;
            } else {
               sum += val;
               sum2 += val*val;
            }
         }
      }
   }
/*
 * normalise so that the product of the data with the weights is the
 * number of counts in the object
 */
   fac = sum/sum2;
   for(y = 0;y < NCOEFF_PSF;y++) {
      for(x = 0; x < NCOEFF_PSF;x++) {
         coeffs_psf[color][y][x] *= fac;
      }
   }
   
   neff_psf[color] = sum*fac;

   return(neff_psf[color]);
}

/*****************************************************************************/
/*
 * Actually calculate the PSF counts from the sinc shifted REGION in the prof
 */
void
phPsfCountsFromSincRegion(OBJC *objc,	/* the OBJC in question */
			  int c,	/* desired colour */
			  const CELL_STATS *prof, /* contains sinced region */
			  const FIELDPARAMS *fiparams) /* describe field */
{
   const FRAMEPARAMS *fparams;		/* == &fiparams->frame[c] */
   OBJECT1 *obj1;			/* object being measured */
   double psfCounts;			/* measured psfCounts */

   shAssert(objc != NULL);
   shAssert(c >= 0 && c < objc->ncolor && objc->color[c] != NULL);
   shAssert(prof != NULL && prof->syncreg != NULL);
   shAssert(prof->syncreg->ncol >= 2*NCOEFF_PSF - 1 &&
            prof->syncreg->nrow >= 2*NCOEFF_PSF - 1);
   shAssert(prof->syncreg->type == TYPE_U16 && prof->syncreg->rows != NULL);
   shAssert(fiparams != NULL && c >= 0 && c < fiparams->ncolor);
   shAssert(coeffs_psf[c] != NULL);

   fparams = &fiparams->frame[c];
   obj1 = objc->color[c];

   obj1->psfCounts = psfCounts = 
			      calc_psf_counts(prof->syncreg, c, fparams->bkgd);
   obj1->psfCountsErr = calc_photom_sigma(psfCounts,
					  neff_psf[c], obj1, fparams);
/*
 * If some of the pixels have been interpolated over, the number of DN
 * actually detected will be smaller than psfCounts, and we have to allow for
 * this in the noise calculation. So find out how many DN were really
 * detected (slightly sloppily, as we don't sinc shift when doing this)
 */
   if(obj1->flags & OBJECT1_INTERP && fparams->data != NULL) {
      float flux_correction = 1;	/* correction for interpolated flux */
      float neff;			/* n_effective for measurement */
      const SPANMASK *sm = (SPANMASK *)fparams->data->mask;
      const CHAIN *objmasks = sm->masks[S_MASK_INTERP];

      if(objmasks != NULL) {
	 float drow, dcol;
	 phOffsetDo(fiparams, objc->rowc, objc->colc, 
		    fiparams->ref_band_index, c,
		    NULL, NULL, &drow, NULL, &dcol, NULL);

	 psfCounts = get_real_counts(fparams->data, objmasks,
				     objc->rowc + drow, objc->colc + dcol,
				     (const float **)coeffs_psf[c],
				     NCOEFF_PSF, &neff);
	 psfCounts -= SOFT_BIAS*neff;	/* subtract away SOFT_BIAS */
	 obj1->psfCountsErr = calc_photom_sigma(psfCounts,neff, obj1, fparams);

	 psfCounts -= fparams->bkgd*neff; /* subtract any non-subtracted sky */
	 if(psfCounts < 1 || obj1->psfCounts < 1) {
	    flux_correction = 0;	/* can't find a reliable correction */
	 } else {
	    flux_correction = obj1->psfCounts/psfCounts;
	 }
      }
/*
 * If the flux correction is too small there's a problem; we cannot simply
 * force the correction to be >= 1 as the sloppy non-shifted calculation
 * of get_real_counts() can be wrong by 100%; this is OK as it only affects
 * the errors, not the counts. The choice of 0.5 is XXX
 */
      if(flux_correction < 0.5 || (obj1->flags2 & OBJECT2_BAD_COUNTS_ERROR)) {
	 if(obj1->psfCountsErr < fabs(obj1->psfCounts)) {
	    obj1->flags2 |= OBJECT2_BAD_COUNTS_ERROR;
	    obj1->psfCountsErr = fabs(obj1->psfCounts);
	 }
      } else {
	 obj1->psfCountsErr *= flux_correction;
	 if(flux_correction > 1.2) {	/* XXX */
	    obj1->flags2 |= OBJECT2_PSF_FLUX_INTERP;
	 }
      }
   }
/*
 * apply aperture correction
 */
   obj1->psfCounts *= fparams->psf_app_correction;
   obj1->psfCountsErr *= fparams->psf_app_correction;

   obj1->psfCountsErr = sqrt(pow(obj1->psfCountsErr, 2) +
		       pow(obj1->psfCounts*fparams->psf_app_correctionErr, 2));

}

/*
 * helper function to evaluate the PSF counts
 */
static double
calc_psf_counts(const REGION *reg,	/* region with sinc-centered object */
		int c,			/* the colour in question */
		float bkgd)		/* non-subtracted background */
{
   int rowc = reg->ncol/2;		/* centre of */
   int colc = reg->nrow/2;		/*        the object in question */
   int i, j;
   U16 **rptr = reg->rows;		/* == reg->rows */
   double sum = 0;			/* desired sum */
   int val;                             /* value of sum of some pixels */
   
   sum += rptr[rowc][colc]*coeffs_psf[c][0][0];
   for(i = 1; i < NCOEFF_PSF;i++) {
      val = rptr[rowc][colc - i] + rptr[rowc][colc + i];
      sum += val*coeffs_psf[c][0][i];
   }
   for(i = 1; i < NCOEFF_PSF;i++) {
      val = rptr[rowc - i][colc] + rptr[rowc + i][colc];
      sum += val*coeffs_psf[c][i][0];
   }
   for(i = 1; i < NCOEFF_PSF;i++) {
      for(j = 1; j < NCOEFF_PSF;j++) {
         val = rptr[rowc - j][colc - i] + rptr[rowc - j][colc + i] +
           rptr[rowc + j][colc - i] + rptr[rowc + j][colc + i];
         sum += val*coeffs_psf[c][j][i];
      }
   }
   sum -= SOFT_BIAS*neff_psf[c];    /* subtract away SOFT_BIAS */

   sum -= bkgd*neff_psf[c];

   return(sum);
}

/*****************************************************************************/
/*
 * wrapper for phPsfCountsFromSincRegion() for the use of PSP.
 *
 * If the central pixel is negative, we assume that the passed in region
 * has already been centred up, and is of size
 *   (SYNC_REG_SIZE + 2*SYNCEXTRA)x(SYNC_REG_SIZE + 2*SYNCEXTRA)
 * such regions are most conveniently obtained by unpacking a CELL_STATS...
 *
 * N.b. bkgd doesn't include the SOFT_BIAS
 */
float
phPsfCountsFromRegion(const REGION *reg, /* region containing object */
		      int c,		/* band that data was taken in */
		      float rowc,	/* central pixel of object, */
		      float colc,	/*   or -ve if region is pre-shifted */
		      float ap_correction, /* aperture correction to apply */
		      float bkgd,	/* reg's sky level (not subtracted) */
		      float gain,	/* gain of amplifier */
		      float dark_variance, /* variance of dark background */
		      float bkgd_sub,	/* pre-subtracted background level */
		      float *err)	/* error in counts, or NULL */
{
   FIELDPARAMS *fiparams = phFieldparamsNew("ugriz");
   OBJC *objc = phObjcNew(NCOLOR);
   OBJECT1 obj1 = { 0 };		/* initialise const int id */
   CELL_STATS prof;
   float psfCounts;			/* the desired answer */
   REGION *shifted = NULL;		/* sinc-centred copy of reg */

   shAssert(reg != NULL && reg->type == TYPE_U16);
   shAssert(c >= 0 && c < NCOLOR);
   objc->color[c] = &obj1;		/* remember not to free it */
/*
 * shift region to desired centre, unless it's already been done for us
 */
   if(rowc < 0) {			/* yes; already centred up */
      shAssert(colc < 0);
      shAssert(reg->nrow == SYNC_REG_SIZE + 2*SYNCEXTRA &&
	       reg->ncol == SYNC_REG_SIZE + 2*SYNCEXTRA);

      prof.syncreg = reg;
   } else {				/* we have to do the work */
      rowc -= reg->row0;
      colc -= reg->col0;
      shifted = phRegIntShift(NULL, reg, NULL, 11,
			      0.5*reg->nrow - rowc, 0.5*reg->ncol - colc);
/*
 * repack inputs. Note that the names sky/bkgd are not quite
 * what frames usually uses, where sky is assumed to have been subtracted.
 */
      prof.syncreg = shSubRegNew("phPsf2*NCOEFF_PSF - 1);CountsFromRegion", shifted,
				 2*NCOEFF_PSF - 1, 2*NCOEFF_PSF - 1,
				 reg->nrow/2 - NCOEFF_PSF + 1,
				 reg->ncol/2 - NCOEFF_PSF + 1, NO_FLAGS);
      shAssert(prof.syncreg != NULL);
   }

   fiparams->frame[c].gain = gain;
   fiparams->frame[c].dark_variance = dark_variance;
   fiparams->frame[c].bkgd = bkgd;
   fiparams->frame[c].psf_app_correction = ap_correction;
   objc->color[c]->sky = bkgd_sub;
   objc->color[c]->skyErr = 0;
/*
 * do the work
 */
   phPsfCountsFromSincRegion(objc, c, &prof, fiparams);
   psfCounts = objc->color[c]->psfCounts;
   if(err != NULL) {
      *err = objc->color[c]->psfCountsErr;
   }
/*
 * cleanup
 */
   if(shifted != NULL) {
      shRegDel((REGION *)prof.syncreg);
      shRegDel(shifted);
   }
   phFieldparamsDel(fiparams);
   objc->color[c] = NULL; phObjcDel(objc, 1);

   return(psfCounts);
}

/*****************************************************************************/
/*
 * apply galaxy aperture corrections; the best fit model is already corrected
 */
static void
galaxy_ap_corrections(OBJC *objc,	/* the object to correct */
		      int c,		/* desired band */
		      const FIELDPARAMS *fiparams)
{
   float ac;				/* {deV,exp}_ap_correction at obj1 */
   OBJECT1 *obj1 = objc->color[c];	/* unpacked */
   const FRAMEPARAMS *fparams = &fiparams->frame[c]; /* unpacked */
   
   if(deV_ap_correction[c] == NULL) {	/* not available */
      shAssert(exp_ap_correction[c] == NULL);
      return;
   }

   ac = phPolynomialEval(deV_ap_correction[c], obj1->rowc, obj1->colc);
   obj1->counts_deV *= ac;
   obj1->counts_deVErr *= ac;
   obj1->counts_deVErr = sqrt(pow(obj1->counts_deVErr, 2) +
			      pow(fparams->deV_ap_correctionErr, 2));
   
   ac = phPolynomialEval(exp_ap_correction[c], obj1->rowc, obj1->colc);
   obj1->counts_exp *= ac;
   obj1->counts_expErr *= ac;
   obj1->counts_expErr = sqrt(pow(obj1->counts_expErr, 2) +
			      pow(fparams->exp_ap_correctionErr, 2));
}

/*****************************************************************************/
/*
 * Given an OBJC will models fitted, return the counts in the best-fit
 * model in the canonical band, fitted to some other band
 */
#define LOCAL_SKY 0

static void
calc_model_counts(OBJC *objc,		/* object to fit */
		  int c,		/* desired colour */
		  const CELL_STATS *cstats, /* extracted profile */
		  const FIELDPARAMS *fiparams) /* gain etc. */
{
   float ab, phi, r;			/* == obj1->r_deV or obj1->r_exp etc.*/
   int cc;				/* canonical colour */
   int class;				/* type of model desired */
   FRAMEPARAMS *fparams = &fiparams->frame[c];
   int i;
   int is_deV;				/* the deV fit is better than the exp*/
   OBJECT1 *obj1 = objc->color[c];	/* object being measured */
   OBJECT1 *obj1_cc;			/* canonical object */
#if LOCAL_SKY
   float sky, *skyp = &sky;		/* sky level */
#else
   float *skyp = NULL;			/* don't fit sky */
#endif
   int sky_noise_only = 1;		/* cell variance == sky noise? */
   int use_difference = 0;		/* include difference of cell pairs
					   in cell variance */
/*
 * find and unpack the canonical band
 */
   cc = -1;
   for(i = 0; i < objc->ncolor; i++) {
      if(objc->color[i] != NULL &&
	 objc->color[i]->flags2 & OBJECT2_CANONICAL_BAND) {
	 cc = i;
	 break;
      }
   }
   shAssert(cc >= 0);

   obj1_cc = objc->color[cc];
   is_deV = (obj1_cc->deV_L > obj1_cc->exp_L) ? 1: 0;

#if !LOCAL_SKY
   if(c == cc) {			/* easy; we've already done the work */
      obj1->counts_model = is_deV ?
	objc->color[cc]->counts_deV : objc->color[cc]->counts_exp;
      obj1->counts_modelErr = is_deV ?
	objc->color[cc]->counts_deVErr : objc->color[cc]->counts_expErr;
   } else
#endif
   {
/*
 * we have to do some work.
 */
      if(is_deV) {
	 class = DEV_MODEL;
	 ab = obj1_cc->ab_deV; phi = obj1_cc->phi_deV; r = obj1_cc->r_deV;
      } else {
	 class = EXP_MODEL;
	 ab = obj1_cc->ab_exp; phi = obj1_cc->phi_exp; r = obj1_cc->r_exp;
      }
      phFitCellAsKnownModel(objc, c, cstats, fiparams, 0, class, ab, phi, r,
			    use_difference, sky_noise_only,
			    &obj1->counts_model, &obj1->counts_modelErr, skyp,
			    NULL);
/*
 * Apply the aperture correction
 */
      if(fiparams->use_galaxy_ap_correction) {
	 float model_ap_correction, model_ap_correctionErr;
      
#if 1
	 if(deV_model_ap_correction[c] == NULL) {
	    shAssert(exp_model_ap_correction[c] == NULL);
	    
	    model_ap_correction = 1.0;
	    model_ap_correctionErr = 0.0;
	 } else {
	    if(is_deV) {
	       model_ap_correction =
		 phPolynomialEval(deV_model_ap_correction[c],
						       obj1->rowc, obj1->colc);
	       model_ap_correctionErr = fparams->deV_model_ap_correctionErr;
	    } else {
	       model_ap_correction =
		 phPolynomialEval(exp_model_ap_correction[c],
						       obj1->rowc, obj1->colc);
	       model_ap_correctionErr = fparams->exp_model_ap_correctionErr;
	    }
	 }
#else
	 if(is_deV) {
	    model_ap_correction = fparams->deV_model_ap_correction;
	    model_ap_correctionErr = fparams->deV_model_ap_correctionErr;
	 } else {
	    model_ap_correction = fparams->exp_model_ap_correction;
	    model_ap_correctionErr = fparams->exp_model_ap_correctionErr;
	 }
#endif
      
	 obj1->counts_model *= model_ap_correction;
	 obj1->counts_modelErr *= model_ap_correction;
	 obj1->counts_modelErr = sqrt(pow(obj1->counts_modelErr, 2) +
				      pow(model_ap_correctionErr, 2));
      }
   }
}

/*****************************************************************************/
/*
 * calculate the "texture" of an object
 */
static void
calc_texture(OBJC *objc,		/* object to measure */
	     int c,			/* desired colour */
	     const CELL_STATS *cstats,	/* extracted profile */
	     const FIELDPARAMS *fiparams) /* gain etc. */
{
   float ab, phi, r;			/* == obj1->r_deV or obj1->r_exp etc.*/
   float chisq_model;			/* chisq for original model fit */
   float chisq_texture;			/* chisq for no-difference-variance
					   fit to profile */
   int class;				/* type of model desired */
   int nu_model, nu_texture;		/* d.o.f. for chisq_... */
   OBJECT1 *const obj1 = objc->color[c]; /* object being measured */

   if(obj1->star_L > obj1->deV_L && obj1->star_L > obj1->exp_L) {
      class = PSF_MODEL;
      ab = phi = r = 0;
      chisq_model = obj1->chisq_star;
      nu_model = obj1->nu_star;
   } else {
      if(obj1->deV_L > obj1->exp_L) {
	 class = DEV_MODEL;
	 ab = obj1->ab_deV; phi = obj1->phi_deV; r = obj1->r_deV;
	 chisq_model = obj1->chisq_deV;
	 nu_model = obj1->nu_deV;
      } else {
	 class = EXP_MODEL;
	 ab = obj1->ab_exp; phi = obj1->phi_exp; r = obj1->r_exp;
	 chisq_model = obj1->chisq_exp;
	 nu_model = obj1->nu_exp;
      }
   }

   chisq_texture =
     phFitCellAsKnownModel(objc, c, cstats, fiparams, 0, class, ab, phi, r,
			   0, 1, NULL, NULL, NULL, &nu_texture);

   obj1->texture = (nu_texture == 0 || nu_model == 0) ? 0 :
			       chisq_texture/nu_texture - chisq_model/nu_model;
}

/*****************************************************************************/
/*
 * Copy the profiles into the OBJECT1, and return the number of "good" points,
 * those with S/N >= SN_min
 */
static int
set_profiles(OBJECT1 *obj1,
	     int color,
	     const CELL_STATS *prof,
	     const FIELDPARAMS *fiparams
	     )
{
   static int clip = 0;			/* how hard to clip */
   const FRAMEPARAMS *fparams = &fiparams->frame[color];
   int i;
   int ngood;				/* number of good points */
   float noise;				/* estimate of sigma for profile */
   static float SN_min = 0;		/* minimum acceptable S/N ratio */

   ngood = -1;
   for(i = 0;i < prof->nannuli;i++) {
      if(i < obj1->nprof) {		/* already measured */
	 noise = obj1->profErr[i];
      } else {
	 obj1->profMed[i] = phProfileMedian(prof, i, clip, 1, &noise);
	 obj1->profMean[i] = phProfileMean(prof, i, clip, 1, NULL);

	 if(i == 0) {			/* central pixel, so we have no
					   knowledge of its noise */
	    noise = calc_photom_sigma(obj1->profMean[0],1,obj1,fparams);
	 }
	 obj1->profErr[i] = noise;
      }
/*
 * if we don't have a measurement, or if noise is a NaN, give up
 */
      if(obj1->profMed[i] < -999 || noise != noise) {
	 if(ngood < 0) {
	    ngood = i;
	 }
      }
      
      if(noise == 0) {			/* no problem there! */
	 continue;
      }
      if(obj1->profMed[i]/noise < SN_min) {
	 if(ngood < 0) {
	    ngood = i;
	 }
      }
   }
   obj1->nprof = prof->nannuli;

   return(ngood >= 0 ? ngood : obj1->nprof);
}

/*****************************************************************************/
/*
 * Truncate an object's radial profile.
 *
 * The profile extends to the larger of the largest radius used for
 * any measurements, namely petro_f2*(estimated r_P), and the number
 * of points declared to be `good' on an S/N criterion in set_profiles()
 */
static void
trim_profiles(OBJC *objc, int color, int ngood, float rmax)
{
   int i;
   int nprof;				/* number of radial bins to keep */
   OBJECT1 *obj1 = objc->color[color];

   for(i = 0;i < obj1->nprof;i++) {
      if(r[i] > rmax && i > ngood) {
	 break;
      }
   }
   nprof = i;
/*
 * we want to keep nprof points, so set the rest to unknown
 */
   for(i = nprof;i < obj1->nprof;i++) {
      obj1->profMean[i] = obj1->profMed[i] = obj1->profErr[i] = -1000;
   }
   obj1->nprof = nprof;

   if(ngood < nprof) {
      obj1->flags |= OBJECT1_BAD_RADIAL;
   }
}

/*****************************************************************************/
/*
 * Calculate Petrosian Quantities
 *
 * Find the cumulative radial profile.
 *
 * If requested, smooth the radial profile allowing for the errors in the
 * points. We explicitly force the smoothed profile to be an even function
 */
static void
get_cumul(const float *prof,		/* profile to cumulate */
	  const float *proferr,		/* error in profile; can be NULL
					   if no smoothing is required */
	  int nprof,			/* number of points */
	  float *cumul,			/* the desired answers */
	  int smooth_profs)		/* smooth? */
{
   float chisq = 0;
   int i;
   float rs[2*(NANN+2) + 1], smoothed[2*(NANN+2) + 1], errs[2*(NANN+2) + 1];
   SPLINE *smoothed_sp;			/* the smoothed spline */
   float val;				/* value of profile, error, or radius*/
      
   if(smooth_profs) {
/*
 * Use a spline to smooth the points to an acceptable chisq.
 *
 * We add two extra points to pin the outer parts of the profile to 0, in the
 * absence of any data. Why two? Because the smoothed spline is a natural
 * spline, i.e. its second derivative vanishes at the last knot, so we
 * want two points to make the resulting linear extrapolation "horizontal"
 */
      static float fake_err = 0.1;	/* error associated with fake points */
      const int nprof_s = nprof + 2;
      
      shAssert(proferr != NULL);
      for(i = 0;i < nprof;i++) {
	 rs[nprof_s + i] = rs[nprof_s - i - 1] = asinh_ph(r[i + 1]);
	 rs[nprof_s - i - 1] = -rs[nprof_s - i - 1];
	 smoothed[nprof_s + i] = smoothed[nprof_s - i - 1] = prof[i];
	 
	 if((val = proferr[i]) == 0.0f) {
	    val = EPSILON_f;
	 }
	 errs[nprof_s + i] = errs[nprof_s - i - 1] = val;
      }
/*
 * add extra points
 */
      rs[nprof_s + i] = rs[nprof_s - i - 1] = -rs[2] + log(2*10.0);
      rs[nprof_s - i - 1] = -rs[nprof_s - i - 1];
      smoothed[nprof_s + i] = smoothed[nprof_s - i - 1] = 0;
      errs[nprof_s + i] = errs[nprof_s - i - 1] = fake_err;

      i++;
      rs[nprof_s + i] = rs[nprof_s - i - 1] = -rs[2] + log(3*10.0);
      rs[nprof_s - i - 1] = -rs[nprof_s - i - 1];
      smoothed[nprof_s + i] = smoothed[nprof_s - i - 1] = 0;
      errs[nprof_s + i] = errs[nprof_s - i - 1] = fake_err;

      smoothed_sp =
	  phSplineFindSmoothed(rs, smoothed, errs, 2*nprof_s, 1, &chisq, errs);
      shAssert(smoothed_sp != NULL);

      cumul[0] = 0;
      for(i = 0;i < nprof;i++) {
	 cumul[i + 1] = cumul[i] + smoothed_sp->coeffs[0][nprof_s + i]*area[i];
      }

      phSplineDel(smoothed_sp);
   } else {				/* no need to smooth */
      cumul[0] = 0;
      for(i = 0;i < nprof;i++) {
	 cumul[i + 1] = cumul[i] + area[i]*prof[i];
      }
   }
}

/*****************************************************************************/
/*
 * Given a cumulative profile, differentiate it to find the surface brightness
 * profile at the points r[], and return it as a spline
 */
static SPLINE *
get_sb_profile(const float *cumul,	/* cumulative profile */
	       const SPLINE *cumul_sp,	/* cumulative profile as spline */
	       int n)			/* number of points desired */
{
   float as_sb[NANN + 1];		/* asinh(sb[]) */
   int i;
   float sb[NANN + 1];			/* surface brightness profile */
   SPLINE *sb_sp;			/*              spline of sb[] */

   phSplineDerivative(cumul_sp, as_r, sb, n);

   for(i = 1;i < n;i++) {
      sb[i] *= sqrt((1 + cumul[i]*cumul[i])/(1 + r[i]*r[i])); /* undo asinh */
      sb[i] /= (2*M_PI*r[i]);		/* undo the cumulation */
      as_sb[i] = asinh_ph(sb[i]);
   }
   sb_sp = phSplineFindTautEven(as_r + 1, as_sb + 1, n - 1, Gamma);
   shAssert(sb_sp != NULL);

   return(sb_sp);
}

/*****************************************************************************/
/*
 * find and return the Petrosian radius, setting appropriate bits in *flags
 */
static float
get_petrosianRad(float pf1,		/* Petrosian ratio */
		 float pf2,		/* minimum surface brightness */
		 const SPLINE *sb_sp,	/* spline of surface brightness */
		 const float *ratio,	/* Petrosian ratio */
		 int ngood,		/* number of good points in this band*/
		 int nratio,		/* dimension of ratio */
		 int *flags)		/* O: flags to set */
{
   float as_rmax;			/* asinh(max radius to look for r_P) */
   int i;
   int n_petro_r;			/* number of Petrosian radii */
   float petro_r[3*NANN];		/* Petrosian radii; we can't have more
					   than 3 in any radial interval */
   float rr[3], vals[3];		/* for calls to phSplineInterpolate */

   *flags = 0;

   if(nratio < 2) {
      n_petro_r = 0;
   } else {
      SPLINE *ratio_sp = phSplineFindTautEven(as_r,ratio,nratio + 1,Gamma);
      shAssert(ratio_sp != NULL);

      as_rmax = as_r[ngood];
      if(phSplineRoots(ratio_sp, 0, 0, as_rmax, rr, 3) != 0) {
	 as_rmax = rr[0];
      }
      n_petro_r = phSplineRoots(ratio_sp, pf1, 0, as_rmax, petro_r, 3*NANN);
      shAssert(n_petro_r >= 0);

      phSplineDel(ratio_sp); 
   }
/*
 * Throw out Petrosian radii with a surface brightness < petro_f2
 */
   for(i = 0;i < n_petro_r;i++) {
      rr[0] = petro_r[i];
      phSplineInterpolate(sb_sp,rr,vals,1);
      vals[0] = sinh_ph(vals[0]);

      if(vals[0] < pf2) {		/* not acceptable */
	 int j;
	 for(j = i;j < n_petro_r - 1;j++) {
	    petro_r[j] = petro_r[j + 1];
	 }
	 i--; n_petro_r--;

	 *flags |= OBJECT1_PETROFAINT;
      }
   }
/*
 * See what we got, and act accordingly
 */
   if(n_petro_r > 0) {			/* we have a Petrosian radius */
      if(n_petro_r > 1) {		/* too many radii; use the largest */
	 *flags |= OBJECT1_MANYPETRO;
      }
      return(sinh_ph(petro_r[n_petro_r - 1])); /* n.b. interpolation
						       was in asinh(r) */
   }
/*
 * We failed to find a Petrosian radius
 */
   *flags |= OBJECT1_NOPETRO;
   if(!(*flags & OBJECT1_PETROFAINT)) {	/* we didn't reject any candidates */
      *flags |= OBJECT1_NOPETRO_BIG;
   }     

   return(0);
}

/*****************************************************************************/
/*
 * Return the error in the radius containing frac of the counts.  This
 * is done naively
 */
static float
petro_rNN_error(const SPLINE *cumul_sp,	/* cumulative profile */
		float rmax,		/* maximum radius to look for radii */
		float frac,		/* desired fraction of counts */
		float counts,		/* Petrosian counts for object */
		float countsErr)	/* error in Petrosian counts */
{
   int nm, np;				/* number of solutions for -+ */
   float rm, rp;			/* radii for counts -+ countsErr */

   nm = phSplineRoots(cumul_sp, asinh_ph(frac*(counts - countsErr)),
		     0, asinh_ph(rmax), &rm, 1);
   if(nm == 0) {
      return(-1000);
   }
   rm = sinh_ph(rm);
      
   np = phSplineRoots(cumul_sp, asinh_ph(frac*(counts + countsErr)),
		     0, asinh_ph(rmax), &rp, 1);
   if(np == 0) {
      return(-1000);
   }
   rp = sinh_ph(rp);

   return(0.5*(rp - rm));
}

/*****************************************************************************/
/*
 * Find the Petrosian radius, flux, and R50/R90
 */
static SPLINE *
calc_petrosian(OBJC *objc,
	       int color,
	       int ngood,
	       const FIELDPARAMS *fiparams)
{
   OBJECT1 *const obj1 = objc->color[color]; /* the OBJECT1 in question */
   float cumul[NANN + 1];		/* cumulative radial profile */
   float as_cumul[NANN + 1];		/* asinh(cumul[]) */
   SPLINE *cumul_sp;			/* spline of as_cumul[] */
   const float gain = fiparams->frame[color].gain; /* gain of CCD */
   int i;
   int ngood_ref;			/* ngood for reference colour */
   int nprof = obj1->nprof;		/* number of radial points */
   int nratio;				/* number of values of Petro. ratio */
   float pf2;				/* petro_f2, in counts/pix^2 */
   const float *const prof = (median_profs ? obj1->profMed : obj1->profMean);
   float r_P;				/* canonical Petrosian radius */
   float r_PErr;			/*                        and error */
   float r_P_p4;			/* r_P*petro_4 */
   float ratio[NANN + 1];		/* Petrosian ratio and */
   float ratioErr[NANN + 1];		/*      Petrosian ratio's s.d. */
   float rr[3], vals[3];		/* for calls to phSplineInterpolate */
   SPLINE *sb_sp;			/* spline of surface brightness prof */
   float skyVar;			/* per-pixel variance due to sky */
   const int smooth_petro_radii = 	/* smooth Petrosian ratio? */
     fiparams->smooth_petro_radii;
   const int petro_gcv_errors =		/* Find R_P errors by GCV? */
     fiparams->petro_gcv_errors;
/*
 * Number of points in the cumulative profile in the reference
 * band
 */
   ngood_ref = (color == fiparams->ref_band_index) ?
			     ngood : objc->color[fiparams->ref_band_index]->nprof;

   if(nprof < 2) {
      obj1->flags |= OBJECT1_NOPROFILE | OBJECT1_NOPETRO;
      return(NULL);
   }
/*
 * Find the (possibly smoothed) radial profile
 */
   get_cumul(prof, obj1->profErr, nprof, cumul, smooth_petro_radii);
/*
 * profile in reference colour may go further than in this band; if so set
 * extrapolation to 0
 */
   for(i = nprof;i < ngood_ref;i++) {
      cumul[i + 1] = cumul[nprof];
   }
   nprof = ngood_ref;
/*
 * Spline the cumulative profile. We use an asinh() to control
 * the dynamic range in the profile.
 *
 * The choice of boundary conditions at the origin's a little tricky. The
 * gradient d(asinh(cumul))/d(asinh(r)) is 0 at the origin --- but only
 * _very_ close to the origin. For a constant surface brightness source,
 * once cumul >> 1, the gradient becomes very large (~1/r for r << 1).
 * Experiment shows that the best results are achieved by not imposing
 * any symmetry (and thus gradient) constraints at r == 0
 */   
   for(i = 0;i <= nprof;i++) {
      as_cumul[i] = asinh_ph(cumul[i]);
   }

   cumul_sp = phSplineFindTaut(as_r, as_cumul, nprof + 1, Gamma);
   shAssert(cumul_sp != NULL);

   sb_sp = get_sb_profile(cumul, cumul_sp, nprof + 1);
/*
 * find Petrosian ratio as a function of radius, and the error in the ratio. 
 *
 * Note that if the object is only marginally detected, and if the
 * sky is overestimated, the cumulative profile can go through zero
 */
   skyVar = obj1->sky/gain + fiparams->frame[color].dark_variance;
   ratio[0] = 1; ratioErr[0] = 0;
   
   nratio = (ngood < ngood_ref) ? ngood_ref : ngood;
   for(i = 1;i <= nratio;i++) {
      float a, b;
      a = 0.8;
      b = (i == nratio) ? 1 : 1.25;	/* 1.25*rr[nratio] would extrapolate */
      rr[0] = asinh_ph(a*r[i]);
      rr[1] = asinh_ph(r[i]);
      rr[2] = asinh_ph(b*r[i]);
      phSplineInterpolate(cumul_sp,rr,vals,3);
      vals[0] = sinh_ph(vals[0]);
      vals[1] = sinh_ph(vals[1]);
      vals[2] = sinh_ph(vals[2]);
      
      if(vals[1] == 0) {
	 ratio[i] = ratio[i - 1];
	 if(!petro_gcv_errors) {
	    ratioErr[i] = 1e10;
	 }
      } else {
	 const float pi_r2 = M_PI*pow(r[i],2);
	 
	 ratio[i] = (vals[2] - vals[0])/(vals[1]*(b*b - a*a));
	 if(!petro_gcv_errors) {
/*
 * Now the error in ratio. First the photon-noise contribution, but if
 * the fractional error in the nearest radial profile point is larger
 * than that, adopt that fractional error for the Petrosian ratio too
 */
	    if(vals[2] - vals[0] == 0) {
	       ratioErr[i] = 1/(vals[1]*(b*b - a*a))*
						sqrt(pi_r2*skyVar*(b*b - a*a));
	    } else {
	       double tmp =
		 (fabs(vals[2] - vals[0])/gain + pi_r2*skyVar*(b*b - a*a))/
		   pow(vals[2] - vals[0], 2) -
		     2*(fabs(vals[1] - vals[0])/gain + pi_r2*skyVar*(1 - a*a))/
					    (vals[1]*fabs(vals[2] - vals[0])) +
			   (fabs(vals[1])/gain + pi_r2*skyVar)/pow(vals[1], 2);

	       if(tmp < 0 && tmp >= -1e-6) { /* allow a little slop */
		  tmp = 0;
	       }
	       
	       if(tmp < 0) {
		  shError("Variance of Petrosian Radius = %g < 0",tmp);
		  ratioErr[i] = 1e10;
	       } else {
		  ratioErr[i] = fabs(ratio[i])*sqrt(tmp);
	       }
	    }
	    
	    if(ratioErr[i] < 1e-10) {	/* avoid possible division by 0 */
	       ratioErr[i] = 1e-10;
	    }
	    
	    if(ratio[i]*obj1->profErr[i - 1] > ratioErr[i]*prof[i - 1]) {
	       if(prof[i - 1] > 1e-3) {
		  ratioErr[i] = ratio[i]*obj1->profErr[i - 1]/prof[i - 1];
	       }
	    }
	 }
      }
   }
/*
 * The interpolation scheme in use for the cumulative profile has problems
 * near the centre; for a constant surface brightness profile, the
 * Petrosian ratio is measured to be:
 *    radius/pixels:   0.564, 1.693, 2.585, 4.406, 7.506
 *    Petrosian ratio: 0.945, 1.016, 1.001, 0.999, 1.000
 * The slightly low 0.56pixel value could confuse the spline fitter, so we
 * silently adjust it; similarly, we lower ratio[2] a little --- note that
 * we are talking about radii of under an arcsecond (rr[2] ~ 0.7arcsec)
 */
   ratio[1] /= 0.945; ratioErr[1] /= 0.945;
   if(ratio[1] > 1.0) {
      ratio[1] = 1.0;
   }
   ratio[2] /= 1.016; ratioErr[2] /= 1.016;
#if 0
/*
 * smooth the Petrosian ratio
 */
   if(smooth_petro_radii) {
      float chisq = 0;			/* desired chi^2 */
      SPLINE *ratio_sp;			/* (smoothed) spline of Petro. ratio */
      
      ratio_sp = phSplineFindSmoothed(as_r, ratio, ratioErr,
				      nratio + 1, 1, &chisq, NULL);
      for(i = 0; i <= nratio; i++) {
	 ratio[i] = ratio_sp->coeffs[0][i];
      }
      
      phSplineDel(ratio_sp);
   }
#endif
   if(petro_gcv_errors) {		/* calculate errors by generalised
					   cross variance; not actually
					   smoothing as var == 0 */
      SPLINE *ratio_sp;			/* spline of Petro. ratio (not used) */
      
      for(i = 0; i <= nratio; i++) {
	 ratioErr[i] = 1;
      }
      ratio_sp = phSplineFindSmoothed(as_r, ratio, ratioErr,
						nratio + 1, 0, NULL, ratioErr);
      phSplineDel(ratio_sp);
   }
/*
 * if desired, save Petrosian ratio in the TEST_INFO
 */
#if SAVE_PETRO
   if(objc->test == NULL) objc->test = phTestInfoNew(objc->ncolor);
   
   objc->test->nPetroRatio[color] = nratio;
   memcpy(objc->test->petroRatio[color], ratio, nratio*sizeof(float));
   memcpy(objc->test->petroRatioErr[color], ratioErr, nratio*sizeof(float));
#endif
/*
 * We have the Petrosian Ratio, but for measuring fluxes we must use
 * an unsmoothed profile. If we smoothed before, get an unsmoothed
 * version now.
 */
   if(smooth_petro_radii) {
      phSplineDel(cumul_sp);

      get_cumul(prof, NULL, obj1->nprof, cumul, 0);
      for(i = obj1->nprof;i < ngood_ref;i++) {
	 cumul[i + 1] = cumul[obj1->nprof];
      }

      for(i = 0;i <= nprof;i++) {
	 as_cumul[i] = asinh_ph(cumul[i]);
      }
      cumul_sp = phSplineFindTaut(as_r,as_cumul,nprof + 1,Gamma);
      shAssert(cumul_sp != NULL);
   }
/*
 * Now spline the Petrosian ratio, and look for all solutions to ratio[] = f_1
 *
 * pf2 is the counts/pixel^2 corresponding to petro_f2 mag/arcsec^2
 */
   pf2 = fiparams->frame[color].sb_counts*
		    pow(10,0.4*(fiparams->frame[color].sb_mag - fiparams->petro_f2));
   
   obj1->petroRad = get_petrosianRad(fiparams->petro_f1, pf2, sb_sp,
				     ratio, ngood, nratio, &i);
   obj1->flags |= i;

   if(!(obj1->flags & OBJECT1_NOPETRO)) { /* success; estimate the error */
      int flags_p, flags_m;		/* flags when finding rP_[pm] */
      float rP_p, rP_m;			/* Petrosian radius +/- 1sigma */
      float ratio_pm[NANN + 1];		/* Petrosian ratio +/- 1sigma */

      for(i = 0;i <= nratio;i++) {
	 ratio_pm[i] = ratio[i] - ratioErr[i];
      }
      rP_m = get_petrosianRad(fiparams->petro_f1, pf2, sb_sp,
			      ratio_pm, ngood, nratio, &flags_m);
      
      for(i = 0;i <= nratio;i++) {
	 ratio_pm[i] = ratio[i] + ratioErr[i];
      }
      rP_p = get_petrosianRad(fiparams->petro_f1, pf2, sb_sp,
			      ratio_pm, ngood, nratio, &flags_p);

      if(flags_m & OBJECT1_NOPETRO) {	/* rP_m failed */
	 if(flags_p & OBJECT1_NOPETRO) {
	    obj1->petroRadErr = -1000;	/* ... and so did rP_p */
	 } else {
	    obj1->petroRadErr = rP_p - obj1->petroRad;
	 }
      } else {
	 if(flags_p & OBJECT1_NOPETRO) { /* rP_p failed */
	    obj1->petroRadErr = obj1->petroRad - rP_m;
	 } else {
	    obj1->petroRadErr = 0.5*fabs(rP_p - rP_m);
	 }
      }
   } else {				/* failed to find petroRad */
      if(obj1->flags & OBJECT1_NOPETRO_BIG) {
	 obj1->petroRad = r[obj1->nprof];
	 obj1->petroRadErr = -1000;
      } else {
	 obj1->petroRad = petro_f5_pix;
	 obj1->petroRadErr = -1000;
      }
   }
/*
 * Find the Petrosian flux, and 50/90% light radii. Note that we arranged
 * matters so that the canonical colour was processed first.
 *
 * The error in the Petrosian flux is compounded of the photon noise
 * within pf4*r_P and the uncertainty in r_P; these are taken to be independent
 * and thus added in quadrature.
 *
 * Because fiparams->petro_f4 is (almost certainly) greater than one,
 * the region used to determine the Petrosian flux may well extend
 * beyond the region for which we extracted a profile. If this is true,
 * we could do one of three things: enlarge the area extracted; extrapolate
 * the profile; or use the flux within the measured part of the profile.
 *
 * Because we've already measured the profile as far out as we think possible,
 * and because the safe way to extrapolate the flux would be to set the surface
 * brightness to zero, the first method is ruled out, and the two latter
 * suggestions are equivalent. We therefore follow the last suggestion
 */
   r_P = objc->color[fiparams->ref_band_index]->petroRad;
   r_PErr = objc->color[fiparams->ref_band_index]->petroRadErr;
   r_P_p4 = r_P*fiparams->petro_f4;

   if(r_P_p4 > sinh_ph(cumul_sp->knots[cumul_sp->nknot - 1])) {
      r_P_p4 = sinh_ph(cumul_sp->knots[cumul_sp->nknot - 1]);
   }
   rr[0] = asinh_ph(r_P_p4);
   phSplineInterpolate(cumul_sp,rr,vals,1);
   obj1->petroCounts = sinh_ph(vals[0]);
   obj1->petroCountsErr =
     calc_photom_sigma(obj1->petroCounts, M_PI*r_P_p4*r_P_p4, obj1,
						      &fiparams->frame[color]);

   if(r_PErr > 0) {			/* include error due to error in r_P */
      rr[0] = asinh_ph(r_P_p4 - fiparams->petro_f4*r_PErr);
      if(rr[0] < 0) {
	 rr[0] = 0;
      }
      phSplineInterpolate(cumul_sp,rr,vals,1);
      vals[0] = sinh_ph(vals[0]);
      
      rr[0] = asinh_ph(fiparams->petro_f4*(r_P + r_PErr));
      if(rr[0] > cumul_sp->knots[cumul_sp->nknot - 1]) {
	 rr[0] = cumul_sp->knots[cumul_sp->nknot - 1];
      }
      phSplineInterpolate(cumul_sp,rr,&vals[1],1);
      vals[1] = sinh_ph(vals[1]);
      
      obj1->petroCountsErr = sqrt(pow(obj1->petroCountsErr,2) +
				  pow(0.5*(vals[1] - vals[0]),2));
   }
/*
 * 50/90% radii next
 */
   if(0.5*obj1->petroCounts < cumul[0]) {
      obj1->petroR50 = 0;
      obj1->petroR50Err = -1000;
   } else {
      i = phSplineRoots(cumul_sp, asinh_ph(0.5*obj1->petroCounts),
			0, asinh_ph(r_P_p4), rr, 3);
      if(i == 0) { rr[0] = 0; i = 1; }	/* XXX */
      shAssert(i != 0);
      if(i != 1) {
	 obj1->flags |= OBJECT1_MANYR50;
#if 0
	 shError("More than one 50%% light radius for OBJECT1 %d",obj1->id);
#endif
      }
      obj1->petroR50 = sinh_ph(rr[0]);
      obj1->petroR50Err = petro_rNN_error(cumul_sp, r_P_p4, 0.50,
					  obj1->petroCounts,
					  obj1->petroCountsErr);
   }
   
   if(0.9*obj1->petroCounts < cumul[0]) {
      obj1->petroR90 = 0;
      obj1->petroR90Err = -1000;
   } else {
      i = phSplineRoots(cumul_sp, asinh_ph(0.9*obj1->petroCounts),
			0, asinh_ph(r_P_p4), rr, 3);
      if(i == 0) { rr[0] = 0; i = 1; }	/* XXX */
      shAssert(i != 0);
      if(i != 1) {
	 obj1->flags |= OBJECT1_MANYR90;
#if 0
	 shError("More than one 90%% light radius for OBJECT1 %d",obj1->id);
#endif
      }
      obj1->petroR90 = sinh_ph(rr[0]);
      obj1->petroR90Err = petro_rNN_error(cumul_sp, r_P_p4, 0.90,
					  obj1->petroCounts,
					  obj1->petroCountsErr);
   }
/*
 * clean up
 */
   phSplineDel(sb_sp);

   return(cumul_sp);
}

/*****************************************************************************/
/*
 * See note above p_phEllipseFit for why these are global
 */
#define NPARAM 5			/* number of parameters to fit */
static MAT *A = NULL;			/* The LSQ problem is */
static VEC *b = NULL;			/*    A x = b */
static MAT *Q = NULL;			/* eigen value */
static VEC *l = NULL;			/*             decomposition of A */

void
p_phInitEllipseFit(void)
{
   A = phMatNew(NPARAM, NPARAM);
   Q = phMatNew(NPARAM, NPARAM);
   b = phVecNew(NPARAM);
   l = phVecNew(NPARAM);

   (void)p_phEllipseFit(NULL, NULL, 1);	/* initialise the trig arrays */
}

void
p_phFiniEllipseFit(void)
{
   phMatDel(A); A = NULL;
   phVecDel(b); b = NULL;
   phMatDel(Q); Q = NULL;
   phVecDel(l); l = NULL;
}

/*
 * <AUTO EXTRACT>
 *
 * This function isn't really meant to be a global at all, except that
 * I wanted to be able to call it seperately from measure objects, in
 * order to quantify the bias introduced by the use of NSEC sectors
 * in estimating the r[] values.
 *
 * Given the radial distances to a curve in the directions i*(2*pi/NSEC),
 * and the errors in those distances, return a VEC whose components are
 * respectively the centre, the major and minor axis, and the p.a., measured
 * in a positive direction from "north".
 */
VEC *
p_phEllipseFit(const float *r,		/* NSEC radii to be fitted */
	       const float *rErr,	/* and their errors NOTUSED */
	       int debias)		/* debias a and b? */

{
   float b0, b1, b2, a1, a2;		/* unpacked Fourier coeffs */
   VEC *coeffs;				/* desired coefficients */
   static float cost[4*NSEC], sint[4*NSEC]; /* == {cos,sin}(theta_i) */
   int i;
   int ia;				/* index into trig arrays */
   static float l_min = 1e-6;		/* minimum acceptable eigen value */

   float phi;				/* position angle of major axis */
   float mcc, mcr, mrr;			/* second moments of fit, in row,col */
   float mxx, myy;			/* second moments in diagonal frame */
   float ri;				/* == r[i] */
   float sum00, sum01, sum02, sum03, sum04, /* sums to fill out A and b */
                sum11, sum12, sum13, sum14,
                       sum22, sum23, sum24,
                              sum33, sum34,
                                     sum44;
/*
 * setup the trig arrays if r == NULL.
 * They are calculated every 360/(2*NSEC) degrees (i.e. every _half_ sector)
 */
   if(r == NULL) {
      float theta;
      for(i = 0; i < 4*NSEC; i++) {
	 theta = M_PI*i/(float)NSEC;
	 cost[i] = cos(theta); sint[i] = sin(theta);
      }
      return(NULL);
   }
/*
 * fill out LHS of normal equations
 */
   shAssert(A != NULL && NPARAM == 5);
   sum00 = sum01 = sum02 = sum03 = sum04 = 0;
   sum11 = sum12 = sum13 = sum14 = 0;
   sum22 = sum23 = sum24 = 0;
   sum33 = sum34 = 0;
   sum44 = 0;
   for(i = 0; i < NSEC; i++) {
      if(r[i] < 0) continue;

      ia = 2*i;				/* index into trig arrays */

      sum00++;
      sum01 += cost[ia]; 
      sum02 += sint[ia];
      sum03 += cost[2*ia]; 
      sum04 += sint[2*ia];

      sum11 += cost[ia]*cost[ia];
      sum12 += cost[ia]*sint[ia];
      sum13 += cost[ia]*cost[2*ia]; 
      sum14 += cost[ia]*sint[2*ia];

      sum22 += sint[ia]*sint[ia];
      sum23 += sint[ia]*cost[2*ia]; 
      sum24 += sint[ia]*sint[2*ia];

      sum33 += cost[2*ia]*cost[2*ia]; 
      sum34 += cost[2*ia]*sint[2*ia];

      sum44 += sint[2*ia]*sint[2*ia];
   }

   if(sum00 < 5) {			/* too few points to fit */
      return(NULL);
   }

   A->me[0][0] = sum00;
   A->me[0][1] = A->me[1][0] = sum01;
   A->me[0][2] = A->me[2][0] = sum02;
   A->me[0][3] = A->me[3][0] = sum03;
   A->me[0][4] = A->me[4][0] = sum04;

   A->me[1][1] = sum11;
   A->me[1][2] = A->me[2][1] = sum12;
   A->me[1][3] = A->me[3][1] = sum13;
   A->me[1][4] = A->me[4][1] = sum14;

   A->me[2][2] = sum22;
   A->me[2][3] = A->me[3][2] = sum23;
   A->me[2][4] = A->me[4][2] = sum24;

   A->me[3][3] = sum33;
   A->me[3][4] = A->me[4][3] = sum34;

   A->me[4][4] = sum44;
/*
 * and now the RHS
 */
   sum00 = sum01 = sum02 = sum03 = sum04 = 0;
   for(i = 0; i < NSEC; i++) {
      if((ri = r[i]) < 0) continue;
      ri = sqrt(ri);			/* make ellipse flatter; after the fit
					   is squared, the radii will be >= 0 */

      ia = 2*i;				/* index into trig arrays */

      sum00 += ri;
      sum01 += ri*cost[ia]; 
      sum02 += ri*sint[ia];
      sum03 += ri*cost[2*ia]; 
      sum04 += ri*sint[2*ia];
   }

   b->ve[0] = sum00;
   b->ve[1] = sum01;
   b->ve[2] = sum02;
   b->ve[3] = sum03;
   b->ve[4] = sum04;
/*
 * OK, time to solve for the coefficients. We'll do this by eigen value
 * decomposition, as the matrix can be singular (e.g. if the points lie
 * on a circle)
 */
   (void)phEigen(A, Q, l);

   for(i = 0; i < NPARAM; i++) {
      if(l->ve[i] < l_min) {
	 l->ve[i] = 0;
      }
   }
   coeffs = phEigenBackSub(Q, l, b);
/*
 * Convert those Fourier coefficients to a description of an ellipse
 */
   b0 = coeffs->ve[0];
   b1 = coeffs->ve[1];
   a1 = coeffs->ve[2];
   b2 = coeffs->ve[3];
   a2 = coeffs->ve[4];
/*
 * find the second moments of e.g. (col - <col>)^2.
 */
   mcc = 0.125*(4*b0*b0 + a1*a1 + b1*b1 + 4*b0*b2 + 2*b2*b2 + 2*a2*a2);
   mcr = 0.5*b0*a2;
   mrr = 0.125*(4*b0*b0 + a1*a1 + b1*b1 - 4*b0*b2 + 2*b2*b2 + 2*a2*a2);
/*
 * Find phi that maximises mxx
 */
   if(mcr == 0 && mrr - mcc == 0) {
      phi = 0;
   } else {
      phi = 0.5*atan2(2*mcr, mrr - mcc);
   }

   mxx = mrr + (mcc - mrr)*pow(cos(phi),2) - mcr*sin(2*phi);
   myy = (mrr + mcc) - mxx;
   if(myy > mxx) {			/* phi is p.a. of minor axis */
      float tmp = mxx; mxx = myy; myy = tmp;
      phi += M_PI/2;
   }
/*
 * convert those values of <r^2> to axial ratios. Note that for an ellipse
 * aligned with the coordinate axes,
 *  <x^2> = a^2b/(a + b),  <y^2> = ab^2/(a + b)
 * so a^2 = <x^2>/<y^2>*(<x^2> + <y^2>)
 *
 * Also remember that we took the square root of the radii, so we don't need a
 * square root in the axis lengths
 */
   coeffs->ve[0] = a1/2;
   coeffs->ve[1] = b1/2;
   coeffs->ve[2] = mxx/myy*(mxx + myy);
   coeffs->ve[3] = myy/mxx*(mxx + myy);
   coeffs->ve[4] = (M_PI/2 - phi)*180/M_PI;

   while(coeffs->ve[4] < -90) {		/* reduce to range [-90,90] */
      coeffs->ve[4] += 180;
   }
   while(coeffs->ve[4] > 90) {
      coeffs->ve[4] -= 180;
   }
/*
 * if desired, correct for the effects of the sector arrays on the measured
 * ellipse parameters
 */
   if(debias) {
      p_phEllipseFitDebias(coeffs);
   }
#if 0
/*
 * It's more robust if we require that the area of the fitted ellipse equals
 * that of the area included in the sectors within distance r[]; so fiddle
 * the axes accordingly
 */
   sum00 = sum01 = 0;
   for(i = 0; i < NSEC; i++) {
      if(r[i] < 0) continue;

      sum00++; sum01 += r[i]*r[i];
   }
   shAssert(sum00 > 0);
   {
      float area = sum01/sum00;		/* area/pi */
      float fac = sqrt(area/(coeffs->ve[2]*coeffs->ve[3])); /* correction */
      coeffs->ve[2] *= fac;
      coeffs->ve[3] *= fac;
   }
#endif

   return(coeffs);
}

/*****************************************************************************/
/*
 * See note above p_phEllipseFit for why this is global
 *
 * Given a measured set of coefficients, correct for the biases introduced
 * by the use of sector arrays. These coefficients were arrived at by looking
 * at the outputs of the stand alone programme, find_iso_bias
 */
void
p_phEllipseFitDebias(VEC *coeffs)
{
   float a0, a1, b0, b1, a2, b2;	/* coeffs for corrections for a and b */
   float arm;				/* measured axis ratio */
   float phim, phit;			/* measured and true position angles */
   float p0;				/* coeffs for phi's correction */

   shAssert(NSEC == 12);		/* these are magic coefficients! */

   arm = coeffs->ve[3]/coeffs->ve[2];

   if(arm < 0.03) {			/* (arm-0.03)^3/4 == -NaN */
      return;				/* cannot debias */
   }

   phim = coeffs->ve[4];
   phim = (phim + 360) - 30*((int)(phim + 360)/30); /* reduce to [0,30) */
   if(arm < 0.4) {
      p0 = -23.8280 + 108.634*sqrt(arm) + arm*(-141.275 + arm*72.6992);
   } else {
      p0 = 0;
   }

   phit = phim + p0*sin(M_PI*(phim/15 - 1));
   shAssert(phit >= -0.001 && phit <= 30.001);
   
   if(arm < 0.3) {
      a0 = 1/(-0.0651703 + arm*(2.95038 + arm*(6.47144 -arm*22.1218)));
   } else {
      a0 = 1 + (1 - arm)*(-0.0767746 + (1 - arm)*(1.17578 +
				     (1 - arm)*(-3.45642 + (1 - arm)*3.76807)));
   }
   if(arm < 0.3) {
      a1 = 1/(0.794716 + arm*(-42.2788 + arm*(497.617 - arm*2067.24)));
   } else if(arm <= 0.7) {
      a1 = -0.180717 + arm*(0.862976 + arm*(-1.37524 + arm*0.730316));
   } else {
      a1 = 0;
   }

   if(arm < 0.125) {
      b0 = 1.64243 + 25.4849*pow(arm - 0.03, 0.75) +
						   arm*(-44.9199 + arm*41.8945);
   } else {
      b0 = 0.999866 + (1 - arm)*(-0.0163803 +
				   (1 - arm)*(0.0110474 + (1 - arm)*0.0772781));
   }
   if(arm <= 0.05) {
      b1 = -0.303833 - 4.41992*pow(arm - 0.03, 0.75) +
						     arm*(6.3125 + arm*28.0469);
   } else if(arm <= 0.3) {
      b1 = -0.234013 + arm*(1.6137 - arm*2.82491);
   } else {
      b1 = 0;
   }

   if(arm <= 0.4) {
      a2 = b2 = -2.578 + arm*(106.662 + arm*(-279.419 + arm*244.074));
   } else {
      a2 = b2 = 11;
   }

   coeffs->ve[2] *= a0 + a1/(1 + pow((phit - 15)/a2, 2));
   coeffs->ve[3] *= b0 + b1/(1 + pow((phit - 15)/b2, 2));
   coeffs->ve[4] += phit - phim;
}

/*****************************************************************************/
/*
 * find the shape of an isophote at a given surface brightness
 */
#define LG_3 0.4771212547		/* log_10(3) */

static void
calc_iso_ellipse(OBJC *objc,		/* object in question */
		 int c,			/* in this band */
		 const CELL_STATS *cellprof, /* extracted profile */
		 const FIELDPARAMS *fiparams) /* gain etc. */
{
   float as_rmax;			/* asinh(max radius for iso_cts) */
   float ar;				/* asinh(r) */
   VEC *coeffs[3];			/* desired coefficients */
   float cumul[NANN + 1];		/* cumulation of prof[] */
   float as_cumul[NANN + 1];		/* asinh(cumul[]) */
   const struct pstats *cell;		/* a cell in cellprof */
   SPLINE *cumul_sp;			/* spline of as_cumul[] */
   const float iso_cts = fiparams->frame[c].sb_counts; /* counts/pixel at
							  desired isophote */
   int i;
   const int nprof = cellprof->nannuli_c; /* number of radial points */
   OBJECT1 *const obj1 = objc->color[c]; /* the OBJECT1 in question */
   float prof[NANN];			/* radial profile in a sector */
   float profErr[NANN];			/* error in prof[] */
   float r_iso[3][NSEC];		/* radii of desired isophotes,
					   {0,1,2} --> iso(1 + {0.5,0,-0.5}) */
   SPLINE *sb_sp;			/* spline of surface brightness prof */
   int sec;				/* counter for sectors */
/*
 * first handle the central cell
 */
   cell = &cellprof->cells[0];
   prof[0] = median_profs ? cell->qt[1] : cell->mean;
   profErr[0] = cell->sig;
   as_rmax = as_r[nprof];

   if(prof[0] < iso_cts) {
      obj1->flags |= OBJECT1_ELLIPFAINT;
      return;
   }

   for(sec = 0; sec < NSEC; sec++) {
/*
 * Find the (possibly smoothed) radial profile in sector
 */
      for(i = 1; i < nprof; i++) {
	 cell = &cellprof->cells[1 + (i - 1)*NSEC + sec];
	 prof[i] = median_profs ? cell->qt[1] : cell->mean;
	 profErr[i] = cell->sig;
      }
      get_cumul(prof, profErr, nprof, cumul, fiparams->smooth_profs);
/*
 * Spline the cumulative profile. See notes in calc_petrosian() for
 * choice of boundary conditions
 */   
      for(i = 0;i <= nprof;i++) {
	 as_cumul[i] = asinh_ph(cumul[i]);
      }

      cumul_sp = phSplineFindTaut(as_r, as_cumul, nprof + 1, Gamma);
      shAssert(cumul_sp != NULL);
      
      sb_sp = get_sb_profile(cumul, cumul_sp, nprof + 1);
/*
 * find the radius corresponding to {1.5, 1, 0.5}iso_cts
 */
      for(i = 0; i <= 2; i++) {
	 if(phSplineRoots(sb_sp, asinh_ph((1 + 0.5*(1 - i))*iso_cts),
						    0, as_rmax, &ar, 1) <= 0) {
	    ar = -1;
	 } else {
	    ar = sinh_ph(ar);
	 }
	 r_iso[i][sec] = ar;
      }
      phSplineDel(cumul_sp); phSplineDel(sb_sp);
   }
/*
 * fit an Fourier series to each of those set of radii, and hence
 * derive the centre, axes, and p.a. of the best-fit ellipse
 */
   for(i = 0; i <= 2; i++) {
      if((coeffs[i] = p_phEllipseFit(r_iso[i], NULL, 1)) == NULL) {
	 obj1->flags |= OBJECT1_ELLIPFAINT;

	 while(--i >= 0) {
	    phVecDel(coeffs[i]);
	 }
	 return;
      }
   }

   obj1->iso_rowc = obj1->rowc + coeffs[1]->ve[0];
   obj1->iso_rowcErr = -1000;
   obj1->iso_colc = obj1->colc + coeffs[1]->ve[1];
   obj1->iso_colcErr = -1000;
   obj1->iso_a = coeffs[1]->ve[2];
   obj1->iso_aErr = -1000;
   obj1->iso_b = coeffs[1]->ve[3];
   obj1->iso_bErr = -1000;
   obj1->iso_phi = coeffs[1]->ve[4];
   obj1->iso_phiErr = -1000;

   obj1->iso_rowcGrad = (coeffs[2]->ve[0] - coeffs[1]->ve[0])/(2.5*LG_3);
   obj1->iso_colcGrad = (coeffs[2]->ve[1] - coeffs[1]->ve[1])/(2.5*LG_3);
   obj1->iso_aGrad = (coeffs[2]->ve[2] - coeffs[1]->ve[2])/(2.5*LG_3);
   obj1->iso_bGrad = (coeffs[2]->ve[3] - coeffs[1]->ve[3])/(2.5*LG_3);
   obj1->iso_phiGrad = (coeffs[2]->ve[4] - coeffs[1]->ve[4])/(2.5*LG_3);

   for(i = 0; i <= 2; i++) {
      phVecDel(coeffs[i]);
   }
}

/*****************************************************************************/
/*
 * Calculate an OBJC's shape parameters, the Stokes parameters Q and U, 
 * and the derived quantities (a-b)/(a+b) and position angle of major axis
 */
static void
calc_shape_old(OBJC *objc,
	       int color,
	       const CELL_STATS *prof,
	       const FIELDPARAMS *fiparams,
	       const SPLINE *cumul_sp)
{
   float petroR, as_petroR;		/* obj1->petroRad, asinh(petroR) */
   float cts;				/* the number of counts in a cell */
   int i;
   OBJECT1 *const obj1 = objc->color[color]; /* the OBJECT1 in question */
   int rowc, colc;			/* centre of object */
   int ret;				/* a return code */
   double sum, sum_Q, sum_U;
   SPLINE *stokes_sp;			/* spline of Q[] and then U[] */
   float val;				/* a returned value */
   const REGION *const syncreg = prof->syncreg;
   const struct pstats *cells = prof->cells; /* unpack the prof */
   const struct cellmod *mgeom = prof->mgeom;
   float Q[NANN], U[NANN];		/* cumulative Stokes parameters */

   rowc = syncreg->nrow/2; colc = syncreg->ncol/2;

   sum = sum_Q = sum_U = 0.0;
   for(i = 0;i < prof->nannuli_c;i++) {
      if(phStokesParamEval(syncreg,i,SOFT_BIAS,rowc,colc,'Q',0,&val) == 0) {
	 sum_Q += val;
	 cts = area[i]*obj1->profMean[i];
	 sum += cts;
	 ret = phStokesParamEval(syncreg,i,SOFT_BIAS,rowc,colc,'U',0,&val);
	 shAssert(ret == 0);
	 sum_U += val;
      } else {
	 int j;
	 float dsum = 0, dQ = 0, dU = 0; /* changes to sums due to this
					    annulus */
	 float sum_med = 0;		/* sum of medians */
	 float med;			/* median level in a cell */

	 for(j = NSEC*(i - 1) + 1;j <= NSEC*i;j++) {
	    if(i == 0) {		/* central pixel */
	       j = 0;
	    }
	    med = cells[j].qt[1];
	    cts = cells[j].nel*med;

	    dsum += cts;
	    sum_med += med;
	    dQ += med*mgeom[j].Q;
	    dU += med*mgeom[j].U;
	    if(i == 0) {		/* only one cell in central "annulus"*/
	       break;
	    }
	 }
	 sum += dsum;
	 if(sum_med == 0) {		/* so are dQ and dU */
	    ;
	 } else {
	    sum_Q += dsum*dQ/sum_med;
	    sum_U += dsum*dU/sum_med;
	 }
      }
      if(sum == 0) {
	 Q[i] = U[i] = 0;
      } else {
	 Q[i] = sum_Q/sum;
	 U[i] = sum_U/sum;
      }
   }
/*
 * interpolate U and Q to the Petrosian Radius
 */
   petroR = objc->color[fiparams->ref_band_index]->petroRad;
   as_petroR = asinh_ph(petroR);

   stokes_sp = phSplineFindTautEven(&as_r[1],Q,prof->nannuli_c,Gamma);
   shAssert(stokes_sp != NULL);
   phSplineInterpolate(stokes_sp,&as_petroR,&obj1->Q,1);
   phSplineDel(stokes_sp);

   stokes_sp = phSplineFindTautEven(&as_r[1],U,prof->nannuli_c,Gamma);
   shAssert(stokes_sp != NULL);
   phSplineInterpolate(stokes_sp,&as_petroR,&obj1->U,1);
   phSplineDel(stokes_sp);
/*
 * now do the error analysis and debiasing. We need to know the
 * signal-to-noise ratio for the flux within petroRad.
 *
 * In calculating the variances, we use the formulae appropriate to the
 * case where the sky noise dominates.
 */
   {
      float majaxis, P;			/* initial estimates of
					   phi and (a-b)/(a+b) */
      float petroFlux;			/* flux within 1*petroRad */
      float petroFluxErr;		/* Error in flux within 1*petroRad */
      float varQ, varU;			/* variances of Q and U in the frame
					   aligned with the principal axes */

      shAssert(cumul_sp != NULL);
      
      phSplineInterpolate(cumul_sp,&as_petroR,&petroFlux,1);
      petroFlux = sinh_ph(petroFlux);
      petroFluxErr = calc_photom_sigma(petroFlux,M_PI*petroR*petroR,obj1,
						      &fiparams->frame[color]);
      majaxis = 0.5*atan2(obj1->U,obj1->Q);	/* initial estimates */
      P = sqrt(pow(obj1->Q,2) + pow(obj1->U,2));
      
      if(petroFlux < 1e-5) {
	 varQ = varU = 1e10;
      } else if(petroFluxErr > 10*petroFlux) { /* don't debias if SN < 0.1 */
	 varU = pow(petroFluxErr/petroFlux,2)/2; /*assume sky noise dominates*/
	 varQ = (1 + 2*P*P)*varU;
      } else {
	 float itau = petroFluxErr/petroFlux;
	 
#if 0					/* debias */
	 obj1->Q -= P*pow(itau,2)*cos(2*majaxis);
	 obj1->U -= P*pow(itau,2)*sin(2*majaxis);
#endif
	 
	 majaxis = 0.5*atan2(obj1->U,obj1->Q); /* debiased estimates */
	 if(majaxis != majaxis) {		/* NaN, i.e. {U,Q} == 0 */
	    obj1->Q = obj1->U = 0;
	    majaxis = 0;
	 }
	 P = sqrt(pow(obj1->Q,2) + pow(obj1->U,2));

	 varU = pow(itau,2)/2;		/* assume that sky noise dominates */
	 varQ = (1 + 2*P*P)*varU;
      }

      if(varQ < 0.99e5) {
	 obj1->QErr = sqrt(pow(cos(2*majaxis),2)*varQ +
						   pow(sin(2*majaxis),2)*varU);
	 obj1->UErr = sqrt(pow(sin(2*majaxis),2)*varQ +
						   pow(cos(2*majaxis),2)*varU);
      }
   
      obj1->majaxis = majaxis;
      obj1->aratio = P;
   }
/*
 * check for NaNs
 */
   if(obj1->Q != obj1->Q || obj1->U != obj1->U) {
      obj1->Q = obj1->U = -1000;
      obj1->QErr = obj1->UErr = -9999;
      obj1->flags |= OBJECT1_NOSTOKES;
   }
}

/*****************************************************************************/
/*
 * Fit the centre of the object with an expansion in powers of the PSF
 */
#if FIT_PSF_EXPANSION
static void
fit_psf_expansion(OBJC *objc,		/* the object under consideration */
		  int c,		/* the colour of interest */
		  const CELL_STATS *prof, /* object's radial profile */
		  const FIELDPARAMS *fiparams) /* describe field */
{
   MAT *Apsf;				/* The LSQ problem is */
   VEC *bpsf;				/*    A x = b */
   MAT *Qpsf;				/* eigen value */
   VEC *lpsf;				/*             decomposition of A */
   VEC *x;				/* desired coefficients */
   const PSF_BASIS *basis = fiparams->frame[c].psfBasis; /* the PSF */
   double dot;				/* an inner products */
   int i, j;
   PSF_REG *regs[4];			/* the PSF and its derivative */
   REGION *sregs[4];			/* centre of each regs[]->reg */
   const REGION *sincreg;		/* object's sinc-shifted region */
   TEST_INFO *test = objc->test;	/* unpacked for convenience */

   if(basis == NULL) {
      return;
   }
/*
 * find the PSF and its second moments at the position of object
 */
   regs[0] = phPsfKLReconstruct(basis,
				objc->color[c]->rowc, objc->color[c]->colc,
								    TYPE_FL32);
   regs[1] = phPsfKLReconstruct(basis->deriv[DBASIS_DROW2],
				 objc->color[c]->rowc, objc->color[c]->colc,
								    TYPE_FL32);
   regs[2] = phPsfKLReconstruct(basis->deriv[DBASIS_DROWDCOL],
				objc->color[c]->rowc, objc->color[c]->colc,
								    TYPE_FL32);
   regs[3] = phPsfKLReconstruct(basis->deriv[DBASIS_DCOL2],
				objc->color[c]->rowc, objc->color[c]->colc,
								    TYPE_FL32);
/*
 * We need a subregion of each of those the same size as the sinc region
 */
   sincreg = prof->syncreg;
   for(i = 0; i < 4; i++) {
      sregs[i] = shSubRegNew("", regs[i]->reg, sincreg->nrow, sincreg->ncol,
			     regs[i]->reg->nrow/2 - sincreg->nrow/2,
			     regs[i]->reg->ncol/2 - sincreg->ncol/2, NO_FLAGS);
   }
/*
 * Calculate required inner products
 */
   Apsf = phMatNew(4, 4);
   Qpsf = phMatNew(4, 4);
   bpsf = phVecNew(4);
   lpsf = phVecNew(4);

   for(i = 0; i < 4; i++) {
      phRegionDotRegion(&dot, sincreg, sregs[i], 0);
      bpsf->ve[i] = dot;
      for(j = i; j < 4; j++) {
	 phRegionDotRegion(&dot, sregs[i], sregs[j], 0);
	 Apsf->me[i][j] = Apsf->me[j][i] = dot;
      }
   }
/*
 * Solve for coefficients
 */
   (void)phEigen(Apsf, Qpsf, lpsf);
   x = phEigenBackSub(Qpsf, lpsf, bpsf);

   for(i = 0; i < 4; i++) {
      test->psf_exp_coeffs[c][i] = x->ve[i];
   }
/*
 * clean up
 */
   phMatDel(Apsf);
   phVecDel(bpsf);
   phMatDel(Qpsf);
   phVecDel(lpsf);
   phVecDel(x);

   for(i = 0; i < 4; i++) {
      shRegDel(sregs[i]);
      phPsfRegDel(regs[i]);
   }
}
#endif

/*
 * calculate the shape parameters U and Q from the pixels in the object's
 * detection mask
 */
#if 0
static void
calc_shape_from_mask(OBJC *objc,
		     int color,
		     const CELL_STATS *prof,
		     const FIELDPARAMS *fiparams)
{
   float bkgd;				/* background inc. soft bias */
   float drow, dcol;			/* offsets from reference colour */
   double dsum;				/* sum of I within a span */
   const REGION *data = fiparams->frame[color].data; /* the pixel data */
   float I;				/* Intensity in a pixel above bkgd */
   int i, j;
   int ix1, ix2, iy;			/* == s[i].{x1,x2,y} */
   int n;				/* number of pixels in a span */
   int nspan;				/* number of spans in ai */
   OBJECT1 *obj1 = objc->color[color];	/* the object in question */
   U16 *pptr;				/* pointer to row in image */
   U16 **rows;				/* == data->rows */
   int r1, r2, c1, c2;			/* corners of sinc-shifted region */
   float rowc, colc;			/* canonical centre in this band */
   float Rad2;				/* == (x^2 + y^2) */
   float Rmax2;				/* (maximum allowable radius)^2 */
   SPAN *s;				/* == objmask->span */
   double sum, Usum, Qsum;		/* running sums */
   double sumIRad2;			/* tmp sum for Q within a span */
   const REGION *syncreg;		/* sinc-shifted region around object */
   float x, y;				/* centre of a pixel */
   double x2, y2;			/* == x^2, y^2 */
/*
 * find canonical centre in this bands' coordinate system
 */
   rowc = objc->rowc; colc = objc->colc;
   phOffsetDo(fiparams, rowc, colc, color, fiparams->ref_band_index,
	      NULL, NULL, &drow, NULL, &dcol, NULL);
   rowc -= drow;
   colc -= dcol;
/*
 * Figure out the coordinate patch that's in the sinc-shifted region. The
 * sync region was set for us when we extracted the radial profile
 */
   {
      const CELL_STATS *prof = phProfileGetLast();
      shAssert(prof->id == obj1->id);	/* it's this object's profile */
      syncreg = prof->syncreg;
   }
   shAssert(syncreg != NULL);
   r1 = syncreg->row0; r2 = r1 + syncreg->nrow - 1;
   c1 = syncreg->col0; c2 = c1 + syncreg->ncol - 1;
/*
 * Don't consider pixels more than the r' Petrosian radius from centre
 */
   Rmax2 = pow(objc->color[fiparams->ref_band_index]->petroRad, 2);
/*
 * Define the macro that actually increments the weighted sums 
 */
#define INCR_UQ_SUMS \
   I = *pptr++ - bkgd; \
   \
   if(Rad2 <= Rmax2 && Rad2 != 0.0) { \
      dsum += I; \
      I /= Rad2; \
      sumIRad2 += I;			/* move I*y^2/r^2 out of inner loop */\
      Usum += x*y*I; \
   } \
   \
   x++;					/* go to next pixel */ \
   Rad2 += 2*x - 1			/* (x + 1)^2 = x^2 + 2*(x+1) - 1 */
/*
 * Calculate the desired moments by going through the detection mask
 * span by span. We use a few tricks to move operations out of the inner loop.
 */
   bkgd = fiparams->frame[color].bkgd + SOFT_BIAS;
   
   nspan = obj1->mask->nspan;
   s = obj1->mask->s;
   rows = data->rows;
   
   sum = Usum = Qsum = 0.0;
   for(i = 0; i < nspan; i++) {
      iy = s[i].y; ix1 = s[i].x1; ix2 = s[i].x2;
      dsum = sumIRad2 = 0.0;

      if(iy < r1 || iy > r2) {		/* no need to worry about syncreg */
	 y = iy  + 0.5 - rowc; y2 = y*y;
	 x = ix1 + 0.5 - colc;
	 n = ix2 - ix1 + 1;
	 
	 if(y2 > Rmax2) {		/* span's too far from centre */
	    continue;
	 }

	 Rad2 = x*x + y2;
	 pptr = &rows[iy][ix1];
	 for(j = 0; j < n; j++) {
	    INCR_UQ_SUMS;
	 }
      } else {
/*
 * Span may intersect the sinc-shifted region.  Process the partial spans
 * to left and right of the sinc-region first, adjusting ix1 and ix2 as
 * we go, and then deal with the central pixels
 */
	 if(ix1 < c1) {			/* to left of syncreg */
	    y = iy  + 0.5 - rowc; y2 = y*y;
	    x = ix1 + 0.5 - colc;
	    if(ix2 < c1) {
	       n = ix2 - ix1 + 1;
	    } else {
	       n = c1 - ix1;
	    }
	    
	    Rad2 = x*x + y2;
	    pptr = &rows[iy][ix1];
	    for(j = 0; j < n; j++) {
	       INCR_UQ_SUMS;
	    }

	    if(ix2 < c1) {		/* span's complete */
	       continue;
	    }
	    ix1 = c1;
	 }

	 if(ix2 > c2) {			/* to right of syncreg */
	    const int xstart = (ix1 > c2) ? ix1 : c2 + 1;
	    y = iy + 0.5 - rowc; y2 = y*y;
	    x = xstart + 0.5 - colc;
	    n = ix2 - xstart + 1;
	    
	    Rad2 = x*x + y2;
	    pptr = &rows[iy][xstart];
	    for(j = 0; j < n; j++) {
	       INCR_UQ_SUMS;
	    }

	    if(ix1 > c2) {		/* span's complete */
	       continue;
	    }
	    ix2 = c2;
	 }

	 n = ix2 - ix1 + 1;
	 if(n > 0) {			/* syncreg itself */
	    iy -= r1; ix1 -= c1; ix2 -= c1; /* syncreg coordinate system */
	    y = iy  - syncreg->nrow/2; y2 = y*y;
	    x = ix1 - syncreg->nrow/2;
	 
	    Rad2 = x*x + y2;
	    pptr = &syncreg->rows[iy][ix1];
	    pptr = &rows[r1 + iy][c1 + ix1]; /* XXX */
	    for(j = 0; j < n; j++) {
	       INCR_UQ_SUMS;
	    }
	 }
      }
      sum += dsum;
      Qsum += dsum - 2*y2*sumIRad2;
   }
   Usum *= 2;				/* U == 2 <xy/r^2> */
/*
 * set the desired variables
 */
   if(sum == 0.0) {
      obj1->Q = obj1->U = -1000;
      obj1->QErr = obj1->UErr = -9999;
      obj1->flags |= OBJECT1_NOSTOKES;
      return;
   }
   
#if 0
   obj1->Q = Qsum/sum;
   obj1->U = Usum/sum;
#else
   obj1->QErr = Qsum/sum;		/* XXXXX */
   obj1->UErr = Usum/sum;
#endif
}
#endif

static void
calc_shape(OBJC *objc,
	   int color,
	   const CELL_STATS *prof,
	   const FIELDPARAMS *fiparams,
	   const SPLINE *cumul_sp)
{
   calc_shape_old(objc, color, prof, fiparams, cumul_sp);
#if 0
   calc_shape_from_mask(objc, color, prof, fiparams);
#endif
}

/*****************************************************************************/
/*
 * Classify an object
 */
static OBJ_TYPE
do_classify(const float *fac,		/* fiddle factors */
	    int nfac,			/* number of fiddle factors */
	    float star_L,		/* NOTUSED */
	    float exp_L, float deV_L,
	    float psfCounts, float psfCountsErr,
	    float counts_deV, float counts_deVErr,
	    float counts_exp, float counts_expErr,
	    float counts_model, float counts_modelErr) /* NOTUSED */
{
   const float gal_counts = (deV_L > exp_L) ? counts_deV : counts_exp;
   const float gal_countsErr = (deV_L > exp_L) ? counts_deVErr : counts_expErr;
   OBJ_TYPE type;

   shAssert(nfac >= 3);
   type = (fac[0]*(gal_counts + fac[1]*gal_countsErr) <
		      psfCounts + fac[2]*psfCountsErr) ? OBJ_STAR : OBJ_GALAXY;

   return(type);
}

/*****************************************************************************/
/*
 * Does object appear to be a cosmic ray?
 */
static int
classify_as_cr(OBJC *objc,
	       const FIELDPARAMS *fiparams,
	       int c)
{
   static float cond3_fac = 1.5;	/* fiddle factor for condition #3 */
   const FRAMEPARAMS *fparams = &fiparams->frame[c];
   const float gain = fparams->gain;
   int is_cr = 0;			/* object appears to be a CR */
   static int ncr_min = 3;		/* min. number of CR pixels to be CR */
   OBJECT1 *obj1 = objc->color[c];	/* the object in question */
   CELL_STATS *prof;			/* extracted profile */
   int size = 11;			/* n{row,col} of region to search */
   REGION *sreg;			/* region around centre of object */
   
#if 0
   static float cond3_fac2 = 0.6;	/* 2nd fiddle factor for condition #3*/
   static float min_sigma = 3;		/* min sig. above sky in pixel for CR*/
   static int min_e = 100;		/* min number of e^- in an CRs */
   int ncr;				/* number of cosmic rays detected */
   const float skysig = sqrt(obj1->sky/gain + fparams->dark_variance);

   sreg = shSubRegNew("", fparams->data, size, size,
		      obj1->rowc - size/2, obj1->colc - size/2, 0);
   
   shAssert(sreg != NULL);
   sreg->mask = (MASK *)phSpanmaskNew(size, size);
   
   ncr = phFindCR(sreg, obj1->sky, fparams->bkgd, skysig, fparams->psf, gain,
		  min_sigma, min_e, cond3_fac, cond3_fac2, 1);
   
   if(ncr > 0 &&
      phRectIntersectMask(((SPANMASK *)sreg->mask)->masks[S_MASK_CR],
			  size/2 - 1, size/2 - 1, size/2 + 1, size/2 + 1)) {
      is_cr = 1;
   }
   
   phSpanmaskDel((SPANMASK *)sreg->mask); sreg->mask = NULL;
   shRegDel(sreg);

   is_cr = 0;				/* XXX */
#endif

#if 1
   prof = phProfileGetLast();		/* set in phProfileExtract */
   shAssert(prof->syncreg != NULL && prof->id == obj1->id);
   
   sreg = (REGION *)prof->syncreg;
   size = sreg->nrow; shAssert(size == sreg->ncol);

   if(KLPsf[c] != NULL &&
      (obj1->flags & OBJECT1_DETECTED) &&
      !(obj1->flags2 & OBJECT2_PSF_FLUX_INTERP)) {
      float delta;			/* ratio of PSF to its central value */
      int i, j;
      float I0;				/* object's central intensity */
      float I0Err;			/* I0's standard deviation */
      float Ip;				/* object's intensity at a point */
      float IpErr;			/* Ip's standard deviation */
      const float P0 = KLPsfReg[c]->rows[size/2][size/2] - SOFT_BIAS;
      static float err_fac = 0.05;	/* min. error == pixel*err_fac */
/*
 * See if the gradients in the object are greater than those in the PSF.
 *
 * This algorithm is identical to that discussed at the top of CR.c,
 * except that we use the known PSF, and only consider pixels around
 * the known centre of the object.
 */
      I0 = sreg->rows[size/2][size/2] - fparams->bkgd - SOFT_BIAS;
      I0Err = sqrt((I0 + obj1->sky)/gain + pow(err_fac*fabs(I0),2));
      
      for(i = size/2 - 1; i <= size/2 + 1; i++) {
	 for(j = size/2 - 1; j <= size/2 + 1; j++) {
	    if(i == size/2 && j == size/2) {
	       continue;
	    }
	    
	    Ip = sreg->rows[i][j] - fparams->bkgd - SOFT_BIAS;
	    IpErr = sqrt((Ip + obj1->sky)/gain + pow(err_fac*fabs(Ip),2));

	    delta = (KLPsfReg[c]->rows[i][j] - SOFT_BIAS)/P0;

	    if((Ip + cond3_fac*IpErr) < delta*(I0 - cond3_fac*I0Err)) {
	       is_cr++;
	    }
	 }
      }
   }
#endif

   return(is_cr >= ncr_min ? 1 : 0);
}

/*****************************************************************************/
/*
 * classify an OBJECT1
 */
static void
classify_obj1(OBJC *objc,
	      const FIELDPARAMS *fiparams,
	      int c)
{
   OBJECT1 *obj1 = objc->color[c];

   if(objc->type == OBJ_SKY) {
      obj1->type = objc->type;
   } else if(objc->type == OBJ_UNK) {
      obj1->type = do_classify(fiparams->sg_classifier, NCLASSIFIER_FIDDLE,
			       obj1->star_L, obj1->exp_L, obj1->deV_L,
			       obj1->psfCounts, obj1->psfCountsErr,
			       obj1->counts_deV, obj1->counts_deVErr,
			       obj1->counts_exp, obj1->counts_expErr,
			       obj1->counts_model, obj1->counts_modelErr);
/*
 * Is object consistent with being a CR?
 */
      if(classify_as_cr(objc, fiparams, c)) {
	 obj1->flags2 |= OBJECT2_MAYBE_CR;
      }
/*
 * Is the symmetrically-placed pixel saturated?
 */
      {
	 const REGION *data = fiparams->frame[c].data; /* the pixel data */
	 const float size = 3;	/* size of area to search for interp*/
	 const int c0 = (data->ncol - obj1->colc) - size/2 + 0.5;
	 const int r0 = obj1->rowc - size/2 + 0.5;
	 const SPANMASK *sm = (SPANMASK *)data->mask;
	 shAssert(sm != NULL && sm->cookie == SPAN_COOKIE);
	 
	 if(phRectIntersectMask(sm->masks[S_MASK_SATUR],
				c0, r0, c0 + size + 0.5, r0 + size + 0.5)) {
	    obj1->flags2 |= OBJECT2_MAYBE_EGHOST;
	 }
      }
   }
}

/*****************************************************************************/
/*
 * Global classification of object
 */
static void
classify_objc(OBJC *objc,
	      const FIELDPARAMS *fiparams) /* fiber color etc. */
{
   int be_picky;			/* should I be picky about flags? */
   float chisq_star, chisq_deV, chisq_exp; /* obj1->chisq_star etc. */
   float psfCounts, psfCountsErr;	/* obj1->psfCounts etc. */
   float counts_deV, counts_deVErr;	/* obj1->counts_deV etc. */
   float counts_exp, counts_expErr;	/* obj1->counts_exp etc. */
   float counts_model, counts_modelErr;	/* obj1->counts_model etc. */
   int i;
   int nu_star, nu_deV, nu_exp;		/* obj1->nu_star */
   float star_L, exp_L, deV_L;		/* obj1->star_L etc. */
   const OBJECT1 *obj1;			/* == objc->color[i] */

   if(objc->type != OBJ_UNK) {
      return;
   }
   
   chisq_star = chisq_deV = chisq_exp = 0;
   nu_star = nu_deV = nu_exp = 0;

   psfCounts = psfCountsErr = 0;
   counts_deV = counts_deVErr = 0;
   counts_exp = counts_expErr = 0;
   counts_model = counts_modelErr = 0;
/*
 * Go through all colours seeing if there are any nice clean bands to
 * classify on; if there aren't we shall have to be more accepting of
 * questionable data
 */
   be_picky = 0;
   for(i = 0;i < objc->ncolor;i++) {
      obj1 = objc->color[i];
      if((obj1->flags & OBJECT1_DETECTED) &&
	 !(obj1->flags2 & (OBJECT2_DEBLEND_NOPEAK | OBJECT2_INTERP_CENTER))) {
	 be_picky = 1;
	 break;				/* a good band */
      }
   }
/*
 * collect information from all the bands
 */
   for(i = 0;i < objc->ncolor;i++) {
      obj1 = objc->color[i];
      if(!(obj1->flags & OBJECT1_DETECTED) ||
	 obj1->flags2 & OBJECT2_DEBLEND_NOPEAK) {
	 continue;
      }
      if(be_picky && (obj1->flags2 & OBJECT2_INTERP_CENTER)) {
	 continue;
      }
      
      psfCounts += obj1->psfCounts;
      psfCountsErr += pow(obj1->psfCountsErr, 2);
      chisq_star += obj1->chisq_star;
      nu_star += obj1->nu_star;
      
      counts_deV += obj1->counts_deV;
      counts_deVErr += pow(obj1->counts_deVErr, 2);
      chisq_deV += obj1->chisq_deV;
      nu_deV += obj1->nu_deV;
      
      counts_exp += obj1->counts_exp;
      counts_expErr += pow(obj1->counts_expErr, 2);
      chisq_exp += obj1->chisq_exp;
      nu_exp += obj1->nu_exp;
      
      counts_model += obj1->counts_model;
      counts_modelErr += pow(obj1->counts_modelErr, 2);
   }
   
   if(nu_star == 0) {			/* no information */
      objc->type = OBJ_UNK;
      return;
   }
/*
 * actually classify
 */
   psfCountsErr = sqrt(psfCountsErr);
   counts_deVErr = sqrt(counts_deVErr);
   counts_expErr = sqrt(counts_expErr);
   counts_modelErr = sqrt(counts_modelErr);

   star_L = phChisqProb(chisq_star, nu_star, L_SOFT);
   deV_L = phChisqProb(chisq_deV, nu_deV, L_SOFT);
   exp_L = phChisqProb(chisq_exp, nu_exp, L_SOFT);

   objc->type = do_classify(fiparams->sg_classifier, NCLASSIFIER_FIDDLE,
			    star_L, exp_L, deV_L,
			    psfCounts, psfCountsErr,
			    counts_deV, counts_deVErr,
			    counts_exp, counts_expErr,
			    counts_model, counts_modelErr);
}

/*****************************************************************************/
/*
 * find the fraction of light in the PSF for an object
 */
static void
find_fracPSF(OBJC *objc,
	     int color)
{
   OBJECT1 *obj1 = objc->color[color];

   obj1->fracPSF = (obj1->type == OBJ_STAR) ? 1 : 0;
}

/*****************************************************************************/
/*
 * Update various fields in the FIELDSTAT from the given OBJC.
 *
 * We cannot update the median colours, but we can accumulate the
 * information
 */
static void
update_fieldstat(OBJC *objc,
		 const FIELDPARAMS *fiparams,
		 FIELDSTAT *fstat)
{
   float fiberCts, psfCts;		/* fiber and psf counts in this and */
   float fiberCtsR, psfCtsR;		/* next redder band, (relative to 20th
					   mag object) */
   const float minCts = pow(10, -0.4);	/* minimum acceptable counts
					   (relative to 20th mag object) */
   int i;
   int n;				/* index of this object in
					   medians arrays */
   
   fstat->nobjects++;
   fstat->nchild += objc->nchild;

   switch (objc->type) {
    case OBJ_STAR:
      fstat->nstars++;
      break;
    case OBJ_GALAXY:
      fstat->ngals++;
      break;
    default:
      break;
   }
/*
 * Only save colours etc. of stars as we are using this for QA, and don't want
 * to worry about shapes of galaxies, or the effects of large scale structure
 * on the star/galaxy ratio (and thus the median colour)
 */
   if(objc->type == OBJ_STAR) {
/*
 * Save the object's colours so that we can return them in the FIELDSTAT
 *
 * It might seem natural to realloc the arrays if there are too many
 * objects, but this confuses the memory allocator --- more specifically,
 * it leads to an allocation that cannot be freed at the end of the field
 */
      if(medians.nobj == medians.size) {
	 medians.nobj--;
      }
      n = medians.nobj;
      
      for(i = 0;i < objc->ncolor;i++) {
	 shAssert(minCts > 0.0);
	 
	 if(i == objc->ncolor - 1) {	/* no colour to compare to */
	    fiberCts = fiberCtsR = psfCts = psfCtsR = minCts + 1;	 
	 } else {
	    fiberCts = objc->color[i]->fiberCounts/fiparams->frame[i].flux20;
	    psfCts = objc->color[i]->psfCounts/fiparams->frame[i].flux20;
	    fiberCtsR =
	      objc->color[i+1]->fiberCounts/fiparams->frame[i+1].flux20;
	    psfCtsR = objc->color[i+1]->psfCounts/fiparams->frame[i+1].flux20;
	 }
	 
	 if(fiberCtsR < minCts || psfCtsR < minCts ||
	    fiberCts < minCts || psfCts < minCts) { /* ignore this band */
	    medians.fiberColors[i][n] = medians.psfColors[i][n] = 1e10;
	    medians.Q[i][n] = medians.U[i][n] = -1000;
	 } else {
	    medians.fiberColors[i][n] = -2.5*log10(fabs(fiberCts/fiberCtsR));
	    medians.psfColors[i][n] = -2.5*log10(fabs(psfCts/psfCtsR));
	    medians.Q[i][n] = objc->color[i]->Q;
	    medians.U[i][n] = objc->color[i]->U;
	 }
      }
      
      medians.nobj++;
   }
}

/*****************************************************************************/

static int
cmp(const void *a, const void *b)
{
   float fa = *(float *)a;
   float fb = *(float *)b;

   return((fa < fb) ? -1 : ((fa == fb) ? 0 : 1));
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 *
 * Fill a FIELDSTAT's median {{fiber,psf}Color,Q,U} fields, ignoring all
 * values over 1e10; also set the n_eff fields if not yet done
 */
void
phFieldstatSetFromMO(FIELDSTAT *fieldstat)
{
   int i, j;
   const int nobj = medians.nobj;
   int n;				/* number in this band */
   
   shAssert(medians.size > 0);

   if(nobj == 0) {
      for(i = 0;i < medians.ncolor;i++) {
	 fieldstat->median_fiberColor[i] = fieldstat->median_psfColor[i] = 0;
      }
      return;
   }

   for(i = 0;i < medians.ncolor;i++) {
      qsort(medians.fiberColors[i],nobj,sizeof(medians.fiberColors[i][0]),cmp);
      for(n = nobj - 1; n >=0 && medians.fiberColors[i][n] >= 1e10; n--) {
	 continue;
      }
      if(n <= 0) {
	 fieldstat->median_fiberColor[i] = 1e10;
      } else {
	 fieldstat->median_fiberColor[i] = ((n%2 == 1) ?
	    medians.fiberColors[i][n/2] :
	    0.5*(medians.fiberColors[i][n/2-1] + medians.fiberColors[i][n/2]));
      }

      qsort(medians.psfColors[i],nobj,sizeof(medians.psfColors[i][0]),cmp);
      for(n = nobj - 1; n >=0 && medians.psfColors[i][n] >= 1e10; n--) {
	 continue;
      }
      if(n <= 0) {
	 fieldstat->median_psfColor[i] = 1e10;
      } else {
	 fieldstat->median_psfColor[i] = ((n%2 == 1) ?
		medians.psfColors[i][n/2] :
		0.5*(medians.psfColors[i][n/2-1] + medians.psfColors[i][n/2]));
      }

      qsort(medians.Q[i],nobj,sizeof(medians.Q[i][0]),cmp);
      qsort(medians.U[i],nobj,sizeof(medians.U[i][0]),cmp);
      for(j = 0; j < nobj && medians.Q[i][j] <= -999; j++) continue;
      n = nobj - j;
      if(n <= 0) {
	 fieldstat->Q[i] = fieldstat->U[i] = 1e10;
      } else {
	 fieldstat->Q[i] = (n%2 == 1) ? medians.Q[i][j + n/2] :
		       0.5*(medians.Q[i][j + n/2 - 1] + medians.Q[i][j + n/2]);
	 fieldstat->U[i] = (n%2 == 1) ? medians.U[i][j + n/2] :
		       0.5*(medians.U[i][j + n/2 - 1] + medians.U[i][j + n/2]);
      }
   }
}

/*****************************************************************************/
/*
 * fit a straight line to a set of points, y = a + b*t
 *
 * Return the value of chi^2 for the fit
 */
static float
fit_1d_velocity(const float *t,		/* "x" values */
		const float *y_i,	/* input "y" values */
		const float *yErr,	/* errors in y */
		int n,			/* number of points to fit */
		float *a,		/* estimated y-value at t == 0 */
		float *aErr,		/* error in a */
		float *covar,		/* covariance of a and b (not sqrt) */
		float *b,		/* value of slope */
		float *bErr)		/* error in b */
{
   double det;				/* determinant needed to find a,b */
   int i;
   float ivar;				/* == 1/yErr^2 */
   double sum, sum_t, sum_tt, sum_y, sum_ty;
   float y[NCOLOR];			/* == y_i[] - ybar */
   double ybar;				/* straight mean of y */

   shAssert(n > 0 && n <= NCOLOR);
/*
 * Start by subtracting off y's mean, to improve numerical stability
 */
   sum_y = 0.0;
   for(i = 0; i < n; i++) {
      sum_y += y_i[i];
   }
   ybar = sum_y/n;
/*
 * Now fit to data
 */
   sum = sum_t = sum_tt = sum_y = sum_ty = 0.0;
   for(i = 0; i < n; i++) {
      y[i] = y_i[i] - ybar;

      ivar = (yErr[i] < 1e-10) ? 1e20 : 1/(yErr[i]*yErr[i]);
      sum += ivar;
      sum_t += ivar*t[i];
      sum_tt += ivar*t[i]*t[i];
      sum_y += ivar*y[i];
      sum_ty += ivar*y[i]*t[i];
   }
   det = (sum*sum_tt - sum_t*sum_t);
   if(fabs(det) < 1e-10) {
      *a = 0; *aErr = -9999;
      *covar = 0;
      *b = 0; *bErr = -9999;
      return(-1);
   }

   *a =  (sum_tt*sum_y - sum_t*sum_ty)/det;
   *aErr = sqrt(sum_tt/det);
   *covar = -sum_t/det;
   *b = (sum*sum_ty - sum_t*sum_y)/det;
   *bErr = sqrt(sum/det);

   if(*aErr > 1e10) { *aErr = -9999; }
   if(*bErr > 1e10) { *bErr = -9999; }
/*
 * Find chi^2
 */
   sum = 0;
   for(i = 0; i < n; i++) {
      ivar = (yErr[i] < 1e-10) ? 1e10 : 1/yErr[i];
      sum += pow(ivar*((*a) + (*b)*t[i] - y[i]), 2);
   }

   *a += ybar;

   return(sum);
}

/*****************************************************************************/
/*
 * <AUTO EXTRACT>
 * 
 * Estimate an object's velocity from its position in the bands it's
 * detected in.
 * 
 * Return the reduced chi^2 for the fit; there are 2*(ndetect - 2) degrees
 * of freedom, so a negative chi^2 is returned for ndetect <= 2
 */
float
phVelocityFind(OBJC *objc,		/* OBJC whose velocity is desired */
	       const FIELDPARAMS *fiparams, /* info about the frame */
	       float *row,		/* fitted row positions, or NULL */
	       float *rowErr,		/* errors in row, or NULL */
	       float *col,		/* fitted column positions, or NULL */
	       float *colErr)		/* errors in col, or NULL */
{
   int c;
   float chisq;				/* chisq for the fit */
   float colc[NCOLOR], colcErr[NCOLOR];	/* column position and error */
   float covar;				/* covariance of pos and velocity */
   float drow, dcol;			/* offsets from reference colour */
   float drowErr, dcolErr;		/* errors in drow, dcol */
   float errMin = 1e-2;			/* minimum possible positional error */
   float psfMags[NCOLOR];		/* PSF magnitudes in all bands */
   float psfMagsErr[NCOLOR];		/* errors in psfMaqs */
   const OBJECT1 *obj1;			/* == objc->colors[c] */
   int ndetect = 0;			/* number of bands with detections */
   float pos, posErr;			/* estimated position at t==0 + error*/
   float rowc[NCOLOR], rowcErr[NCOLOR];	/* row position and error */
   float t[NCOLOR];			/* "time" of detection, measured
					   in frames not seconds */
   
   shAssert(objc != NULL && fiparams != NULL);
/*
 * estimate magnitudes in each band to get the astrometric errors right
 */
   for(c = 0; c < objc->ncolor; c++) {
      const float cts = objc->color[c]->psfCounts;
      if(cts <= 0) {
	 psfMags[c] = 30; psfMagsErr[c] = 1e10;
      } else {
	 psfMags[c] = 20 - 2.5*log10(cts/fiparams->frame[c].flux20);
	 psfMagsErr[c] = 2.5/log(10)*objc->color[c]->psfCountsErr/cts;
      }
   }
/*
 * estimate position for each band in canonical band's coordinate system
 */
   for(c = 0; c < objc->ncolor; c++) {
      if((obj1 = objc->color[c]) != NULL && obj1->flags & OBJECT1_DETECTED) {
	 if(obj1->flags & OBJECT1_CANONICAL_CENTER) {
	    continue;
	 }

	 if((obj1->flags  & OBJECT1_PEAKCENTER) ||
	    (obj1->flags2 & OBJECT2_INTERP_CENTER)) { /* don't use position */
	    continue;
	 }
	 
	 t[ndetect] = fiparams->frame[c].dframe;
	 phOffsetDo(fiparams, obj1->rowc, obj1->colc, 
		    c, fiparams->ref_band_index,
		    psfMags, psfMagsErr, &drow, &drowErr, &dcol, &dcolErr);

	 rowc[ndetect] = obj1->rowc + drow;
	 rowcErr[ndetect] = sqrt(pow(obj1->rowcErr,2) + pow(drowErr,2) +
				 pow(errMin,2));
	 colc[ndetect] = obj1->colc + dcol;
	 colcErr[ndetect] = sqrt(pow(obj1->colcErr,2) + pow(dcolErr,2) +
				 pow(errMin,2));
	 ndetect++;
      }
   }

   if(ndetect == 0) {
      return(-1);
   }
/*
 * Fit row data
 */
   chisq = fit_1d_velocity(t, rowc, rowcErr, ndetect, &pos, &posErr, &covar,
			   &objc->rowv, &objc->rowvErr);
/*
 * Estimate row, rowErr, etc.
 */
   for(c = 0; c < objc->ncolor; c++) {
      const int dframe = fiparams->frame[c].dframe;
      if(row != NULL) {
	 row[c] = pos + objc->rowv*dframe;
      }
      if(rowErr != NULL) {
	 if(objc->rowvErr < 0) {
	    rowErr[c] = -9999;
	 } else {
	    rowErr[c] = sqrt(pow(posErr,2) + 2*dframe*covar +
						 pow(dframe*objc->rowvErr, 2));
	 }
      }
   }
/*
 * and now column
 */
   chisq += fit_1d_velocity(t, colc, colcErr, ndetect, &pos, &posErr, &covar,
			    &objc->colv, &objc->colvErr);
/*
 * Estimate col, colErr, etc.
 */
   for(c = 0; c < objc->ncolor; c++) {
      const int dframe = fiparams->frame[c].dframe;
      if(col != NULL) {
	 col[c] = pos + objc->colv*dframe;
      }
      if(colErr != NULL) {
	 if(objc->colvErr < 0) {
	    colErr[c] = -9999;
	 } else {
	    colErr[c] = sqrt(pow(posErr,2) + 2*dframe*covar +
						 pow(dframe*objc->colvErr, 2));
	 }
      }
   }
/*
 * Divide chi^2 by number of degrees of freedom
 */
   if(ndetect <= 2) {
      if(chisq > 5e-4) {
	 shError("phVelocityFind: chisq = %.2g (OBJC id = %d)",
		 chisq, objc->id);
      }
      chisq = -1;
   } else {
      chisq /= 2*(ndetect - 2);
   }

   return(chisq);
}

/*****************************************************************************/
/*
 * 
 */
#if DO_PHIL

#if 1
static int
calcmom(float xcen, float ycen,		/* centre of object */
	int *ix1, int *ix2, int *iy1, int *iy2, /* bounding box to consider */
	float bkgd,			/* data's background level */
	int interpflag,			/* interpolate within pixels? */
	double w11, double w22, double w12, /* weights */
	float detw,
	double *w1, double *w2, double *ww12,
	double *sumx, double *sumy,	/* desired */
	double *sumxx, double *sumyy,	/* desired */
	double *sumxy, double *sum,	/*       sums */
	double *sum1, double *sum2,
	const REGION *data)		/* the data */
{   
   int i,j;
   U16 *row;				/* pointer to a row of data */
   float xx,xx2,yy,yy2;
   float xl,xh,yl,yh;
   float expon;
   float tmod,ymod;
   float xxx,yyy;
   float weight;
   float t1;
   float grad;
   double sumt,sumxt,sumyt,sumxxt,sumyyt,sumxyt;

   grad = 4.*sqrt( ((w11 > w22) ? w11 : w22));

   *ix1 = xcen - grad - 0.5;
   *ix1 = (*ix1 < 0) ? 0 : *ix1;
   *iy1 = ycen - grad - 0.5;
   *iy1 = (*iy1 < 0) ? 0 : *iy1;
     
   *ix2 = xcen + grad + 0.5;
   if(*ix2 >= data->ncol) {
     *ix2 = data->ncol - 1;
   }
   *iy2 = ycen + grad + 0.5;
   if(*iy2 >= data->nrow) {
     *iy2 = data->nrow - 1;
   }
   
   
   *w1 = w11/detw;
   *w2 = w22/detw;
   *ww12 = w12/detw;

   sumt = sumxt = sumyt = sumxxt = sumyyt = sumxyt = 0;
   
   for(i = *iy1;i <= *iy2;i++) {
      row = data->rows[i];
      yy = i-ycen;
      yy2 = yy*yy;
      yl = yy - 0.375;
      yh = yy + 0.375;
      for(j = *ix1;j <= *ix2;j++) {
	 xx = j-xcen;
	 if(interpflag) {
	   
	    xl = xx - 0.375;
	    xh = xx + 0.375;
	    expon = xl*xl**w2 + yl*yl**w1 - 2.*xl*yl**ww12;
	    t1 = xh*xh**w2 + yh*yh**w1 - 2.*xh*yh**ww12;
	    expon = (expon > t1) ? expon : t1;
	    t1 = xl*xl**w2 + yh*yh**w1 - 2.*xl*yh**ww12;
	    expon = (expon > t1) ? expon : t1;
	    t1 = xh*xh**w2 + yl*yl**w1 - 2.*xh*yl**ww12;
	    expon = (expon > t1) ? expon : t1;
	    if(expon <= 9.0) {
	       tmod = row[j] - bkgd;
	       for(yyy = yl;yyy <= yh;yyy+=0.25) {
		  yy2 = yyy*yyy;
		  for(xxx = xl;xxx <= xh;xxx+=0.25) {
		     xx2 = xxx*xxx;
		     expon = xx2**w2 + yy2**w1 - 2.*xxx*yyy**ww12;
		     weight = exp(-0.5*expon);
		     ymod = tmod*weight;
		     /* *sum1+=(xx2+yy2)*ymod;*/
		     /* *sum2+=(xx2-yy2)*ymod;*/
		     sumxt+=ymod*(xxx+xcen);
		     sumyt+=ymod*(yyy+ycen);
		     sumxxt+=xx2*ymod;
		     sumyyt+=yy2*ymod;
		     sumxyt+=xxx*yyy*ymod;
		     sumt+=ymod;
		  }
	       }
	    }
	 } else {
	    xx2 = xx*xx;
	    expon = xx2**w2 + yy2**w1 - 2.*xx*yy**ww12;

	    if(expon <= 9.0) {
	       weight = exp(-0.5*expon);
	       ymod = (row[j] - bkgd)*weight;
	       /*	       *sum1+=(xx2+yy2)*ymod;*/
	       /*	       *sum2+=(xx2-yy2)*ymod;*/
	       sumxt+=ymod*j;
	       sumyt+=ymod*i;
	       sumxxt += xx2*ymod;
	       sumyyt += yy2*ymod;
	       sumxyt += xx*yy*ymod;
	       sumt += ymod;
	    }
	 }
      }
   }

   *sum=sumt; *sumx=sumxt; *sumy=sumyt; *sumxx=sumxxt; *sumyy=sumyyt ; *sumxy=sumxyt;
#if 0
   *sum1=sum1t; *sum2=sum2t;
#else
   *sum1 = *sum2 = 0;
#endif

   if(*sum <= 0) {
      return(-1);
   } else {
      return(0);
   }
}
#endif

static void
calcerr(float xcen, float ycen,		/* centre of object */
	int ix1, int ix2, int iy1, int iy2, /* bounding box to consider */
	float bkgd,			/* data's background level */
	int interpflag,			/* interpolate within pixels? */
	float w1, float w2, float ww12,	/* weights */
	float sumxx, float sumyy, float sumxy, /* quadratic sums */
	double sum1, double sum2,
	float *errxx, float *erryy, float *errxy, /* errors in sums, or NULL */
	double *sums4,			/* sum x^4 */
	double *s1, double *s2,
	const REGION *data)		/* the data */
{   
   int i,j;
   U16 *row;				/* pointer to a row of data */
   double sxx = 0, sxy = 0, syy = 0, s4 = 0; /* unalias err{xx,xy,yy}, sums4 */
   float xx,xx2,yy,yy2;
   float xl,xh,yl,yh;
   float tmp;
   float tmod,ymod;
   float xxx,yyy;
   float weight;
   float expon;
   double sum3,sum4;
      
   sum3=sum4=0;
   for(i = iy1;i <= iy2;i++) {
      row = data->rows[i];
      yy = i-ycen;
      yy2 = yy*yy;
      if(interpflag) {
	 yl = yy - 0.375;
	 yh = yy + 0.375;
      }
      for(j = ix1;j <= ix2;j++) {
	 xx = j-xcen;
	 if(interpflag) {
	    xl = xx - 0.375;
	    xh = xx + 0.375;

	    expon = xl*xl*w2 + yl*yl*w1 - 2.*xl*yl*ww12;
	    tmp = xh*xh*w2 + yh*yh*w1 - 2.*xh*yh*ww12;
	    if(tmp > expon) {
	       expon = tmp;
	    }
	    
	    tmp = xl*xl*w2 + yh*yh*w1 - 2.*xl*yh*ww12;
	    if(tmp > expon) {
	       expon = tmp;
	    }
	    
	    tmp = xh*xh*w2 + yl*yl*w1 - 2.*xh*yl*ww12;
	    if(tmp > expon) {
	       expon = tmp;
	    }

	    if(expon <= 9.0) {
	       tmod = row[j] - bkgd;
	       for(yyy = yl;yyy <= yh;yyy+=0.25) {
		  yy2 = yyy*yyy;
		  for(xxx = xl; xxx <= xh; xxx += 0.25) {
		     xx2 = xxx*xxx;
		     expon = xx2*w2 + yy2*w1 - 2*xxx*yyy*ww12;
		     weight = exp(-0.5*expon);
		     ymod = tmod*weight;
		     /*		     sum3+=pow(weight*(xx2+yy2 - sum1),2);*/
		     /*		     sum4+=pow(weight*(xx2-yy2 - sum2),2);*/
		     sum3+=pow(weight*(xx2-yy2 - sum1*(xx2+yy2)),2);
		     sum4+=pow(weight*(xxx*yyy - sum2*(xx2+yy2)),2);
		     sxx += pow(weight*(xx2 - sumxx),2);
		     syy += pow(weight*(yy2 - sumyy),2);
		     sxy += pow(weight*(xxx*yyy - sumxy),2);
		     s4 += expon*expon*ymod;
		  }
	       }
	    }
	 } else {
	    xx2 = xx*xx;
	    expon = xx2*w2 + yy2*w1 - 2.*xx*yy*ww12;
	    if(expon <= 9.0) {
	       weight = exp(-0.5*expon);
	       ymod = (row[j] - bkgd)*weight;
	       sum3 += pow(weight*(xx2-yy2 - sum1*(xx2+yy2)),2);
	       sum4 += pow(weight*(xx*yy - sum2*(xx2+yy2)),2);
	       sxx += pow(weight*(xx2 - sumxx),2);
	       syy += pow(weight*(yy2 - sumyy),2);
	       sxy += pow(weight*(xx*yy - sumxy),2);
	       s4 += expon*expon*ymod;
	    }
	 }	       
      }
   }
/*
 * Pack up desired results
 */
   if(errxx != NULL) {
      *errxx = sxx;
   }
   if(errxy != NULL) {
      *errxy = sxy;
   }
   if(erryy != NULL) {
      *erryy = syy;
   }
   *sums4 = s4; *s1 = sum3; *s2 = sum4;
}

#define MAXIT 100
#define XINTERP 0.0
#define XINTERP2 0.0
#define TOL1 0.001
#define TOL2 0.01
/*
 * Calculate the adaptive moments.
 *
 * Note: this routine assumes that the _centre_ of a pixel is (0.0, 0.0)
 * which is not the SDSS convention. Caveat Lector.
 */
static void
calc_adaptive_moments(OBJC *objc,
		      int color,
		      const FIELDPARAMS *fparams )
{
   const float bkgd = SOFT_BIAS + fparams->frame[color].bkgd;
   const REGION *data = fparams->frame[color].data;
   float d;				/* weighted size of object */
   float detw;				/* determinant of matrix containing
					   moments of weighting function */
   float detm;				/* determinant of m[12][12] matrix */
   float detn;				/* determinant of n[12][12] matrix */
   float e1, e2;			/* current and old values of */
   float e1old = 1.e6, e2old =1.e6;	/*        shape parameters e1 and e2 */
   int imom;                            /* iteration number */
   int interpflag;			/* interpolate finer than a pixel? */
   int ix1,ix2,iy1,iy2;                 /* corners of the box being analyzed */
   float m11,m22,m12;			/* weighted Quad. moments of object */
   float m11old = 1.e6;			/* previous version of m11 (Qxx) */
   float n11, n22, n12;			/* matrix elements used to determine
					   next guess at weighting fcn */
   OBJECT1 *obj1 = objc->color[color];
   float shiftmax;                      /* max allowed centroid shift */
   float sigsky;                        /* sky standard deviation */
   double sum2t;			       
   double sum;				/* sum of intensity*weight */
   double sums4,sums4p;			 /* higher order moment used for
					   PSF correction */
   double sumx, sumy;			/* sum ((int)[xy])*intensity*weight */
   double sumxx, sumxy, sumyy;		/* sum {x^2,xy,y^2}*intensity*weight */
   double w1,w2,ww12;			/* w11 etc. divided by determinant */ 
   double w11,w22,w12;			/* moments of the weighting fcn */
   float xcen = obj1->colc - 0.5;	/* centre of */
   float ycen = obj1->rowc - 0.5;	/*           object */
   /*
    * default values
    */
   REGION *reg;
#if 0
   double summ;
#endif
   double sum1,sum2;
   double s1,s2,s1p,s2p,sum1p,sum2p;

   sums4 = -9999;

   shiftmax = 2*obj1->petroRad;
   if(shiftmax < 2) {
      shiftmax = 2;
   } else if(shiftmax > 10) {
      shiftmax = 10;
   }

   w11 = w22 = 1.5;
   w12 = 0;
   
   for(imom = 0; imom < MAXIT; imom++) {
#if 0
      printf("XXX colour %d  id %d imom %d",
	     objc->color[color]->id, color, imom);
#endif
      detw = w11*w22 - w12*w12;
      if(w11 < XINTERP || w22 < XINTERP || detw < XINTERP2) {
	 interpflag = 1;
      } else {
	 interpflag = 0;
      }
      shAssert(detw > 0.0);

      if(calcmom(xcen,ycen,&ix1,&ix2,&iy1,&iy2,bkgd,interpflag,w11,w22,w12,detw,
	      &w1,&w2,&ww12,&sumx,&sumy,&sumxx,&sumyy,&sumxy,&sum,&sum1,&sum2,data) < 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }
             
      if(sum <= 0) {			/* too faint to process */
	obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	return;
      }
#if 0
/*
* Find new centre
*/

      xcen = sumx/sum;
      ycen = sumy/sum;
#endif
      if(fabs(sumx/sum - (obj1->colc - 0.5)) > shiftmax || 
	 fabs(sumy/sum - (obj1->rowc - 0.5)) > shiftmax) {
	 obj1->Mrr = sumx/sum - (obj1->rowc - 0.5);
	 obj1->Mcc = sumy/sum - (obj1->colc - 0.5);
	 obj1->flags2 |= OBJECT2_AMOMENT_SHIFT;
	 return;
      }
      /*
       * OK, we have the centre. Proceed to find the second moments.
       * This is only needed if we update the centre; if we use the
       * obj1->{row,col}c we've already done the work
       */

      sum1/=sum;
      sum2/=sum;
      m11 = sumxx/sum;
      m22 = sumyy/sum;
      m12 = sumxy/sum;

      if(m11 <= 0 || m22 <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }

      d = m11 + m22;
      e1 = (m11 - m22)/d;
      e2 = 2.*m12/d;
      if(fabs(e1 - e1old) < TOL1 && fabs(e2 - e2old) < TOL1 &&
	 fabs(m11/m11old - 1.) < TOL2 ) {
#if 0
	summ=sum;
#endif
/* convergence criteria met */
	break;				
      }
 /*
 * Didn't converge, calculate new values for weighting function
 */
      e1old = e1;
      e2old = e2;
      m11old = m11;
      detm = m11*m22 - m12*m12;
      if(detm <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }
#if 0
      printf(" detm %g m11 %g m12 %g m22 %g", detm, m11, m12, m22);
#endif

      detm = 1./detm;
      detw = 1./detw;
      n11 = m22*detm - w22*detw;
      n22 = m11*detm - w11*detw;
      n12 = -m12*detm + w12*detw;
      detn = n11*n22 - n12*n12;
      if(detn <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }

#if 0
      printf(" detn %g n11 %g n12 %g n22 %g\n", detn, n11, n12, n22);
#endif
      detn = 1./detn;
      w11 = n22*detn;
      w22 = n11*detn;
      w12 = -n12*detn;

      if(w11 <= 0 || w22 <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }
   }

   if(imom == MAXIT) {
      obj1->flags2 |= OBJECT2_AMOMENT_MAXITER;
      return;
   }       
/*
 * calculate uncertainties
 */
   if(sumxx + sumyy == 0.0) {
      obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
      return;
   }
   
   sum1 = (sumxx - sumyy)/(sumxx + sumyy);
   sum2 = sumxy/(sumxx + sumyy);
   
   calcerr(xcen, ycen, ix1, ix2, iy1, iy2, bkgd, interpflag,
	   w1, w2, ww12, m11, m22, m12,sum1,sum2,
	   &obj1->MccErr, &obj1->MrrErr, &obj1->McrErr, &sums4, &s1, &s2, data);

   sigsky = sqrt(obj1->sky/fparams->frame[color].gain+
		 fparams->frame[color].dark_variance);
   sum2t = sigsky/sum;
   if(interpflag) {
     sum2t *= 4;
   }
   s1=sqrt(s1)*sigsky/(sumxx+sumyy);
   s2=2.*sqrt(s2)*sigsky/(sumxx+sumyy);

   obj1->Mcc = w11; obj1->MccErr = 2*sum2t*sqrt(obj1->MccErr);
   obj1->Mcr = w12; obj1->McrErr = 2*sum2t*sqrt(obj1->McrErr);
   obj1->Mrr = w22; obj1->MrrErr = 2*sum2t*sqrt(obj1->MrrErr);
   sums4 /= sum;

   obj1->Mcr4 = sums4;
/*
 * Now determine moments of PSF
 */
   if(fparams->frame[color].psfBasis == NULL) {
      shError("No PSF_BASIS for band %d", color);
      return;				/* XXX needs a diagnostic bit? */
   }
   
   shAssert(KLPsf[color] != NULL);
   reg = KLPsf[color]->reg;

   w11 = w22 = 1.5;
   w12 = 0;

   e1old=e2old=m11old=1.e6;

   xcen=reg->ncol/2-0.5;
   ycen=reg->nrow/2-0.5;
   for(imom = 0; imom < MAXIT; imom++) {
      detw = w11*w22 - w12*w12;
      if(w11 < XINTERP || w22 < XINTERP || detw < XINTERP2) {
	 interpflag = 1;
      } else {
	 interpflag = 0;
      }
      shAssert(detw > 0.0);

      if(calcmom(xcen,ycen,&ix1,&ix2,&iy1,&iy2,bkgd,interpflag,w11,w22,w12,detw,
	      &w1,&w2,&ww12,&sumx,&sumy,&sumxx,&sumyy,&sumxy,&sum,&sum1p,&sum2p,reg) < 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }
             
      if(sum <= 0) {			/* too faint to process */
	obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	return;
      }

      xcen=sumx/sum;
      ycen=sumy/sum;
      m11 = sumxx/sum;
      m22 = sumyy/sum;
      m12 = sumxy/sum;

      if(m11 <= 0 || m22 <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }

      d = m11 + m22;
      e1 = (m11 - m22)/d;
      e2 = 2.*m12/d;
      if(fabs(e1 - e1old) < TOL1 && fabs(e2 - e2old) < TOL1 &&
	 fabs(m11/m11old - 1.) < TOL2 ) {
/* convergence criteria met */
	break;				
      }
 /*
 * Didn't converge, calculate new values for weighting function
 */
      e1old = e1;
      e2old = e2;
      m11old = m11;
      detm = m11*m22 - m12*m12;
      if(detm <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }

      detm = 1./detm;
      detw = 1./detw;
      n11 = m22*detm - w22*detw;
      n22 = m11*detm - w11*detw;
      n12 = -m12*detm + w12*detw;
      detn = n11*n22 - n12*n12;
      if(detn <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }

      detn = 1./detn;
      w11 = n22*detn;
      w22 = n11*detn;
      w12 = -n12*detn;

      if(w11 <= 0 || w22 <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 return;
      }
   }
/*
 * pack up results
 */
   obj1->Mcc_psf = w11;
   obj1->Mcr_psf = w12;
   obj1->Mrr_psf = w22;

   if(imom == MAXIT) {
      obj1->flags2 |= OBJECT2_AMOMENT_MAXITER;
      return;
   }       

   calcerr(xcen, ycen, ix1, ix2, iy1, iy2, bkgd, interpflag,
	   w1, w2, ww12, m11, m22, m12,sum1,sum2,
	   NULL, NULL, NULL, &sums4p, &s1p, &s2p, reg);
   
   sums4p /= sum;

   obj1->Mcr4_psf = sums4p;
}
#endif
