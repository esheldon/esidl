
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
	float w11, float w22, float w12,	/* weights */
	float detw,
	float *w1,float *w2, float *ww12,
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
	float *errxx, float *erryy, float *errxy, /* errors in sums */
	double *sums4,			/* ?? */
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
	       /*	       sum3+=pow(weight*(xx2+yy2 - sum1),2);*/
	       /*	       sum4+=pow(weight*(xx2-yy2 - sum2),2);*/
	       sum3+=pow(weight*(xx2-yy2 - sum1*(xx2+yy2)),2);
	       sum4+=pow(weight*(xx*yy - sum2*(xx2+yy2)),2);
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
   *errxx = sxx; *erryy = syy; *errxy = sxy; *sums4 = s4; *s1=sum3; *s2=sum4;
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
   /*   float expon;*/                         /* exponent in the weighting fcn */ 
   /*   float grad;	*/			/* radius of box used for analysis */
   int imom;                            /* iteration number */
   int interpflag;			/* interpolate finer than a pixel? */
   int ix1,ix2,iy1,iy2;                 /* corners of the box being analyzed */
   float m11,m22,m12;			/* weighted Quad. moments of object */
   float m11old = 1.e6;			/* previous version of m11 (Qxx) */
   float n11, n22, n12;			/* matrix elements used to determine
					   next guess at weighting fcn */
   OBJECT1 *obj1 = objc->color[color];
   /*   float qxxerr,qyyerr,qxyerr;*/		/* errors in qxx etc. */
   /*U16 *row;*/				/* pointer to a row of data */
   float shiftmax;                      /* max allowed centroid shift */
   float sigsky;                        /* sky standard deviation */
   double sum2t;			       
   double sum;				/* sum of intensity*weight */
   double sums4,sums4p;			 /* higher order moment used for
					   PSF correction */
   double sumx, sumy;			/* sum ((int)[xy])*intensity*weight */
   double sumxx, sumxy, sumyy;		/* sum {x^2,xy,y^2}*intensity*weight */
   /*float tmp;*/                            /* temporary variable */
   /*float tmod;*/                          /* background subtracted flux value*/
   float w1,w2,ww12;                    /* w11 etc. divided by determinant */ 
   float w11,w22,w12;                   /* moments of the weighting fcn */
   /*  float weight;	*/		/* value of the weighting function */
   /*float ymod;*/				/* weight*(pixel value - bkgd) */
   float xcen = obj1->colc - 0.5;	/* centre of */
   float ycen = obj1->rowc - 0.5;	/*           object */
   /*   float xl,xh,yl,yh;*/                   /* max and min subsampling points for
					   interpolating weight fcn*/
   /*   float xxx,yyy;  */                     /* relative coordinates where moment
					  is being evaluated */
   /*   float yy,yy2,xx,xx2; */                 /* relative coords and their squares
					   where moments are evaluated*/
   /*
    * default values
    */
   PSF_REG *psf_reg;			/* reconstructed PSF region etc. */
   REGION *reg;
   float PccErr, PrrErr, PcrErr;
   float e1err,e2err;
   /*   float M1err,M2err;*/                 /*output uncertainties on PSF corrected e1 and e2*/
   float smear;                       /*smear polarizability*/
   float psfsize,objsize;            /*size of psf and object*/
   double summ;
   /*   float tmp;*/
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
	summ=sum;
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

   if(sums4 == 0.0) {
      obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
      return;
   }

   sigsky = sqrt(obj1->sky/fparams->frame[color].gain+
		 fparams->frame[color].dark_variance);
   sum2t = sigsky/sum;
   if(interpflag) {
     sum2t *= 4;
   }
   /*   s1*=4*sum2t*sum2t;*/
   /*   s2*=4*sum2t*sum2t;*/
   s1=sqrt(s1)*sigsky/(sumxx+sumyy);
   s2=2.*sqrt(s2)*sigsky/(sumxx+sumyy);

   obj1->Mcc = w11; obj1->MccErr = 2*sum2t*sqrt(obj1->MccErr);
   obj1->Mcr = w12; obj1->McrErr = 2*sum2t*sqrt(obj1->McrErr);
   obj1->Mrr = w22; obj1->MrrErr = 2*sum2t*sqrt(obj1->MrrErr);
#if 0 
   M1err=4*sqrt(M_PI)*sigsky*pow((obj1->Mcc*obj1->Mrr - obj1->Mcr*obj1->Mcr),0.25)/(4*summ-sums4);
   if(interpflag)M1err*=16;
#endif
   sums4 /= sum;
   
   /******************************
    * Now determine moments of PSF
    ******************************/

   if(fparams->frame[color].psfBasis == NULL) {
      shError("No PSF_BASIS for band %d", color);
      return;				/* XXX needs a diagnostic bit? */
   }
   
   /* calculate moments for PSF*/
   psf_reg = phPsfKLReconstruct(fparams->frame[color].psfBasis, ycen, xcen, TYPE_U16);


   reg=psf_reg->reg;

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
	 phPsfRegDel(psf_reg);
	 return;
      }
             
      if(sum <= 0) {			/* too faint to process */
	obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	phPsfRegDel(psf_reg);
	return;
      }

      xcen=sumx/sum;
      ycen=sumy/sum;
      m11 = sumxx/sum;
      m22 = sumyy/sum;
      m12 = sumxy/sum;

      if(m11 <= 0 || m22 <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 phPsfRegDel(psf_reg);
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
	 phPsfRegDel(psf_reg);
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
	 phPsfRegDel(psf_reg);
	 return;
      }

      detn = 1./detn;
      w11 = n22*detn;
      w22 = n11*detn;
      w12 = -n12*detn;

      if(w11 <= 0 || w22 <= 0) {
	 obj1->flags2 |= OBJECT2_AMOMENT_FAINT;
	 phPsfRegDel(psf_reg);
	 return;
      }
   }

   if(imom == MAXIT) {
      obj1->flags2 |= OBJECT2_AMOMENT_MAXITER;
      phPsfRegDel(psf_reg);
      return;
   }       

   calcerr(xcen, ycen, ix1, ix2, iy1, iy2, bkgd, interpflag,
	   w1, w2, ww12, m11, m22, m12,sum1,sum2,
	   &PccErr, &PrrErr, &PcrErr, &sums4p, &s1p, &s2p, reg);
   
   sums4p /= sum;

   psfsize=w11+w22;
   objsize=obj1->Mcc+obj1->Mrr;
   smear=(psfsize)/(objsize)*(4/sums4p-1)/(4/sums4-1);
   if(smear<=0.8){
     e1old=(obj1->Mcc-obj1->Mrr)/objsize;
     e2old=2.*obj1->Mcr/objsize;
     e1=(obj1->Mcc-obj1->Mrr)/objsize-smear*(w11-w22)/psfsize;
     e2=2*(obj1->Mcr/objsize-smear*w12/psfsize);
     e1/=(1-smear);
     e2/=(1-smear);
#if 0
     M1err/=(1-smear);
     tmp=s1/pow(sum1,4);
     e1err=sqrt(s2/(sum1*sum1)+sum2*sum2*tmp);
     e2err=2*sqrt(obj1->McrErr*obj1->McrErr/(sum1*sum1)+obj1->Mcr*obj1->Mcr*tmp);
#endif
     e1err=s1;
     e2err=s2;
     e2err/=(1-smear);
     e1err/=(1-smear);


#if 0
     e2err = 2/objsize*sqrt(obj1->McrErr*obj1->McrErr +
			  obj1->MccErr*obj1->MccErr*tmp +
			  obj1->MrrErr*obj1->MrrErr*tmp);
     e1err = 2/(objsize*objsize)*sqrt(pow(obj1->Mcc*obj1->MrrErr, 2) +
				      pow(obj1->Mrr*obj1->MccErr, 2));
#endif
     
#if 0
     fprintf(stdout,"%f %f %f %f %f %f %f %f\n",obj1->petroCounts,e1,e2,M1err,
	     e1err,e2err,obj1->Mrr,obj1->Mcc);

     fprintf(stdout,"%f %f %f %f %f %f %f\n",obj1->petroCounts,obj1->Mrr,
	     obj1->Mcc,obj1->Mcr,obj1->MrrErr,obj1->MccErr,obj1->McrErr);
#endif

   }


   obj1->Me1 = e1; obj1->Me1Err = e1err;
   obj1->Me2 = e2; obj1->Me2Err = e2err;
   
#if 0
   fprintf(stderr,"%d %f %f %f %f %f %f %f %f %f %d %f %f \n",obj1->id,
	   xcen,ycen,obj1->Mcc,obj1->Mrr,obj1->Mcr,obj1->MccErr,
	   obj1->MrrErr,obj1->McrErr,sums4,imom,obj1->colc,obj1->rowc);
#endif

   phPsfRegDel(psf_reg);

}
#endif

