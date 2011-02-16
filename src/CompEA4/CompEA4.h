#define DEFVAL     -9999.0  /* Default value for e1,e2,Rout */
#define MAXRVAL        0.8  /* maximum for which we apply dilution correction*/

#define FLAG_BADE      0x1  /* 2^0 Failed first sqrt test in CompEA4 */
#define FLAG_BADEMULT1 0x2  /* 2^1 Failed in first shearmult */
#define FLAG_BADEMULT2 0x4  /* 2^2 Failed in second shearmult */
#define FLAG_BADERED1  0x8  /* 2^3 Failed sqrt e1red test */
#define FLAG_BADERED2  0x10 /* 2^4 Failed sqrt e1red test */
#define FLAG_BADRVAL   0x20 /* 2^5 R==0.0 */

typedef struct {

  IDL_MEMINT Ndata;

  float *e1; 
  float *e2; 
  float *m_rr_cc; 
  float *rho4; 


  float *e1_psf; 
  float *e2_psf; 
  float *m_rr_cc_psf;
  float *rho4_psf;

} CompEA4InStruct;

typedef struct {

  IDL_MEMINT Ndata;

  double *e1out;
  IDL_VPTR e1outVptr;

  double *e2out;
  IDL_VPTR e2outVptr;

  double *Rout;
  IDL_VPTR RoutVptr;


  IDL_MEMINT *flags;
  IDL_VPTR flagsVptr;

} CompEA4OutStruct;



int shearmult(double e1a, double e2a, double e1b, double e2b, double
	      *e1out, double *e2out);

int CompEA4(double Tratio, double e1p, double e2p, double a4p, double
	    e1o, double e2o, double a4o,
	    double *e1, double *e2, double *Rout);

float* 
CompEA4GetFloatPtr(IDL_VPTR vptr, IDL_MEMINT *n, IDL_MEMINT nexpected);
CompEA4InStruct * CompEA4GetInputs(IDL_VPTR argv[]);
CompEA4OutStruct * CompEA4GetOutputs(IDL_MEMINT n);
