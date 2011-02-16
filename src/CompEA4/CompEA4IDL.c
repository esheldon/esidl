#include <stdlib.h>
#include <stdio.h>
#include "idl_export.h"
#include "CompEA4.h"

void compea4(int argc, IDL_VPTR argv[])
{

	IDL_VPTR e1out, e2out, Rout, flags;
	CompEA4InStruct *inStruct;
	CompEA4OutStruct *outStruct;

	double 
		Tratio, a4, a4_psf;

	IDL_MEMINT i;

	if (argc != 12)
	{
		IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
				"Syntax: compea4, e1, e2, rrcc, rho4, e1psf, e2psf, rrccpsf, rho4psf, e1out, e2out, Rout, flags");
	}

	e1out = argv[8];
	e2out = argv[9];
	Rout  = argv[10];
	flags = argv[11];

	inStruct = CompEA4GetInputs(argv);
	if (inStruct == NULL)
		return;

	outStruct = CompEA4GetOutputs(inStruct->Ndata);
	if (outStruct == NULL)
		return;

	for(i=0;i<inStruct->Ndata;++i)
	{

		Tratio = (double) inStruct->m_rr_cc_psf[i]/inStruct->m_rr_cc[i];

		a4 = (double) inStruct->rho4[i]/2. - 1.;
		a4_psf = (double) inStruct->rho4_psf[i]/2. - 1.;

		// these floats get implicitly converted to doubles when passed
		outStruct->flags[i] = 
			(IDL_MEMINT) CompEA4( Tratio, 
					(double) inStruct->e1_psf[i], 
					(double) inStruct->e2_psf[i], 
					a4_psf, 
					(double) inStruct->e1[i], 
					(double) inStruct->e2[i], 
					a4, 
					&(outStruct->e1out[i]), 
					&(outStruct->e2out[i]), 
					&(outStruct->Rout[i]) );

		// Note I used to return 1-R
		//outStruct->Rout[i] = 1.0 - outStruct->Rout[i];

	}

	// Output Variables. No copy is done since they are
	// temporary variables.
	IDL_VarCopy(outStruct->e1outVptr, e1out);
	IDL_VarCopy(outStruct->e2outVptr, e2out);
	IDL_VarCopy(outStruct->RoutVptr, Rout);
	IDL_VarCopy(outStruct->flagsVptr, flags);

	free(inStruct);
	free(outStruct);

}

float *CompEA4GetFloatPtr(IDL_VPTR vptr, IDL_MEMINT *n, IDL_MEMINT nexpected)
{
	float *data;

	if (vptr->type != IDL_TYP_FLOAT)
	{
		IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, 
				"All inputs must be type FLOAT");
		return(NULL);
	}

	IDL_VarGetData(vptr, n, (char **) &data, FALSE);

	if (nexpected != -1)
	{
		if (*n != nexpected)
		{
			IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, 
					"All inputs must be same length");
			return(NULL);
		}
	}


	return(data);
}
CompEA4InStruct * CompEA4GetInputs(IDL_VPTR argv[])
{

	IDL_VPTR tmpVptr;
	IDL_MEMINT n;
	float *tptr;

	CompEA4InStruct *inStruct;
	inStruct = (CompEA4InStruct *) calloc(1, sizeof(CompEA4InStruct));

	// e1
	inStruct->e1 = CompEA4GetFloatPtr(argv[0], &n, -1);

	if (inStruct->e1 == NULL)
	{      
		free(inStruct);
		return(NULL);
	}

	inStruct->Ndata = n;


	// e2
	inStruct->e2 = CompEA4GetFloatPtr(argv[1], &n, inStruct->Ndata);
	if (inStruct->e2 == NULL)
	{      
		free(inStruct);
		return(NULL);
	}

	// m_rr_cc
	inStruct->m_rr_cc = CompEA4GetFloatPtr(argv[2], &n, inStruct->Ndata);
	if (inStruct->m_rr_cc == NULL)
	{      
		free(inStruct);
		return(NULL);
	}

	// rho4
	inStruct->rho4 = CompEA4GetFloatPtr(argv[3], &n, inStruct->Ndata);

	if (inStruct->rho4 == NULL)
	{      
		free(inStruct);
		return(NULL);
	}




	// e1_psf
	inStruct->e1_psf = CompEA4GetFloatPtr(argv[4], &n, inStruct->Ndata);

	if (inStruct->e1_psf == NULL)
	{      
		free(inStruct);
		return(NULL);
	}
	// e2_psf
	inStruct->e2_psf = CompEA4GetFloatPtr(argv[5], &n, inStruct->Ndata);

	if (inStruct->e2_psf == NULL)
	{      
		free(inStruct);
		return(NULL);
	}

	// m_rr_cc_psf
	inStruct->m_rr_cc_psf = CompEA4GetFloatPtr(argv[6], &n, inStruct->Ndata);

	if (inStruct->m_rr_cc_psf == NULL)
	{      
		free(inStruct);
		return(NULL);
	}

	// rho4_psf
	inStruct->rho4_psf = CompEA4GetFloatPtr(argv[7], &n, inStruct->Ndata);

	if (inStruct->rho4_psf == NULL)
	{      
		free(inStruct);
		return(NULL);
	}

	return(inStruct);

}

CompEA4OutStruct * CompEA4GetOutputs(IDL_MEMINT n)
{

	IDL_ARRAY_DIM dim;
	CompEA4OutStruct * outStruct;

	dim[0] = n;

	outStruct = calloc(1,sizeof(CompEA4OutStruct));
	outStruct->Ndata = n;

	outStruct->e1out = 
		(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, 1, dim, 
				IDL_ARR_INI_NOP, &(outStruct->e1outVptr));
	outStruct->e2out = 
		(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, 1, dim, 
				IDL_ARR_INI_NOP, &(outStruct->e2outVptr));

	outStruct->Rout = 
		(double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, 1, dim, 
				IDL_ARR_INI_NOP, &(outStruct->RoutVptr));

	outStruct->flags = 
		(IDL_MEMINT *) IDL_MakeTempArray(IDL_TYP_MEMINT, 1, dim, 
				IDL_ARR_INI_NOP, &(outStruct->flagsVptr));

	return(outStruct);


}

#define ARRLEN(arr) (sizeof(arr)/sizeof(arr[0]))

int IDL_Load(void)
{

	/* This must be static. It is a struct. */
	/* The name in strings is the name by which it will be called from IDL and
	   MUST BE CAPITALIZED 
	   5th parameter will say if it accepts keywords and some other flags 
	   For more info see page 325 of external dev. guide */
	static IDL_SYSFUN_DEF2 proc_addr[] = {
		{ (IDL_SYSRTN_GENERIC) compea4, "COMPEA4", 0, IDL_MAXPARAMS, 0, 0},
	};

	/* False means it is not a function */
	return IDL_SysRtnAdd(proc_addr, IDL_FALSE, ARRLEN(proc_addr));

}
