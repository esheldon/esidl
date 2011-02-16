
int nParams(int argc);

typedef struct {
  IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */

  IDL_VPTR nrows;
  int nrows_there;

} KW_RESULT;

static IDL_KW_PAR kw_pars[] = {

  IDL_KW_FAST_SCAN, 

  /* These must be in alphabetical order */

  /* 
     IDL_KW_VIN: can be just an input, like var=3
     IDL_KW_OUT: can be output variable, var=var
     IDL_KW_ZERO: If not there, variable is set to 0

     If the VIN is set, then the _there variables will be
     0 when an undefined name variable is sent: no good for
     returning variables!
  */

  {"NROWS", IDL_TYP_UNDEF, 1, IDL_KW_OUT | IDL_KW_ZERO, 
   IDL_KW_OFFSETOF(nrows_there), IDL_KW_OFFSETOF(nrows) },

  { NULL }
};
