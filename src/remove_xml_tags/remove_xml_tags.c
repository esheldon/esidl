/*
 *
 * remove_xml_tags
 *
 * Remove xml tags strings <...> from the input string
 * Keep track of <TR> tags and return optionally through the nrows keyword
 *
 * From IDL:  newstring = remove_xml_tags(string, nrows=nrows)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "idl_export.h"

#include "remove_xml_tags.h"

KW_RESULT kw;

IDL_VPTR remove_xml_tags(int argc, IDL_VPTR *argv, char *argk)
{

  /* The input string */
  IDL_VPTR stringVptr;
  char *string;
  int slen;

  /* The output string */
  IDL_VPTR newStringVptr;
  char *newstring;

  /* Byte values for characters of interest */
  static char char_left_bracket = 60;
  static char char_right_bracket = 62;
  static char char_space = 32;
  static char char_R = 82;
  static char char_r = 114;
  static char char_T = 84;
  static char char_t = 116;

  /* Are we within a tag? */
  char inTag=0;

  /* Keep track of first two tags letters */
  char tagChar1;
  char tagChar2;
  short TRCheck=0;

  /* How many rows */
  int nrows=0;

  int i;
  
  /* Get the keywords */
  (void) IDL_KWProcessByOffset(argc, argv, argk, kw_pars, 
			       (IDL_VPTR *) 0, 1, &kw);

  if (nParams(argc) != 1)
    {
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		  "Syntax: newstring = remove_xml_tags(oldstring, nrows=)");
    }

  stringVptr = argv[0];
  if (stringVptr->type != IDL_TYP_STRING)
    {
      IDL_MESSAGE(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		  "Input parameter must be a string");
    }
  slen = stringVptr->value.str.slen;
  string = stringVptr->value.str.s;

  /*
    Copy to output variable.  Equivalent to working in-place
    if the input variable is temporary
  */

  if (stringVptr->flags & IDL_V_TEMP)
    newStringVptr = stringVptr;
  else 
    {
      newStringVptr = IDL_Gettmp();
      IDL_VarCopy(stringVptr, newStringVptr);
    }
  newstring = newStringVptr->value.str.s;

  /* Note, we need to use newstring here, just in case
     the input was temporary, since we have just copied
     the *pointer* and the other variable no longer exists! */

  for (i=0; i<slen; i++)
    {

      /* inTag will also contain the index into the tag, i.e. the
	 number of characters in the tag minus 1 */

      if (inTag != 0) 
	{

	  /* We have reached the end of the tag.  */
	  if (newstring[i] == char_right_bracket)
	    {
	      newstring[i] = char_space;

	      if ( (tagChar1 == char_T && tagChar2 == char_R) ||
		   (tagChar1 == char_t && tagChar2 == char_r) )
		nrows += 1;

	      tagChar1=0;
	      tagChar2=0;

	      inTag = 0;
	    }
	  /* We are still in the tag */
	  else
	    {
	      if (inTag == 1) tagChar1=newstring[i];
	      else if (inTag == 2) tagChar2=newstring[i];

	      newstring[i] = char_space;
	      inTag += 1;
	    }
	}
      /* Not in a tag */
      if (newstring[i] == char_left_bracket) 
	{
	  newstring[i] = char_space;
	  inTag = 1;
	}

    }

  if (kw.nrows_there) {
    IDL_StoreScalarZero(kw.nrows, IDL_TYP_LONG);
    kw.nrows->value.l = nrows;
  }

  return(newStringVptr);

}

/*===========================================================================
 *
 * nParams
 *
 * How many positional arguments were sent?
 *
 *===========================================================================*/

int nParams(int argc)
{

  int nKeywords;

  nKeywords =     
    kw.nrows_there;

  return 
    argc - nKeywords;

}



#define ARRLEN(arr) (sizeof(arr)/sizeof(arr[0]))

int IDL_Load(void)
{

  /* This must be static. It is a struct. */
  /* The name in strings is the name by which it will be called from IDL and
     MUST BE CAPITALIZED 
     5th parameter will say if it accepts keywords and some other flags 
     For more info see page 325 of external dev. guide */
  static IDL_SYSFUN_DEF2 function_addr[] = {
    { (IDL_SYSRTN_GENERIC) remove_xml_tags, "REMOVE_XML_TAGS", 
      0, IDL_MAXPARAMS, 
      IDL_SYSFUN_DEF_F_KEYWORDS, 0},
  };

  /* False means it is not a function */
  return IDL_SysRtnAdd(function_addr, IDL_TRUE, ARRLEN(function_addr));

}
