FUNCTION testsub::init, maxbcg_sample

  retval = self->maxbcg_correlate::init(maxbcg_sample)
  return, retval

END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; An example sub-sample
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; bin definitions
PRO testsub::ngals_bins, nbin, lowlim, highlim
  on_error, 2
  IF n_params() LT 3 THEN BEGIN 
      message,'-Syntax: ts->ngals_bins, nbin, lowlim, highlim'
  ENDIF 
  CASE nbin OF
      2: BEGIN 
          lowlim = [10, 20]
          highlim = [20, 300]
      END 
      ELSE: message,"Unsupported nbin = "+ntostr(nbin)
  ENDCASE 
END 

; Generate the were string
FUNCTION testsub::ngals_where_string, nbin
  on_error, 2
  IF n_elements(nbin) EQ 0 THEN BEGIN 
      message,'-Syntax: nbin=ts->ngals_where_string(subtype)'
  ENDIF 
  self->ngals_bins, nbin, lowlim, highlim  
  FOR i=0L, nbin-1 DO BEGIN 

      tstring = $
        'struct.ngals GE '+ntostr(lowlim[i])+' AND '+$
        'struct.ngals LT '+ntostr(highlim[i])

      add_arrval, tstring, where_string

  ENDFOR 
  return,where_string

END 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This code called by the correlation class. Just add 
; to the lists when you add new sub-samples
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Return number of bins for a subsample label.
FUNCTION testsub::subtype_nbin, subtype
  on_error, 2
  IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: nbin=ts->nbin(subtype)'
  ENDIF 
  CASE subtype OF
      'ngals2': return,2
      ELSE: message,'Unknown subtype: '+strn(subtype)
  ENDCASE 
END 

; Return where string for subsample label.
FUNCTION testsub::where_string, subtype, nbin=nbin
    IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: ws=mb->where_string(subtype)'
  ENDIF 

  nbin = self->subtype_nbin(subtype)

  cat = self->catalog()
  CASE subtype OF 
      'ngals2': BEGIN 
          return, self->ngals_where_string(nbin)
      END 
      ELSE: message,'Unknown type: '+ntostr(subtype)
  ENDCASE 
END 


PRO testsub__define
  struct = {                             $
             testsub,                    $
             some_dummy_variable: 0,     $
             INHERITS maxbcg_correlate  $
           }
END 
