FUNCTION spectra_name, stripes, indir=indir, spec=spec, lrg=lrg, lss=lss, letter=letter, post=post, sample=sample, all=all, southern=southern
  
  np=n_params()
  IF np LT 1 THEN BEGIN 
      print,'-Syntax: filename = spectra_name, stripes, indir=indir, /spec, /lrg, /lss, letter=letter, post=post, sample=sample, /all, /southern'
      return,''
  ENDIF 

  ;; The directory
  IF n_elements(indir) EQ 0 THEN BEGIN 
      IF keyword_set(spec) THEN BEGIN 
          combdir = sdssidl_config('SHAPECORR_DIR') + 'combined/'
          indir = combdir + 'spec/'
      ENDIF ELSE BEGIN 
          indir = vagc_lensinput_dir(lss=lss)
      ENDELSE 
  ENDIF 

  ;; Various samples
  IF keyword_set(spec) THEN BEGIN 
      catname = 'spec'
  ENDIF ELSE BEGIN 
      catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)
  ENDELSE 

  IF keyword_set(lrg) THEN BEGIN 
      catname = 'lrg_'+catname
  ENDIF ELSE IF keyword_set(southern) THEN BEGIN 
      catname = 'southern_'+catname
  ENDIF 

  IF keyword_set(all) THEN BEGIN 
      stripe_string = '09_10_11_12_13_14_15_16_26_27_28_29_30_31_32_33_34_35_36_37_42_43_44_76_82_86'
  ENDIF ELSE BEGIN 
      ;; Convert possible array of stripes to a string
      stripe_string = stripearr2string(stripes)
  ENDELSE 

  ;; The file
  name = indir + 'stripe'+stripe_string+'_'+catname+'.fit'

  return, name
END 
