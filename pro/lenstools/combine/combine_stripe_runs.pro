PRO combine_stripe_runs, stripe, clrs, nomake_scat=nomake_scat, hirata=hirata, hudson=hudson


  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: combine_stripe_runs, stripe, clrs, nomake_scat=nomake_scat, hirata=hirata, hudson=hudson'
      print,'Runs get_combie_stripe_runs to choose maximum rerun for each unique run'
      return
  ENDIF 

  IF keyword_set(hirata) THEN BEGIN
      hirstr = '_h' 
      hdr_hir = 'yes'
  ENDIF ELSE BEGIN 
      hirstr = ''
      hdr_hir = 'no'
  ENDELSE 

  get_combine_stripe_runs, stripe, runs, reruns, /silent
  clrstr = clrarr2string(clrs)

  print
  print,'#############################################################'
  print,'Doing stripe: '+ntostr(stripe)+' colors: '+clrstr

  dir = '/net/cheops1/data0/corrected/'
  nrun = n_elements(runs)
  FOR i=0L, nrun-1 DO BEGIN 
      
      rstr = ntostr(runs[i])
      rrstr = ntostr(reruns[i])
      file = dir + 'corr'+rstr+'/'+rrstr+'/combined/'+$
        'run'+rstr+'_'+rrstr+'_srcgal_'+clrstr+hirstr+'.fit'

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; we want to make the file if it doesn't exist
      ;; even if nomake_scat is set
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF ( NOT fexist(file) ) OR (NOT keyword_set(nomake_scat)) THEN BEGIN 

          make_scat53_meane, runs[i], reruns[i], clrs, hirata=hirata

      ENDIF 

  ENDFOR 

  combine_stripe, runs, reruns, clrs, hirata=hirata, hudson=hudson

END 
