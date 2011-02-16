PRO run_test_error_formulae, theta, outfile, maxerr=maxerr, etot=etot

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: run_test_error_formulae, theta_array'
      return
  ENDIF 

  IF n_elements(outfile) EQ 0 THEN BEGIN 
      outfile='/sdss2/data0/sdss/tmp/test_error_form_N1.fit'
      WHILE fexist(outfile) DO outfile=newname(outfile)
  ENDIF 
  print
  print,'Outfile = ',outfile
  print
  ntrial = 100L
  
  aratio = .5

  nobjx=10L
  nobjy=10L
  xsize=500L
  ysize=500L

  sigma=4.0                     ;pixels, width of gaussian object
  counts0 = 5000.
  countsstep = 1000.
  ncounts = 20

  counts=fltarr(ncounts)
  FOR i=0L, ncounts-1 DO counts[i]=counts0+i*countsstep

  ntheta=n_elements(theta)

  arrval=fltarr(ntrial*nobjx*nobjy)
  st = create_struct('inpute1',0.0,$
                     'inpute2',0.0,$
                     'theta',0.0,$
                     's2n', 0.0, $
                     'e1',arrval,$
                     'e2',arrval, $
                     'momerr',arrval,$
                     $
                     'meane1', 0.0,$
                     'meane2', 0.0,$
                     $
                     'mease1err', 0.0,$
                     'mease2err',0.0,$
                     'meanmomerr', 0.0, $
                     'meanmomerrerr', 0.0)

  struct=replicate(st, ncounts*ntheta)

  i=0L
  FOR ic=0L, ncounts-1 DO BEGIN 

      FOR it=0L, ntheta-1 DO BEGIN 

          IF 0 THEN BEGIN 
              test_error_formulae, ntrial, $
                nobjx, nobjy, xsize, ysize, $
                sigma, aratio, theta[it], counts[ic], $
                e1, e2, $
                momerr, inpute1, inpute2, $
                s2n
          ENDIF ELSE BEGIN 
              test_error_formulae_rdiff, ntrial, $
                nobjx, nobjy, xsize, ysize, $
                sigma, aratio, theta[it], counts[ic], $
                e1, e2, $
                momerr, inpute1, inpute2, $
                s2n, maxerr=maxerr, etot=etot
          ENDELSE 

          w=where(e1 NE 0.0,nw)
          
          
          ;; want to measure the scatter around 
          ;; true e, so do separate measure of noise that way
          sigmae=sqrt(0.32^2 + momerr[w]^2)
          wmom, e1[w], sigmae, meane1, mease1err
          sigmae=sqrt(0.32^2 + momerr[w]^2)
          wmom, e2[w], sigmae, meane2, mease2err

          ;; mean uncertainties, errors
          sig=replicate(1.0, nw)
          wmom, momerr[w], sig, meanmomerr, meanmomerrerr
          
          print
          print,'s/n',s2n
          print,'mean momerr',meanmomerr
          print,'meas e1err',mease1err
          print,'meas e2err',mease2err
          print,'meas mean e1: ',meane1
          print,'meas mean e2: ',meane2
          print,'input e1: ',inpute1
          print,'input e2: ',inpute2
          print
          
          struct[i].inpute1=inpute1
          struct[i].inpute2=inpute2
          struct[i].theta = theta[it]
          struct[i].s2n=s2n
          struct[i].e1=e1
          struct[i].e2=e2
          struct[i].momerr = momerr
          
          struct[i].meane1 = meane1
          struct[i].meane2 = meane2
          struct[i].mease1err=mease1err
          struct[i].mease2err=mease2err
          struct[i].meanmomerr=meanmomerr
          struct[i].meanmomerrerr=meanmomerrerr
          
          i=i+1L
      ENDFOR 

  ENDFOR 
      
  mwrfits, struct, outfile, /create


return
END 
