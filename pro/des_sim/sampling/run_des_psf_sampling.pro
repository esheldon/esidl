PRO run_des_psf_sampling_outfile, seeing, arcperpix, setnum, stFile, psFile

  front = 'sampling'+strn(setnum, length=2, padchar='0')

  IF sdssidl_config('hostname') EQ 'treebeard' THEN BEGIN 
      outDir = '/data0/des_sim/sampling_errors/'
  ENDIF ELSE BEGIN 
      outDir = '/net/cheops1/data0/esheldon/des_sim/sampling_errors/'
  ENDELSE 

  ;; ucent: uniform noise added to centroid
  seeing_string = ntostr(seeing, 4, /round)
  arcperpix_string = ntostr(arcperpix, 5, /round)
  
  stFile = $
    outDir + front + $
    '_psf'+seeing_string+'_pixel'+arcperpix_string+'.st'

  psFile = $
    outDir + front + $
    '_psf'+seeing_string+'_pixel'+arcperpix_string+'.eps'

END 

PRO run_des_psf_sampling, seeing_set, apix_set, setnum

  tm = systime(1)

  nstarx = 400L
  
  des_psf_sampling_getsets, seeing_set, apix_set, seeing, apix
  print,apix


  ;; Shape

  ;; Using multiple of these is important for checking the *bias* 
  ;; as a function of scale
;  etot = [0.02, 0.05, 0.1, 0.15, 0.2]
;  posAngle = [-reverse(findgen(7)), findgen(7)]*15.*!pi/180.
;  etot = 0.1
;  posAngle = 22.5               ; degrees
;  posAngle = [-90.0, -67.5, -45.0, -22.5, $
;              0.0, $
;               22.5,  45.0,  67.5,  90.0]

  etot = 0.1
  posAngle = -90. + 11.25*findgen(17)

  napix   = n_elements(apix)


  FOR ip=0L, napix-1 DO BEGIN 

      run_des_psf_sampling_outfile, $
        seeing, apix[ip], setnum, $
        stFile, psFile
      
      print
      print,'Output st file: ',stFile
      
      
      begplot,name=psFile,/color,/encap
      des_psf_sampling, seeing, apix[ip], etot, posAngle, $
        nstarx, struct, /doplot
      endplot,/trim_bbox
      
      print
      print,'Writing to file: ',stFile
      write_idlstruct, struct, stFile
      
  ENDFOR 


  ptime,systime(1)-tm

END 
