
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    SDSS_SURVEY_ROT
;       
; PURPOSE:
;    Find field-by-field rotation angle between CCD coordinates (row,col)
;    and SDSS corrected survey coord. (clambda,ceta) for each bandpass.
;
; CALLING SEQUENCE:
;    sdss_survey_rot, run, rerun, camcol [,rotstruct,outdir=outdir]
;
; INPUTS: 
;    run: The SDSS run in integer form
;    rerun: The SDSS rerun in integer form
;    camcol: The SDSS camcol in integer form
;
; OPTIONAL INPUTS:
;    NONE (see KEYWORD PARAMETERS)
;
; KEYWORD PARAMETERS:
;    outdir: The output directory (full path with / on the end). Default is 
;            CWD.
;       
; OUTPUTS: 
;    fits file containing the rotation structure (see optional outputs)
;
; OPTIONAL OUTPUTS:
;    rotstruct: A structure containing the rotation angle for each field.
;      IDL> help,rotstruct,/str
;      ** Structure <40421c88>, 4 tags, length=64, refs=1:
;         FIELD           LONG                11
;         HANGLE          FLOAT     Array[5]
;         VANGLE          FLOAT     Array[5]
;         ANGLE           FLOAT     Array[5]
;    ANGLE is the one used. It is the mean of angles calculated in the 
;    vertical and horizontal directions (VANGLE, HANGLE), which are slightly 
;    different due to the small deviation from perfect rotation.
;
; CALLED ROUTINES:
;    sdss_read()
;    ROWCOL2MUNU
;    GC2SURVEY
;    RD2XY
;    MWRFITS
;    APLOT
;    LEGEND
;    SIMPCTABLE
; 
; PROCEDURE: 
;    Convert a set of (row,col) to (lambda,eta) and find rotation
;    angle between the two coord. systems.
;	
;
; REVISION HISTORY:
;    Creation date: 23-OCT-2000 Erin Scott Sheldon
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO sdss_survey_rot,run,rerun,camcol, indir=indir, outdir=outdir, status=status

  ;; status is 1 unless we reach end
  status=1

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: sdss_survey_rot, run, rerun, camcol, status=status'
      print,''
      print,'Use doc_library,"sdss_survey_rot"  for more help.'  
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output file  names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  run_status = sdss_runstatus()
  w=where(run_status.run EQ run AND run_status.stripe NE -1,nw)

  IF n_elements(outdir) EQ 0 THEN BEGIN 
      fetch_dir, run, camcol, rerun, dir, corratldir=outdir
  ENDIF 

  outfile = outdir + 'surveyrot_'+run2string(run)+'_'+ntostr(camcol)+'_'+$
    ntostr(rerun)+'.st'
  asciifile = outdir + 'surveyrot_'+run2string(run)+'_'+ntostr(camcol)+'_'+$
    ntostr(rerun)+'.tab'

  psname = outdir + 'surveyrot_'+run2string(run)+'_'+ntostr(camcol)+'_'+$
    ntostr(rerun)+'_N1.ps'

  WHILE exist(psname) DO BEGIN
      psname = newname(psname)
  ENDWHILE 

  print
  print,'Output File: ',outfile
  print,'PS File: ',psname
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; start outputting postscript data
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gratio = 1.36603
  begplot, name=psname, /color
  ;simpctable
 
  clr = [0,1,2,3,4]
  plotclr = intarr(5)
  plotclr[0] = !cyan
  plotclr[1] = !green
  plotclr[2] = !magenta
  plotclr[3] = !red
  plotclr[4] = !p.color

  colors = ['u','g','r','i','z']

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; astrometry structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  radeg = 180./!dpi
  val = .4d/3600.
  astr={cd: double(identity(2)),$
        cdelt: [val,val], $
        crpix: dblarr(2), $
        crval: [0., 0.]}
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define set of pixels to rotate
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lrow = 1489-1
  lcol = 2048-1
  nnr=20.
  vert=dblarr(2,nnr)
  horz=dblarr(2,nnr)
  ccol = 1023.5
  crow = 744.0
  FOR i=0L, nnr-1 DO BEGIN
      vert[0,i] = ccol          ;columns
      vert[1,i] = lrow-i        ;rows
      IF i EQ 10 THEN vert[1,i] = 1000
      horz[0,i] = lcol - i      ;columns
      horz[1,i] = crow          ;rows
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nclr = n_elements(clr)
  arrval = fltarr(nclr)
  rotst=create_struct('fieldid', 0ULL,$
                      'run',0L,$
                      'rerun',0,$
                      'camcol',0b,$
                      'field',0,$
                      'hangle', arrval, $
                      'vangle', arrval, $
                      'angle', arrval)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; loop over bandpasses
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR ic=0L, nclr-1 DO BEGIN 

      ;; read in asTrans file
      trans=sdss_read('astrans',run,camcol,band=clr[ic],rerun=rerun, $
                      node=node,inc=inc,status=astatus)
      IF if astatus ne 0 THEN BEGIN
          endplot,/noprint
          return
      ENDIF 

      ;; Make output structure
      IF ic EQ 0 THEN BEGIN 
          fields = trans.field
          nf = n_elements(fields)
          rotstruct = replicate(rotst, nf)

          rotstruct.run = run
          rotstruct.rerun = rerun
          rotstruct.camcol = camcol
          rotstruct.field = fields

          rotstruct.fieldid = sdss_photoid(rotstruct.run,   $
                                           rotstruct.rerun, $
                                           rotstruct.camcol,$
                                           rotstruct.field)
      ENDIF 

      mvert = fltarr(nf)
      mhorz = fltarr(nf)
      angle = fltarr(nf)
      FOR fi=0L, nf-1 DO BEGIN 
          field = fields[fi]
          
          ;; convert center (row,col) to (mu,nu) great circle coords
          rowcol2munu, trans, field, crow, ccol, cmu, cnu
          ;; convert (mu,nu) to (clambda,ceta) corrected survey coords.
          gc2csurvey, cmu, cnu, node, inc, clam, ceta

          ;; same for the vertical/horizontal points
          rowcol2munu, trans, field, vert[1,*], vert[0,*], vmu, vnu
          rowcol2munu, trans, field, horz[1,*], horz[0,*], hmu, hnu
          gc2csurvey, vmu, vnu, node, inc, vlam, veta
          gc2csurvey, hmu, hnu, node, inc, hlam, heta

          ;; tangent project
          astr.crval = [ceta, clam]
          rd2xy,ceta,clam,astr,cxx,cyy
          rd2xy, veta, vlam, astr, vxx, vyy
          rd2xy, heta, hlam, astr, hxx, hyy

          ;; angle of points in (lambda,eta)
          lametaaang = atan(vyy, vxx)
          lametabang = atan(hyy, hxx)

          mvert[fi] = mean(!dpi/2. - lametaaang)
          mhorz[fi] = mean(0. - lametabang)
          angle[fi] = mean( [!dpi/2. - lametaaang, 0. - lametabang] )

      ENDFOR 

      ;; begin plotting
      IF ic EQ 0 THEN BEGIN 
          xtitle = 'Field'
          ytitle = 'Rotation Angle (degrees)'
          title = 'Run '+run2string(run)+' Camcol '+ntostr(camcol)
          yrange=[-3,3]
          aplot,gratio,fields,fltarr(nf),yrange=yrange,ystyle=1,$
            xtitle=xtitle,ytitle=ytitle,title=title
          oplot,fields,radeg*angle,color=plotclr[ic]
      ENDIF ELSE BEGIN 
          oplot,fields,radeg*angle,color=plotclr[ic]
      ENDELSE 
      
      rotstruct.hangle[ic] = mhorz
      rotstruct.vangle[ic] = mvert
      rotstruct.angle[ic] = angle
  ENDFOR 

  legend,colors,linestyle=[0,0,0,0,0],color=plotclr,/center,/top, $
    thick=[5,5,5,5,5],charsize=1

  endplot,/noprint

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output fits file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Writing output file: ',outfile
  write_idlstruct, rotstruct, outfile
  print
  print,'Writing ascii file: ',asciifile
  ascii_write, rotstruct, asciifile
;  mwrfits, rotstruct, outfile, /create

  status=0

  return
END 









