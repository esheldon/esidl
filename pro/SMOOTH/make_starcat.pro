PRO make_starcat, run, rerun, clr, iorder, remove=remove

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: make_shape_cat
;       
; PURPOSE: 
;   Make huge dilution corrected catalogs of stars and galaxies 
;   (split into separate files by color) 
;	
;
; CALLING SEQUENCE:
;      make_shape_cat, run, rerun
;                 
;
; INPUTS: run: the run in integer form
;         rerun: rerun in integer form
;
; OUTPUTS:
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	Erin Scott Sheldon UofMichigan 3/5/2000
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 3 THEN BEGIN 
      print,'-Syntax: make_starcat, run, rerun, clr, iorder'
      print,''
      print,'Use doc_library,"dilutecorr"  for more help.'  
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)

  IF NOT keyword_set(remove) THEN remove=0
  IF n_elements(iorder) EQ 0 THEN iorder=2
  newfront = 'adatc'
  selectclr = 2.                ; select based on red
  edef = 1.e10                  ; default value for e1 and e2 if not meas.

  rr = ntostr(run)
                                           
  nend = '.fit'
                                
  colmin = 1                    ; What columns to use
  colmax = 6
                                
  colors = ['u','g','r','i','z']

  nsig = 4.

  maxmag = 22.
  minmag = 18.

  typ = 'blah3'+ntostr(long(systime(1)))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Catalog of corrected shapes and positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  struct = create_struct(name=typ, $
                         'e1',0., $
                         'e2',0., $
                         'uncert', 0., $
                         'ra', double(0.), $
                         'dec', double(0.), $
                         'petrocounts', 0.)  

  sdssidl_setup
  outdir=!sdss_shapecorr_dir+'corr'+rr+'/'+ntostr(rerun)+'/combined/'
  FOR col = colmin, colmax DO BEGIN

      cstr = ntostr(col)
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Setup the column
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tmp = 0
      fetch_dir, run, col, rerun, dir, atldir, $
        corrdir=corrdir, corratldir=fitdir

      fetch_file_list,corrdir, files, fnums,front=newfront, $
        fieldmin=fieldmin, fieldmax=fieldmax
      nfields = n_elements(fnums)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; read in corrected fields
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      print,'reading from ',corrdir
      print,'outputting to ',outdir

      numlist=lonarr(nfields)
      ptrlist=ptrarr(nfields)   ;the array of pointers
      ntotal=0L

      outhdr = ['END']
      FOR fi=0, nfields-1 DO BEGIN 

          field = fnums[fi]
          ind = fi*4L

          tmp=0
          openr, unit, files[fi], /get_lun
          IF fi EQ 0 THEN BEGIN 
              tmp = mrdfits3(unit, 1, 0, /silent)
          ENDIF ELSE BEGIN
              tmp = mrdfits3(unit, 1, 0, /silent, /deja_vu)
          ENDELSE 
          free_lun, unit

          IF fi EQ 0 THEN BEGIN
              firstra = min(tmp.ra)
              sxaddpar, outhdr, 'FIRSTRA', firstra
          ENDIF ELSE IF fi EQ nfields-1 THEN BEGIN 
              lastra  = max(tmp.ra)
              sxaddpar, outhdr, 'LASTRA', lastra
          ENDIF 

          ;; Select objects by r-band magnitude.  This is probably not a good
          ;; idea for clusters because you will get a lot of contamination.

          ;; Extract stars


          stind = where( tmp.petrocounts[selectclr] LT maxmag AND $
                         tmp.petrocounts[selectclr] GE minmag AND $
                         tmp.starflag[selectclr] EQ 1, nstar)
                        

          IF nstar NE 0 THEN BEGIN 

              ;; Stars for this color
              w=where(tmp[stind].e1[clr] NE edef AND $
                      tmp[stind].e2[clr] NE edef, nstar)
              IF nstar NE 0 THEN BEGIN 
                  stind = stind[w]
                  w=where(tmp[stind].r[clr] LT .8 AND $
                          tmp[stind].r[clr] GE 0., nstar)
                  IF nstar NE 0 THEN BEGIN 
                      stind=stind[w]
          
                      startmp = replicate(struct, nstar)
                      corr = 1./(1.- tmp[stind].r[clr])
                      startmp.e1 = tmp[stind].e1[clr]*corr
                      startmp.e2 = tmp[stind].e2[clr]*corr
                      startmp.uncert = tmp[stind].momerr[clr]*corr

                      startmp.petrocounts = tmp[stind].petrocounts[clr]
                      startmp.ra  = tmp[stind].ra
                      startmp.dec = tmp[stind].dec
              
                      numlist[fi] = nstar
                      ntotal = ntotal+nstar
                      ptrlist[fi] = ptr_new(startmp)
                  ENDIF
              ENDIF 
          ENDIF 
          IF fi MOD 20 EQ 0 THEN print,'Field: ',ntostr(field),'  ',$
            colors[clr]+$
            '-band Stars: '+ntostr(nstar)
      ENDFOR 
      
      print,'Camcol: ',cstr,' total ',colors[clr],'-band galaxies: ',$
        ntostr(ntotal)
      stars = replicate(struct, ntotal)
      beg=0
      FOR i=0,nfields-1 DO BEGIN 
          IF numlist[i] NE 0 THEN BEGIN 
              stars(beg:beg+numlist(i)-1)=*ptrlist(i)
          ENDIF 
          ptr_free,ptrlist(i)   ;kill the heap variables associated 
                                ;with the pointer
          beg=beg+numlist(i)

      ENDFOR 
      
      outname = outdir + 'run'+rr+'_star_'+cstr+colors[clr]+nend
      print,'Output file: ',outname
      mwrfits, temporary(stars), outname, outhdr, /create

      print
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Combine all columns together
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Combining columns into larger structures'
  print
  print,'Stars'
  print


  firstra=20000.
  lastra=-20000.

  ptrlist = ptrarr(6)
  FOR col = colmin, colmax DO BEGIN
      cstr = ntostr(col)
      name = outdir + 'run'+rr+'_star_'+cstr+colors[clr]+nend
      tmpclr = mrdfits(name,1,hdr, structyp=typ, /silent)
      ;; get first, last ra of all columns
      fra = sxpar(hdr, 'FIRSTRA')
      lra = sxpar(hdr, 'LASTRA')

      IF fra LT firstra THEN firstra = fra
      IF lra GT lastra THEN lastra = lra
      IF remove THEN spawn,'rm '+name

      ptrlist(col-colmin) = ptr_new(tmpclr)
      n=n_elements(tmpclr)
      IF col EQ 1 THEN numlist = n ELSE numlist = [numlist, n]
  ENDFOR 

  nt = total(numlist)
  tmp = replicate(struct, nt)

  beg = 0
  FOR col = colmin, colmax DO BEGIN
      num = numlist[col-colmin]
      tmp[ beg:beg+num-1 ] = *ptrlist[col-colmin]
      ptr_free,ptrlist[col-colmin]
      beg = beg + num
  ENDFOR 

  outhdr = ['END']
  sxaddpar, outhdr, 'FIRSTRA', firstra
  sxaddpar, outhdr, 'LASTRA', lastra
  outname = outdir + 'run'+rr+'_star_'+colors[clr]+nend
  mwrfits,tmp, outname, outhdr, /create


  ptime,systime(1)-time

  return
END 
