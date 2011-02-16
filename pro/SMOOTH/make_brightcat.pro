PRO make_brightcat, run, rerun, clr, iorder, remove=remove

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
      print,'-Syntax: make_fcat, run, rerun, clr, iorder'
      print,''
      print,'Use doc_library,"dilutecorr"  for more help.'  
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)

  galtype = 3

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

  maxmag = 16.
  minmag = 15.

  typ = 'blah2'+ntostr(long(systime(1)))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Catalog of corrected shapes and positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  brightstruct = create_struct(name=typ, $
                             'ra', double(0.), $
                             'dec',double(0.), $
                             'petrocounts', fltarr(5), $
                             'petrorad', 0., $
                             'e1', 0., $
                             'e2', 0., $
                             'uncert', 0., $
                             'r', 0., $
                             'photoz', 0.)

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
      ;; read in psf info
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      momname, run, col, selectclr, iorder, psffile
      psffile=fitdir+psffile
      print
      print,'Reading in moment file: ',psffile
      readcol,psffile, sjunk, junk,aa,bb,cc,dd,ee,ff,$
        junk,junk,junk,junk,junk,junk,junk,junk,junk,rms,$
        format='A,I,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D'          

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

          ;; Reddening correct
          tmp.petrocounts = tmp.petrocounts-tmp.reddening

          ;; Extract galaxies: nsig greater than psfsize

          x = tmp.colc[selectclr]
          y = tmp.rowc[selectclr]
          size = tmp.ixx[selectclr] + tmp.iyy[selectclr]
          psfsize = aa[ind] + bb[ind]*x   $
                            + cc[ind]*y   $
                            + dd[ind]*x*y $
                            + ee[ind]*x*x $
                            + ff[ind]*y*y
          psfrms = rms[ind]
          gind = where( tmp.petrocounts[selectclr] LT maxmag AND $
                        tmp.petrocounts[selectclr] GE minmag AND $
                        (size GT psfsize+nsig*psfrms), ngal)
                        
          x=0 & y=0 & size=0 & psfsize=0 & psfrms=0

          IF ngal NE 0 THEN BEGIN 
              ;; Bright Galaxies for this color
              w=where(tmp[gind].e1[clr] NE edef AND $
                      tmp[gind].e2[clr] NE edef, ngal)
              IF ngal NE 0 THEN BEGIN 
                  gind = gind[w]
                  w=where(tmp[gind].r[clr] LT .8 AND $
                          tmp[gind].r[clr] GE 0. AND $
                          tmp[gind].objc_type EQ galtype, ngal)
                  IF ngal NE 0 THEN BEGIN 
                      gind=gind[w]
          
                      brighttmp = replicate(brightstruct, ngal)
                      corr = 1./(1.- tmp[gind].r[clr])
                      brighttmp.e1 = tmp[gind].e1[clr]*corr
                      brighttmp.e2 = tmp[gind].e2[clr]*corr
                      brighttmp.uncert = tmp[gind].momerr[clr]*corr

                      brighttmp.r = tmp[gind].r[clr]
                      brighttmp.petrocounts = tmp[gind].petrocounts
                      brighttmp.petrorad    = tmp[gind].petrorad[clr]
                      brighttmp.ra  = tmp[gind].ra
                      brighttmp.dec = tmp[gind].dec
                      brighttmp.photoz = tmp[gind].photoz
              
                      numlist[fi] = ngal
                      ntotal = ntotal+ngal
                      ptrlist[fi] = ptr_new(brighttmp)
                  ENDIF
              ENDIF 
          ENDIF 
          IF fi MOD 20 EQ 0 THEN print,'Field: ',ntostr(field),'  ',$
            colors[clr]+$
            '-band Bright Galaxies: '+ntostr(ngal)
      ENDFOR 
      
      print,'Camcol: ',cstr,' total ',colors[clr],'-band galaxies: ',$
        ntostr(ntotal)
      brights = replicate(brightstruct, ntotal)
      beg=0
      FOR i=0,nfields-1 DO BEGIN 
          IF numlist[i] NE 0 THEN BEGIN 
              brights(beg:beg+numlist(i)-1)=*ptrlist(i)
          ENDIF 
          ptr_free,ptrlist(i)   ;kill the heap variables associated 
                                ;with the pointer
          beg=beg+numlist(i)

      ENDFOR 
      
      outname = outdir + 'run'+rr+'_brightgal_'+cstr+colors[clr]+nend
      print,'Output file: ',outname
      mwrfits, temporary(brights), outname, outhdr, /create

      print
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Combine all columns together
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Combining columns into larger structures'
  print
  print,'Bright Galaxies'
  print


  firstra=20000.
  lastra=-20000.

  ptrlist = ptrarr(6)
  FOR col = colmin, colmax DO BEGIN
      cstr = ntostr(col)
      name = outdir + 'run'+rr+'_brightgal_'+cstr+colors[clr]+nend
      tmpclr = mrdfits(name,1,hdr,structyp=typ, /silent)
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
  tmp = replicate(brightstruct, nt)

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
  outname = outdir + 'run'+rr+'_brightgal_'+colors[clr]+nend
  mwrfits,tmp, outname, outhdr, /create


  ptime,systime(1)-time

  return
END 
