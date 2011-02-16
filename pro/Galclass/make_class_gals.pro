PRO make_class_gals, run, rerun, clr, iorder, remove=remove

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: make_class_gals
;       
; PURPOSE: 
;   Make huge dilution corrected catalogs of stars and galaxies 
;   (split into separate files by color) 
;	
;
; CALLING SEQUENCE:
;      make_class_gals, run, rerun
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
      print,'-Syntax: make_class_gals, run, rerun, clr [, iorder, remove=remove]'
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

  maxmag = 18.
  minmag = 16.

  typ = 'blah1'+ntostr(long(systime(1)))
  typ2 = 'blah2'+ntostr(long(systime(1)))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Catalog of corrected shapes and positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  stmp=create_struct(name='typ', $
                     'petror50', fltarr(5), $
                     'petror90', fltarr(5), $
                     'nprof', fltarr(5), $
                     'profmean', fltarr(15,5), $
                     'classification',0)
  

  sdssidl_setup
  outdir=!sdss_shapecorr_dir+'corr'+rr+'/'+ntostr(rerun)+'/combined/'

  ;; Judith's directory
;  outdir = '/sdss3/data4/jracusin/'

GOTO,jump
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
      fetch_file_list,dir,tsobjfiles
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
          tsobj=0
          openr, unit, tsobjfiles[fi], /get_lun

          IF fi EQ 0 THEN BEGIN 
              tsobj = mrdfits3(unit, 1, 0, /silent)
          ENDIF ELSE BEGIN
              tsobj = mrdfits3(unit, 1, 0, /silent, /deja_vu)
          ENDELSE 
          free_lun, unit

          tmp = mrdfits(files[fi], 1, /silent)

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
                        (size GT psfsize+nsig*psfrms) AND $
                        tmp.objc_type EQ galtype, ngal)
                        
          x=0 & y=0 & size=0 & psfsize=0 & psfrms=0

          IF ngal NE 0 THEN BEGIN 
              ;; Lens Galaxies for this color
              w=where(tmp[gind].e1[clr] NE edef AND $
                      tmp[gind].e2[clr] NE edef, ngal)
              IF ngal NE 0 THEN BEGIN 
                  gind = gind[w]
                  w=where(tmp[gind].r[clr] LT .8 AND $
                          tmp[gind].r[clr] GE 0., ngal)
                  IF ngal NE 0 THEN BEGIN 
                      gind=gind[w]

                      photo_match, tmp[gind].run, tmp[gind].rerun, $
                                   tmp[gind].camcol, tmp[gind].field, $
                                   tmp[gind].id, $
                                   tsobj.run, tsobj.rerun, tsobj.camcol, $
                                   tsobj.field, tsobj.id, $
                                   m1, gind2
                                   
                      crap = replicate(stmp, ngal)
                      copy_struct, tsobj[gind2], crap
                      tmp2 = replicate(tmp[0], ngal)
                      copy_struct, tmp[gind], tmp2

                      tmpst=create_struct(name=typ2,tmp2[0], crap[0])
                      lenstmp = replicate(tmpst, ngal)
                      copy_struct, tmp2, lenstmp
                      copy_struct, crap, lenstmp

                      corr = 1./(1.- lenstmp.r[clr])
                      lenstmp.e1[clr] = lenstmp.e1[clr]*corr
                      lenstmp.e2[clr] = lenstmp.e2[clr]*corr
                      lenstmp.momerr[clr] = lenstmp.momerr[clr]*corr
              
                      numlist[fi] = ngal
                      ntotal = ntotal+ngal
                      ptrlist[fi] = ptr_new(lenstmp)

                      tmp2=0
                      crap=0

                  ENDIF
              ENDIF 
          ENDIF 
          IF fi MOD 20 EQ 0 THEN print,'Field: ',ntostr(field),'  ',$
            colors[clr]+$
            '-band Lens Galaxies: '+ntostr(ngal)
      ENDFOR 
      
      print,'Camcol: ',cstr,' total ',colors[clr],'-band galaxies: ',$
        ntostr(ntotal)
      lenses = replicate(tmpst, ntotal)
      beg=0
      FOR i=0,nfields-1 DO BEGIN 
          IF numlist[i] NE 0 THEN BEGIN 
              lenses(beg:beg+numlist(i)-1)=*ptrlist(i)
          ENDIF 
          ptr_free,ptrlist(i)   ;kill the heap variables associated 
                                ;with the pointer
          beg=beg+numlist(i)

      ENDFOR 
      
      outname = outdir + 'run'+rr+'_classgal_'+cstr+colors[clr]+nend
      print,'Output file: ',outname
      mwrfits, temporary(lenses), outname, outhdr, /create

      print
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Combine all columns together
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Combining columns into larger structures'
  print
  print,'Lens Galaxies'
  print
jump:

  firstra=20000.
  lastra=-20000.

  ptrlist = ptrarr(6)
  FOR col = colmin, colmax DO BEGIN
      cstr = ntostr(col)
      name = outdir + 'run'+rr+'_classgal_'+cstr+colors[clr]+nend
      tmpclr = mrdfits(name,1,hdr,structyp=typ2, /silent)
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
  tmp = replicate(tmpst, nt)

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
  outname = outdir + 'run'+rr+'_classgal_'+colors[clr]+nend
  mwrfits,tmp, outname, outhdr, /create


  ptime,systime(1)-time

  return
END 
