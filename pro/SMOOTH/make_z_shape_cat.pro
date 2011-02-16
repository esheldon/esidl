PRO make_z_shape_cat, run

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: make_z_shape_cat
;       
; PURPOSE: 
;   Make huge dilution corrected catalogs of stars and galaxies 
;   (split into separate files by color) as picked by the
;   routines extract_stars and extract_lensgalaxies
;	
;
; CALLING SEQUENCE:
;      make_z_shape_cat, run
;                 
;
; INPUTS: run: the run in integer form
;
; OUTPUTS: files for each color containing every object in the run.  Tags
;          are ra, dec, e1, e2, momerr
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	Erin Scott Sheldon UofMichigan 7/19/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
      print,'-Syntax: make_z_shape_cat, run'
      print,''
      print,'Use doc_library,"dilutecorr"  for more help.'  
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  selectclr = 2.                ; select based on red
  edef = 1.e10                  ; default value for e1 and e2 if not meas.

  rr = ntostr(run)
                                           
;  dir = '/sdss3/usrdevel/philf/run'+rr+'/'        ; read directory 
   dir = '/sdss3/usrdevel/philf/master/'
  outdir = '/sdss4/data1/esheldon/CORRECTED/'     ; write directory
  psfile=outdir+'_out.ps'
  nend = '.fit'

  typ2 = 'haha'+ ntostr(long(systime(1)))
  typ3 = 'hehe'+ntostr(long(systime(1)))

                                
  colmin = 1                    ; What columns to use
  colmax = 6
                                
  clrmin = 1                    ; What filters to use
  clrmax = 3

  colors = ['u','g','r','i','z']

  ;; output to postscript
  makeps, psfile

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; We know something about these runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cstr = ntostr(colmin)
  f = dir+'adat'+cstr+'cz.fit'
  fits_info, f, /silent, n_ext = maxf

  CASE run OF
    756: BEGIN
           start = 187
           nframes = 607
         END 
    752: BEGIN
           start = 2
           nframes = 607
         END
    ELSE: BEGIN
           start = 1
           nframes = maxf
          END 
  ENDCASE 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Catalog of corrected shapes and positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  struct = create_struct(name = typ2, 'e1',0., $
                         'e2',0., $
                         'uncert', 0., $
                         'ra', double(0.), $
                         'dec', double(0.), $
                         'petrocounts', 0., $
                         'photoz', 0.)
  
  lensstruct = create_struct(name=typ3, $
                             'ra', double(0.), $
                             'dec',double(0.), $
                             'petrocounts', fltarr(5), $
                             'petrorad', 0., $
                             'e1', 0., $
                             'e2', 0., $
                             'uncert', 0., $
                             'r', 0., $
                             'photoz', 0.)

  ; Just get the stuff we need
  tl = ['OBJC_TYPE','OBJC_FLAGS','FLAGS','PETROCOUNTS',$
        'PETRORAD','RA','DEC','E1','E2','R','MOMERR', 'PHOTOZ']

  nlgtot = 0L
  nsgtot = 0L
  nstot = 0L
  FOR col = colmin, colmax DO BEGIN

      tmp = 0.
      cstr = ntostr(col)
      typ = 'blah'+ntostr(long(systime(1)))+cstr
      file=dir+'adat'+cstr+'cz.fit'

      read_photo_col, file, tmp, $
                      taglist=tl, $
                      start=start, $
                      nframes=nframes, $
                      maxf=maxf, $
                      struct_typ = typ

      print
      print,'---------------------- Stars -----------------------'
      ; Extract the stars
      extract_stars, tmp, selectclr, ind1    
  
      print
      print,'---------------------- Src Gal -----------------------'
      ; Extract sources
      extract_lensgalaxies, tmp, selectclr, ind2, $ 
                      max_mag=22., $
                      min_mag=18.
      wind2 = where(tmp[ind2].photoz NE -10., n2)
      ind2 = ind2[wind2]

      print
      print,'---------------------- Lens Gal -----------------------'
      ; Extract lenses
      extract_lensgalaxies, tmp, selectclr, ind3, $ 
                      max_mag=18.0, $
                      min_mag=16.0
      wind3 = where(tmp[ind3].photoz NE -10., n3)
      ind3 = ind3[wind3]

      FOR clr = clrmin, clrmax DO BEGIN 

          ;; stars
          w=where(tmp[ind1].e1[clr] NE edef AND $
                  tmp[ind1].e2[clr] NE edef, n1)
          cind1=ind1[w]

          tmpclr = 0.
          tmpclr = replicate(struct, n1)

          tmpclr.e1 = tmp[cind1].e1[clr]
          tmpclr.e2 = tmp[cind1].e2[clr]
          tmpclr.uncert = tmp[cind1].momerr[clr]
          tmpclr.ra = tmp[cind1].ra
          tmpclr.dec = tmp[cind1].dec

          outname = outdir + 'run'+rr+'_zstar_'+cstr+colors[clr]+nend
          mwrfits, tmpclr, outname, /create
      
          ;; Source Galaxies
          IF clr NE selectclr THEN BEGIN 
              w=where(tmp[ind2].e1[clr] NE edef AND $
                      tmp[ind2].e2[clr] NE edef AND $
                      tmp[ind2].r[clr] LT .8, n2)
              cind2=ind2[w]
          ENDIF ELSE BEGIN
              n2=n_elements(ind2)
              cind2=ind2
          ENDELSE 

          tmpclr = 0.
          tmpclr = replicate(struct, n2)
          corr = 1./(1. - tmp[cind2].r[clr])
          tmpclr.e1 = tmp[cind2].e1[clr]*corr
          tmpclr.e2 = tmp[cind2].e2[clr]*corr
          tmpclr.uncert = tmp[cind2].momerr[clr]*corr

          tmpclr.petrocounts    = tmp[cind2].petrocounts[clr]

          tmpclr.ra = tmp[cind2].ra
          tmpclr.dec = tmp[cind2].dec
          tmpclr.photoz = tmp[cind2].photoz

          outname = outdir + 'run'+rr+'_zsrcgal_'+cstr+colors[clr]+nend
          mwrfits, tmpclr, outname, /create
        
          ;; Lens Galaxies
          IF clr NE selectclr THEN BEGIN 
              w=where(tmp[ind3].e1[clr] NE edef AND $
                      tmp[ind3].e2[clr] NE edef AND $
                      tmp[ind3].r[clr] LT .8, n3)
              cind3=ind3[w]
          ENDIF ELSE BEGIN
              n3=n_elements(ind3)
              cind3=ind3
          ENDELSE 

          tmpclr=0.
          tmpclr = replicate(lensstruct, n3)
          corr = 1./(1.- tmp[cind3].r[clr])
          tmpclr.e1 = tmp[cind3].e1[clr]*corr
          tmpclr.e2 = tmp[cind3].e2[clr]*corr
          tmpclr.uncert = tmp[cind3].momerr[clr]*corr

          tmpclr.r = tmp[cind3].r[clr]
          tmpclr.petrocounts = tmp[cind3].petrocounts
          tmpclr.petrorad    = tmp[cind3].petrorad[clr]
          tmpclr.ra  = tmp[cind3].ra
          tmpclr.dec = tmp[cind3].dec
          tmpclr.photoz = tmp[cind3].photoz

          outname = outdir + 'run'+rr+'_zlensgal_'+cstr+colors[clr]+nend
          mwrfits, tmpclr, outname, /create

          IF clr EQ 2 THEN BEGIN 
              nstot = nstot+n1
              nsgtot = nsgtot+n2
              nlgtot = nlgtot+n3
          ENDIF 

      ENDFOR 
  ENDFOR 

  print
  print,ntostr(nstot)+' stars in r band of run '+rr
  print,ntostr(nsgtot)+' source galaxies in r band of run '+rr
  print,ntostr(nlgtot)+' lens galaxies in r band of run '+rr

  ;; Source Galaxies
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Combine all columns together for each color 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Combining columns into larger structures'
  print
  print,'Source Galaxies'
  print

  FOR clr = clrmin, clrmax DO BEGIN 

      ptrlist = ptrarr(6)
      FOR col = colmin, colmax DO BEGIN
          cstr = ntostr(col)
          name = outdir + 'run'+rr+'_zsrcgal_'+cstr+colors[clr]+nend
          tmpclr = mrdfits(name,1,hdr, structyp=typ2,/silent)
          spawn,'rm '+name

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
      outname = outdir + 'run'+rr+'_zsrcgal_'+colors[clr]+nend
      mwrfits,tmp, outname, /create
  ENDFOR 

  ;; Lens Galaxies

  print
  print,'Lens Galaxies'
  print

  FOR clr = clrmin, clrmax DO BEGIN 

      ptrlist = ptrarr(6)
      FOR col = colmin, colmax DO BEGIN
          cstr = ntostr(col)
          name = outdir + 'run'+rr+'_zlensgal_'+cstr+colors[clr]+nend
          tmpclr = mrdfits(name,1,hdr, structyp=typ3,/silent)
          spawn,'rm '+name

          ptrlist(col-colmin) = ptr_new(tmpclr)
          n=n_elements(tmpclr)
          IF col EQ 1 THEN numlist = n ELSE numlist = [numlist, n]
      ENDFOR 

      nt = total(numlist)
      tmp = replicate(lensstruct, nt)

      beg = 0
      FOR col = colmin, colmax DO BEGIN
          num = numlist[col-colmin]
          tmp[ beg:beg+num-1 ] = *ptrlist[col-colmin]
          ptr_free,ptrlist[col-colmin]
          beg = beg + num
      ENDFOR 
      outname = outdir + 'run'+rr+'_zlensgal_'+colors[clr]+nend
      mwrfits,tmp, outname, /create
  ENDFOR 

  ;; Stars
  ;; Combine all columns together for each color 

  print,'stars'
  print

  FOR clr = clrmin, clrmax DO BEGIN 

      ptrlist = ptrarr(6)
      FOR col = colmin, colmax DO BEGIN
          cstr = ntostr(col)
          name = outdir + 'run'+rr+'_zstar_'+cstr+colors[clr]+nend
          tmpclr = mrdfits(name,1,hdr, structyp=typ2,/silent)
          spawn,'rm '+name
          
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
      outname = outdir + 'run'+rr+'_zstar_'+colors[clr]+nend
      mwrfits,tmp, outname, /create
  ENDFOR 

  ;; close postscript file
  ep

  return
END 
