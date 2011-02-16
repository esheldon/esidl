PRO test_phototype, run, rerun, clr, iorder, nsig=nsig, remove=remove

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: make_shape_cat
;       
; PURPOSE: 
;   Make huge dilution corrected catalogs of source galaxies
;	
; CALLING SEQUENCE:
;      make_shape_cat, run, rerun, clr, iorder, typecut=typecut, 
;           remove=remove
;                 
;
; INPUTS: run: the run in integer form
;         rerun: rerun in integer form
;         clr: the color index from which to make the catalog
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
      print,'-Syntax: make_scat, run, rerun, clr, iorder'
      print,''
      print,'Use doc_library,"dilutecorr"  for more help.'  
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)

  galtype = 3

  IF NOT keyword_set(typecut) THEN typecut=0
  IF typecut THEN BEGIN
      print
      print,'Will cut on photo type'
      print
  ENDIF 
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

  IF n_elements(nsig) EQ 0 THEN nsig = 4.

  maxmag = 22.
  minmag = 18.

  typ = 'blah1'+ntostr(long(systime(1)))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Catalog of corrected shapes and positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  struct = create_struct(name=typ,$
                         'rho4',fltarr(5), $
                         'momerr', 0., $
                         'e1', 0., $
                         'petrocounts', fltarr(5), $
                         'objc_type', 0)
  sdssidl_setup
  outdir=!sdss_shapecorr_dir+'corr'+rr+'/'+ntostr(rerun)+'/combined/'
  FOR col = 3, 3 DO BEGIN

      IF run EQ 1140 THEN BEGIN 
          IF col EQ 2 THEN GOTO,jump
      ENDIF 

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

      numlist1=lonarr(nfields)
      ptrlist1=ptrarr(nfields)   ;the array of pointers
      ntotal1=0L

      numlist2=lonarr(nfields)
      ptrlist2=ptrarr(nfields)   ;the array of pointers
      ntotal2=0L

      outhdr = ['END']

fistart=100
nfields = 600
      FOR fi=fistart, nfields-1 DO BEGIN 

          field = fnums[fi]
          ind = fi*4L

          tmp=0
          openr, unit, files[fi], /get_lun
          IF fi EQ fistart THEN BEGIN 
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
          tmp.counts_model = tmp.counts_model-tmp.reddening

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

              ;; Source Galaxies for this color
              w=where(tmp[gind].e1[clr] NE edef AND $
                      tmp[gind].e2[clr] NE edef, ngal)
              IF ngal NE 0 THEN BEGIN 
                  gind = gind[w]
                  w=where(tmp[gind].r[clr] LT .8 AND $
                          tmp[gind].r[clr] GE 0., ngal)
                  IF ngal NE 0 THEN BEGIN 
                      gind = gind[w]
                      corr = 1./(1.- tmp[gind].r[clr])
                      w=where(tmp[gind].momerr[clr]*corr LT 0.64,ngal)
                      corr=0
                      IF ngal NE 0 THEN BEGIN 
                          gind=gind[w]
                          
                          srctmp1 = replicate(struct, ngal)
                          srctmp1.rho4 = tmp[gind].rho4
                          
                          corr = 1./(1.- tmp[gind].r[clr])
                          srctmp1.momerr = tmp[gind].momerr[clr]*corr
                          srctmp1.e1 = tmp[gind].e1[clr]*corr
                          srctmp1.petrocounts = tmp[gind].petrocounts
                          numlist1[fi] = ngal
                          ntotal1 = ntotal1+ngal
                          ptrlist1[fi] = ptr_new(srctmp1)

                          ;;;;;;; LOOK FOR STUFF THAT ARE PHOTO GAL
                          w2=where( (tmp[gind].type[1] EQ galtype AND $
                                     tmp[gind].type[2] EQ galtype) OR $
                                    (tmp[gind].type[2] EQ galtype AND $
                                     tmp[gind].type[3] EQ galtype), ngal2)
                          
                          IF (ngal2 NE 0) AND (ngal2 NE ngal) THEN BEGIN 
                              gout = gind
                              remove, w2, gout
                              nout = n_elements(gout)
                              
                              srctmp2 = replicate(struct, nout)
                              corr2 = 1./(1.- tmp[gout].r[clr])
                              srctmp2.rho4 = tmp[gout].rho4
                              srctmp2.momerr = tmp[gout].momerr[clr]*corr2
                              srctmp2.e1 = tmp[gout].e1[clr]*corr2
                              srctmp2.petrocounts = tmp[gout].petrocounts
                              srctmp2.objc_type = tmp[gout].objc_type
                              numlist2[fi] = nout
                              ntotal2 = ntotal2+nout
                              ptrlist2[fi] = ptr_new(srctmp2)
                          ENDIF 
                      ENDIF 
                  ENDIF
              ENDIF 
          ENDIF 
          IF fi MOD 20 EQ 0 THEN print,'Field: ',ntostr(field),'  ',$
            colors[clr]+$
            '-band Src Galaxies: '+ntostr(ngal)
      ENDFOR 

      print,'All: ',ntotal1
      print,'OUT: ',ntotal2

      sources1 = replicate(struct, ntotal1)
      beg=0
      FOR i=0,nfields-1 DO BEGIN 
          IF numlist1[i] NE 0 THEN BEGIN 
              sources1(beg:beg+numlist1(i)-1)=*ptrlist1(i)
          ENDIF 
          ptr_free,ptrlist1(i)   ;kill the heap variables associated 
                                ;with the pointer
          beg=beg+numlist1(i)

      ENDFOR 
      
      sources2 = replicate(struct, ntotal2)
      beg=0
      FOR i=0,nfields-1 DO BEGIN 
          IF numlist2[i] NE 0 THEN BEGIN 
              sources2(beg:beg+numlist2(i)-1)=*ptrlist2(i)
          ENDIF 
          ptr_free,ptrlist2(i)   ;kill the heap variables associated 
                                ;with the pointer
          beg=beg+numlist2(i)

      ENDFOR 


      hist=histogram(sources1.petrocounts[clr],bin=.2,$
                     min=18., max=22., reverse_indices=rev_ind)
      nhist = n_elements(hist) ;last will be garbage
      magarray = fltarr(nhist)
      rho4array = fltarr(nhist)
      rho4_sigarray = fltarr(nhist)
      rho4_errarray = fltarr(nhist)
      momerr = fltarr(nhist)
      momsig = fltarr(nhist)

      FOR i=0,nhist-1 DO BEGIN 
          w=rev_ind[ rev_ind[i]:rev_ind[i+1]-1 ]
          nw = n_elements(w)
          print,nw
          mag_m = mean(sources1[w].petrocounts[clr])
          rho4_m = moment(sources1[w].rho4[clr], maxmoment=2)
          
          magarray[i] = mag_m
          rho4array[i] = rho4_m[0]
          rho4_sigarray[i] = sqrt(rho4_m[1])
          rho4_errarray[i] = rho4_sigarray[i]/sqrt(nw)

          momerr[i] = mean( sources1[w].momerr )
          momsig[i] = sdev( sources1[w].e1 )
      ENDFOR 

      xmarg = !x.margin
      !x.margin = [10,10]   

      simpctable, rct, gct, bct

      psname = '/sdss3/usrdevel/esheldon/PLOTS/run'+$
        ntostr(run)+'_camcol'+ntostr(col)+'_rho4_mag_.ps'
      begplot, name=psname, /color, /land

      !p.background = !white
      !p.color = !black
      erase & multiplot, [1,2]

      xrange=[18,22]

      plot, [0], /nodata, yrange=[0,4], xrange=xrange, background=!white, color=!black, $
        ytitle='!7q!3!U4!N',charsize=1.3
      oplot, sources1.petrocounts[clr], sources1.rho4[clr], psym=3, color=!black
      oplot, magarray, rho4array, color=!magenta
      oplot, magarray, rho4array+rho4_sigarray, color=!magenta, line=2
      oplot, magarray, rho4array-rho4_sigarray, color=!magenta, line=2
      legend,['All: '+ntostr(ntotal1), 'Mean(All)'],/right,bcolor=!black, $
        psym=[3,0], $
        textcolor=[!black,!black], colors=[!black, !magenta]

      multiplot
      plot, [0], /nodata, yrange=[0,4], xrange=xrange, background=!white, color=!black, $
        ytitle='!7q!3!U4!N',xtitle="r'",charsize=1.3
      oplot, sources2.petrocounts[clr], sources2.rho4[clr], psym=3, color=!black
      oplot, magarray, rho4array, color=!magenta
      oplot, magarray, rho4array+rho4_sigarray, color=!magenta, line=2
      oplot, magarray, rho4array-rho4_sigarray, color=!magenta, line=2
      legend,['OUT: '+ntostr(ntotal2), 'Mean(All)'], bcolor=!black, $
        textcolor=[!black, !black], /right, colors=[!black, !magenta],psym=[3,0]

      multiplot, /reset

      plot, [0], /nodata, yrange=[0,4], xrange=xrange, background=!white, color=!black, $
        ytitle='!7q!3!U4!N',xtitle="r'",charsize=1.3
      oplot, sources1.petrocounts[clr], sources1.rho4[clr], psym=3, color=!black
      oplot, sources2.petrocounts[clr], sources2.rho4[clr], psym=3, color=!blue
      oplot, magarray, rho4array, color=!magenta
      oplot, magarray, rho4array+rho4_sigarray, color=!magenta, line=2
      oplot, magarray, rho4array-rho4_sigarray, color=!magenta, line=2
      legend,['All: '+ntostr(ntotal1), 'OUT: '+ntostr(ntotal2), 'Mean(All)'],$
        /right,bcolor=!black, psym=[3, 3, 0], $
        textcolor=[!black,!blue, !magenta], colors=[!black, !blue, !magenta]
      
      erase & multiplot, [1,2]

      plot, [0], /nodata, yrange=[0,4], xrange=xrange, background=!white, color=!black, $
        charsize=1.3, ytitle='!S!7q!3!U4!N!R!A-!N'
      oploterror, magarray, rho4array, rho4_errarray, psym=1, color=!black

      multiplot
      plot, [0], /nodata, background=!white, color=!black, $
        xrange=xrange, yrange=[0,.6], $
        xtitle="r'",charsize=1.3;, ytitle='!7r!3(!7q!3!U4!N)'
      oplot, magarray, rho4_sigarray, color=!black
      oplot, magarray, momerr, line=2, color=!black
      oplot, magarray, momsig, line=3, color=!black
      legend,['!7r!3(!7q!3!U4!N)','!7r!3(e!D1!N)','err(e!D1!N)'],line=[0,3,2], bcolor=!black, $
        textcolor=[!black, !black, !black]

      multiplot, /reset

      !x.margin = xmarg

      endplot, /noprint
      pslandfix, psname, xw=60

;      gifname = '/sdss3/usrdevel/esheldon/PLOTS/run'+ntostr(run)+'_camcol'+ntostr(col)+'_rho4_mag_.gif'
;      write_gif, gifname, tvrd(), rct, gct, bct

      print

      jump:
  ENDFOR 

  ptime,systime(1)-time

  return
END 
