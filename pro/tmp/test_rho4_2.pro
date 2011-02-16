PRO test_rho4_2, run, rerun, camcol, iorder, correct1, correct2, objsize, petrosize, start=start, nframes=nframes

  IF n_params() eq 0 then begin
      print,'-syntax test_rho4_2, run, rerun, camcol, iorder, correct1, correct2, objsize, petrosize, start=start, nframes=nframes'
      return
  ENDIF

  time=systime(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  newfront = 'adat'             ;The front of files with adaptive moments
  newfront='adatc'
  outfront = 'adatc'

  colors=['u','g','r','i','z']
;  bands = [1,2,3]               ; we will process g,r,i
  bands = [2]
  nband = n_elements(bands)
  
  rstr = ntostr(run)
  rrstr = ntostr(rerun)
  cstr = ntostr(camcol)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Setup the column
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir,$
    corratldir = fitdir
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes, $
                   fieldmin=fieldmin, fieldmax=fieldmax
  ntot = fieldmax-fieldmin+1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up new filenames. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rename_tsobj, files, corrdir, newfront, adatfiles, nchar
  rename_tsobj, files, corrdir, outfront, outfiles, outnchar

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; setup moment arrays for each color
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Reading in moment arrays'
  FOR iband=0, nband-1 DO BEGIN
      cindex = bands[iband]
      momname, run, camcol, cindex, iorder, psffile
      psffile=fitdir+psffile
      readcol,psffile, sjunk, junk,a,b,c,d,e,f,format='A,I,D,D,D,D,D,D,D'

      IF iband EQ 0 THEN BEGIN 
          aa = fltarr(nband, n_elements(a))
          bb = aa
          cc = aa
          dd = aa
          ee = aa
          ff = aa
      ENDIF 

      aa[iband,*]=a[*]
      bb[iband,*]=b[*]
      cc[iband,*]=c[*]
      dd[iband,*]=d[*]
      ee[iband,*]=e[*]
      ff[iband,*]=f[*]
  ENDFOR 

  print,'-----------------------------------------------------'
  print,' Run: ',rstr,' Rerun: ',rrstr,' Camcol: ',cstr
  print,'-----------------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields and find corrected e1, e2 and R
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Outputs files are '+outfront+'-*'
  print

  numlist=lonarr(nframes)
  ptrlist1=ptrarr(nframes)
  ptrlist2=ptrarr(nframes)
  ptrlist3=ptrarr(nframes)
  ptrlist4=ptrarr(nframes)
  ntotal=0L

  FOR ic = 0L, nframes-1 DO BEGIN 
      infile = adatfiles[ic]
      outf   = outfiles[ic]

;      IF outf EQ infile THEN BEGIN
;          print,'Error: infile = outfile'
;          print,infile,outf
;          return
;      ENDIF 

      field = fnums[ic]
      fstr = ntostr(field)
      
      ind  = ic*4L
      ind1 = ind+1
      ind2 = ind+2
      ind3 = ind+3

      ;; Fields have already been trimmed of overlap.  Don't need to
      ;; use read_tsobjmin
      openr, lun, infile, /get_lun
      IF ic EQ 0 THEN BEGIN 
          pstruct = mrdfits3(lun, 1, 0, /silent)
      ENDIF ELSE BEGIN
          pstruct = mrdfits3(lun, 1, 0, /silent, /deja_vu)
      ENDELSE 
      free_lun, lun

                                ; Initialize arrays for this field

      nps = n_elements(pstruct)
      good = lindgen(nps)       ; indices of good measurements
      ngood = nps
      x = fltarr(nps) & y = x
      psfsize = x & psfrho4 = x & objsize = x 
      cor1 = x & cor2 = x
      psfe1 = x & psfe2 = x & obje1 = x & obje2 = x

      print,'Run: ',rstr,' Camcol: ',cstr,' Field: ',fstr,' Objects: ',ntostr(ngood)

      FOR iband = 0, nband-1 DO BEGIN 
         
          cindex = bands[iband]

          x[*] = pstruct.colc[cindex]
          y[*] = pstruct.rowc[cindex]
      
          ;; find psf size in pixels at the position of each object
          psfsize[*] = aa[iband,ind] + bb[iband,ind]*x   $
                                     + cc[iband,ind]*y   $
                                     + dd[iband,ind]*x*y $
                                     + ee[iband,ind]*x*x $
                                     + ff[iband,ind]*y*y

          ;; higher order moment of psf
          psfrho4[*] = aa[iband,ind3] + bb[iband,ind3]*x   $
                                      + cc[iband,ind3]*y   $
                                      + dd[iband,ind3]*x*y $
                                      + ee[iband,ind3]*x*x $
                                      + ff[iband,ind3]*y*y
      
          ;; size of objects in pixels
          objsize[*] = pstruct.ixx[cindex] + pstruct.iyy[cindex]

          ;; Bad measurements
          ww = where(objsize eq 0., nww)
          vv = where(psfsize eq 0., nvv)
uu = where(pstruct.petror90err[2] GT 0.,nuu)
          IF (nww NE 0) OR (nvv NE 0) OR (uu NE 0)THEN BEGIN
              ww=[ww,vv,uu]
              wminone=where(ww NE -1)
              ww=ww[wminone]
              ww=ww[rem_dup(ww)]

              remove, ww, good
              
                                ;Set bad to defaults
              pstruct[ww].e1[cindex] = 1.e10
              pstruct[ww].e2[cindex] = 1.e10
              pstruct[ww].R[cindex] = -1.
              ngood = n_elements(good)

          ENDIF 

          wtmp=where(pstruct[good].petrocounts[2]-pstruct[good].reddening[2] LE 22.0 AND $
                     pstruct[good].momerr[2] LT 0.64)

          good = good[wtmp]
          cor1[good] = psfsize[good]/objsize[good]* $
               (4./psfrho4[good]-1.)/(4./pstruct[good].rho4[cindex]-1.)
          cor2[good] = (4./psfrho4[good]-1.)/(4./pstruct[good].rho4[cindex]-1.)

          corr1=cor1[good]
          corr2=cor2[good]
          petrosize=pstruct[good].petror90[2]
          osize=objsize[good]

          numlist[ic] = n_elements(good)
          ptrlist1[ic] = ptr_new(corr1)
          ptrlist2[ic] = ptr_new(corr2)
          ptrlist3[ic] = ptr_new(osize)
          ptrlist4[ic] = ptr_new(petrosize)
          petrosize=0
          osize=0
          corr1=0
          corr2=0
;          psfe1[good] = (aa[iband,ind1] + bb[iband,ind1]*x[good]   $
;                    + cc[iband,ind1]*y[good]   $
;                    + dd[iband,ind1]*x[good]*y[good] $
;                    + ee[iband,ind1]*x[good]*x[good] $
;                    + ff[iband,ind1]*y[good]*y[good] )/psfsize[good]
;
;          psfe2[good] = 2.*(aa[iband,ind2] + bb[iband,ind2]*x[good]   $
;                    + cc[iband,ind2]*y[good]   $
;                    + dd[iband,ind2]*x[good]*y[good] $
;                    + ee[iband,ind2]*x[good]*x[good] $
;                    + ff[iband,ind2]*y[good]*y[good] )/psfsize[good]
;      
;          obje1[good] = ( pstruct[good].ixx[cindex] $
;                            -pstruct[good].iyy[cindex] )/objsize[good]
;          obje2[good] = 2.*pstruct[good].ixy[cindex]/objsize[good]
;      
;          pstruct[good].e1[cindex] = obje1[good] - $
;                                        psfe1[good]*cor1[good] 
;          pstruct[good].e2[cindex] = obje2[good] - $
;                                        psfe2[good]*cor1[good]
;          pstruct[good].R[cindex]  = cor1[good]
          
;          print,'Band: ',colors[cindex],' Ngood: ',ntostr(ngood)

          IF iband NE nband-1 THEN BEGIN
              good = 0 & good = lindgen(nps)
              ngood = nps
          ENDIF           

      ENDFOR                    ;End loop over bandpasses

;      mwrfits, pstruct, outf, /create

                                ;Free up the memory
      pstruct = 0
      good = 0
      x = 0 & y = 0
      psfsize = 0 & psfrho4 = 0
      objsize = 0
      cor1 = 0
      cor2 = 0
      corr1 = 0
      corr2 = 0
      psfe1 = 0 & psfe2 = 0
      obje1 = 0 & obje2 = 0
      
  ENDFOR                        ;End loop over fields
   
  correct1 = fltarr(total(numlist))
  correct2 = correct1
  petrosize = correct1
  objsize = correct1

  beg=0L
  FOR i=0L,nframes-1 DO BEGIN 
      correct1[beg:beg+numlist[i]-1]=*ptrlist1[i]
      correct2[beg:beg+numlist[i]-1]=*ptrlist2[i]
      objsize[beg:beg+numlist[i]-1]=*ptrlist3[i]
      petrosize[beg:beg+numlist[i]-1]=*ptrlist4[i]
      ptr_free,ptrlist1[i]       
      ptr_free,ptrlist2[i]      
      ptr_free,ptrlist3[i]
      ptr_free,ptrlist4[i]
                                ;kill the heap variables associated 
                                ;with the pointer
      beg=beg+numlist[i]
  ENDFOR 

  ptime,systime(1)-time

  return
END 

