PRO combine_classcat, run1, rerun1, run2, rerun2, clr, zcat=zcat

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: combine_classcat, run1, rerun1, run2, rerun2, clr, zcat=zcat'
      return
  ENDIF 

  time = systime(1)

  r1str = ntostr(run1)
  rr1str = ntostr(rerun1)
  r2str = ntostr(run2)
  rr2str = ntostr(rerun2)

  sdssidl_setup
  
;  indir1 = !sdss_shapecorr_dir+'corr'+r1str+'/'+rr1str+'/combined/'
;  indir2 = !sdss_shapecorr_dir+'corr'+r2str+'/'+rr2str+'/combined/'

  indir1 = '/sdss3/data4/jracusin/'
  indir2 = indir1

  outdir = indir1

  uncert_cut = 0.64

  nnm = 'classgal'
  IF n_elements(zcat) NE 0 THEN nnm = 'z'+nnm
  typ = 'blah3'+ntostr(long(systime(1)))
  colors=['u','g','r','i','z']
  nclr = n_elements(clr)

  r1=0
  r2=0
  struct=0
  name = indir1 + 'run'+ r1str +'_'+ nnm +'_'+colors[clr]+'_done.fit'
  print
  print,'---------------------------------------------------------'
  print,'Reading ',name
  r1 = mrdfits(name, 1, hdr1, structyp = typ, /silent)
  fra1 = sxpar(hdr1, 'FIRSTRA')
  lra1 = sxpar(hdr1, 'LASTRA')

  name = indir2 + 'run'+ r2str +'_'+ nnm +'_'+colors[clr]+'_done.fit'
  print,'Reading ',name
  r2 = mrdfits(name, 1, hdr2, structyp = typ, /silent)
  fra2 = sxpar(hdr2, 'FIRSTRA')
  lra2 = sxpar(hdr2, 'LASTRA')

  keep1=lindgen( n_elements(r1) )
  keep2=lindgen( n_elements(r2) )
  ;; cut on uncertainty
      
;  w1=where(r1[keep1].uncert LT uncert_cut, n1)
;  w2=where(r2[keep2].uncert LT uncert_cut, n2)

  
         
;  keep1 = keep1[w1]
;  keep2 = keep2[w2]
  

  ;; Check if runs go across ra=0
  
  IF (fra1 EQ 0) AND (fra2 EQ 0) THEN BEGIN 
      fra1 = min(r1.ra) & lra1 = max(r1.ra)
      fra2 = min(r2.ra) & lra2 = max(r2.ra)
  ENDIF 
  IF fra1 GT lra1 THEN flag1=1 ELSE flag1=0
  IF fra2 GT lra2 THEN flag2=1 ELSE flag2=0

  IF flag1 AND flag2 THEN BEGIN ;Both runs go across ra=0.
      
      firstra = max([fra1, fra2])
      lastra  = min([lra1, lra2])

      w1 = where(r1[keep1].ra GE firstra OR r1[keep1].ra LE lastra, n1)
      w2 = where(r2[keep2].ra GE firstra OR r2[keep2].ra LE lastra, n2)
      
      bigflag = 1
  ENDIF ELSE BEGIN
      IF flag1 AND (NOT flag2) THEN BEGIN 
          firstra = fra2        ;Only r1 goes across ra=0
          lastra = min([lra1, lra2])
      ENDIF ELSE IF (NOT flag1) AND flag2 THEN BEGIN 
          firstra = fra1        ;Only r2 goes across ra=0
          lastra = min([lra1, lra2])
      ENDIF ELSE BEGIN          ;Neither goes across ra=0
          firstra = max([fra1, fra2])
          lastra  = min([lra1, lra2])
      ENDELSE 
      
      w1 = where(r1[keep1].ra GE firstra AND r1[keep1].ra LE lastra, n1)
      w2 = where(r2[keep2].ra GE firstra AND r2[keep2].ra LE lastra, n2)
      bigflag = 0
  ENDELSE 
  
  keep1=keep1[w1]
  keep2=keep2[w2]

  ;; set up new struct and free some memory
  ntot = n1+n2
  struct = replicate(r1[0], ntot)
  IF bigflag EQ 0 THEN BEGIN 

      struct[0:n1-1] = r1[keep1]
      r1=0

      struct[n1:ntot-1] = r2[keep2]
      r2=0

      print,'Sorting by objc_rowc'
      s_ind = sort( struct.ra )

  ENDIF ELSE BEGIN 
      w11 = where( r1[keep1].ra GE firstra, n11 )
      w21 = where( r2[keep2].ra GE firstra, n21 )
      ntot1 = n11+n21

      w12 = where( r1[keep1].ra LE lastra,  n12 )
      w22 = where( r2[keep2].ra LE lastra,  n22 )
      ntot2 = n12+n22
          
      IF ntot NE (ntot1 + ntot2) THEN BEGIN
          print,'What'
          return
      ENDIF 

      stemp1 = replicate(r1[0], ntot1)
      stemp2 = replicate(r1[0], ntot2)

      stemp1[0:n11-1] = r1[ keep1[w11] ]
      stemp1[n11:n11+n21-1] = r2[ keep2[w21] ]
      s_ind1 = sort( stemp1.ra )

      stemp2[0:n12-1] = r1[ keep1[w12] ]
      stemp2[n12:n12+n22-1] = r2[ keep2[w22] ]
      s_ind2 = sort( stemp2.ra )
          
      struct[0:ntot1-1] = stemp1[s_ind1]
      stemp1 = 0
      struct[ntot1:ntot-1] = stemp2[s_ind2]
      stemp2 = 0
      s_ind = lindgen(ntot)
          
      s_ind1 = 0
      s_ind2 = 0
  ENDELSE 
  keep1=0
  keep2=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get rid of duplicates.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Throwing out duplicates'
  rmclose_radec, struct[s_ind], keep
  s_ind = s_ind[keep]
  ntot = n_elements(keep)

  print,'Number in rectangular overlap region: ',ntostr(ntot)
      
  sname=outdir+'run'+r1str+'_'+r2str+'_'+nnm+'_'+colors[clr]+'_overlap.fit'
  outhdr = ['END']
  sxaddpar, outhdr, 'FIRSTRA', firstra
  sxaddpar, outhdr, 'LASTRA', lastra
  print,'Writing combined file: ',sname
  mwrfits, struct[s_ind], sname, outhdr, /create

  s_ind=0
  keep=0


  ptime,systime(1)-time

  return
END 
