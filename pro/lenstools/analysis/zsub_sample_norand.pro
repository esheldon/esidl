;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    SUB_SAMPLE
;       
; PURPOSE:
;    Take the outputs of zobjshear and zrandshear and find
;    subsamples based on an input boolean expression. Then
;    make all 4 outputs: the lensum, the sum, the shear, 
;    and the redshift files.
;
; CALLING SEQUENCE:
;    sub_sample, lensumfile, randsumfile, wstr, outdir=outdir, 
;            addstr=addstr, indices=indices, uselens=uselens, 
;            userand=userand, addweight=addweight
;
; INPUTS: 
;    lensumfile: the lensum file output by zobjshear
;    randsumfile: the lensum file output by zrandshear
;    wstr: the boolean expression in string form. It should use the
;          structure lensum. (see example)  You can leave off this
;          input if you use the indices input (see below)
;
; OPTIONAL INPUTS:
;    outdir: the output directory. The default is the same directory
;            as the input file.
;    addstr: a string too add to the front of the name. Default is ''
;            This will not overwrite the lensumfile. That file is output
;            with a N1.fit or N2.fit, etc., ending and the program increments
;            that number by 1 until the file is not found.
;    indices: instead of a boolean expression, the indices can be fed in, 
;             (see examples)
;    uselens: Send in the lensum structure, don't read from disk. This can
;             save time if you are doing many subsamples. NOTE the lensumfile
;             and randsumfile must still be entered so the output name can
;             be created.
;    userand: same as uselens
;    addweight: an additional weighting to add to the lensing outputs. Must
;               be same size as lensum.
;
; KEYWORD PARAMETERS:
;    NONE
;       
; OUTPUTS: 
;    the new lensum,sum,shear,and redshift files.
;
; OPTIONAL OUTPUTS:
;    NONE.
;
; EXAMPLES:
;    ;; EXAMPLE 1
;    dir='/sdss4/data1/esheldon/GAL_GAL/spectra/'
;    outdir='~/myfiles/'
;    lensumfile=dir+'zgal_gal_stripe82_g_lensum_N1.fit'
;    randsumfile=dir+'zrand_stripe82_g_lensum_N1.fit'
;    wstr = 'lensum.z1d lt 0.4'
;    addstr = 'subz'
;    sub_sample, lensumfile, randsumfile, wstr, addstr=addstr, outdir=outdir
;
;    ;; EXAMPLE2 Its hard to remember all the flag numbers, so
;    ;; do the choosing outside of SUB_SAMPLE but let it make
;    ;; all the new structures.
;
;    lensum=mrdfits(lensumfile, 1)
;    make_flag_struct,fs & fs.bright='N' & fs.satur='N' & .....etc
;    flag_select, lensum, fs, 2, si
;    addstr='flagcut'
;    sub_sample, lensumfile, randsumfile, $
;         indices=si, uselens=lensum, addstr=addstr, outdir=outdir
;
; CALLED ROUTINES:
;    MRDFITS
;    HEADFITS
;    SXPAR 
;    GETZTAG
;    MATCH2RAND
;    COMBINE_ZLENSUM
;    DIRSEP
;    REPSTR
;    EXIST
;    NEWNAME
;    ZSSHHDR
;    ZSUMHDR
;    MWRFITS
;
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    Creation ??-JUNE-2000  Erin Scott Sheldon UofMich
;                                       
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO zsub_sample_norand, lensumfile, wstr, outdir=outdir, addstr=addstr_in, indices=indices, uselens=uselens, addweight=addweight, overwrite=overwrite, addhdr=addhdr

  nparam = n_params()
  IF nparam LT 1 THEN BEGIN 
      print,'-Syntax: zsub_sample, lensumfile, wstr, outdir=outdir, addstr=addstr,indices=indices,uselens=uselens, addweight=addweight, overwrite=overwrite'
      return
  ENDIF 
  IF nparam LT 2 THEN BEGIN 
      IF n_elements(indices) EQ 0 THEN BEGIN
          print,'-Syntax:zsub_sample, lensumfile, wstr, outdir=outdir, addstr=addstr,indices=indices,uselens=uselens, addweight=addweight, overwrite=overwrite'
          return
      ENDIF 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; makes subsamples
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time=systime(1)
  
  IF (n_elements(addstr_in) EQ 0) THEN BEGIN 
      addstr = '' 
  ENDIF ELSE BEGIN
      IF addstr_in EQ '' THEN BEGIN
          addstr='' 
      ENDIF ELSE BEGIN 
          len = strlen(addstr_in)
          qq=strmid(addstr_in, len-1)
          IF qq EQ '_' THEN addstr = addstr_in ELSE addstr = addstr_in+'_'
      ENDELSE 
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get structs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF fexist(lensumfile) THEN BEGIN 

      print
      print,'----------------------------------------------'
      IF n_elements(uselens) EQ 0 THEN BEGIN 
          print,'Reading lens sum file: ',lensumfile
          lensum = mrdfits(lensumfile, 1, hdr, /silent)
      ENDIF ELSE BEGIN 
          print,'Using input lensum'
          lensum=uselens
          hdr = headfits(lensumfile, exten=1)
      ENDELSE 

      hval = sxpar(hdr, 'H')
      rminkpc = sxpar(hdr, 'RMINKPC')
      rmaxkpc = sxpar(hdr, 'RMAXKPC')
      binsize = sxpar(hdr, 'BINWIDTH')

  ENDIF ELSE BEGIN
      print
      print,'Files not found'
      return
  ENDELSE 

  wz = getztag(lensum)
  IF wz EQ -1 THEN message,'No redshift tag found!'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Add weighting if requested
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  norig = n_elements(lensum)
  nadd = n_elements(addweight)
  nbin = n_elements(lensum[0].rsum)

  IF nadd NE 0 THEN BEGIN 
      IF nadd NE norig THEN message,'addweight must be same size as lensum'
      print
      print,'Adding input weighting'
      print
      addweight2 = addweight^2
      FOR i=0L, nbin-1 DO BEGIN 

          lensum.etansum[i] = lensum.etansum[i]*addweight
          lensum.etanerrsum[i] = lensum.etanerrsum[i]*addweight2
          lensum.eradsum[i] = lensum.eradsum[i]*addweight
          lensum.eraderrsum[i] = lensum.eraderrsum[i]*addweight2

          lensum.tansigsum[i] = lensum.tansigsum[i]*addweight
          lensum.tansigerrsum[i] = lensum.tansigerrsum[i]*addweight2
          lensum.radsigsum[i] = lensum.radsigsum[i]*addweight
          lensum.radsigerrsum[i] = lensum.radsigerrsum[i]*addweight2

          lensum.wsum[i] = lensum.wsum[i]*addweight

      ENDFOR 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; where command selects subsamples
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(indices) EQ 0 THEN BEGIN 
      command = 'w = where('+wstr+',nlens)'
      IF NOT execute(command) THEN return
  ENDIF ELSE BEGIN 
      nlens = n_elements(indices)
      w=indices
  ENDELSE 

  print
  IF nlens EQ 0 THEN BEGIN 
      print,'No lenses passed the cut!'
      return
  ENDIF ELSE BEGIN 
      print,ntostr(nlens),' lenses passed cut'
  ENDELSE 

  lensum = temporary(lensum[w])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make shear and sum structs for
  ;; this subsample 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Creating sum files and shear files for subsample'

  combine_zlensum, lensum, binsize, rminkpc, rmaxkpc, hval, shstruct

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make redshift structures
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zs = create_struct('z', 0., $
                     'scritinv', 0.)


  zs=create_struct(zs, 'lambda',0d,'eta',0d)

  zstruct = replicate(zs, nlens)
  IF tag_exist(lensum,'lambda') AND tag_exist(lensum,'eta') THEN BEGIN 
      zstruct.lambda = lensum.lambda
      zstruct.eta = lensum.eta
  ENDIF ELSE BEGIN 
      eq2survey, lensum.ra,lensum.dec,lam,eta
      zstruct.lambda = lam
      zstruct.eta = eta
  ENDELSE 

  zstruct.z = lensum.(wz[0])
  zstruct.scritinv = lensum.scritinv
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make output file names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dirsep, lensumfile, indir, tlf

  IF n_elements(outdir) EQ 0 THEN outdir = indir

  lf = outdir + addstr + tlf
  lshf  = outdir + addstr+ repstr(tlf, 'lensum_', '')
  zf    = outdir + addstr+ repstr(tlf, 'lensum', 'z')

  IF NOT keyword_set(overwrite) THEN BEGIN 
      WHILE fexist(lf) OR fexist(rf) DO BEGIN 
          lf = newname(lf)
          lshf = newname(lshf)
          zf = newname(zf)
      ENDWHILE 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make output headers
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; shear header for lenses
  shhdr = zshhdr(shstruct)

  ;; Sum file for lenses
  sumhdr = zsumhdr(shstruct)
  lensumhdr = sumhdr


  print
  print,'Outputting files:'
  print

  print,lshf
  mwrfits2, shstruct, lshf, shhdr, /create, /destroy
  print,zf
  mwrfits2, zstruct, zf, /create, /destroy
  print,lf
  mwrfits2, lensum, lf, lensumhdr, /create, /destroy

  print
  print,'----------------------------------------------'
  ptime,systime(1)-time


  return
END 
