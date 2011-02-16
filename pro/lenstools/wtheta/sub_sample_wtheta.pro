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


PRO sub_sample_wtheta, lensumfile, randsumfile, wstr, outdir=outdir, addstr=addstr_in, indices=indices, uselens=uselens, userand=userand, addweight=addweight, overwrite=overwrite, norandlensum=norandlensum;, lumw=lumw

  nparam = n_params()
  IF nparam LT 2 THEN BEGIN 
      print,'-Syntax: sub_sample, lensumfile, randsumfile, wstr, outdir=outdir, addstr=addstr,indices=indices,uselens=uselens, userand=userand, addweight=addweight, lumw=lumw, overwrite=overwrite, norandlensum=norandlensum'
      return
  ENDIF 
  IF nparam LT 3 THEN BEGIN 
      IF n_elements(indices) EQ 0 THEN BEGIN
          print,'-Syntax: sub_sample, lensumfile, randsumfile, wstr, outdir=outdir, addstr=addstr,indices=indices,uselens=uselens, userand=userand, addweight=addweight, lumw=lumw, overwrite=overwrite, norandlensum=norandlensum'
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

  IF fexist(lensumfile) AND fexist(randsumfile) THEN BEGIN 

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
      IF n_elements(userand) EQ 0 THEN BEGIN 
          print,'Reading rand sum file: ',randsumfile
          randsum = mrdfits(randsumfile, 1, /silent)
      ENDIF ELSE BEGIN 
          print,'Using input randsum'
          randsum=userand
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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; where command selects subsamples
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(indices) EQ 0 THEN BEGIN 
      command = 'w = where('+wstr+',nlens)'
      tmp = execute(command)
      IF tmp NE 1 THEN return
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
  ;; match to random catalog (must get the
  ;; redshift distribution right)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Matching to random points'

  match2rand, lensum, randsum, wrand
  randsum = temporary(randsum[wrand])
  nrand = n_elements(wrand)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make shear and sum structs for
  ;; this subsample and randoms
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Creating sum files and shear files for subsample'

;  IF keyword_set(lumw) THEN BEGIN 
      trepstr = 'lumlensum'
      combine_wthetalumw_lensum, lensum,  binsize, rminkpc, rmaxkpc, hval, sumstruct, wstruct
      combine_wthetalumw_lensum, randsum,binsize,rminkpc,rmaxkpc,hval,rndsumstruct,rndwstruct
;  ENDIF ELSE BEGIN 
;      trepstr = 'lensum'
;      combine_wtheta_lensum, lensum,  binsize, rminkpc, rmaxkpc, hval, sumstruct, wstruct
;      combine_wtheta_lensum, randsum, binsize, rminkpc, rmaxkpc, hval, rndsumstruct, rndwstruct
;  ENDELSE 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make redshift structures
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zs = create_struct('z', 0.)
  zstruct = replicate(zs, nlens)
  zstruct.z = lensum.(wz[0])

  rzstruct = replicate(zs, nrand)
  rzstruct.z = randsum.z

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make output file names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dirsep, lensumfile, indir, tlf
  dirsep, randsumfile, ddd, trf
  IF n_elements(outdir) EQ 0 THEN outdir = indir

  lf = outdir + addstr + tlf
  rf = outdir + addstr + trf

  lsumf = outdir + addstr+ repstr(tlf, trepstr, 'sum')
  rsumf = outdir + addstr+ repstr(trf, trepstr, 'sum')
  
  lshf  = outdir + addstr+ repstr(tlf, trepstr+'_', '')
  rshf  = outdir + addstr+ repstr(trf, trepstr+'_', '')

  zf    = outdir + addstr+ repstr(tlf, trepstr, 'z')
  rzf   = outdir + addstr+ repstr(trf, trepstr, 'z')

  IF NOT keyword_set(overwrite) THEN BEGIN 
      WHILE exist(lf) OR exist(rf) DO BEGIN 
          lf = newname(lf)
          lsumf = newname(lsumf)
          lshf = newname(lshf)
          rf = newname(rf)
          rsumf = newname(rsumf)
          rshf = newname(rshf)
          zf = newname(zf)
          rzf = newname(rzf)
      ENDWHILE 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make output headers
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; header for lenses
  hdr = whdr(wstruct)
  sumhdr=hdr
  lensumhdr=hdr

  ;; header for random
  rhdr = whdr(rndwstruct)
  rsumhdr=rhdr
  rlensumhdr=rhdr

  print
  print,'Outputting files:'
  print

  print,lshf
  mwrfits2, wstruct, lshf, hdr, /create, /destroy
  print,lsumf
  mwrfits2, sumstruct, lsumf, sumhdr, /create, /destroy
  print,zf
  mwrfits2, zstruct, zf, /create, /destroy
  print,lf
  mwrfits2, lensum, lf, lensumhdr, /create, /destroy

  print
  print,rshf
  mwrfits2, rndwstruct, rshf, rhdr, /create, /destroy
  print,rsumf
  mwrfits2, rndsumstruct, rsumf, rsumhdr, /create, /destroy
  print,rzf
  mwrfits2, rzstruct, rzf, /create, /destroy
  IF NOT keyword_set(norandlensum) THEN BEGIN 
      print,rf
      mwrfits2, randsum, rf, rlensumhdr, /create, /destroy
  ENDIF 

  print
  print,'----------------------------------------------'
  ptime,systime(1)-time

  return
END 
