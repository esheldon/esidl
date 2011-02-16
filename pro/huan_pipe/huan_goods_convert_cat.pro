PRO huan_goods_convert_cat

  dir = "~/Huan/goods/"
  catfile = dir + 'cat/h_goods_si_r1.0z_cat_mod.txt'
  rmsfile = dir + 'images/'

  defval = -9999.
  defval2 = -9999d
  defval3 = -9999L
  arrval = replicate(defval, 3)
  arrval2 = replicate(defval, 11)
  ;; 'id_iau',         "J033155.16-274540.0", $
  ss = create_struct('ALPHA_J2000',    defval2, $
                     'DELTA_J2000',    defval2, $
                     'SECT_REFNUM',    defval3, $
                     'X_SECT',         defval, $
                     'Y_SECT',         defval, $
                     'X_MOSAIC',       defval, $
                     'Y_MOSAIC',       defval, $
                     'XPEAK_MOSAIC',   defval3, $
                     'YPEAK_MOSAIC',   defval3, $
                     'XPEAK_WORLD',    defval, $
                     'YPEAK_WORLD',    defval, $
                     'XMIN_MOSAIC',    defval3, $
                     'YMIN_MOSAIC',    defval3, $
                     'XMAX_MOSAIC',    defval3, $
                     'YMAX_MOSAIC',    defval3, $
                     'ISOAREA_IMAGE',  defval3, $ ;17
                     'THETA_IMAGE',    defval, $
                     'ELLIPTICITY',    defval, $
                     'ELONGATION',     defval, $
                     'ERRTHETA_IMAGE', defval, $
                     'KRON_RADIUS',    defval, $
                     'FLUX_RADIUS',    arrval, $
                     'FWHM_IMAGE',     defval, $
                     'CLASS_STAR',     defval, $
                     'FLAGS',          defval3, $ ;28
                     'IMAFLAGS_ISO',   defval3, $
                     'NIMAFLAGS_ISO',  defval3, $ ; 30
                     'BACKGROUND',     defval, $
                     'FLUX_MAX',       defval, $
                     'MAG_ISO',        defval, $
                     'MAGERR_ISO',     defval, $
                     'FLUX_ISO',       defval, $
                     'FLUXERR_ISO',    defval, $ ; 36
                     'MAG_ISOCOR',     defval, $
                     'MAGERR_ISOCOR',  defval, $
                     'FLUX_ISOCOR',    defval, $
                     'FLUXERR_ISOCOR', defval, $
                     'MAG_AUTO',       defval, $
                     'MAGERR_AUTO',    defval, $
                     'FLUX_AUTO',      defval, $
                     'FLUXERR_AUTO',   defval, $
                     'MAG_BEST',       defval, $
                     'MAGERR_BEST',    defval, $
                     'FLUX_BEST',      defval, $
                     'FLUXERR_BEST',   defval, $ ; 48
                     'MAG_APER',       arrval2, $ ; 49-59
                     'MAGERR_APER',    arrval2, $ ; 60-70
                     'FLUX_APER',      arrval2, $ ; 71-81
                     'FLUXERR_APER',   arrval2, $ ; 82-92
                     'X2_IMAGE',       defval, $
                     'Y2_IMAGE',       defval, $
                     'XY_IMAGE',       defval, $
                     'ERRX2_IMAGE',    defval, $
                     'ERRY2_IMAGE',    defval, $
                     'ERRXY_IMAGE',    defval, $
                     'A_IMAGE',        defval, $
                     'B_IMAGE',        defval, $
                     'ERRA_IMAGE',     defval, $
                     'ERRB_IMAGE',     defval, $
                     'ID_MOSAIC',      defval3)


  ncomments = 61
  nl = numlines(catfile)
  nobj = nl - 61

  struct = replicate(ss, nobj)

  openr, lun, catfile, /get_lun

  ln = ""
  FOR i=0L,ncomments-1 DO BEGIN 
      readf, lun, ln
      ;;print,ln
  ENDFOR 

  ;; Now read the struct
  print,'Reading catfile: ',catfile
  readf, lun, struct

  free_lun, lun

; For checking the fucked up sections
;  oplot,struct.x_sect-1,struct.y_sect-1,color=!yellow,psym=8,symsize=0.25
;  rmd = rem_dup(struct.sect_refnum)
;  usec = struct[rmd].sect_refnum
;  nn = n_elements(rmd)
;  FOR i=0L, nn-1 DO BEGIN 
;      sec = usec[i]
;      print,'Section ',sec
;      w=where(struct.sect_refnum EQ sec)

;      oplot,struct[w].x_sect-1,struct[w].y_sect-1,color=!red,psym=8,symsize=0.25
;      key=prompt_kbrd('hit a key')
;  ENDFOR 

;stop

  ;; Add sky rms to the file
  print
  print,'Adding tags'
  add_tags, temporary(struct), ['ivar', 'rms_sky'], ['0.0','0.0'], struct

  ;; Get the ones in our sections and read in the
  ;; weight map (1/var)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The catalog is fucked up.  Objects with sect_refnum = 32 -> 23 image
  ;;                                                       43 -> 34 image
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  image_sections = [22, 23, 33, 34]
  cat_sections   = [22, 32, 33, 43]
  nsec = n_elements(image_sections)

  print
  print,'Adding rms_sky'

  FOR i=0L, nsec-1 DO BEGIN 

      imsec  = image_sections[i]
      catsec = cat_sections[i]
      ivarfile = dir + 'images/h_si_sect'+ntostr(imsec)+'_v1.0_wht_img.fits'
      print
      print,'Reading ivar file: ',ivarfile
      ivarmap = mrdfits(ivarfile)

      w=where( struct.sect_refnum EQ catsec, nw)

      tstruct = struct[w]

      x = tstruct.x_sect-1
      y = tstruct.y_sect-1

      ivar = ivarmap[ x, y ]

      tstruct.ivar = ivar
      w2 = where(tstruct.ivar NE 0, nw2, comp=wcomp, ncomp=nwcomp)

      tstruct[w2].rms_sky = sqrt(1./ivar[w2])
      maxsky = max(tstruct[w2].rms_sky)

      ;; If no sky there, then set to max
      IF nwcomp NE 0 THEN tstruct[wcomp].rms_sky = maxsky*1000.

      outfile = dir + 'cat/h_goods_si_sect'+ntostr(imsec)+'_r1.0z_cat_mod.fit'
      print
      print,'Writing cat file for this section: ',outfile
      mwrfits2, tstruct, outfile, /create, /destroy

  ENDFOR 


END 
