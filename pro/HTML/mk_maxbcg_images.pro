PRO mk_maxbcg_images, bcg, big=big, small=small

  mb = obj_new('maxbcg','dr406')
  IF n_elements(bcg) EQ 0 THEN BEGIN 
      bcg = mb->get()
  ENDIF 

  ;; images of the top 100
  ndo = 1000
  ndo_fchart = 1000


  ;; make image with size of 1.5*search radius
  circle_size = 1.0             ; Mpc
  circle_size_med = 0.5
  circle_size_small = 0.1       ; 100 kpc
  box_size = 1.0*1.15            ; half size of side

  ;;IF NOT tag_exist(bcg,'tngals', index=ngals_tag) THEN message,'tngals does not exist'
  IF NOT tag_exist(bcg,'ngals200', index=ngals_tag) THEN message,'ngals200 does not exist'

  wpaper = where(bcg.paper)


  s = sort(bcg[wpaper].(ngals_tag))
  s = wpaper[s]
  w = reverse(s)
  outDir = '~/plots/maxbcg/fchart/'
  linkDir = 'plots/MaxBCG/fchart/'

  IF keyword_set(big) THEN BEGIN 
      dataFile = outDir + 'image_info_big.dat'
      ext = '_big.jpg'
  ENDIF ELSE BEGIN 
      dataFile = outDir + 'image_info.dat'
      ext = '.jpg'
  ENDELSE 

  openw, lun, dataFile,/get_lun

  markstruct = {ra:0d, dec:0d, type:'circle', radius:0.0, linestyle:0}
  markstruct = replicate(markstruct, 3)

  FOR ii=0L, ndo-1 DO BEGIN 

      i=w[ii]

      print
      print,'-----------------------------------------------'
      print,'ii = '+ntostr(ii+1)+'/'+ntostr(ndo)
      print,'Redshift = '+ntostr(bcg[i].z)
      print,'Ngals    = '+ntostr(bcg[i].(ngals_tag))
      print,'-----------------------------------------------'

      rad = box_size/angdist_lambda(bcg[i].z) ; radians
      rad = rad*180d/!pi*60.    ; arcminutes
      circrad = circle_size/angdist_lambda(bcg[i].z) ; radians
      circrad = circrad*180d/!pi*60.0 ; arcminutes
      circrad_med   = circrad*circle_size_med/circle_size
      circrad_small = circrad*circle_size_small/circle_size

      ;; all centered on same spot
      markstruct[*].ra = bcg[i].ra
      markstruct[*].dec = bcg[i].dec

      markstruct[0].radius = circrad
      markstruct[0].linestyle = 0
      markstruct[1].radius = circrad_med
      markstruct[1].linestyle = 2
      markstruct[2].radius = circrad_small
      markstruct[2].linestyle = 1


      idstr = strn(bcg[i].id,padchar='0',len=6)
      ngalstr = strn(bcg[i].(ngals_tag),padchar='0',len=3)
      jpegFchartName = 'maxBCG_fchart_ngals'+ngalstr+'_id'+idstr+ext
      jpegFchart = outDir + jpegFchartName

      bcg_rabs = calc_sdss_lumsolar(bcg[i].bcg_rlum, clr=2, /inverse)
      rabs200 = calc_sdss_lumsolar(bcg[i].rlum200, clr=2, /inverse)
      bcg_iabs = calc_sdss_lumsolar(bcg[i].bcg_ilum, clr=3, /inverse)
      iabs200 = calc_sdss_lumsolar(bcg[i].ilum200, clr=3, /inverse)

      bcg_iabsstr = ntostr(bcg_iabs, 5, /round)
      addTitle = 'id = '+idstr+'  ngals = '+strn(bcg[i].(ngals_tag),len=3)+'  z = '+ntostr(bcg[i].z, 4, /round);+$
;        ' M!S!Di!N!RBCG!N = '+bcg_iabsstr

      ;; Don't regenerate
      IF NOT fexist(jpegFchart) AND ii NE 71 AND (ii+1) LE ndo_fchart THEN BEGIN 
          rgbfchart, bcg[i].ra, bcg[i].dec, $
            radius = rad, $
            extra_markstruct=markstruct, $
            addTitle=addTitle, $
            jpegFchart=jpegFchart ;, /nodisplay, z_resolution=zres
      ENDIF 

      ;; Factors of two account for halving the resolution
      radpix = rad*60.0/0.4
      scale = 0.4

      radpix = radpix/2.0
      scale = scale*2.0
      width = 2*radpix
      height = 2*radpix

      ;; Sky server limits to 2048x2048
      IF width GT 2048 THEN BEGIN 
          fac = width/2048.0
          scale = scale*fac

          width  = 2048
          height = 2048
      ENDIF 

      URL = sdss_fchart_url(bcg[i].ra, bcg[i].dec, $
                            scale = scale, $
                            width=width, height=height)
      naviURL = sdss_fchart_url(bcg[i].ra, bcg[i].dec, $
                                scale = scale, $
                                width=width, height=height, /navigator, $
                                opt='GS')


;      printf, lun, $
;        bcg[i].id, $
;        bcg[i].(ngals_tag)


      printf, lun, $
        jpegFchartName+'   ',$
        bcg[i].id, $
        bcg[i].(ngals_tag), $
        bcg_rabs, $
        rabs200, $
        bcg_iabs, $
        iabs200, $
        bcg[i].z, $
        bcg[i].bcg_spec_z,$
        bcg[i].ra, $
        bcg[i].dec, $
        '  '+URL, $
        '  '+naviURL, $
        format='(a,I,I,g,g,g,g,g,g,g,g,a,a)'

      flush, lun
  ENDFOR 

  free_lun, lun

  obj_destroy,mb

END 
