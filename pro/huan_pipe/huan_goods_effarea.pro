PRO huan_goods_effarea, type

  IF n_elements(type) EQ 0 THEN typestr = '' ELSE typestr=type+'_'

;  maxwt = 350000.

  dir = '/net/cheops1/data0/esheldon/Huan/goods/images/'
  
  f22 = dir+typestr+'h_si_sect22_v1.0_wht_img.fits'
  print,'Reading ',f22
  wt22 = mrdfits(f22)
  maxwt = max(wt22)
  a22 = long( total(wt22)/maxwt )

  f23 = dir+typestr+'h_si_sect23_v1.0_wht_img.fits'
  wt23 = mrdfits(f23)
  maxwt = max(wt23)
  a23 = long( total(wt23)/maxwt )

  f33 = dir+typestr+'h_si_sect33_v1.0_wht_img.fits'
  wt33 = mrdfits(f33)
  maxwt = max(wt33)
  a33 = long( total(wt33)/maxwt )

  f34 = dir+typestr+'h_si_sect34_v1.0_wht_img.fits'
  wt34 = mrdfits(f34)
  maxwt = max(wt34)
  a34 = long( total(wt34)/maxwt )


  arcperpix = 0.03
  warea = a22+a23+a33+a34
  warea = warea*arcperpix^2/60.0^2

  area = (8192.*8192.)*arcperpix^2/60.0^2*4

  print,'Area = ',area
  print,'wArea = ',warea

END 
