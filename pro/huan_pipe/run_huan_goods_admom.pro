PRO run_huan_goods_admom

  type = 'des5yr'

  dir = '~/Huan/goods/'

  secstr = ['22','23','33','34']

  seeing = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, $
            0.50, 0.55, 0.60, 0.65, $
            0.675, $
            0.70, 0.76, $
            0.78, $
            0.80, 0.85, $
            0.87, $
            0.90, 0.95, 1.00, 1.05, $
            1.10, 1.15, 1.20]

  nseeing=n_elements(seeing)
  
  FOR i=0L, nseeing-1 DO BEGIN 

      seestr = ntostr(seeing[i], 4, /round)

      catlist = dir + $
        'cat/'+type+'_'+seestr+'_h_goods_si_sect'+secstr+'_r1.0z_cat_mod.fit'
      imlist  = dir + $
        'images/'+type+'_'+seestr+'_h_si_sect'+secstr+'_v1.0_drz_img.fits'
      
      huan_goods_admom, imlist, catlist

  ENDFOR 

END 
