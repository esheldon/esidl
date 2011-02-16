PRO mk_andreas_page, struct

  dir = '/net/cheops2/data1/lensinputs/AndreasClusters/'
  file = dir + 'Mr18.groups.dat'

  structdef = {groupid:0L, $
               ra:0.0, $
               dec:0.0, $
               z:0.0, $
               ngals: 0L, $
               Mr_tot: 0.0, $
               gmr_tot: 0.0, $
               sigma_v: 0.0, $
               r_perp_rms: 0.0, $
               r_edge: 0.0}

  IF n_elements(struct) EQ 0 THEN BEGIN 
      read_struct, file, structdef, struct, nlines=nclust
  ENDIF 

  outfile = '/net/cheops1/home/www/html/weaklens/AndreasClusters/cluster.dat'
  openw, lun, outfile, /get_lun

  rmin = 0.5                    ; minimum radius in Mpc

  nclust = n_elements(struct)
  FOR i=0L, nclust-1 DO BEGIN 

      ;; box of size 1.5*rms radius
      radius = struct[i].r_perp_rms*1.5 > 0.5 ; Mpc/h
;      radius = radius/0.7       ; h = 0.7

      rad = radius/angdist_lambda(struct[i].z)

      ;; now in pixels
      radpix = rad*180d/!dpi*60d*60d/0.4
      scale = 0.4

      ;; by default, half the resolution
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

      URL = sdss_fchart_url(struct[i].ra, struct[i].dec, $
                            scale = scale, $
                            width=width, height=height)
      naviURL = sdss_fchart_url(struct[i].ra, struct[i].dec, $
                                scale = scale, $
                                width=width, height=height,$
                                /navigator,opt='GS')


      printf, lun, $
        struct[i].groupid, $
        struct[i].ra, $
        struct[i].dec, $
        struct[i].z, $
        struct[i].ngals, $
        struct[i].mr_tot, $
        struct[i].gmr_tot, $
        struct[i].sigma_v, $
        struct[i].r_perp_rms, $
        struct[i].r_edge, $
        '  '+URL, $
        '  '+naviURL, $
        format = '(I,g,g,g,I,g,g,g,g,g,a,a)'
        
  ENDFOR 

  free_lun, lun

END 
