PRO maketotwtheta

  pold=!p.multi
  !p.multi=[0,1,3]
  charsize = 1.5
  yrange=[-.02, .12]

  indir = '/sdss4/data1/esheldon/GAL_GAL/DENSITY/'
  run1str='752'
  run2str='756'
  front1 = 'gal_gal_'
  front2  = 'for'+front1
  clr=[1,2,3]
  nclr = n_elements(clr)
  colors=['u','g','r','i','z']
  rmin = 10./60.                ;arcmin

  ;;get shear profiles and foreground overdensities

  psfile = indir+'wtheta_plots.ps'
  begplot, name=psfile
  xtitle='Projected Radius (arcsec)'

  FOR i=0, nclr-1 DO BEGIN

      shfile_dens = indir+front1+run1str+'_'+run2str+'_'+$
        colors[clr[i]]+'_N1.dat'
      shfile_low  = indir+front1+run1str+'_'+run2str+'_'+$
        colors[clr[i]]+'_N2.dat'

      forfile_dens = indir+front2+run1str+'_'+run2str+'_'+$
        colors[clr[i]]+'_N1.dat'
      randfile_dens = indir+front2+'rand_'+run1str+'_'+run2str+'_'+$
        colors[clr[i]]+'_N1.dat'
      forfile_low  = indir+front2+run1str+'_'+run2str+'_'+$
        colors[clr[i]]+'_N2.dat'
      randfile_low  = indir+front2+'rand_'+run1str+'_'+run2str+'_'+$
        colors[clr[i]]+'_N2.dat'

      rdobjshear, shfile_dens, shstr_dens, shnlens_dens, /silent
      rdobjshear, shfile_low, shstr_low, shnlens_low, /silent

      rdobjshear, forfile_dens, forstr_dens, fornlens_dens, tm, binsize,/silent
      rdobjshear, randfile_dens, randstr_dens, randnlens_dens, /silent
      rdobjshear, forfile_low, forstr_low, fornlens_low, /silent
      rdobjshear, randfile_low, randstr_low, randnlens_low, /silent

      ;;; Calculate area in each bin
      nbin = n_elements(forstr_dens)
      area = fltarr(nbin)
      bsize = binsize/60.   ;arcminutes
      FOR ii=0, nbin-1 DO BEGIN
          r1 = rmin + ii*bsize
          r2 = r1+bsize
          area[ii] = !pi*(r2^2 - r1^2)
      ENDFOR 
      
      nshear = n_elements(shstr_dens)
      dens_shear = fltarr(nbin) & dens_shearerr=fltarr(nbin)
      low_shear = dens_shear & low_shearerr=dens_shear
      dens_shear[0:nshear-1] = shstr_dens.shear
      dens_shearerr[0:nshear-1] = shstr_dens.shearerr
      low_shear[0:nshear-1] = shstr_low.shear
      low_shearerr[0:nshear-1] = shstr_low.shearerr


      message='        meanr           area         shear          shearerr      real gal/lens     error      rand gal/lens      error      overdens/area     error'
      
      ;;;;;;;;;;;;; output for high density regions
      print,colors[clr[i]]+' output for high density regions'
      forperlens_dens = forstr_dens.npair/fornlens_dens
      forperlenserr_dens = sqrt(forstr_dens.npair)/fornlens_dens

      randperlens_dens = randstr_dens.npair/randnlens_dens
      randperlenserr_dens = sqrt(randstr_dens.npair)/randnlens_dens

      over_dens = (forperlens_dens-randperlens_dens)/area
      overerr_dens = sqrt(forperlenserr_dens^2 + randperlenserr_dens^2)/area

      outfile = indir+'wtheta_dense_'+colors[clr[i]]+'.dat'
      openw, lun, outfile, /get_lun
      !textunit = lun

      printf,lun,message
      format='10F15.5'
      forprint, forstr_dens.meanr, area, dens_shear, $
                dens_shearerr, $
                forperlens_dens, forperlenserr_dens, $
                randperlens_dens, randperlenserr_dens, $
                over_dens, overerr_dens, $
                format=format,TEXT=5, /silent

      close, lun
      free_lun, lun

      ;; ;;;;;;;;;;;;;;;;output for less dense regions
      print,colors[clr[i]]+' output for low density regions'
      forperlens_low = forstr_low.npair/fornlens_low
      forperlenserr_low = sqrt(forstr_low.npair)/fornlens_low

      randperlens_low = randstr_low.npair/randnlens_low
      randperlenserr_low = sqrt(randstr_low.npair)/randnlens_low

      over_low  = (forperlens_low-randperlens_low)/area
      overerr_low = sqrt(forperlenserr_low^2 + randperlenserr_low^2)/area

      outfile = indir+'wtheta_low_'+colors[clr[i]]+'.dat'
      openw, lun, outfile, /get_lun
      !textunit = lun

      printf,lun,message
      forprint, forstr_low.meanr, area, $
                low_shear, low_shearerr, $
                forperlens_low, forperlenserr_low, $
                randperlens_low, randperlenserr_low, $
                over_low, overerr_low, $
                format=format,TEXT=5, /silent

      close, lun
      free_lun, lun

      n_dens = n_elements(forstr_dens.meanr)
      
      title='Excess of FG gals    '+colors[clr[i]]+' band'
      ploterr,forstr_dens.meanr,over_dens, overerr_dens, linestyle=0, $
        charsize=charsize, yrange=yrange, ystyle=1, title=title,xtitle=xtitle
      oploterr,forstr_low.meanr,over_low , overerr_low, linestyle=2, $
        charsize=charsize, yrange=yrange, ystyle=1
      oplot,[0,10000],[0,0]
      legend,['More Dense','Less Dense'],linestyle=[0,2],position=[800,.11]


  ENDFOR 
  ep
  !p.multi=pold
  
return
END 
