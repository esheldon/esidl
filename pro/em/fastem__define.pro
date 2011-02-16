FUNCTION fastem::init
  self.allgauss = ptr_new(/alloc)
  self.reduced_allgauss = ptr_new(/alloc)
  return,1
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; General tools for dealing with fastem outputs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Look at the outputs from fastem
FUNCTION fastem::gauss_struct, ndim, extra=extra

  IF n_elements(ndim) EQ 0 THEN BEGIN 
      message,'-Syntax: gauss = em->gauss_struct(ndim, /extra)'
  ENDIF 

  struct = $
    { number: 0, $
      prob: 0.0d, $
      cen: dblarr(ndim), $
      cov: dblarr(ndim, ndim) }

  IF keyword_set(extra) THEN BEGIN 
      struct = $
        create_struct(struct, $
                      'detcov', 0d, $
                      'invcov', struct.cov, $
                      'invstatus', 0, $
                      'detinvcov', 0d )
  ENDIF 

  return,struct

END 


PRO fastem::read_gauss_header, lun, ngauss, ndim

  ;; Read the number of gaussians
  line = ''
  readf, lun, line
  tsp = strsplit(line, /extract)
  ngauss = long(tsp[1])

  ;; Read the number of dimensions
  readf, lun, line

  tsp = strsplit(line, /extract)
  ndim = long(tsp[1])


END 


FUNCTION fastem::read_gauss, gaussfile, ndim, ngauss

  IF n_params() LT 1 THEN BEGIN
      message,'-Syntax: gauss = fastem->read_gauss(gaussfile [, ndim, ngauss])'
  ENDIF 

  openr, lun, gaussfile, /get_lun
  self->read_gauss_header, lun, ngauss, ndim

  ;; read into standard gaussian struct
  struct = replicate( self->gauss_struct(ndim), ngauss)
  readf, lun, struct

  free_lun, lun

  return,struct

END 







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Create red galaxy spectra color inputs for the EM code fastem
;; Bin in z and plot colors, write colors out to disk.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; holds both inputs and outputs
FUNCTION fastem::redspec_datadir
  dir = '/net/cheops2/home/esheldon/work/red_galaxy_colors/'
  return,dir
END 

PRO fastem::bin_color_byz, z, color, bz, bcolor, bcolor_sdev, index_array

;  zbin = 0.02
;  binner, z, color, zbin, bz, bcolor, bcolor_err

  nperbin = 6400
  binner_bynum, z, color, nperbin, bz, bcolor, bcolor_err, $
    ybinned_sdev=bcolor_sdev, index_array=index_array

END 


PRO fastem::redspec_colors, sp, cs

  IF n_elements(sp) EQ 0 THEN BEGIN 
      query = 'select modelflux,modelflux_ivar,absmodelmag,extinction,primtarget,sectarget,z from specgal'
      print,'Query: ',query
      sp = pgsql_query_async(query)
  ENDIF 

  IF n_elements(cs) EQ 0 THEN BEGIN 

      ;; observed
      u = 22.5 - 2.5*alog10(sp.modelflux[0]) - sp.extinction[0]
      g = 22.5 - 2.5*alog10(sp.modelflux[1]) - sp.extinction[1]

      cs = {umg:u-g}
      u = u

      r = 22.5 - 2.5*alog10(sp.modelflux[2]) - sp.extinction[2]
      cs = create_struct(cs, 'gmr', g-r)
      g = 0

      i = 22.5 - 2.5*alog10(sp.modelflux[3]) - sp.extinction[3]
      cs = create_struct(cs, 'rmi', r-i)
      i=0

      ;; absolute
      mumg = sp.absmodelmag[0] - sp.absmodelmag[1]
      cs = create_struct(cs, 'mumg', mumg)
      mumg = 0

      mgmr = sp.absmodelmag[1] - sp.absmodelmag[2]
      cs = create_struct(cs, 'mgmr', mgmr)
      mgmr = 0
      mrmi = sp.absmodelmag[2] - sp.absmodelmag[3]
      cs = create_struct(cs, 'mrmi', mrmi)
      mrmi = 0

  ENDIF 


  !p.multi=[0,2,2]
  wred = $
    where( (sp.primtarget and sdss_flag('target','galaxy_red')) ne 0 $
           AND $
           cs.gmr gt 0 and cs.gmr lt 3 and cs.rmi gt 0 and cs.rmi lt 1.5)
  wred2 = $
    where( (sp.primtarget and sdss_flag('target','galaxy_red')) ne 0 $
           and $
           cs.gmr gt 0 and cs.gmr lt 3 and cs.rmi gt 0 and cs.rmi lt 1.5 AND $
           sp.z GT 0.15)


  ;; The zmean and index_array should be same
  self->bin_color_byz, sp[wred].z, cs.gmr[wred], zmean, gmrmean, gmrerr
  self->bin_color_byz, $
    sp[wred].z, cs.rmi[wred], zmean, rmimean, rmierr, index_array

  print,'zmean = ',zmean
  ;; Absolute mag colors
  xtitle = 'abs umg'
  ytitle = 'abs gmr'
  xrange = [0.5, 3.0]
  yrange = [0.5, 1.5]
  ploth, cs.mumg[wred], cs.mgmr[wred], xrange=xrange, yrange=yrange,$
    xtitle=xtitle,ytitle=ytitle
  ploth, cs.mumg[wred2], cs.mgmr[wred2], xrange=xrange, yrange=yrange,$
    xtitle=xtitle,ytitle=ytitle, $
    title='zcut'

  xtitle = 'abs gmr'
  ytitle = 'abs rmi'
  xrange = [0.8, 1.2]
  yrange = [0.32, 0.55]
  ploth, cs.mgmr[wred], cs.mrmi[wred], xrange=xrange, yrange=yrange,$
    xtitle=xtitle,ytitle=ytitle
  ploth, cs.mgmr[wred2], cs.mrmi[wred2], xrange=xrange, yrange=yrange,$
    xtitle=xtitle,ytitle=ytitle, $
    title='zcut'

  IF 0 THEN BEGIN 
      wp = where(cs.mgmr[wred] GT 0.8 AND cs.mgmr[wred] LT 1.2 AND $
                 cs.mrmi[wred] GT 0.32 AND cs.mrmi[wred] LT 0.55)
      colprint, $
        file='~/tmp/redgal_gmr_rmi.dat', $
        cs.mgmr[wred[wp]], cs.mrmi[wred[wp]]
      
      wp = where(cs.mgmr[wred2] GT 0.8 AND cs.mgmr[wred2] LT 1.2 AND $
                 cs.mrmi[wred2] GT 0.32 AND cs.mrmi[wred2] LT 0.55)
      colprint, $
        file='~/tmp/redgal_gmr_rmi_zcut.dat', $
        cs.mgmr[wred2[wp]], cs.mrmi[wred2[wp]]
  ENDIF 



  ;; Observed colors
  xtitle = 'observed umg'
  ytitle = 'observed gmr'
  xrange = [1.0, 2.5]
  yrange = [0.4, 2.0]
  ploth, cs.umg[wred], cs.gmr[wred], xrange=xrange, yrange=yrange
  ploth, cs.umg[wred2], cs.gmr[wred2], xrange=xrange, yrange=yrange, $
    title='zcut'

  xtitle = 'observed gmr'
  ytitle = 'observed rmi'
  xrange = [-0, 2.5]
  yrange = [0.2, 1.3]
  ploth, cs.gmr[wred], cs.rmi[wred], xrange=xrange, yrange=yrange


;  nbin = n_elements(zmean)
;  FOR i=0L, nbin-1 DO oplot, gmrmean
  print,zmean
;  oplot, gmrmean, rmimean,color=!red
 
 
  plotz = [0.02, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, $
           0.40, 0.44, 0.48, 0.52, 0.56, 0.60]

  interp_gmr = interpol(gmrmean, zmean, plotz)
  interp_rmi = interpol(rmimean, zmean, plotz)
  interp_gmrerr = interpol(gmrerr, zmean, plotz)
  interp_rmierr = interpol(rmierr, zmean, plotz)

  nz = n_elements(plotz)
  clr = rainbow_rgb( arrscl(findgen(nz), 225.0, 360.0) )

  FOR i=0L, nz-1 DO BEGIN 
      oploterror, $
        [interp_gmr[i]], [interp_rmi[i]], $
        [interp_gmrerr[i]], [interp_rmierr[i]], $
        psym=8, symsize=1, color=clr[i], errc=clr[i], hat=0

      IF i NE nz-1 THEN BEGIN 
          oplot, $
            [interp_gmr[i], interp_gmr[i+1]], $
            [interp_rmi[i], interp_rmi[i+1]], $
            color=clr[i], thick=1.5
      ENDIF 

  ENDFOR 

;  oploterror, interp_gmr, interp_rmi, interp_gmrerr, interp_rmierr, $
;    psym=8, color=!red, errc=!red, hat=0

  ploth, cs.gmr[wred2], cs.rmi[wred2], xrange=xrange, yrange=yrange, $
    title='zcut'

  !p.multi = 0

  ;; Write colors in z bins out to files
  nz = n_elements(zmean)

  outdir = self->redspec_datadir()
  dstruct = {gmr:0.0, rmi:0.0}
  FOR i=0L, nz-1 DO BEGIN 
          
      w = *index_array[i]
      w = wred[w]
          
      nw = n_elements(w)
      desc_struct = $
        {minz: min(sp[w].z), $
         maxz: max(sp[w].z), $
         meanz: float(zmean[i]) $
        }
          
      dat_struct = replicate(dstruct, nw) 
      dat_struct.gmr = cs.gmr[w]
      dat_struct.rmi = cs.rmi[w]
      
      desc_file = self->redspec_descfile(i+1)
      print
      print,'Writing description file: ',desc_file
      write_idlstruct, desc_struct, desc_file, /ascii

      data_file = self->redspec_datafile(i+1)
      print
      print,'Writing data file: ',data_file
      ascii_write, dat_struct, data_file


  ENDFOR 


END 

FUNCTION fastem::redspec_datafile, number

  IF n_elements(number) EQ 0 THEN BEGIN 
      message,'-Syntax: em->redspec_datfile(number)'
  ENDIF 

  dir = self->redspec_datadir()
  nn = n_elements(number)
  IF nn EQ 1 THEN files = '' ELSE files = strarr(nn)

  FOR i=0L, nn-1 DO BEGIN 
      istr = strn(long(number[i]), len=2, padchar='0')
      files[i] = dir + 'red_galaxy_colors'+istr+'.dat'
  ENDFOR 

  return,files

END 
FUNCTION fastem::redspec_descfile, number

  IF n_elements(number) EQ 0 THEN BEGIN 
      message,'-Syntax: em->redspec_descfile(number)'
  ENDIF 

  dir = self->redspec_datadir()
  nn = n_elements(number)
  IF nn EQ 1 THEN files = '' ELSE files = strarr(nn)

  FOR i=0L, nn-1 DO BEGIN 
      istr = strn(long(number[i]), len=2, padchar='0')
      files[i] = dir + 'red_galaxy_desc'+istr+'.st'
  ENDFOR
  return,files
END 
FUNCTION fastem::redspec_gaussfile, number

  IF n_elements(number) EQ 0 THEN BEGIN 
      message,'-Syntax: em->redspec_gaussfile(number)'
  ENDIF 

  dir = self->redspec_datadir()
  nn = n_elements(number)
  IF nn EQ 1 THEN files = '' ELSE files = strarr(nn)

  FOR i=0L, nn-1 DO BEGIN 
      istr = strn(long(number[0]), len=2, padchar='0')
      files[i] = dir + 'red_galaxy_colors'+istr+'.cen'
  ENDFOR
  return,files

END 


FUNCTION fastem::redspec_read_datafile, number

  IF n_elements(number) EQ 0 THEN BEGIN 
      message,'-Syntax: struct = em->redspec_read_datfile(number)'
  ENDIF 

  file = self->redspec_datafile(number[0])

  tstruct = {gmr:0.0, rmi:0.0}

  print
  read_struct, file, tstruct, struct
  return,struct

END 
FUNCTION fastem::redspec_read_all_datafiles

  nn = self->redspec_n_zbins()

  data = ptrarr(nn)
  number = 1+lindgen(nn)

  FOR i=0L, nn-1 DO BEGIN 
      data[i] = ptr_new( self->redspec_read_datafile(number[i]), $
                         /no_copy)
  ENDFOR 
  return,data

END 

FUNCTION fastem::redspec_read_descfile, number

  IF n_elements(number) EQ 0 THEN BEGIN 
      message,'-Syntax: desc = em->redspec_read_descfile(number)'
  ENDIF 

  file = self->redspec_descfile(number)
;  print
;  print,'Reading description file: ',file
  struct = read_idlstruct_multi(file)
  return,struct

END 

;; Write config scripts for running fastem
PRO fastem::redspec_write_config

  dir = self->bin_color_dir()
  files = file_search(dir + 'red_galaxy_colors*.dat')

  config_files = repstr(files, '.dat', '.cfg')
  center_files = repstr(files, '.dat', '.cen')
  ps_files = repstr(files, '.dat', '.ps')
  nf = n_elements(config_files)

  FOR i=0L, nf-1 DO BEGIN 
      print,'Writing to config file: ',config_files[i]
      openw, lun, config_files[i], /get_lun

      dirsep, files[i], tdir, tfile
      shortname = repstr(tfile, '.dat','')
      printf,lun,'BEGIN_DATASETS'
      printf,lun,shortname+'  '+files[i]
      printf,lun,'END_DATASETS'
      printf,lun
      printf,lun,'in  '+shortname
      printf,lun,'nsteps  20000'
      printf,lun,'nsecs  -1'
      printf,lun,'savecenter  '+center_files[i]
      printf,lun,'saveps  '+ps_files[i]
      printf,lun,'nclusters 2'
      printf,lun,'varyclusters false'
      printf,lun,'draw false'

      free_lun, lun
  ENDFOR 

END 



FUNCTION fastem::redspec_multigauss_struct, ngauss, ndim

  tstruct = self->gauss_struct(ndim, /extra)
  struct = {ngauss: ngauss, $
            ndim: ndim, $
            minz: 0.0, $
            maxx: 0.0, $
            meanz: 0.0, $
            gauss: replicate(tstruct, ngauss)}
  return,struct
  
END 

FUNCTION fastem::redspec_read_gauss, number

  IF n_elements(number) EQ 0 THEN BEGIN 
      message,'-Syntax: gauss = em->redspec_read_gauss(number)'
  ENDIF 


  gaussfile = self->redspec_gaussfile(number)
  descfile = self->redspec_descfile(number)
  print
  print,'Reading gaussian center file: ',gaussfile

  tst = self->read_gauss(gaussfile, ndim, ngauss)

  struct = self->redspec_multigauss_struct(ngauss, ndim)

  ;; copy in extra info
  FOR i=0L, ngauss-1 DO BEGIN 

      ;; struct with extra info.  Copy in what we already know
      exstruct = self->gauss_struct(ndim, /extra)
      copy_struct, tst[i], exstruct

      ;; new stuff
      detcov = determ(tst[i].cov, /double, /check)
      invcov = invert(tst[i].cov, status)
      detinvcov = determ(invcov, /double, /check)

      exstruct.detcov = detcov
      exstruct.invcov = invcov
      exstruct.detinvcov = detinvcov

      struct.gauss[i] = exstruct

  ENDFOR 


  ;; The description file
  dst = self->redspec_read_descfile(number)

  struct.minz = dst.minz
  struct.maxx = dst.maxz
  struct.meanz = dst.meanz
  return,struct

END 

FUNCTION fastem::redspec_n_zbins

  dir = self->redspec_datadir()
  files = file_search(dir+'red_galaxy_colors*.dat', count=n_zbins)

  return,n_zbins

END 

;; read all the output gaussian files
FUNCTION fastem::redspec_read_allgauss

  n_zbins = self->redspec_n_zbins()
  bin_numbers = 1+lindgen(n_zbins)
  
  struct = {$
             nbin: n_zbins, $
             gaussians: ptrarr(n_zbins) $
           }
  
  FOR i=0L, n_zbins-1 DO BEGIN 

      struct.gaussians[i] = $
        ptr_new( self->redspec_read_gauss(bin_numbers[i]), /no_copy)

  ENDFOR 

  return,struct


END 













FUNCTION fastem::redspec_allgauss
  IF n_elements(*self.allgauss) EQ 0 THEN BEGIN 
      allgauss = self->redspec_read_allgauss()
      self.allgauss = ptr_new(allgauss, /no_copy)
  ENDIF 

  allgauss = *self.allgauss
  return,allgauss
END 
FUNCTION fastem::redspec_reduced_allgauss


  IF n_elements(*self.reduced_allgauss) EQ 0 THEN BEGIN 

      allgauss = self->redspec_allgauss()
      ndim = (*allgauss.gaussians[0]).ndim

      ;; just keep track of most important gaussian right now
      
      tst = self->gauss_struct(ndim, /extra)
      tst = create_struct('meanz', 0.0, tst)
      
      reduced_allgauss = replicate(tst, allgauss.nbin)

      FOR i=0L, allgauss.nbin-1 DO BEGIN 

          gauss_st = *allgauss.gaussians[i]
          reduced_allgauss[i].meanz = gauss_st.meanz

          gauss = gauss_st.gauss

          maxprob = max(gauss.prob, maxid)

          tst = self->gauss_struct(ndim, /extra)
          copy_struct, gauss[maxid], tst
          tst = create_struct('meanz', gauss_st.meanz, tst)

          reduced_allgauss[i] = tst

          ;; only using one gauss here
          reduced_allgauss[i].prob = 1.0
      ENDFOR 
      self.reduced_allgauss = ptr_new(reduced_allgauss, /no_copy)
  ENDIF 

  reduced_allgauss = *self.reduced_allgauss
  return,reduced_allgauss

END 




;; For 2-d color data
PRO fastem::redspec_display_gauss

  simpctable, colorlist=colorlist
  w = where(colorlist NE !red AND colorlist NE !dodgerBlue)
  colorlist = colorlist[w]

  gmr_range = [0.0, 2.5]
  rmi_range = [0.2, 1.3]
  xtitle = 'g-r'
  ytitle = 'r-i'

  fem = self->redspec_allgauss()
  data = self->redspec_read_all_datafiles()
  n_z_bin = fem.nbin

  cencolor = [!red, !dodgerBlue]

  FOR i=0L, n_z_bin-1 DO BEGIN 

      cstruct = *data[i]
      IF i EQ 0 THEN BEGIN 
          plot, $
            cstruct.gmr, cstruct.rmi, $
            psym=3, $
            xrange = gmr_range, xstyle=1+1, yrange=rmi_range, ystyle=1+2, $
            xtitle=xtitle, ytitle=ytitle, /iso
      ENDIF ELSE BEGIN 
          oplot, $
            cstruct.gmr, cstruct.rmi, $
            psym=3;, color=colorlist[i]
      ENDELSE 
  ENDFOR 

  FOR i=0L, n_z_bin-1 DO BEGIN 

      print,'-------------------------------------------------'

      cstruct = *data[i]
      gstruct = *fem.gaussians[i]
      
;      IF i EQ 0 THEN BEGIN 
;          plot, $
;            cstruct.gmr, cstruct.rmi, $
;            psym=3, $
;            xrange = gmr_range, xstyle=1+1, yrange=rmi_range, ystyle=1+2, $
;            xtitle=xtitle, ytitle=ytitle ;, /iso
;      ENDIF ELSE BEGIN 
;          oplot, $
;            cstruct.gmr, cstruct.rmi, $
;            psym=3
;      ENDELSE 
      s = reverse( sort(gstruct.gauss.prob) )
      print,'zmean = '+ntostr(gstruct.meanz)
      FOR gi=0L, gstruct.ngauss-1 DO BEGIN 

          ii = s[gi]

          gauss = gstruct.gauss[ii]

          print,'prob = '+ntostr(gauss.prob)

          IF gauss.prob GT 0.01 THEN BEGIN 
              oplot, [gauss.cen[0]], [gauss.cen[1]], $
                     psym=7, color=cencolor[gi], symsize=1.5

              self->calc_ellipse, gauss.cov, smax, smin, theta
              tvellipse, smax, smin, gauss.cen[0], gauss.cen[1], theta,$
                         color=cencolor[gi], /data
          ENDIF ELSE BEGIN 
              print,'** Skipping'
          ENDELSE 
          
      ENDFOR 

      key = prompt_kbrd()
      IF key EQ 'q' THEN BEGIN 
          self->redspec_free, fem
          return 
      ENDIF 
  ENDFOR 



  self->redspec_free, fem

END 

PRO fastem::calc_ellipse, cov, smax, smin, theta

  
  invcov = invert( cov, status, /double)
  IF status NE 0 THEN message,'Could not invert'


  ;; Eigenvalues


  ;; This was my way, it works too but their routines are probably
  ;; more stable

  IF 0 THEN BEGIN 
      lam11 = 0.5*( invcov[0,0]+invcov[1,1] + $
                    sqrt(4.0*invcov[0,1]^2 + (invcov[0,0] - invcov[1,1] )^2 ) )
      lam22 = 0.5*( invcov[0,0]+invcov[1,1] - $
                    sqrt(4.0*invcov[0,1]^2 + (invcov[0,0] - invcov[1,1] )^2 ) )
      
      cos2 = (invcov[1,1] - (lam11/lam22)*invcov[0,0])/(lam22-lam11^2/lam22)
      sin2 = 1.0 - cos2
      
      
      ;; minus because I screwed up the transform
      sgn = sign(sin2/cos2)
      theta = ( atan( sqrt(sin2),sqrt(cos2)) )*180.0/!pi
      
      lammax = max([lam11, lam22], min=lammin)
      smin = sqrt(1.0/lammax)
      smax = sqrt(1.0/lammin)
      
  ENDIF ELSE BEGIN 

      A = invcov
      TRIRED, A, D, E
      TRIQL, D, E, A

;      print,'Evals = ',D
;      print,'Eigenvectors: ',A
;      print,'Eigenvec angle1: ',atan(A[0,1],A[0,0])*180.0/!pi
;      print,'Eigenvec angle2: ',atan(A[1,0],A[1,1])*180.0/!pi

      lammax = max(D, di, min=lammin)
      smin = sqrt(1.0/lammax)
      smax = sqrt(1.0/lammin)

      theta = atan(A[di,0],A[di,1])*180.0/!pi

  ENDELSE 

  print,'Smax: '+ntostr(smax)+'  Smin: '+ntostr(smin)+' theta: '+ntostr(theta)

END 



FUNCTION fastem::redspec_interpolate, z

  reduced_allgauss = self->redspec_reduced_allgauss()
  ndim = n_elements(reduced_allgauss[0].cen)

  outstruct = self->gauss_struct(ndim,/extra)

  nz = n_elements(z) 
  outstruct = replicate(outstruct, nz)

  ;; Interpolate the results
  outstruct.cen[0] = $
    interpol(reduced_allgauss.cen[0], reduced_allgauss.meanz, z)
  outstruct.cen[1] = $
    interpol(reduced_allgauss.cen[1], reduced_allgauss.meanz, z)

  outstruct.cov[0,0] = $
    interpol(reduced_allgauss.cov[0,0], reduced_allgauss.meanz, z)
  outstruct.cov[0,1] = $
    interpol(reduced_allgauss.cov[0,1], reduced_allgauss.meanz, z)
  outstruct.cov[1,0] = outstruct.cov[0,1]
  outstruct.cov[1,1] = $
    interpol(reduced_allgauss.cov[1,1], reduced_allgauss.meanz, z)

  outstruct.detcov = $
    interpol(reduced_allgauss.detcov, reduced_allgauss.meanz, z)

  outstruct.invcov[0,0] = $
    interpol(reduced_allgauss.invcov[0,0], reduced_allgauss.meanz, z)
  outstruct.invcov[0,1] = $
    interpol(reduced_allgauss.invcov[0,1], reduced_allgauss.meanz, z)
  outstruct.invcov[1,0] = outstruct.invcov[0,1]
  outstruct.invcov[1,1] = $
    interpol(reduced_allgauss.invcov[1,1], reduced_allgauss.meanz, z)

  outstruct.detinvcov = $
    interpol(reduced_allgauss.detinvcov, reduced_allgauss.meanz, z)

  return,outstruct





  invcov = dblarr(2,2,nz)
  detcov = dblarr(nz)
  FOR i=0L, nz-1 DO BEGIN 
      detcov[i] = determ(outstruct[i].cov)
      invcov[*,*,i] = invert(outstruct[i].cov)
;      print,z[i],outstruct[i].detcov, det, outstruct[i].invcov[0,0], invcov[0,0]
  ENDFOR 

  !p.multi = [0,2,2]
  plot, z, outstruct.detcov
  oplot,z, detcov, color=!red

  plot, z, outstruct.invcov[0,0]
  oplot,z,invcov[0,0,*], color=!red

  plot, z, outstruct.invcov[0,1]
  oplot,z,invcov[0,1,*], color=!red

  plot, z, outstruct.invcov[1,1]
  oplot,z,invcov[1,1,*], color=!red
  !p.multi = 0

  return,outstruct

END 


;; free the pointer array and struct
PRO fastem::redspec_free, struct

  IF n_elements(struct) EQ 0 THEN BEGIN 
      print,'-Syntax: em->redspec_free, struct'
      print,'Struct comes from em->redspec_read_allgauss()'
      return
  ENDIF 
  ptr_free, struct.gaussians
  delvarx, struct

END 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; z-spec characerization for photozs.  For now this only has simulated 
;; data
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; The various files 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION fastem::zphot_zspec_dir
  return,'/net/cheops2/home/esheldon/work/zphot_error/'
END 

;; read in the photoz-specz pairs
FUNCTION fastem::dessim_dir
  return,self->zphot_zspec_dir()+'des_simdata/'
END 
FUNCTION fastem::dessim_zphot_zspec_read

  dir = self->dessim_dir()
  file = dir + 'zpzstrain.dat'
  tstruct = {zphot: 0.0, specz: 0.0}
  read_struct, file, tstruct, struct

  return,struct
END 


FUNCTION fastem::dessim_input_file, binnum, equal=equal
  IF n_elements(binnum) EQ 0 THEN BEGIN 
      print,'-Syntax: file = em->dessim_input_file(binnum, /equal)'
  ENDIF 
  dir = self->dessim_dir()
  file = dir + 'pzphot'
  IF keyword_set(equal) THEN file = file+'_equal'
  file = file + '_bin'+strn(binnum,len=2,padchar='0')+'.dat'

  return,file
END
FUNCTION fastem::dessim_input_read, binnum, equal=equal
  file = self->dessim_input_file(binnum, equal=equal)
  rdfloat, file, zphot
  return,zphot
END  
FUNCTION fastem::dessim_desc_file, equal=equal
  dir = self->dessim_dir()
  file = dir + 'pzphot'
  IF keyword_set(equal) THEN file = file+'_equal'
  file = file + '_desc.dat'
  return,file
END 
FUNCTION fastem::dessim_desc_read, equal=equal
  file = self->dessim_desc_file(equal=equal)
  struct = {bin:0, minz:0.0, maxz:0.0, meanz:0.0, number:0L}
  nlines = numlines(file)
  nread =  nlines-2
  struct = replicate(struct, nread)
  
  openr,lun,file,/get_lun
  line=''
  readf, lun, line
  readf, lun, line

  readf, lun, struct

  free_lun, lun
  return,struct
END 
FUNCTION fastem::dessim_config_file, binnum, equal=equal
  file = self->dessim_input_file(binnum,equal=equal)
  file = repstr(file, '.dat', '.cfg')
  return,file
END 
FUNCTION fastem::dessim_gaussian_file, binnum, equal=equal
  file = self->dessim_input_file(binnum,equal=equal)
  file = repstr(file, '.dat', '.cen')
  return,file
END 
FUNCTION fastem::dessim_gaussian_read, binnum, equal=equal, status=status
  status = 1
  file = self->dessim_gaussian_file(binnum, equal=equal)
  IF fexist(file) THEN BEGIN 
      status = 0
      gauss = self->gauss_read_1d(file)
  ENDIF ELSE BEGIN 
      print,'File does not exist: ',file
      gauss = -1
  ENDELSE 
  return,gauss
END 
FUNCTION fastem::dessim_script_file, equal=equal
  file = self->dessim_desc_file(equal=equal)
  file = repstr(file, '_desc.dat', '_run.sh')
  return,file
END 



;; make the input files for fastem

PRO fastem::dessim_write_config, binnum, equal=equal

  datafile = self->dessim_input_file(binnum,equal=equal)
  conf = self->dessim_config_file(binnum,equal=equal)
  gaussfile = self->dessim_gaussian_file(binnum,equal=equal)
  dirsep, datafile, tdir, tfile
  name = repstr(tfile, '.dat')

  print,'Writing config file: ',conf
  openw, lun, conf, /get_lun
  
  printf,lun,'BEGIN_DATASETS'
  printf,lun,name+'  '+datafile
  printf,lun,'END_DATASETS'
  printf,lun
  printf,lun,'in  '+name
  printf,lun,'savecenter  '+gaussfile
  printf,lun,'nsecs  -1'
  printf,lun,'nsteps  10000'

  free_lun, lun
END 

PRO fastem::dessim_write_input

  struct = self->dessim_zphot_zspec_read()
  nst = n_elements(struct)
  outst = {zphot: 0.0}


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; binned with equal sized z bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  
  print,'Evenly binned'
  print,'----------------------------------------'
  minz = min(struct.specz, max=maxz)
  nbin = 30
  zbin = (maxz-minz)/nbin
  h = histogram(struct.specz, min=minz, max=maxz, binsize = zbin, rev=rev)

  nh = n_elements(h)

  desc_file = self->dessim_desc_file()
  script_file = self->dessim_script_file()

  print
  print,'Opening descriptioin file: ',desc_file
  openw, desc_lun, desc_file, /get_lun
  printf,desc_lun,$
         'bin','minz','maxz','meanz','number', $
         format = '(5A10)'
  printf,desc_lun,strjoin(replicate('-',5*10))

  print,'Opening script file: ',script_file
  openw, script_lun, script_file, /get_lun
  printf,script_lun, '#!/bin/sh'
  printf,script_lun




  For i=0L, nh-1 DO BEGIN 
      IF rev[i] NE rev[i+1] THEN BEGIN 
          w = rev[ rev[i]:rev[i+1]-1 ]
          nn = n_elements(w)
          minz = min(struct[w].specz, max=maxz)
          meanz = mean(struct[w].specz)


          outstruct = replicate(outst, nn)
          outstruct.zphot = struct[w].zphot

          file = self->dessim_input_file(i)
          print,'Writing to file: ',file
          ascii_write, outstruct, file


          printf,desc_lun,$
                 i,minz,maxz,meanz,nn, $
                 format = '(I10,F10,F10,F10,I10)'

          printf,script_lun,$
                 '~/src/fastem_linux_gui/fastem_wx config '+$
                 self->dessim_config_file(i)
          printf,script_lun,'sleep 5'



          self->dessim_write_config, i


      ENDIF ELSE BEGIN 
          nn = 0.0
          meanz = -1.0
      ENDELSE 

;      print,i,meanz,nn
  ENDFOR 

  free_lun, desc_lun
  free_lun, script_lun

  spawn,['chmod','755',script_file],/noshell


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; binned with equal number per bin
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  print,'Equal number per bin'
  print,'----------------------------------------'

  nperbin = long(float(nst)/nbin)
  s = sort(struct.specz)
  l = lindgen(nst)
  h = histogram(l, bin=nperbin, rev=rev)
  
  desc_file = self->dessim_desc_file(/equal)
  script_file = self->dessim_script_file(/equal)


  print
  print,'Opening descriptioin file: ',desc_file
  openw, desc_lun, desc_file, /get_lun
  printf,desc_lun,$
         'bin','minz','maxz','meanz','number', $
         format = '(5A10)'
  printf,desc_lun,strjoin(replicate('-',5*10))

  print,'Opening script file: ',script_file
  openw, script_lun, script_file, /get_lun
  printf,script_lun, '#!/bin/sh'
  printf,script_lun

  FOR i=0L, nh-1 DO BEGIN 
      IF rev[i] NE rev[i+1] THEN BEGIN 
          w = rev[ rev[i]:rev[i+1]-1 ]
          w = s[w]
          nn = n_elements(w)

          minz = min(struct[w].specz, max=maxz)
          meanz = mean(struct[w].specz)


          outstruct = replicate(outst, nn)
          outstruct.zphot = struct[w].zphot

          file = self->dessim_input_file(i, /equal)
          print
          print,'Writing to file: ',file
          ascii_write, outstruct, file



          printf,desc_lun,$
                 i,minz,maxz,meanz,nn, $
                 format = '(I10,F10,F10,F10,I10)'

          printf,script_lun,$
                 '~/src/fastem_linux_gui/fastem_wx config '+$
                 self->dessim_config_file(i, /equal)
          printf,script_lun,'sleep 5'


          self->dessim_write_config, i, /equal

      ENDIF ELSE BEGIN 
          nn = 0.0
          meanz = -1.0
      ENDELSE 

;      print,i,meanz,nn
  ENDFOR 

  free_lun, desc_lun
  free_lun, script_lun

  spawn,['chmod','755',script_file],/noshell

END 



PRO fastem::dessim_plot_gaussfit, binnums, equal=equal, dops=dops, nolegend=nolegend

  IF keyword_set(dops) THEN BEGIN 
      dir = '/net/cheops2/home/esheldon/plots/photoz/uchicago/em_fits/'
      file = dir + 'pzphot'
      IF keyword_set(equal) THEN file = file+'_equal'
      file = file+'_fits.ps'
      begplot,name=file,/color

      !p.multi = [0,0,2]
  ENDIF 

  desc = self->dessim_desc_read(equal=equal)

  zbin = 0.01
  nz = 1000

  nplot = n_elements(binnums)

  FOR i=0L, nplot-1 DO BEGIN 

      num = binnums[i]

      data = self->dessim_input_read(num, equal=equal)
      gauss = self->dessim_gaussian_read(num, equal=equal, status=status)

      IF status EQ 0 THEN BEGIN 
      
          minz = min(data, max=maxz)
          
          gmod = self->gaussmodel_1d(gauss, minz, maxz, nz)
          
          plothist, data, bin=zbin, /norm, $
                    xtitle='zphot',ytitle='P(zphot | zspec)'
          self->gmix_plot_model_1d, gmod, /overplot, nolegend=nolegend
          
          leg = ['minzspec: '+ntostr(desc[num].minz), $
                 'maxzspec: '+ntostr(desc[num].maxz), $
                 'meanzspec: '+ntostr(desc[num].meanz), $
                 'nobj: '+ntostr(desc[num].number)]
          legend,leg,/right, box=0, charsize=1
          
          
          key = prompt_kbrd()
          IF key EQ 'q' THEN return
      ENDIF 
  ENDFOR 

  IF keyword_set(dops) THEN BEGIN 
      endplot
      !p.multi=0
  ENDIF 

END 







;; The test file Carlos gave me
FUNCTION fastem::fastem_zdiff_read
  file = '/net/cheops2/home/esheldon/work/zphot_error/zphot_m_zspec.dat'
  rdfloat, file, zdiff
  return,zdiff
END 



;; read the 1-d files

FUNCTION fastem::gauss_read_1d, file

  IF n_elements(file) EQ 0 THEN BEGIN 
      message,'-Syntax: gauss = em->gauss_read_1d(file)'
  ENDIF 
  nlines = numlines(file)
  ngauss = nlines-2

  struct = {number: 0,  $
            prob:   0d, $
            cen:    0d, $
            cov:    0d}
            

  struct = replicate(struct, ngauss)
  openr, lun, file, /get_lun

  ;; skip header lines
  line = ''
  readf, lun, line
  print,line
  readf, lun, line
  print,line

  readf, lun, struct

  free_lun, lun

  return,struct

END 


FUNCTION fastem::gaussmodel_1d, gstr, minx, maxx, nx

  ngauss = n_elements(gstr)

  probs = gstr.prob
  means = gstr.cen
  sigmas = sqrt(gstr.cov)


  xi = arrscl( findgen(nx), minx, maxx )

  model = fltarr(nx)
  
  FOR i=0L, ngauss-1 DO BEGIN 
      tmodel = probs[i]*gaussprob(xi, means[i], sigmas[i])
      
      model[*] = model[*] + tmodel[*]
  ENDFOR 

  norm = qgauss(model, xi, 200)

  model = model/norm
  model_struct = {ngauss:ngauss, $
                  x: xi, $
                  model: model}
  FOR i=0L, ngauss-1 DO BEGIN 

      name = 'model'+strn(i+1)

      tmodel = probs[i]*gaussprob(xi, means[i], sigmas[i])
      tmodel = tmodel/norm

      model_struct = $
        create_struct(model_struct, $
                      name, tmodel)

  ENDFOR 

  return,model_struct

END 
























PRO fastem__define

  struct = { fastem, $
             allgauss: ptr_new(), $
             reduced_allgauss: ptr_new() $
           }

END 
