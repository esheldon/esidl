FUNCTION sdss_shapelet::init
  return,1
END 

FUNCTION sdss_shapelet::idstring, strOrRun, rerun, camcol, field, id
  
  IF size(strOrRun,/tname) EQ 'STRUCT' THEN BEGIN 

      run    = strOrRun.run
      rerun  = strOrRun.rerun
      camcol = strOrRun.camcol
      field  = strOrRun.field
      id     = strOrRun.id

  ENDIF ELSE IF n_params() EQ 5 THEN BEGIN 

      run = strOrRun

  ENDIF ELSE BEGIN
      print,'-Syntax: idString = sh->idstring(struct)'
      print,'                -- OR --'
      print,'         idString = sh->idstring(run,rerun,camcol,field,id)'
      return,'ERROR'
  ENDELSE 

  runstr = run2string(run)
  camcolstr = ntostr(long(camcol))
  rerunstr = ntostr(long(rerun))
  fieldstr = field2string(field)
  idstr = id2string(id)

  idstring = runstr+'-'+rerunstr+'-'+camcolstr+'-'+fieldstr+'-'+idstr

  return,idstring

END 







FUNCTION sdss_shapelet::dir
  outDir = expand_tilde('~/shapelet_outputs')
  return,outdir
END 
FUNCTION sdss_shapelet::image_dir, stripe, clr, createdir=createdir
  on_error, 2
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: dir=ssh->image_dir(stripe, clr, /createdir)'
      print
      message,'Halting'
  ENDIF 
  dir = self->dir()
  sstr = stripe2string(stripe)
  cstr = !colors[clr]
  imageDir = concat_dir(dir, 'stripe'+sstr+'-'+cstr+'-images')

  IF keyword_set(createdir) AND NOT fexist(imageDir) THEN BEGIN 
      file_mkdir, imageDir
  ENDIF 

  return,imageDir
END 
FUNCTION sdss_shapelet::output_file, stripe, clr, dtype_in, nmax
  on_error, 2
  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: file=ssh->output_file(stripe, clr, dtype, nmax)'
      print
      message,'Halting'
  ENDIF 
  dir = self->dir()
  sstr = stripe2string(stripe)
  cstr = !colors[clr]
  nstr = strn(nmax,len=2,padchar='0')

  dtype = strlowcase( strtrim(dtype_in,2) )

  outfile = $
    'stripe'+sstr+'-vagc-shapelets-'+cstr+'-'+dtype+'-nmax'+nstr+'.st'
  outFile = concat_dir(dir, outfile)
  return, outfile

END 
FUNCTION sdss_shapelet::output_read, stripe, clr, dtype, nmax, status=status
  on_error, 2
  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: file=ssh->output_read(stripe, clr, dtype, nmax)'
      print
      message,'Halting'
  ENDIF 
  file = self->output_file(stripe,clr,dtype,nmax)
  print,'Reading file: ',file
  struct = read_idlstruct(file, status=status)
  return,struct
END 



FUNCTION sdss_shapelet::combined_stripe_file, stripe, dtype, nmax
  on_error, 2
  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: file=ssh->combined_stripe_file(stripe, dtype, nmax)'
      print
      message,'Halting'
  ENDIF 
  tfile = self->output_file(stripe, 2, dtype, nmax)
  file = repstr(tfile, '-r', '')
  return,file
END 
FUNCTION sdss_shapelet::combined_stripe_read, stripe, dtype, nmax, status=status
  on_error, 2
  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: file=ssh->combined_stripe_read(stripe, dtype, nmax)'
      print
      message,'Halting'
  ENDIF 
  file = self->combined_stripe_file(stripe,dtype,nmax)
  print,'Reading file: ',file
  struct = read_idlstruct_multi(file, status=status)
  return,struct
END 

FUNCTION sdss_shapelet::combined_file, dtype, nmax
  on_error, 2
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: file=ssh->combined_file(dtype, nmax)'
      print
      message,'Halting'
  ENDIF 
  tfile = self->combined_stripe_file(16, dtype, nmax)
  file = repstr(tfile, 'stripe16-', '')
  return,file
END 
FUNCTION sdss_shapelet::combined_read, dtype, nmax, status=status
  on_error, 2
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: file=ssh->combined_read(dtype, nmax)'
      print
      message,'Halting'
  ENDIF 
  file = self->combined_file(dtype,nmax)
  print,'Reading file: ',file
  struct = read_idlstruct(file, status=status)
  return,struct
END 



FUNCTION sdss_shapelet::png_file, str, clr, createdir=createdir, nodir=nodir
  on_error, 2
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: file=ssh->png_file(stuct, clr, /createdir)'
      print
      message,'Halting'
  ENDIF 
  dir = self->image_dir(str.stripe, clr, createdir=createdir)
  idString = self->idstring(str)
  cstr = !colors[clr]

  nstr = strn(str.nmax, len=2, padchar='0')
  name = $
    'shapeletDecomp-'+cstr+'-'+str.dtype+'-'+'nmax'+nstr+'-' + idString + '.png'
  
  IF NOT keyword_set(nodir) THEN BEGIN 
      file = concat_dir(dir, name)
      return,file
  ENDIF ELSE BEGIN 
      return,name
  ENDELSE 

END 















PRO sdss_shapelet::scripts, set, typestr

  on_error, 2
  dir = '/net/cheops2/home/esheldon/idlscripts/shapelet/'

  setstr = ntostr(set)
  CASE set OF
      1: BEGIN 
          ;; Already did 9
;          stripes = [9,10,11,12,13,14,15,16]
          stripes = [10,11,12,13,14,15,16]
      END 
      2: BEGIN 
          stripes=[26,27,28,29,30,31,32,33]
      END 
      3: BEGIN 
          stripes=[34,35,36,37,42,43,44,76,82,86]
      END 
      ELSE: message,'Unknown set #: '+setstr
  ENDCASE 

  file = dir + 'run_shapelet_decomp_'+typestr+'_'+setstr+'.sh'
  print,'Writing to script: ',file
  openw, lun, file, /get_lun

  printf, lun, '!/bin/sh'
  printf, lun

  ;; Note order or importance!
  colors = [1,2,3,4,0]

  nstripe = n_elements(stripes)

  FOR ci=0,4 DO BEGIN 
      clr = colors[ci]
      
      FOR i=0L, nstripe-1 DO BEGIN 
          
          sstr = ntostr(stripes[i])
          
          printf, lun, 'nice -19 idl<<EOF'
          printf, lun, '  stripe='+sstr
          printf, lun, '  type="'+typestr+'"'
          printf, lun, '  nmax=15'
          printf, lun, '  clr='+ntostr(clr)
          printf, lun, '  run_shapelet_decomp, stripe, type, nmax, clr, $'
          printf, lun, '          /doplot, /zbuff, /png, status=status'
          printf, lun, '  if status ne 0 then exit, status=45'
          printf, lun, 'EOF'
          printf, lun, 'status=$?'
          printf, lun, 'if [ $status -ne 0 ]'
          printf, lun, 'then'
          printf, lun, '    echo Error in `basename $0` for stripe='+sstr+' clr='+!colors[clr]
          printf, lun, '    err="Shapelet Error for stripe='+sstr+' clr='+!colors[clr]+'"'
          printf, lun, '    echo "$err" | mail esheldon@cfcp.uchicago.edu -s "$err"'
          printf, lun, 'fi'
          printf, lun

          printf,lun,'dt=`date`'
          printf,lun,'message="Finished shapelets stripe='+sstr+' clr='+!colors[clr]+'"'
          printf,lun,'echo "$message:  $dt" | mail esheldon@cfcp.uchicago.edu -s "$message"'


      ENDFOR 
  ENDFOR 



  free_lun, lun

  spawn,['chmod','755',file],/noshell

END 


PRO sdss_shapelet::st2fits, pattern=pattern

  dir = '~/shapelet_outputs/'
  IF n_elements(pattern) EQ 0 THEN BEGIN 
      files = $
        file_search(dir+'stripe[0-9][0-9]-vagc-shapelets-*.st', count=count)
  ENDIF ELSE BEGIN
      files = $
        file_search(dir+pattern, count=count)
  ENDELSE 

  FOR i=0L, count-1 DO BEGIN 
      print,files[i]
  ENDFOR 

END 


;; Check which ones have odd values of negative flux in the reconstruction
;; box or low flux within a scale length
PRO sdss_shapelet::flux_tests, struct, band, flux, flux_neg, flux_pos, flux_scale, $
                 w=w, view=view, pngdir=pngdir


  IF n_elements(pngdir) NE 0 THEN BEGIN 
      setupplot,'Z'
      device, set_resolution=[900,900]
      
      loadct, 0
      
      c_colors = [!p.color,!p.color,!p.color,100,100]
  ENDIF 

  num = n_elements(struct)
  flux = replicate(-9999d, num)
  flux_neg = flux
  flux_pos = flux
  flux_scale = flux
  FOR i=0L, num-1 DO BEGIN 

      coeff = reform( struct[i].coeff[band, *] )
      scale = struct[i].scale[band]
      imsize = reform( struct[i].size[band, *] )
      center = reform( struct[i].center[band, *] )
      nmax = struct[i].nmax[band]
      

      recon = self->shapelet::recon(coeff, scale, imsize, nmax, center=center,$
                                    status=status)
      IF status EQ 0 THEN BEGIN 

          flux[i] = total(recon)


          wneg = where(recon LT 0.0, nneg, comp=wpos, ncomp=npos)
          IF nneg GT 0 THEN BEGIN 
              flux_neg[i] = abs( total(recon[wneg],/double) )
          ENDIF 
          IF npos GT 0 THEN BEGIN 
              flux_pos[i] = total(recon[wpos],/double)
          ENDIF 

          ;; create a "radius" array.
          nx = long(imsize[0])
          ny = imsize[1]

          index = lindgen(nx*ny)
          x = index MOD nx
          y = index/nx

          x = x-center[0]
          y = y-center[1]

          r = sqrt( x^2 + y^2 )
          ws = where(r LE scale,  nws)

          IF nws NE 0 THEN BEGIN 
              flux_scale[i] = total( recon[ws] )
          ENDIF 

          IF keyword_set(view) THEN BEGIN 

              frat = ntostr(flux_neg[i]/flux_pos[i])

              IF frat LT 0.01 THEN BEGIN 
                  sfrat = 'fneg/fpos < 0.01'
              ENDIF ELSE BEGIN 
                  sfrat = 'fneg/fpos = '+string(frat, format='(F5.2)')
              ENDELSE 
              sfscale = 'fscale = '+string(flux_scale[i], format='(F5.2)')
              smag = !colors[band]+' = '+string(struct[i].cmodel_counts[band], format='(F5.2)')
              slev = 'levels = '+!csym.plusminus+'[3,5,10,20,50]'+!csym.times+!csym.sigma+'!Dsky!N'

              self->view_recon, struct[i], band, leg=[slev,smag, sfrat, sfscale]


              IF n_elements(pngdir) NE 0 THEN BEGIN 
                  pngfile=self->png_file(struct[i], band, /nodir)
                  pngfile=concat_dir(pngdir, pngfile)
                  print,pngfile
                  write_png, pngfile, tvrd()
              ENDIF 

              
              key = prompt_kbrd('hit a key (q to quit)')
              IF strlowcase(key) EQ 'q' THEN return

          ENDIF 


      ENDIF 
  ENDFOR 


  setupplot,'X'

END 
PRO sdss_shapelet::test_bright, st, type, nmax, band

  IF n_elements(st) EQ 0 THEN BEGIN 
      st = self->combined_stripe_read(34, type, nmax)
  ENDIF 
  
  IF band NE 2 THEN message,'Only support r-band for now'

;  mlim = 15.42, -9999.0

  ;; This gives 100 for stripe 34
;  mlow = 15.42
;  mhigh = 15.694

  mlow = 15.70
  mhigh = 15.85

  ;; choose from randomized sample
  nst = n_elements(st)
  s = sort( randomu(seed, nst) )

  w = where(st[s].shapelet_flags[band] EQ 0 AND $
            st[s].cmodel_counts[band] GT mlow AND $
            st[s].cmodel_counts[band] LT mhigh)

  w = s[w]

  pngdir = '~/shapelet_outputs/test_bright/'+!colors[band]+'_15_70_15_85'
  IF NOT fexist(pngdir) THEN file_mkdir, pngdir
  self->flux_tests, st[w], band, /view, pngdir=pngdir

END 




;; For reconstructions that are de-convolved, you can add the
;; /psf keywords to re-convolve the reconstruction.
PRO sdss_shapelet::view_recon, struct, band, slack=slack, leg=leg, psf=psf

  nst = n_elements(struct)
  nband = n_elements(band)
  IF nst EQ 0 OR nband EQ 0 THEN BEGIN 
      print,'-Syntax: ssh->view_recon, struct, bandpass, /psf'
      print
      message,'Halting'
  ENDIF 

  IF n_elements(slack) EQ 0 THEN slack = 2

  IF nst GT 1 THEN message,'Send only a single struct'

  psfstruct = !sdss->read('psfield', struct.run, struct.camcol, $
                          field=struct.field)


  read_atlas, struct, clr=band, image=atlas, col0=col0, row0=row0, /silent
  skysig = psfstruct.skysig[band]

  ;; These are the parameters currently used in atlas_decomp
  nskysig = 2
  atlas = self->shapelet::trim_image(atlas, 0.0, skysig, nskysig, $
                                     smoothing_window=4, $
                                     slack=slack, minx=minx, maxx=maxx)
                                     

  ;; Assume this the combined decomp file
  sc = size(struct.coeff)
  IF sc[0] EQ 2 THEN BEGIN 
      coeff = reform( struct.coeff[band, *] )
      scale = struct.scale[band]
      imsize = reform( struct.size[band, *] )
      center = reform( struct.center[band, *] )
      nmax = struct.nmax[band]
      flux = struct.flux[band]
  ENDIF ELSE BEGIN 
      coeff = struct.coeff
      scale = struct.scale
      imsize = struct.size
      center = struct.center
      nmax = struct.nmax
      flux = struct.flux
  ENDELSE 

  recon = self->shapelet::recon(coeff, scale, imsize, nmax, center=center)

  IF keyword_set(psf) THEN BEGIN 
      psp = !sdss->psfield_read(struct.run, struct.camcol, struct.field)
      psf_rec, psp, struct.rowc, struct.colc, psfp

      psfimage = *psfp[band]
      psfimage = psfimage/total(psfimage)
      ptr_free, psp, psfp

      recsize = size(recon, /dim)
      psfsize = size(psfimage, /dim)
      psfcen = (psfsize-1.0)/2.0      

      ;; Only use within 3 scale lengths
      psf_fwhm = psfstruct.psf_width[band]/0.4 ; pixels
      psf_sigma = 0.42466091*psf_fwhm

      ;; Trim the psf to 4 sigma
      minx = psfcen[0] - 4.0*psf_sigma > 0
      maxx = psfcen[0] + 4.0*psf_sigma < (psfsize[0]-1)

      miny = psfcen[1] - 4.0*psf_sigma > 0
      maxy = psfcen[1] + 4.0*psf_sigma < (psfsize[1]-1)

      psfimage = psfimage[ minx:maxx, miny:maxy ]

      psfsize = size(psfimage, /dim)

;      makegauss, psfimage, psfsize, fwhm=psf_fwhm


      print,'Seeing = ',psfstruct.psf_width[band]
      print,'FWHM pixels = ', psf_fwhm
      print,'sigma = ',psf_sigma
      print,'total(psfimage) = ',total(psfimage)
      help,psfimage

      IF min(psfsize) GT min(recsize) THEN BEGIN 
          ;; assuming center is at psfsize[0]/2, psfsize[1]/2
          smin = min(recsize)/2.0 - 1.0 > 1.0
          psfcen = psfsize/2.0
          minx = psfcen[0]-smin > 0
          maxx = psfcen[0]+smin < (psfsize[0]-1)
          miny = psfcen[1]-smin > 0
          maxy = psfcen[1]+smin < (psfsize[1]-1)

          psfimage = psfimage[minx:maxx, miny:maxy]

      ENDIF 

      psfimage = psfimage/total(psfimage)
      recon = convol(temporary(recon), psfimage, edge_truncate=1, edge_wrap=0)


  ENDIF 

  recon = recon*flux
  nrecon = imsize[0]*imsize[1]
  recon[*] = recon[*] + skysig*randomu(seed, nrecon, /normal)


  IF !d.name EQ 'X' THEN BEGIN 
      c_colors = [!white, !green, !red, !blue, !darkGreen]
  ENDIF ELSE BEGIN 
      loadct, 0
      c_colors = [!p.color,!p.color,!p.color,100,100]
  ENDELSE 

;  siglevels = [3,5,10,20,50]
;  levels = skysig*siglevels

  ;; Also include negative levels
  siglevels   =    [-50,-20,-10,-5,-3, 3, 5, 10, 20, 50]
  level_lines =    [  2,  2,  2, 2, 2, 0, 0,  0,  0,  0]

  IF !d.n_colors EQ 256 THEN ncolor=100 ELSE ncolor=!grey50
;  level_colors = [replicate(!p.color,5), replicate(ncolor,5)]
  level_colors = [replicate(ncolor,5),replicate(!p.color,5)]
  levels = skysig*siglevels


  !p.multi = [0,2,2]

  acenx = struct.colc[band] - col0[band] - 0.5
  aceny = struct.rowc[band] - row0[band] - 0.5


  image_contour, atlas, $
    levels=levels, c_linestyle=level_lines, c_color=level_colors
  legend, !colors[band]+'-band',/right, box=0, charsize=1, /clear

  image_contour, recon, $
    levels=levels, c_linestyle=level_lines, c_color=level_colors
  legend, 'recon', /right, box=0, charsize=1, /clear

  imdiff = atlas
  imdiff[*] = atlas[*] - recon[*]


  ;; multiply by sqrt(2) due to subtraction
  levels = sqrt(2.0)*skysig*siglevels  

  image_contour, imdiff, $
    levels=levels, c_linestyle=level_lines, c_color=level_colors
  legend,'image-recon',/right,box=0,charsize=1,/clear

  IF n_elements(leg) NE 0 THEN BEGIN 
      plot,[0],/nodata, xstyle=4, ystyle=4
      legend, leg, box=0, charsize=1
  ENDIF 


  !p.multi = 0
END 

;; run the atlas decomposition code and view the result.
PRO sdss_shapelet::view_atlas_recon, struct, dtype, nmax, band, $
                 slack=slack, leg=leg, correct_psf=correct_psf

  num=n_elements(struct)
  FOR i=0L, num-1 DO BEGIN 
      print,'Running atlas_decomp'
      decomp = self->atlas_decomp(struct[i], dtype, nmax, band, $
                                  /trim_atlas, correct_psf=correct_psf, $
                                  slack=slack)
      print,'Running view_recon'
      self->view_recon, decomp, band, leg=leg, slack=slack

      IF i LT (num-1) THEN BEGIN 
          key = prompt_kbrd('hit a key (q to quit)')
          IF key EQ 'q' THEN return 
      ENDIF 

  ENDFOR 
END 



PRO sdss_shapelet::compare_methods, struct, slack=slack

  nmax = 15
  band = 2

  window,0
  window,1
;  window,2
  n = n_elements(struct)
  FOR i=0L, n-1 DO BEGIN 

      print
      print,'-------------------------------------------------------------'
      print
      wset, 0
      
      self->view_atlas_recon, struct[i], 'pcr', nmax, band, leg='PCR', $
        slack=slack

      wset, 1

      self->view_atlas_recon, struct[i], 'approx', nmax, band, leg='APPROX',$
        slack=slack
      
;      wset, 2

;      self->view_atlas_recon, struct[i], 'nss', nmax, band, leg='NSS'


      key = prompt_kbrd('hit a key (q to quit)')
      IF key EQ 'q' THEN return 

  ENDFOR 

END 


PRO sdss_shapelet::compare_slack, struct, slack1, slack2

  nmax = 15
  band = 2

  window,0
  window,1

  s1 = 'slack = '+ntostr(slack1)
  s2 = 'slack = '+ntostr(slack2)

;  window,2
  n = n_elements(struct)
  FOR i=0L, n-1 DO BEGIN 

      print
      print,'-------------------------------------------------------------'
      print
      wset, 0
      
      self->view_atlas_recon, struct[i], 'pcr', nmax, band, leg=s1, $
        slack=slack1

      wset, 1

      self->view_atlas_recon, struct[i], 'pcr', nmax, band, leg=s2,$
        slack=slack2
      
;      wset, 2

;      self->view_atlas_recon, struct[i], 'nss', nmax, band, leg='NSS'


      key = prompt_kbrd('hit a key (q to quit)')
      IF key EQ 'q' THEN return 

  ENDFOR 

END 


PRO sdss_shapelet::compare_nmax, struct, nmax1, nmax2

  type = 'pcr'
  band = 2

  window,0
  window,1

  s1 = 'snmax = '+ntostr(nmax1)
  s2 = 'nmax = '+ntostr(nmax2)

;  window,2
  n = n_elements(struct)
  FOR i=0L, n-1 DO BEGIN 

      print
      print,'-------------------------------------------------------------'
      print
      wset, 0
      
      self->view_atlas_recon, struct[i], type, nmax1, band, leg=s1, $
        slack=10

      wset, 1

      self->view_atlas_recon, struct[i], type, nmax1, band, leg=s2,$
        slack=10
      
;      wset, 2

;      self->view_atlas_recon, struct[i], 'nss', nmax, band, leg='NSS'


      key = prompt_kbrd('hit a key (q to quit)')
      IF key EQ 'q' THEN return 

  ENDFOR 

END 

















PRO sdss_shapelet::open_png, struct, clr

  ;; open the png files output by the shapelet
  ;; wrapper code
  nstruct = n_elements(struct)
  message = "Hit a key: space: next  p: previous  q: quit"
  i=0
  WHILE i LE nstruct-1 DO BEGIN 

      file = self->png_file(struct[i], clr)
    
      im=read_image(file)

      ;; display very slow over network
      erase
      tv,im

      IF nstruct GT 1 THEN BEGIN 

          key=prompt_kbrd(message)
          
          key = strlowcase(key)
          CASE key OF
              'q': return
              'p': i= (i-1) > 0
              ELSE: i=i+1
          ENDCASE 
      ENDIF ELSE i=i+1
  ENDWHILE 


END 

FUNCTION sdss_shapelet::array3_intersect, arr1, arr2, arr3

  intersect12 = array_intersect(arr1, arr2)

  ;; Now of this and i
  intersect123 = array_intersect(intersect12, arr3)
  return, intersect123

END 
PRO sdss_shapelet::match3, arr1, arr2, arr3, m1, m2, m3
  intersect123 = self->array3_intersect(arr1, arr2, arr3)
  match, intersect123, arr1, mi, m1, /sort
  match, intersect123, arr2, mi, m2, /sort
  match, intersect123, arr3, mi, m3, /sort
END 

;; flags for shapelet measurements
FUNCTION sdss_shapelet::flagbit, flagname
  on_error, 2
  CASE strlowcase(flagname) OF 
      'no_measurement': return,0
      'bad_c_p': return,1
      'bad_ortho': return,2
      'bad_flux': return,3
      'neg_flux': return,4
      'scale_flux': return,5
      ELSE: message,'Unknown flag name: '+flagname
  END 
END 
FUNCTION sdss_shapelet::flagval, flagname
  bit = self->flagbit(flagname)
  return, 2^bit
END 

FUNCTION sdss_shapelet::combine_struct, st, num=num

  ;; Now make a new structure with everything combined
  ;; Not doing u,z for now, but still make arrays length 5
  tags = tag_names(st)
  w = where(tags NE 'CLR', ntags)

  ncoeff = st[0].ncoeff
  psf_ncoeff = st[0].psf_ncoeff

  ival = -9999
  iarr = replicate(ival, 5)
  icarr = replicate(ival, 5, ncoeff)
  ipcarr = replicate(ival, 5, psf_ncoeff)
  i2arr = replicate(ival, 5, 2)

  lval = -9999L
  larr = replicate(lval, 5)
  lcarr = replicate(lval, 5, ncoeff)
  lpcarr = replicate(lval, 5, psf_ncoeff)
  l2arr = replicate(lval, 5, 2)



  fval = -9999.0
  farr = replicate(fval, 5)

  dval = -9999d
  darr = replicate(dval, 5)
  dcarr = replicate(dval, 5, ncoeff)
  dpcarr = replicate(dval, 5, psf_ncoeff)
  d2arr = replicate(dval, 5, 2)

  ;; build up structure definition, making
  ;; all measured size 3
  FOR ii=0L, ntags-1 DO BEGIN 

      i = w[ii]

      tag = strlowcase( tags[i] )

      CASE tag OF
          ;; The outputs of decom.  Now arrays of 5
          ;; image
          'scale': addval = darr
          'nmax': addval = larr
          'ncoeff': addval = iarr
          'coeff': addval = dcarr
          'n1': addval = icarr
          'n2': addval = icarr
          'flux': addval = darr
          'center': addval = d2arr
          'size': addval = l2arr
          'c_p': addval = darr
          ;; psf
          'psf_nmax': addval = larr
          'psf_ncoeff': addval = iarr
          'psf_scale': addval = darr
          'psf_coeff': addval = dpcarr
          'psf_n1': addval = ipcarr
          'psf_n2': addval = ipcarr
          ;; Just in pcr
          'npc': addval = iarr
          ;; Just in approx
          'ortho': addval = darr
          ;; Just copied straight from catalog as is
          ELSE: addval = st[0].(i)
      ENDCASE 

      IF n_elements(newst) EQ 0 THEN BEGIN 
          newst = create_struct(tag, addval)
      ENDIF ELSE BEGIN 
          newst = create_struct(newst, tag, addval)
      ENDELSE 
      

  ENDFOR 

  ;; Add flags
  flagarr = replicate(self->flagval('no_measurement'), 5)
  newst = create_struct(newst, 'shapelet_flags', flagarr)

  ;; unmatched tags will go at the end
  otags = $
    ['vagc_index',$
     'run','rerun','camcol','field','id','stripe',$
     'dtype',$
     'center','size', $
     'scale','nmax','ncoeff','coeff','n1','n2','flux',$
     'c_p', 'npc', $
     'psf_scale', 'psf_nmax', 'psf_ncoeff', 'psf_coeff', 'psf_n1', 'psf_n2',$
     'shapelet_flags', $
     'ra', 'dec', 'z', $
     'cmodel_counts', $
     'rowc', 'colc', $
     'm_e1', 'm_e2', 'm_rr_cc', $
     'm_e1_psf', 'm_e2_psf', 'm_rr_cc_psf']
  newst = reorder_tags(newst, otags)

  IF n_elements(num) NE 0 THEN newst = replicate(newst, num)
  return, newst

END 



PRO sdss_shapelet::copy_bandpass, inst, outst, ind, clr

  ;; Now make a new structure with everything combined
  ;; Not doing u,z for now, but still make arrays length 5
  intags = tag_names(inst)
  outtags = tag_names(outst)

  match, intags, outtags, m_in, m_out
  ntags = n_elements(m_in)


  ninst = n_elements(inst)
  ncoeff = inst[0].ncoeff
  psf_ncoeff = inst[0].psf_ncoeff

  ;; build up structure definition, making
  ;; all measured size 3
  FOR ii=0L, ntags-1 DO BEGIN 

      i_in  = m_in[ii]
      i_out = m_out[ii]

      tag = strlowcase( intags[i_in] )

      CASE tag OF
          ;; The outputs of decom.  Now arrays of 5
          ;; image
          'scale': outst[ind].scale[clr] = inst.scale
          'nmax': outst[ind].nmax[clr] = inst.nmax
          'ncoeff': outst[ind].ncoeff[clr] = inst.ncoeff
          'coeff': outst[ind].coeff[clr,*] = $
            reform(inst.coeff, 1, ncoeff, ninst)
          'n1': outst[ind].n1[clr,*] = reform(inst.n1, 1, ncoeff, ninst)
          'n2': outst[ind].n2[clr,*] = reform(inst.n2, 1, ncoeff, ninst)
          'flux': outst[ind].flux[clr] = inst.flux
          'center': outst[ind].center[clr,*] = reform(inst.center, 1, 2, ninst)
          'size': outst[ind].size[clr,*] = reform(inst.size, 1, 2, ninst)
          'c_p': outst[ind].c_p[clr] = inst.c_p
          ;; psf
          'psf_nmax': outst[ind].psf_nmax[clr] = inst.psf_nmax
          'psf_ncoeff': outst[ind].psf_ncoeff[clr] = inst.psf_ncoeff
          'psf_scale': outst[ind].psf_scale[clr] = inst.psf_scale
          'psf_coeff': outst[ind].psf_coeff[clr,*] = $
            reform(inst.psf_coeff, 1, psf_ncoeff, ninst)
          'psf_n1': outst[ind].psf_n1[clr,*] = $
            reform(inst.psf_n1, 1, psf_ncoeff, ninst)
          'psf_n2': outst[ind].psf_n1[clr,*] = $
            reform(inst.psf_n1, 1, psf_ncoeff, ninst)
          ;; Just in pcr
          'npc': outst[ind].npc[clr] = inst.npc
          ;; Just in approx
          'ortho': outst[ind].ortho[clr] = inst.ortho
          ;; Just copied straight from catalog as is
          ELSE: outst[ind].(i_out) = inst.(i_in)
      ENDCASE 

  ENDFOR 

END 


; First time only copied in scalars for some things
PRO sdss_shapelet::fix_outputs, clr

  on_error, 2

  st = self->output_read(9, clr, 'pcr', 15)

  file = self->output_file(9, clr, 'pcr', 15)
  newfile = file + '.new'
  

  newst = rename_tags(st, ['coef','psf_coef'],['coeff','psf_coeff'])

  print,'Writing to new file: ',newfile
  write_idlstruct, newst, newfile





  return
  query = $
    'select run,rerun,camcol,field,id,stripe,'+$
                'ra, dec, rowc, colc, m_rr_cc, m_rr_cc_psf, m_e1, m_e2,'+$
                'm_e1_psf, m_e2_psf, z, cmodel_counts, absmodelmag, '+$
                'match_rerun, match_id, flags, primtarget '+$
    'from specgal_main where stripe = 9'
  print,query

  str = pgsql_query(query)


  st = self->output_read(9, clr, 'pcr', 15)
  
  tags = ['cmodel_counts', 'm_e1', 'm_e2', 'm_e1_psf', 'm_e2_psf']
  desc = replicate('fltarr(5)', n_elements(tags))
  tst = remove_tags(st, tags)
  add_tags, temporary(tst), tags, desc, newst
  
  vagc_index = newst.vagc_index
  newst.cmodel_counts = str[vagc_index].cmodel_counts
  newst.m_e1 = str[vagc_index].m_e1
  newst.m_e2 = str[vagc_index].m_e2
  newst.m_e1_psf = str[vagc_index].m_e1_psf
  newst.m_e2_psf = str[vagc_index].m_e2_psf

  file = self->output_file(9, clr, 'pcr', 15)
  newfile = file + '.new'
  
  IF fexist(newfile) THEN message,'File exists: '+newfile
  print,'Writing file: ',newfile
  write_idlstruct, newst, newfile


END 


PRO sdss_shapelet::rem_dup_byid, gst, rst, ist

  gphotoid = sdss_photoid(gst)
  rphotoid = sdss_photoid(rst)
  iphotoid = sdss_photoid(ist)

  rmd = rem_dup(gphotoid)
  help,gst,rmd
  gst = gst[rmd]

  rmd = rem_dup(rphotoid)
  rst = rst[rmd]

  rmd = rem_dup(iphotoid)
  ist = ist[rmd]



END 

FUNCTION sdss_shapelet::flux_limits, band
  CASE band OF
      0: return, [100.0, 2.e4]
      1: return, [4000.0, 2.e5]
      2: return, [4000.0, 2.e5]
      3: return, [4000.0, 2.e5]
      4: return, [100.0, 6.e4]
      ELSE: message,'Band band: '+ntostr(band)
  ENDCASE 
END 
FUNCTION sdss_shapelet::c_p_max
  return, 1.0
END 
PRO sdss_shapelet::setflags, st

  on_error, 2
  IF n_elements(st) EQ 0 THEN BEGIN 
      print,'-Syntax: sh->setflags, struct'
      print
      message,'Halting'
  ENDIF 

  IF st[0].dtype NE 'pcr' THEN BEGIN 
      message,'No support for dtype = "'+st[0].dtype+'" yet'
  ENDIF 

  c_p_max = self->c_p_max()
  c_p_flag = self->flagval('bad_c_p')
  flux_flag = self->flagval('bad_flux')
  neg_flux_flag = self->flagval('neg_flux')
  scale_flux_flag = self->flagval('scale_flux')

  nomeas_flag = self->flagval('no_measurement')
  FOR band=0,4 DO BEGIN 
      
      ;; If no measurement, don't set these other flags
      w = where( (st.shapelet_flags[band] AND nomeas_flag) EQ 0, nw)

      IF nw NE 0 THEN BEGIN 

          fluxlim = self->flux_limits(band)
          
          bad = where(st[w].flux[band] LT fluxlim[0] OR $
                      st[w].flux[band] GT fluxlim[1], nbad)

          IF nbad NE 0 THEN BEGIN 
              bad = w[bad]
              st[bad].shapelet_flags[band] = $
                st[bad].shapelet_flags[band] + flux_flag
          ENDIF 

          IF st[0].dtype EQ 'pcr' THEN BEGIN 
              bad = where(st[w].c_p[band] GT c_p_max, nbad)
          
              IF nbad NE 0 THEN BEGIN 
                  bad = w[bad]
                  st[bad].shapelet_flags[band] = $
                    st[bad].shapelet_flags[band] + c_p_flag
              ENDIF 
          ENDIF 

          ;; Do Huan's flux tests
          print,'Doing flux tests in '+!colors[band]
          self->flux_tests, st[w], band, f, fn, fp, fs
          wg = where(fn NE -9999 AND fp GT 0.0, ng)
          IF ng NE 0 THEN BEGIN 
              ratio = fn[wg]/fp[wg]
              wbad = where(ratio GT 0.1, nbad)
              IF nbad NE 0 THEN BEGIN 
                  wbad = w[wg[wbad]]
                  st[wbad].shapelet_flags[band] = $
                    st[wbad].shapelet_flags[band] + neg_flux_flag
              ENDIF 

              wbad = where(fs LT 0.1, nbad)
              IF nbad NE 0 THEN BEGIN 
                  wbad = w[wg[wbad]]
                  st[wbad].shapelet_flags[band] = $
                    st[wbad].shapelet_flags[band] + scale_flux_flag
              ENDIF 
              
          ENDIF 
          

      ENDIF 
  ENDFOR 


END 




PRO sdss_shapelet::combine_bands, stripe, dtype, nmax, overwrite=overwrite

  on_error, 2
  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: sh->combine_bands, stripe, dtype, nmax, /overwrite'
      print
      message,'Halting'
  ENDIF 

  cfile = self->combined_stripe_file(stripe, dtype, nmax)
  IF fexist(cfile) AND NOT keyword_set(overwrite) THEN BEGIN  
      key=prompt_kbrd('File '+cfile+' exists, overwrite? (y/n)')
      key = strlowcase(key)
      IF key EQ 'n' THEN return
  ENDIF 

  IF strlowcase(dtype) NE 'pcr' THEN $
    message,'Do not yet support dtype = '+strlowcase(dtype)


  FOR band=0,4 DO BEGIN 

      clr = !colors[band]

      st = self->output_read(stripe, band, dtype, nmax, status=status)

      IF status NE 0 THEN BEGIN 
          ;; No file for this bandpass
          print, 'Outputs for bandpass '+clr+' does not exist'
          key = prompt_kbrd('Continue? (y/n)')
          IF strlowcase(key) EQ 'n' THEN return
      ENDIF ELSE BEGIN 

          vindex = st.vagc_index
          tptr = ptr_new(st, /no_copy)

          add_arrval, band, usebands
          add_arrval, tptr, ptrlist
          add_arrval, vindex, allid

      ENDELSE 
  ENDFOR 

  nuse = n_elements(usebands)
  IF nuse EQ 0 THEN message,'No files found'

  allid = allid[ rem_dup(allid) ]
  ntot = n_elements(allid)

  tst = (*ptrlist[0])[0]
  outstruct = self->combine_struct(tst, num=ntot)

  print,'copying bandpasses into output struct'
  FOR i=0L, nuse-1 DO BEGIN 
      st = *ptrlist[i]
      band = usebands[i]

      vindex = st.vagc_index
      match, vindex, allid, mst, mall, /sort
      self->copy_bandpass, st[mst], outstruct, mall, band
      outstruct[mall].shapelet_flags[band] = 0
  ENDFOR 

  ;; Now find bad measurements
  print,'Setting flags'
  self->setflags, outstruct

  ptr_free, ptrlist

  print,'Writing to .st file: ',cfile
  write_idlstruct, outstruct, cfile
  fitfile = repstr(cfile,'.st','.fit')
  print,'Writing to .fit file: ',fitfile
  mwrfits, outstruct, fitfile, /create


END 

FUNCTION sdss_shapelet::stripes, nstripes
  stripes = $
    [9,10,11,12,13,14,15,16,$
     26,27,28,29,30,31,32,33,$
     34,35,36,37,42,43,44,76,82,86]
  nstripes = n_elements(stripes)
  return,stripes
END 
PRO sdss_shapelet::run_combine_bands, dtype, nmax, overwrite=overwrite
  stripes = self->stripes(nst)
  FOR i=0L, nst-1 DO BEGIN 
      stripe = stripes[i]
      self->combine_bands, stripe, dtype, nmax, overwrite=overwrite
      print,'----------------------------------------------------------'
  ENDFOR 
END 
PRO sdss_shapelet::combine_stripes, dtype, nmax
  outfile = self->combined_file(dtype,nmax)
  print
  print,'Will write to file: ',outfile
  IF fexist(outfile) THEN BEGIN 
      key = ''
      read,'File exists, overwrite? (y/n) ',key
      IF strlowcase(key) NE 'y' THEN return
  ENDIF 
  stripes = self->stripes()
  struct = self->combined_stripe_read(stripes,dtype,nmax)
  
  print
  print,'Writing file: ',outfile
  write_idlstruct, struct, outfile

  fitsfile = repstr(outfile, '.st', '.fit')
  print
  print,'Writing fits file: ',fitsfile
  mwrfits2, struct, fitsfile, /create, /destroy


END 


PRO sdss_shapelet::compare_gri, stripe, dtype, nmax, gst, rst, ist

  IF n_elements(gst) EQ 0 THEN BEGIN 
      gst = self->output_read(stripe, 1, dtype, nmax)
      rst = self->output_read(stripe, 2, dtype, nmax)
      ist = self->output_read(stripe, 3, dtype, nmax)
  ENDIF 


  ggood = where(gst.c_p LT 1)
  rgood = where(rst.c_p LT 1)
  igood = where(ist.c_p LT 1)

  gphotoid = sdss_photoid(gst[ggood])
  rphotoid = sdss_photoid(rst[rgood])
  iphotoid = sdss_photoid(ist[igood])


  self->match3, gphotoid, rphotoid, iphotoid, mg, mr, mi

  mg = ggood[mg]
  mr = rgood[mr]
  mi = igood[mi]

  help,mg,mr,mi

  nmatch = n_elements(mg)
  FOR i=0L, nmatch-1 DO BEGIN 

      print,'Displaying g'
      self->open_png, gst[mg[i]], 1
      key = prompt_kbrd('hit a key')
      IF key EQ 'q' THEN return


      print,'Displaying r'
      self->open_png, rst[mr[i]], 2
      key = prompt_kbrd('hit a key')
      IF key EQ 'q' THEN return


      print,'Displaying i'
      self->open_png, ist[mi[i]], 3
      key = prompt_kbrd('hit a key')
      IF key EQ 'q' THEN return

  ENDFOR 

END 














PRO sdss_shapelet::compare_bckelly_massey_plot, original, recon_bckelly, recon_massey, skysig, idstring, dtype

  !p.multi=[0,2,2]

  nrk=n_elements(recon_bckelly)
  recon_bckelly[*] = recon_bckelly[*] + skysig*randomu(seed,nrk,/normal)

  nrm=n_elements(recon_massey)
  recon_massey[*] = recon_massey[*] + skysig*randomu(seed,nrm,/normal)


  siglevels = [3,5,10,20,50]
  c_lines=[0,0,0,0,0]

  IF (!d.name EQ 'Z') OR (!d.name EQ 'PS') THEN BEGIN 

      IF !d.name EQ 'Z' THEN BEGIN 
          setupplot,'Z'
          device, set_resolution=[770,890]
      ENDIF 
      loadct, 0
      c_colors = [!p.color,!p.color,!p.color,100,100]

  ENDIF ELSE IF !d.name EQ 'X' THEN BEGIN 
      setupplot,'X'
      c_colors = [!white, !green, !red, !blue, !darkGreen]
  ENDIF 

  levels = siglevels*skysig

  image_contour, original, levels=levels, c_colors=c_colors, $
    title=idstring, sky=0.0

  image_contour, recon_bckelly, levels=levels, c_colors=c_colors, $
    title='Kelly Reconstruction "'+dtype+'"', sky=0.0
  image_contour, recon_massey, levels=levels, c_colors=c_colors, $
    title='Massey Reconstruction "ls"', sky=0.0

  legend,ntostr(siglevels)+' '+!csym.sigma,$
    /right,box=0,charsize=0.7, lines=c_lines, colors=c_colors

  

END 


PRO sdss_shapelet::compare_bckelly_massey, str, psFieldInfo=psFieldInfo

  on_error, 2
  IF n_elements(str) EQ 0 THEN BEGIN 
      file = '~/shapelet_outputs/stripe16-vagc-shapelets-approx-nmax15.fit'
      str = mrdfits(file,1)
  ENDIF 

  message,'Need to convert to new stuff'

  nmax = 15
;  dtype = 'approx'
  dtype = 'pcr'

  window,0,xsize=800,ysize=890
;  window,1,xsize=800,ysize=512

  outDir = '~/shapelet_outputs/compare_bckelly_massey/'
  bcKellyDir = outDir + 'bckelly/'
  masseyDir = outDir + 'massey/'
  compareDir = outDir + 'compare/'

  outFileFront = 'stripe16-nmax' + ntostr(nmax) + '-'

  bcKellyTime = 0d
  MasseyTime = 0d

  nstr=n_elements(str)
  

  FOR nn=0L, nstr-1 DO BEGIN 

      !p.multi=[0,2,2]
;      wset,0

      idString = shapelet_idstring(str[nn])

      BCKellyFile = $
        bcKellyDir + outFileFront + idString+'-bckelly-'+dtype+'.png'
      MasseyFile  = $
        masseyDir  + outFileFront + idString+'-massey-ls.png'
      compareFile = $
        compareDir + outFileFront + idString+'-compare-ls-'+dtype+'.png'

      read_atlas, str[nn], imr=imr, col0=col0, row0=row0

      print
      print,'Calling Brandons code'
      print,BCKellyFile

      tm = systime(1)
      decomp = $
        self->shapelet::atlas_decomp(str[nn], dtype, nmax, 2, $
                                     recon=recon, $
                                     psFieldInfo=psFieldInfo)
      tm = systime(1)-tm
      print,'BCKelly time'
      ptime,tm
      bcKellyTime = bcKellyTime + tm

      self->view_recon, imr, recon, psFieldInfo, /addnoise
      write_png, BCKellyFile, tvrd(/true)
      
      w=where(imr EQ 0, nw)

      imrcopy = imr
      imrcopy[w] = psFieldInfo.skysig*randomu(seed, nw, /normal)
      beta = sqrt(str[nn].m_rr_cc[2]/2.0)
      x0 = [str[nn].colc[2],str[nn].rowc[2]] - float([col0[2], row0[2]])

      print
      print,'Calling Masseys code'
      print,MasseyFile
      tm = systime(1)
      shapelets_decomp, imrcopy, beta, nmax, decomp2, recomp2, $
        noise=psFieldInfo.skysig, x0=x0
      tm = systime(1)-tm
      print,'Massey time'
      ptime,tm
      MasseyTime = MasseyTime + tm

      self->view_recon, imr, recomp2, psFieldInfo, /addnoise
      write_png, MasseyFile, tvrd(/true)

;      wset,1
      self->compare_bckelly_massey_plot, $
        imr, recon, recomp2, psFieldInfo.skysig, idstring, dtype
      write_png, compareFile, tvrd(/true)

      key=get_kbrd(1)

  ENDFOR 
  print,'Total BCKelly time'
  ptime,BCKellyTime
  print,'Total Massey time'
  ptime,MasseyTime

END 





















;; Need to add clr
PRO sdss_shapelet::makechris_files, dtype

  on_error, 2
  dtype = strlowcase(dtype)

  dir = '~/shapelet_outputs/'
  outfile = dir + 'vagc-shapelets-r-'+dtype+'-nmax15.fit'
  clr = 2

  inStruct = read_idlstruct_multi(file_search(dir + 'stripe*'+dtype+'*.st'))

  toutst = inStruct[0]

  outst = remove_tags(toutst, 'cmodel_counts')

  outst = create_struct($
                         outst, $
                         'cmodel_counts', fltarr(5), $
                         'cmodelflux_ivar', fltarr(5), $
                         'petrocounts', fltarr(5), $
                         'petroflux_ivar', fltarr(5), $
                         'abspetromag', fltarr(5) $
                       )

  outst = replicate(outst, n_elements(inStruct))

  struct_assign, inStruct, outst, /verbose

  ;; read extra vagc info and copy in. 
  columns=['run','match_rerun','camcol','field','match_id', $
           'cmodel_counts', 'cmodelflux_ivar', $
           'petroflux',   'petroflux_ivar', $
           'abspetromag']

  ;; Note: big file only contains a few tags, so we read in the individual
  ;;       files
  dir = sdssidl_config('shapecorr_dir')+'combined/vagc/'
  files = file_search(dir + 'stripe*matchlocal.fit')
  vagc = mrdfits_multi(files, columns=columns)


  ;; match them up
  print,'Matching'
  outid = photoid(outst.run, outst.rerun, outst.camcol, outst.field, outst.id)
  vagcid = photoid(vagc.run, vagc.match_rerun, vagc.camcol, vagc.field, vagc.match_id)

  match, outid, vagcid, mout, mvagc, /sort
  help,outst, mout, mvagc
  IF n_elements(mout) NE n_elements(outst) THEN message,'not all matched!'


  print,'Copying'
  outst[mout].cmodel_counts = vagc[mvagc].cmodel_counts
  outst[mout].cmodelflux_ivar = vagc[mvagc].cmodelflux_ivar
  outst[mout].petrocounts = 22.5 - alog10(vagc[mvagc].petroflux)
  outst[mout].petroflux_ivar = vagc[mvagc].petroflux_ivar

  outst[mout].abspetromag[0] = vagc[mvagc].abspetromag[0]
  outst[mout].abspetromag[1] = vagc[mvagc].abspetromag[1]
  outst[mout].abspetromag[2] = vagc[mvagc].abspetromag[2]
  outst[mout].abspetromag[3] = vagc[mvagc].abspetromag[3]
  outst[mout].abspetromag[4] = vagc[mvagc].abspetromag[4]

  outst[mout].vagc_index = mvagc

  print
  print,'Writing file: ',outfile
  mwrfits, outst, outfile,/create

END 

















FUNCTION sdss_shapelet::cleanup
  return,1
END 
PRO sdss_shapelet__define

  struct = {$
             sdss_shapelet, $
             sdss_shapelets_dummy: 0, $
             INHERITS shapelet, $
             INHERITS sdss_util $
           }

END 
