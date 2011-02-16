;; These will be used for lensing and intrinsic alignments


FUNCTION percolate::init, sample

  IF n_elements(sample) NE 1 THEN BEGIN 
      message,'You must send a scalar integer for sample'
  ENDIF 

  rmin = 20.0
  logbin = 1

  sample = fix(sample)

  par = $
    { $
      maskfile: '', $
      sample: sample, $
      source_sample: 1, $
      catalog: '',$
      zbuffer: 0.0 $
    }

  type = 'percolate'
  CASE sample OF 
      1: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton re-gaussianization
          par.catalog = 'Mr20'    ; lrg catalog
      END 
      2: BEGIN 
          ;; Intrinsic alignment sample
          rmax = 2800.0
          nbin = 14
          sigmacrit_style = 4   ; intrinsic alignments with redshift cut
          par.catalog = 'Mr20'
          par.source_sample = 1 ; intrinsic align sample # 1
          shape_correction_style = 4 ; means doing intrinsic align. kludge
      END 
      ELSE: message,'Unsupported sample: '+strn(sample)
  ENDCASE 

  self.sample = sample

  print,'type   =                 ',type
  print,'sample =                 ',sample,format='(A,I7)'
  print,'source_sample =          ',par.source_sample,format='(A,I7)'
  print,'rmin =                   ',rmin
  print,'rmax =                   ',rmax
  print,'nbin =                   ',nbin
  print,'sigmacrit_style =        ',sigmacrit_style,format='(A,I7)'
  print,'shape_correction_style = ',shape_correction_style,format='(A,I7)'
  print,'zbuffer =                ',par.zbuffer

  objshear_retval = $
    self->objshear::init(type, rmin, rmax, nbin, $
                         sigmacrit_style, shape_correction_style, $
                         logbin=logbin, par_struct=par)

  return, objshear_retval


END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Original catalogs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION percolate::catalog
  p = self->par_struct()
  return,p.catalog
END 

FUNCTION percolate::catdir
  return,concat_dir(esheldon_config('lensinput_dir'),'percolate/catalog')
END 
FUNCTION percolate::catfile, randnum=randnum
  IF n_elements(randnum) NE 0 THEN BEGIN 
      return,self->randfile(randnum)
  ENDIF 
  cat = self->catalog()
  CASE cat OF
      'Mr18': BEGIN 
          file = 'Mr18.groups.dat'
      END 
      'Mr19': BEGIN 
          file = 'Mr19.groups.dat'
      END 
      'Mr20': BEGIN 
          file = 'Mr20.groups.dat'
      END 
      ELSE: message,'Unknown catalog: '+strn(cat)
  END 
  dir = self->catdir()
  file = concat_dir(dir, file)
  return,file
END 
; For now only have one random catalog
FUNCTION percolate::randfile, randnum
  cat = self->catalog()
  CASE cat OF
      'Mr18': BEGIN 
          file = 'vollim3.random'
      END 
      'Mr19': BEGIN 
          file = 'vollim2.random'
      END 
      'Mr20': BEGIN 
          file = 'vollim1.random'
      END 
      ELSE: message,'Unknown catalog: '+strn(cat)
  END 
  dir = self->dir()
  file = concat_dir(dir, file)
  return,file
END 


FUNCTION percolate::get,status=status
  file = self->catfile()
  structdef = {$
                id:0L,$
                ra:0d, dec:0d, $
                z: 0.0, $
                ngals:0,$
                Mr_tot: 0.0, $
                gmr_tot: 0.0, $
                sigmav: 0.0, $
                r_perp_rms: 0.0, $
                r_edge: 0.0 $
              }
  read_struct, file, structdef, struct
  return, struct
END 


;; override objshear::lenscuts for some special
;; cases
FUNCTION percolate::angcut, lensum, ngoodz

  par = self->par_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Redshift cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  max_allowed_angle = par.max_allowed_angle
  maxz = par.maxz

  minz = 0.0
  print,'Cutting maxz = ',maxz
  print,'Cutting max_allowed_angle = '+ntostr(max_allowed_angle)+' degrees'


  nlens_init = n_elements(lensum)
  wlens = where(lensum.z GT 0.0 AND $
                lensum.z LT maxz AND $
                lensum.angmax LT max_allowed_angle, ngoodz)
  return,wlens
END 

FUNCTION percolate::lenscuts, lstruct, nkeep

  sample = self->sample()
  CASE sample OF 
      1: wlens = self->objshear::lenscuts(lstruct, nkeep)
      2: BEGIN 
          nkeep = n_elements(lstruct)
          wlens = lindgen(nkeep)
      END 
      ELSE: message,"Don't know about sample "+strn(sample)+" yet"
  ENDCASE 
  return,wlens
END 
PRO percolate::setuplens

  tm = systime(1)

  par_struct = self->objshear::par_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file name(s)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(random) THEN randnum=0
  outfile = self->lensinput_file()

  print
  print,'Setting up file: ',outfile

  ;; get the data.  Note: zindex will be the same for all samples
  ;; created from this lens catalog
  cat = self->get()
  ncat = n_elements(cat)
  zindex = lindgen(ncat)

  nlens_init = n_elements(cat)
  
  print,'Creating lens struct'
  lstruct = self->objshear::lensinput_structdef()
  lstruct = replicate(lstruct, nlens_init)
  
  print,'Copying....'
  copy_struct, cat, lstruct
  
  eq2csurvey, cat.ra, cat.dec, clam, ceta
  lstruct.clambda = clam
  lstruct.ceta = ceta
  clam = 0 & ceta = 0
  
  ;; Keep track of the objects and their redshifts
  lstruct.zindex = zindex
  
  ;; Calculate some cosmology-dependent stuff and copy into struct
  self->objshear::calc_cosmo, lstruct
  
  ;; Make some generic lens cuts.  Note, for the intrinsic alignments
  ;; we are not edge checking so use local lenscuts
  
  wlens = self->lenscuts(lstruct, nkeep)
  
  ;; remove the unwanted lenses
  lstruct = lstruct[wlens]
  
  print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'
  
  print
  print,'Writing to file: ',outfile
  write_idlstruct, lstruct, outfile

  ptime,systime(1)-tm

END 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Binning
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO percolate::lum_bins, nbin, lowlim, highlim

  IF n_params() LT 3 THEN BEGIN 
      message,'-Syntax: p->lum_bins, nbin, lowlim, highlim'
  ENDIF 

  CASE nbin[0] OF 
      2: BEGIN 
          ;; This is the median
          lowlim =  [-26.0, -22.11]
          highlim = [-22.11, -21.0]
      END 
      ELSE: BEGIN 
          message,'Unknown nbin: '+strn(nbin)
      END 
  ENDCASE 

END 

FUNCTION percolate::lum_where_string, nbin

  IF n_elements(nbin) EQ 0 THEN BEGIN 
      message,'ws = p->lum_where_string(nbin)'
  ENDIF 
  self->lum_bins, nbin, lowlim, highlim

  IF nbin EQ 1 THEN where_string = '' $
  ELSE where_string = strarr(nbin)
  
  FOR i=0L, nbin-1 DO BEGIN 
      
      where_string[i] = $
        'struct.mr_tot GE '+ntostr(lowlim[i])+' AND '+$
        'struct.mr_tot LT '+ntostr(highlim[i])
      
  ENDFOR 

  return,where_string
END 

FUNCTION percolate::nbin, subtype
  CASE subtype OF
      'lum2': return,2
      ELSE: message,'Unknown subtyep: '+strn(subtype)
  ENDCASE 
END 

FUNCTION percolate::where_string, subtype, nbin=nbin

  IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: ws = p->where_string(subtype, nbin=)'
  ENDIF 
  CASE subtype OF
      'lum2': BEGIN 
          nbin = self->nbin(subtype)
          return,self->lum_where_string(nbin)
      END 
      ELSE: message,'Unknown select subtype: '+ntostr(subtype)
  ENDCASE 

END 


PRO percolate::sub, subtype, randnum=randnum

  average_tags = $
    ['z', 'ngals', 'mr_tot', 'gmr_tot', 'sigmav', 'r_perp_rms','r_edge']

  self->objshear::sub_sample, subtype, randnum=randnum, $
    average_tags=average_tags, zhist=zhist

END 












; Plot the sigma assuming it is actually a plot of the 
; intrinsic tangential shear

PRO percolate::plot_intrinsic, cumulative=cumulative, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      dir = self->plot_dir(/createdir)
      sstr = self->sample_string()

      file = 'percolate_'+sstr+'_intrinsic'
      IF keyword_set(cumulative) THEN file = file + '_cumulative'
      file = file + '.eps'

      file = concat_dir(dir,file)
      begplot, name=file, /encapsulated
  ENDIF 

  sh = self->lensread('combined')

  ytitle='Intrinsic '+!csym.gamma+'!DT!N'


  IF keyword_set(cumulative) THEN BEGIN 
      ytitle = 'Cumulative ' + ytitle
      rad = sh.tmeanr/1000
      shear = sh.tsigma
      shearerr = sh.tsigmaerr
  ENDIF ELSE BEGIN 
      rad = sh.meanr/1000
      shear = sh.sigma
      shearerr = sh.sigmaerr
  ENDELSE 

  yrange = prange(shear, shearerr, /sym,/slack)

  pplot, rad, shear, yerr=shearerr, aspect=!gratio, $
    psym=8,/xlog, ytitle=ytitle, xtitle=!mpcxtitle2, $
    xtickf='loglabels', xticklen=0.04, $
    yrange=yrange

  oplot, [0.001, 100], [0,0]

  ;; print limits
  print,['radius (Mpc)','gamma','gammaerr','tmeanr','tgamma','tgammaerr','lower (95%)','upper (95%)'], $
    format = '( 8(A15) )'
  forprint,$
    rad, sh.sigma, sh.sigmaerr, sh.tmeanr/1000, sh.tsigma, sh.tsigmaerr, sh.tsigma-2*sh.tsigmaerr, sh.tsigma+2*sh.tsigmaerr, $
    format = '( 8(F15) )'

  IF keyword_set(dops) THEN endplot,/trim_bbox

END 

PRO percolate__define

  struct = {$
             percolate, $
             sample: 0, $
             INHERITS objshear $
           }

END 
