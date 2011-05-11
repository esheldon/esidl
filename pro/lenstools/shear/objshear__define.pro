;+
; NAME:
;  OBJSHEAR (class: objshear__define)
;
;
; PURPOSE:
;  Some general object-shear methods
;
;
; CATEGORY:
;  lensing
;
; METHODS:
;  Initialization:
;    os = obj_new('objshear', ....)
;
;    SIGMACRIT_STYLES:
;       1: treat photozs as truth
;       2: use the mean inverse critical density calculated from the 
;          deconvolved overall distribution. 
;       3: Use integral of photoz times prior.
;
;    SHAPE_CORRECTION_STYLE
;       1: old method
;       2: princeton analytic
;       3: princeton re-gaussianization
;       4: we are using the spectra for intrinsic alignments
;
; Steps to run a sample through:
;   - setuplens
;   - write_parfile
;   - write_script
;   - possibly alter scripts to use directories on another machine.
;
; General steps done after you have run a sample through:
; Usually inherit from these in e.g. maxbcg/maxbcg_lensing
;
;   - Do some sub-samples. This automatically does the "combine" stage.
;      sub_sample, subtype
;      sub_sample, subtype, randnum=
;   - combine_allrand, subtype=
;   - corr, subtype=subtype
;   - jackknife, subtype=subtype
;
; To combine two samples, say 10 and 30 Mpc:
;   - combine_samples, sample1, sample2, dtype, subtype=
;       where dtype is 'corr', 'jackknife', etc.
;
;  Plotting routines:
;   - you should inherit from this class and implement
;       ::plot_profile.
;   See maxbcg_lensing class for an example
;
;
; MODIFICATION HISTORY:
;
;-


FUNCTION objshear::init, $
	type_in, $
	rmin, rmax, nbin_OR_binsize, sigmacrit_style, shape_correction_style, $
	par_struct_input = par_struct_input, $
	logbin=logbin

  IF n_params() LT 6 THEN BEGIN 
      message,$
        'You must initialize at least type, rmin, rmax, nbin_or_binsize, sigmacrit_style, shape_correction_style: ',/inf
      message,$
        "   obj=obj_new('objshear', rmin, rmax, nbin_or_binsize, sigmacrit_style, shapee_correction_style, /logbin, par_struct=),",/inf
      return,0
  ENDIF 

  type = strlowcase(type_in)

  par_struct = $
    create_struct('type', type, $
                  'rmin', float(rmin), $
                  'rmax', float(rmax), $
                  'sigmacrit_style', sigmacrit_style, $
                  'shape_correction_style', shape_correction_style)

  IF keyword_set(logbin) THEN BEGIN 
      par_struct = $
        create_struct(par_struct, $
                      'logbin', 1, $
                      'nbin', nbin_OR_binsize, $
                      'binsize', -1.0)              
  ENDIF ELSE BEGIN 
      binsize = nbin_OR_binsize
      par_struct = $
        create_struct(par_struct, $
                      'logbin', 0, $
                      'nbin', long( (rmax - rmin)/binsize ), $
                      'binsize', binsize)              
  ENDELSE 

  IF n_elements(par_struct_input) NE 0 THEN BEGIN 
      par_struct = create_struct(par_struct, par_struct_input)
  ENDIF 

  new_par_struct = self->par_struct_copy(par_struct)

  self.par_structPtr = ptr_new(new_par_struct, /no_copy)

  return,1

END 

;; copy input struct into a default par_struct
FUNCTION objshear::par_struct_copy, input_par_struct

  ;; Set the default paramater structure
  par_struct = create_struct('type',                  '', $
                             'h',                   1.0, $
                             'omega_m',            0.27, $
                             $
                             'sigmacrit_style',       1, $ 
                             'shape_correction_style',3, $ 
                             $
                             'logbin',              1,  $ 
                             $
                             'nbin',               18,  $ 
                             $
                             'binsize',          -1.0,  $ 
                             $
                             'rmin',             20.0,  $
                             'rmax',          11500.0,  $
                             $
                             'comoving',            0,  $
                             'dopairs',             0,  $
                             $
                             'depth',              10,  $
                             'zbuffer',           0.0,  $
                             'maxz',              0.8,  $
                             'max_allowed_angle', 6.0,  $ ; degrees
                             'maskfile',           '',  $
                             'edgecuts',            1,  $
                             'sample',             -1,  $
                             'source_sample',      -1,  $
                             'random_sample',      -1,  $
                             'scinv_sample',       -1,  $
                             'catalog',            ''   $
                            )

  copy_struct, input_par_struct, par_struct
  return,par_struct

END 

; Return a copy of the parameter structure
FUNCTION objshear::par_struct

  IF ptr_valid(self.par_structPtr) THEN BEGIN 
      IF n_elements(*self.par_structPtr) NE 0 THEN BEGIN 
          return,*self.par_structPtr
      ENDIF ELSE BEGIN 
          message,'par_struct is not defined',/inf
          return,-1
      ENDELSE 
  ENDIF ELSE BEGIN 
      message,'The parameter structure is pointing at nothing',/inf
      return,-1
  ENDELSE 

END 








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Some utility programs
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION objshear::type
  par = self->par_struct()
  return,par.type
END 

FUNCTION objshear::catalog
  par = self->par_struct()
  return, par.catalog
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The mask describing the source distribution
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO objshear::maskfile, maskfile, edgecuts=edgecuts

	par = *self.par_structPtr
	if par.maskfile ne '' then begin 
		maskfile = par.maskfile
	endif else if par.shape_correction_style eq 3 $
			and not keyword_set(edgecuts) then begin
		; Edgecuts are the pure 1Mpc edgecut that Ben uses for his cluster
		; finder.
		; Catalogs are drawn from big mask, not for example princeton mask,
		; so just use default in /edgecuts case.
		maskfile = sdssidl_config('pixel_mask_princeton_basic', exist=exist)
		if not exist or not file_test(maskfile) then begin
			message,'Princeton mask not found'
		endif
	endif 

;  IF n_elements(maskfile) NE 0 THEN BEGIN 
;      print,'Using maskfile: ',maskfile
;  ENDIF 

END 


; old random string
FUNCTION objshear::randstring, randnum

  nrand = n_elements(randnum)
  IF nrand EQ 1 THEN randstr = '' $
  ELSE randstr=strarr(nrand)

  FOR i=0L, nrand-1 DO BEGIN 
      randstr[i] = strn(randnum[i],len=2,padchar='0')
  ENDFOR 

  return,randstr

END 

; Strings describing the sample or random sample
FUNCTION objshear::sample
  par = self->par_struct()
  IF tag_exist(par, 'sample') THEN BEGIN 
      return,par.sample
  ENDIF ELSE BEGIN 
      return,-1
  ENDELSE 
END 

FUNCTION objshear::sample_string, sample=sample
    if n_elements(sample) eq 0 then sample=self->sample()
    sstr = ntostr(sample, format='(I02)')
    ; multiple samples get joined by dashes
    sstr = strjoin(sstr, '-')
    return,'sample'+sstr
END 


FUNCTION objshear::random_sample
  par = self->par_struct()
  IF tag_exist(par, 'random_sample') THEN BEGIN 
      return,par.random_sample
  ENDIF ELSE BEGIN 
      return,-1
  ENDELSE 
END 

FUNCTION objshear::random_sample_string
    return, self->sample_string(self->random_sample())
END 


; String for old file names describing rmin/rmax
FUNCTION objshear::rminmax_string
  par = self->par_struct()
  rminstr = ntostr( round( par.rmin ) )
  rmaxstr = ntostr( round( par.rmax ) )
  ext = 'rmin'+rminstr+'_rmax'+rmaxstr

  return,ext
END 

; generate string based on high and low limits.  input the
; limits and the variable name for lowlim <= name <= highlim
FUNCTION objshear::highlow_string, lowlim, highlim, name

  IF n_elements(lowlim) GT 1 THEN message,'only one limit at a time is allowed'
  
  ;; This will correctly do integers and numbers
  ;; with a single value after decimal point

  lstr = ntostr(lowlim)
  hstr = ntostr(highlim)
  
  tt = strsplit(lstr, '.', /extract)
  
  IF n_elements(tt) EQ 2 THEN BEGIN 
      lstr = tt[0] + '.' + strmid(tt[1], 0, 1)
  ENDIF 
  tt = strsplit(hstr, '.', /extract)
  
  IF n_elements(tt) EQ 2 THEN BEGIN 
      hstr = tt[0] + '.' + strmid(tt[1], 0, 1)
  ENDIF 
  
  str = lstr +' <= '+name+' <= '+hstr

  return,str
END 






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; General directory and file name stuff
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Info each directory type
FUNCTION objshear::lens_dirtype_info

  info = {dirtype: '', basetype: '', dirflags: 0}

  dirtypes  = $
    ['input','par',  'pbs',  'output','pairs','combined','corr',  'corr_invert','jackknife','jackknife_invert','random','plot']
  basetypes = $
    ['input','input','input','output','output','output',  'output','output',     'output',   'output',          'output','plot']
  ; 0 if cannot be part of a subtype 1 if can
  dirflags = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]

  ntypes = n_elements(dirtypes)
  info = replicate(info, ntypes)
  info.dirtype = dirtypes
  info.basetype = basetypes
  info.dirflags = dirflags
  return, info
END 

FUNCTION objshear::lens_dirtype_exists, dtype, info=info

  IF n_elements(dtype) EQ 0 THEN BEGIN 
      message,$
        '-syntax: if obj->lens_dirtype_exists(dtype, basetype=) then ...'
  ENDIF 

  info = self->lens_dirtype_info()
  dirtype = strlowcase(dtype)
  w=where(info.dirtype EQ dirtype, nw)
  
  IF nw NE 0 THEN BEGIN 
      info = info[w[0]]
      return,1
  ENDIF ELSE BEGIN 
      basetype = 'Unknown'
      return,0
  ENDELSE 

END 

; base directory
FUNCTION objshear::lensdir_base, dtype, info=info
  on_error, 2
  status = 1
  IF n_elements(dtype) EQ 0 THEN BEGIN 
      message,'-Syntax: basedir = obj->lensdir_base(dirtype)',/inf
      message,'    Will be the lensinput or lensoutput dir'
  ENDIF 

  IF NOT self->lens_dirtype_exists(dtype, info=info) THEN BEGIN 
      message,'Unmatched dirtype: '+strlowcase(dtype),/inf
      message,'Try one of these: '+strjoin(info.dirtype, ', ')
  ENDIF 

  status = 0  

  CASE info.basetype OF
      'input': return,expand_tilde(esheldon_config('lensinput_dir'))
      'output': return,expand_tilde(esheldon_config('lensout_dir'))
      'plot': return,expand_tilde(esheldon_config('plot_dir'))
      ELSE: message,'Unknown basetype: '+info.basetype
  ENDCASE 

END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The workhorses for generating file names
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION objshear::lensdir, dtypein, subtype=subtype, sample=sample, createdir=createdir, base=base

  on_error, 2
  IF n_elements(dtypein) EQ 0 THEN BEGIN 
      message,'-Syntax: dir = obj->lensdir(dirtype, subtype=, sample=, /createdir)',/inf
      message,'         dirtype = (input|output|combined|corr|jackknife|random)'
  ENDIF 
  dtype = strlowcase(dtypein)
  basedir = self->lensdir_base(dtype, info=info)

  ;; this type of lens, e.g. maxbcg
  type = self->type()
  dir = concat_dir(basedir, type)

  IF keyword_set(base) THEN return,dir

  ;; the sample. Note multiple samples get joined by a dash
  sstr = self->sample_string(sample=sample)
  dir = concat_dir(dir, sstr)

  ;; The subtype
  IF info.dirflags EQ 1 AND n_elements(subtype) NE 0 THEN BEGIN 
      IF dtype NE 'plot' THEN dir = concat_dir(dir, 'sub')
      dir = concat_dir(dir, strlowcase(subtype))
  ENDIF

  ;; finally the dirtype if not 'plot'
  IF dtype NE 'plot' THEN BEGIN 
      dir = concat_dir(dir, dtype)
  ENDIF 

  IF keyword_set(createdir) THEN BEGIN 
      IF NOT file_test(dir, /dir) THEN file_mkdir, dir
  ENDIF 

  return,dir
END 



;; Requires that the nbin function be overridden in the
;; descendent
FUNCTION objshear::subtype_nbin, subtype
  message,'You must override this function in the descendent'
END 
FUNCTION objshear::lensfile, dtype, subtype=subtype, bin=bin, sample=sample, encapsulated=encapsulated, nodir=nodir, createdir=createdir, ext=ext

  on_error, 2
  IF n_elements(dtype) EQ 0 THEN BEGIN 
      message,'-Syntax: files = obj->lensfile(dirtype, subtype=, bin=, sample=, /encapsulated, /nodir, /createdir, ext=)',/inf
      message,'         dirtype = (input|output|pairs|combined|corr|jackknife|random)'
  ENDIF 

  dirtype = strlowcase(dtype)
  IF NOT self->lens_dirtype_exists(dtype, info=info) THEN BEGIN 
      message,'Unknown dirtype: '+dirtype,/inf
      message,'Try one of these: '+strjoin(info.dirtype,', ')
  ENDIF 

  type = self->type()
  sstr = self->sample_string(sample=sample)
 
  file = type + '_' + sstr 

  IF info.dirflags EQ 1 AND n_elements(subtype) NE 0 THEN BEGIN 

      file = file + '_' + strlowcase(subtype)

      nbin = n_elements(bin)
      IF n_elements(bin) EQ 0 THEN BEGIN 
          nbin = self->subtype_nbin(subtype)
          bin = lindgen(nbin)
      ENDIF 

      FOR i=0L, nbin-1 DO BEGIN 
          tbinstr = strn(bin[i], length=2, padchar='0')
          add_arrval, tbinstr, binstr
      ENDFOR 

      file = file + '_' + binstr

  ENDIF 

  CASE dirtype OF
      'par': ext='.conf'
      'pairs': ext='.bin'
      'pbs': ext='.pbs'
      ELSE: ext='.st'
  ENDCASE 

  file = file + '_' + dirtype + ext

  IF NOT keyword_set(nodir) THEN BEGIN 
      dir = self->lensdir(dirtype, subtype=subtype, sample=sample, createdir=createdir)
      file = concat_dir(dir, file)
  ENDIF 

  status = 0
  return, file
END 


FUNCTION objshear::_getfiles, files, hdr=hdr, columns=columns, status=status

  IF n_elements(files) EQ 1 THEN print,'Reading file: ',files
  struct = read_idlstruct_multi(files, hdr=hdr, columns=columns, $
                                status=status)

  return,struct

END 

FUNCTION objshear::lensread, dtype, subtype=subtype, bin=bin, sample=sample, columns=columns, hdr=hdr, count=count, nrows=nrows

    on_error, 2
    status = 1
    if n_elements(dtype) eq 0 then begin 
        message,$
            '-Syntax: struct = obj->lensread(dirtype, subtype=, bin=, sample=, columns=, hdr=, count=)',/inf
        message,$
            '         dirtype = (input|output|pairs|combined|corr|jackknife|random)'
    endif 

    files = self->lensfile(dtype, subtype=subtype, bin=bin, sample=sample)
    if dtype eq 'pairs' then begin
        struct = self->read_pairs(files[0], nrows=nrows, status=status)
    endif else begin
        struct = self->_getfiles(files, columns=columns, hdr=hdr, status=status)
    endelse

    if status ne 0 then message,'Could not read file'
    count=n_elements(struct)
    return,struct

end 


; Convert my .st files to .fits
PRO objshear::lensfile_convert, dtype, subtype=subtype, sample=sample, overwrite=overwrite

  file = self->lensfile(dtype, subtype=subtype, sample=sample)

  nfile = n_elements(file)

  FOR i=0L, nfile-1 DO BEGIN 

      print,'-------------------------------------------------------'
      newfile = repstr(file[i], '.st', '.fits')
      
      print
      print,'Converting file: '
      print,file[i]
      print,newfile
      st = read_idlstruct(file[i])
      if fexist(newfile) and not keyword_set(overwrite) then begin
          message,'File exists; send /overwrite'
      endif
      mwrfits, st, newfile, /create
  ENDFOR 

END 











; This is for working with the individual random files that correspond
; the the lens files.  These correspond to the sample directory because
; the randoms have been matched to the lenses and outputs made.  To 
; work with the random inputs and outputs you should see randlensing__define
FUNCTION objshear::randnumfile, randnum, $
                 subtype=subtype, bin=bin, createdir=createdir, nodir=nodir

  on_error, 2

  nrand = n_elements(randnum)
  IF nrand EQ 0 THEN BEGIN 
      message,'-Syntax: file = obj->randnumfile(randnum, subtype=, bin=, /nodir, /createdir)'
  ENDIF 

  IF n_elements(subtype) NE 0 THEN BEGIN 
      nbin = n_elements(bin)
      IF nbin EQ 0 THEN BEGIN 
          nbin = self->subtype_nbin(subtype)
          bin = lindgen(nbin)
      ENDIF 

      ;; use recursion
      IF nbin GT 1 THEN BEGIN 
          IF nrand GT 1 THEN BEGIN 
              message,'Cannot work with multiple rand numbers and '+$
                'bin numbers     ; '+$
                'send a scalar bin= or scalar randnum'
          ENDIF 
          files = strarr(nbin)
          FOR i=0L, nbin-1 DO BEGIN 
              files[i] = self->randnumfile(randnum, subtype=subtype, $
                                           bin=bin[i], createdir=createdir)
          ENDFOR 
          return,files
      ENDIF 
  ENDIF 

  type = self->type()
  sstr = self->sample_string()
  
  file = type + '_' + sstr

  IF n_elements(subtype) NE 0 THEN BEGIN 

      file = file + '_' + strlowcase(subtype)
      binstr = strn(bin[0], length=2, padchar='0')
      file = file + '_' + binstr

  ENDIF 

  file = file + '_random'

  FOR i=0L, nrand-1 DO  BEGIN 
      rstr = '_'+strn(randnum[i],length=2,padchar='0')
      add_arrval, rstr, addrand
  ENDFOR 
  file = file + addrand

  file = file + '.st'

  IF NOT keyword_set(nodir) THEN BEGIN 
      dir = self->lensdir('random',subtype=subtype, createdir=createdir)
      file = concat_dir(dir, file)
  ENDIF 


  return, file

END 

FUNCTION objshear::randnumread, randnum, subtype=subtype, bin=bin, columns=columns, hdr=hdr, status=status

  on_error, 2
  status = 1
  IF n_elements(randnum) EQ 0 THEN BEGIN 
      message,$
        '-Syntax: struct = obj->randnumread(randnum, subtype= ,bin=, columns=, hdr=hdr, status=)'
  ENDIF 

  files = self->randnumfile(randnum, subtype=subtype, bin=bin)
  return,self->_getfiles(files, columns=columns, hdr=hdr, status=status)

END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read random "input" and "output" files 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO objshear::read_random_files, randnum, rlensin, rlensum, hdr, nozhist=nozhist
  IF not keyword_set(nozhist) THEN BEGIN 
      ;; This will change in the future
      rl = obj_new('randlensing', self->random_sample())
      rlensin = rl->randread('input',randnum)
      rlensum = rl->randread('output',randnum,hdr=hdr)
      hdr = idlstruct_hclean(hdr)
      obj_destroy, rl
  ENDIF ELSE BEGIN 
      message,'This must be re-written for /nozhist'
      rlensin = self->lensinput_get(randnum=randnum)
      rlensum = self->lensoutput_get(randnum=randnum,hdr=hdr)
      hdr = idlstruct_hclean(hdr)
  ENDELSE 
END 






; Plots
FUNCTION objshear::plot_dir, subtype=subtype, sample=sample, base=base, createdir=createdir
  dir = self->lensdir('plot', sample=sample, base=base)

  IF n_elements(subtype) NE 0 AND NOT keyword_set(base) THEN BEGIN 
      dir = concat_dir(dir, strlowcase(subtype))
  ENDIF 

  IF keyword_set(createdir) THEN BEGIN 
      IF NOT file_test(dir, /dir) THEN file_mkdir, dir
  ENDIF 

  return,dir
END 

; This is generic for the delta sigma plots of standard output types 
; "combined", "corr", "jackknife".  Currently other types must be 
;  generated separately
FUNCTION objshear::plotfile, dtypein, subtype=subtype, bin=bin, sample=sample, color=color, encapsulated=encapsulated, ratio=ratio, nodir=nodir, createdir=createdir
  on_error, 2
  IF n_elements(dtypein) EQ 0 THEN BEGIN 
      message,'-Syntax: files = obj->plotfile(dirtype, subtype=, bin=, color=, encapsulated=, /ratio, /nodir, /createdir)'
  ENDIF 
  
  dtype = strlowcase(dtypein)

  tfile = self->lensfile(dtype, subtype=subtype, bin=0, sample=sample, $
                         /nodir)
  
  IF keyword_set(encapsulated) THEN ext='.eps' ELSE ext='.ps'
  IF keyword_set(color) THEN cstr = '_color' ELSE cstr=''
  IF keyword_set(ratio) THEN rstr = '_ratio' ELSE rstr=''

  ext = rstr + cstr + ext

  file = repstr(tfile, '.st', ext)

  IF n_elements(subtype) NE 0 THEN BEGIN 
      IF n_elements(bin) NE 0 THEN BEGIN 
          ;; plot for specific bin
          bstr = strn(bin[0], len=2, padchar='0')
          file = repstr(file, '00_'+dtype, bstr+'_'+dtype)
      ENDIF ELSE BEGIN 
          file = repstr(file, '00_'+dtype, dtype)
      ENDELSE 
  ENDIF 

  IF NOT keyword_set(nodir) THEN BEGIN 
      dir = self->lensdir('plot', subtype=subtype, sample=sample, createdir=createdir)
      file = concat_dir(dir, file)
  ENDIF 
  return, file
END 















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;; Tools for setting up lens inputs
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO objshear::calc_cosmo, lensum

	; Calculate some cosmology-dependent quantities and
	; copy into the structure

	par_struct = self->par_struct()

	nlens = n_elements(lensum)
	angmax = replicate(10000d, nlens)
	DL = replicate(1.e-10, nlens) ; small DL are thrown out, so default small

	wtmp = where(lensum.z GT 0)

	; make sure it's kpc since rmax is in kpc
	DL[wtmp] = $
		1000*angdist(0.0, lensum[wtmp].z, $
				     h=par_struct.h, omega_m=par_struct.omega_m)

	angmax = par_struct.rmax/DL*180.0/!pi ; angle in degrees

	lensum.DL = DL
	lensum.angmax = angmax

END 


FUNCTION objshear::lensinput_structdef

  struct = { zindex: 0L, $
             ra: 0d, $
             dec: 0d, $
             clambda: 0d, $
             ceta: 0d, $
             z: 0.0, $
             angmax: 0.0, $
             DL: 0.0, $
             mean_scinv: 0.0, $
             pixelMaskFlags: 0 $
           }

  return,struct

END 


;; This applies the photometric masks, without edge cuts, to a set of input points
FUNCTION objshear::apply_photometric_masks, clambda, ceta
  
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: maskflags = mb->apply_masks(clambda, ceta)'
      return,-1
  ENDIF 

  nlam = n_elements(clambda) & neta = n_elements(ceta)
  IF nlam NE neta THEN message,'clambda and ceta must be same size'
  maskFlags = intarr(nlam)

  print
  print,'Applying simple mask'
  apply_pixel_mask, clambda, ceta, simpleMasked, simpleUnmasked, /simple
  IF simpleMasked[0] NE -1 THEN BEGIN 
      maskFlags[simpleMasked] = maskFlags[simpleMasked] + !FLAGS_MASKED_SIMPLE
  ENDIF 

  ;; These masks aren't ready yet
  IF 0 THEN BEGIN 
      print,'Applying bound mask'
      apply_pixel_mask, $
        clambda, ceta, boundMasked, boundUnmasked, /bound
      IF boundMasked[0] NE -1 THEN BEGIN 
          maskFlags[boundMasked] = maskFlags[boundMasked] + !FLAGS_MASKED_BOUND
      ENDIF 

      print,'Applying combined mask'
      apply_pixel_mask, $
        clambda, ceta, combinedMasked, combinedUnmasked, /combined
      IF combinedMasked[0] NE -1 THEN BEGIN 
          maskFlags[combinedMasked] = $
            maskFlags[combinedMasked] + !FLAGS_MASKED_COMBINED
      ENDIF 

  ENDIF 

  return,maskFlags

END 







function objshear::redshift_cuts, str, nkeep, index=index

	if n_elements(index) eq 0 then begin
		index = lindgen(n_elements(str))
	endif
	ninput=n_elements(index)

	pars = self->par_struct()
	print,'Cutting maxz = ',pars.maxz
	w=where(str[index].z gt 0.0 and str[index].z lt pars.maxz, nkeep)
	if nkeep gt 0 then begin
		w=index[w]
	endif
	print,'Threw out: ',ninput-nkeep
	return, w
end

function objshear::angle_cuts, str, nkeep, index=index

	if n_elements(index) eq 0 then begin
		index = lindgen(n_elements(str))
	endif
	ninput=n_elements(index)


	pars = self->par_struct()
	print,'Cutting max_allowed_angle = ',pars.max_allowed_angle,' degrees'

	w=where(str[index].angmax lt pars.max_allowed_angle, nkeep)
	if nkeep gt 0 then begin
		w=index[w]
	endif
	print,'Threw out: ',ninput-nkeep
	return, w
end

function objshear::pixel_mask_cuts, str, nkeep, index=index, $
		maskflags=maskflags

	pars=self->par_struct()

	if n_elements(index) eq 0 then begin
		index = lindgen(n_elements(str))
	endif
	ninput=n_elements(index)

	; This will return a defined maskfile in certain situations
	self->maskfile, maskfile

	if n_elements(maskfile) ne 0 then begin 
		print
		message,'Using maskfile = '+ntostr(maskfile),/inf
	endif 

	print
	if pars.edgecuts then begin 
		print,'Applying pixel mask with edge cuts'
		apply_pixel_mask, $
			str[index].clambda, str[index].ceta, $
			masked, unmasked, maskFlags, $
			maxangle=str[index].angmax, /twoquad, maskfile=maskfile
	endif else begin 
		print,'Applying pixel mask (no edge cuts)'
		apply_pixel_mask, $
			str[index].clambda, str[index].ceta, $
			masked, unmasked, maskFlags, $
			maskfile=maskfile
	endelse 

	if unmasked[0] eq -1 then begin
		nkeep=0
		w=-1
	endif else begin
		w=index[unmasked]
		nkeep=n_elements(w)
	endelse
	print,'Threw out: ',ninput-nkeep
	return, w
end


function objshear::lenscuts, lensum, nkeep, doplot=doplot

	nlens = n_elements(lensum)
	wlens = self->redshift_cuts(lensum, nkeep)
	if nkeep eq 0 then return, -1

	wlens = self->angle_cuts(lensum, nkeep, index=wlens)
	if nkeep eq 0 then return, -1

	wlens = self->pixel_mask_cuts(lensum, nkeep, index=wlens, $
		maskflags=maskflags)
	if nkeep eq 0 then return, -1
	lensum.pixelmaskflags = maskflags

	return, wlens
end

function objshear::lenscuts_old, lensum, nkeep, doplot=doplot

	;; Make generic lens cuts. Cut on: 
	;;   max_allowed_angle (effectively low-z cut)
	;;   maxz
	;;   z > 0
	;;   photometric mask: basic with edge checking

	par = self->par_struct()

	;; Redshift cuts

	max_allowed_angle = par.max_allowed_angle
	maxz = par.maxz

	minz = 0.0
	print,'Cutting maxz = ',maxz
	print,'Cutting max_allowed_angle = '+ntostr(max_allowed_angle)+' degrees'


	nlens_init = n_elements(lensum)
	wlens = where($
		lensum.z GT minz AND $
		lensum.z LT maxz AND $
		lensum.angmax LT max_allowed_angle, $
		ngoodz, comp=wbadz, ncomp=nbadz)

	if ngoodz EQ 0 then begin 
		print,'***************************************'
		print,'No objects passed redshift cuts'
		nkeep = 0
		return,-1
	endif 

	print,'---------------------------------------------------------'
	print,'Threw out '+ntostr(nlens_init-ngoodz)+' in redshift cuts'
	print,'---------------------------------------------------------'  
	nkeep = ngoodz

	; now apply the pixel mask.  This is now all about how the lenses
	; are related to the distribution of sources

	;; This will return a defined maskfile in certain situations
	self->maskfile, maskfile

	IF n_elements(maskfile) NE 0 THEN BEGIN 
		print
		message,'Using maskfile = '+ntostr(maskfile),/inf
	ENDIF 

	print
	IF par.edgecuts THEN BEGIN 
		print,'Applying pixel mask with edge cuts'
		apply_pixel_mask, $
			lensum[wlens].clambda, lensum[wlens].ceta, $
			masked, unmasked, maskFlags, $
			maxangle=lensum[wlens].angmax, /twoquad, maskfile=maskfile
	ENDIF ELSE BEGIN 
		print,'Applying pixel mask (no edge cuts)'
		apply_pixel_mask, $
			lensum[wlens].clambda, lensum[wlens].ceta, $
			masked, unmasked, maskFlags, $
			maskfile=maskfile
	ENDELSE 

	IF keyword_set(doplot) THEN BEGIN 

		dir = self->lensdir('plot', /createdir)
		pngfile = self->type()+'_'+self->sample_string()+'_maskplot.png'
		pngfile = concat_dir(dir, pngfile)

		setupplot, 'z'
		simpctable, rct, gct, bcg
		device, set_resolution=[1280,1024]
		charsize = 1.5

		plot,lensum.clambda,lensum.ceta,psym=3,$
			xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c'), $
			charsize=charsize

		IF nbadz GT 0 THEN BEGIN 
			clam = [lensum[wbadz].clambda]
			ceta = [lensum[wbadz].ceta]
			oplot,clam,ceta, $
				psym=8, symsize=0.5, color=c2i('red')
		ENDIF 

		IF masked[0] NE -1 THEN BEGIN 
			w=wlens[masked]
			clam = [lensum[w].clambda]
			ceta = [lensum[w].ceta]
			oplot,clam,ceta,psym=3,color=c2i('darkgreen')
		ENDIF 

		deg = ntostr( fix( max_allowed_angle) )
		mess=[!csym.theta+' > '+deg+!csym.degrees, 'masked']

		legend,mess,psym=[8,8],color=c2i(['red','darkgreen']),$
			/right,box=0,charsize=charsize*3.0/4.0

		IF par.edgecuts THEN BEGIN 
			mess = 'Edgecuts applied'
		ENDIF ELSE BEGIN 
			mess = 'No edgecuts'
		ENDELSE 

		mess = ['Sample '+ntostr(self->sample()), $
			'Rmax = '+ntostr((par.rmax/1000.0),4,/round)+' Mpc',$
			mess]
		legend, mess, /left, box=0, charsize=charsize*3.0/4.0

		print
		print,'Plotting to file: ',pngfile
		write_png, pngfile, tvrd(), rct, gct, bcg

		setupplot, 'X'
	ENDIF 

	IF unmasked[0] EQ -1 THEN BEGIN 
		print,  '*************************************'
		print,'No objects passed masks and edgecuts'
		nkeep = 0
		return,-1
	ENDIF 
	lensum[wlens].pixelMaskFlags = maskFlags

	nunmasked = n_elements(unmasked)
	print,'Threw out '+ntostr(nkeep-nunmasked)+' in pixel masks and edge cut'

	wlens = wlens[unmasked]
	nkeep = n_elements(unmasked)

	return,wlens
END 


; This is the most general setup. Should be over-ridden if needed
PRO objshear::setuplens
  
  tm = systime(1)

  ;; Parameters
  par_struct = self->par_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file name(s)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outfile = self->lensfile('input', /createdir)

  ;; make sure dir exists
  tfile = self->lensfile('output', /createdir)
  print
  print,'Setting up file: ',outfile

  ;; get the data.  The get() method must be defined in the descendant.
  ;;  Note: zindex will be the same for all samples created from this lens
  ;;        catalog
  cat = self->get()
  ncat = n_elements(cat)
  zindex = lindgen(ncat)

  nlens_init = n_elements(cat)      

  print,'Creating lens struct'
  lstruct = self->lensinput_structdef()
  lstruct = replicate(lstruct, nlens_init)
  
  print,'Copying....'
  copy_struct, cat, lstruct

  IF NOT tag_exist(cat, 'clambda') THEN BEGIN 
      eq2csurvey, cat.ra, cat.dec, clam, ceta
      lstruct.clambda = clam
      lstruct.ceta = ceta
  ENDIF 
  
  ;; Keep track of the objects and their redshifts
  lstruct.zindex = zindex
  
  ;; Calculate some cosmology-dependent stuff and copy into struct
  self->calc_cosmo, lstruct
  
  ;; Make some generic lens cuts
  wlens = self->lenscuts(lstruct, nkeep, /doplot)
  
  ;; remove the unwanted lenses
  lstruct = lstruct[wlens]
  
  print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'
    
  print
  print,'Writing to file: ',outfile
  write_idlstruct, lstruct, outfile
  

  ptime,systime(1)-tm

END 






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generate random lenses in the source photometric mask.  These can
;; be run through different selections later according to each lens
;; sample
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; number of randoms
FUNCTION objshear::numrand
  sample = self->random_sample()
  rl = obj_new('randlensing', sample)
  numrand = rl->numrand()
  obj_destroy,rl
  return,numrand
END 
; random file numbers
FUNCTION objshear::randnum
  return,lindgen(self->numrand())
END 

; number of randoms per file
FUNCTION objshear::nperfile
  return,500000
END 

; stripes in source catalog.  Used to generate the random
; points
FUNCTION objshear::source_stripes
  ms = obj_new('make_scat')
  stripes = ms->stripes_use()
  obj_destroy, ms
  return,stripes
END 

; redshift range to use in generating the randoms.  Note,
; currenlty the max_allowed_angle cuts off at z=0.038 for
; rmax = 10Mpc

FUNCTION objshear::randzrange
  return,[0.02, 0.35]
END 


; generate random clambda,ceta points based on mask and source
; galaxy stripes
FUNCTION objshear::genrand, nrand

  par = self->par_struct()

  ;; What mask file to use?
  ;; This will return a defined maskfile in certain situations.
  self->maskfile, maskfile

  IF n_elements(maskfile) NE 0 THEN BEGIN 
      print
      message,'Using maskfile = '+ntostr(maskfile),/inf
  ENDIF 

  stripes  = self->source_stripes()
  ;; Output struct
  ss = create_struct('clambda', 0d, $
                     'ceta', 0d, $
                     'maskflags', 0)
  outstruct = replicate(ss, nrand)


  ;; We will just generate them from the photometric mask "basic"
  ;; But keep track of the other masks too

  sdss_genrand, stripes, nrand, rlam, reta, $
                phot_maskfile=maskfile
      
  ;; Don't have other masks for princeton
  IF par.shape_correction_style NE 3 THEN BEGIN 
      maskflags = self->objshear::apply_photometric_masks(rlam, reta)
      outstruct.maskflags = maskflags
  ENDIF 

  outstruct.clambda = rlam
  outstruct.ceta = reta

  return,outstruct

END 

; Generate random redshifts in volume limited sample
FUNCTION objshear::genrandz, nrand

  par = self->par_struct()
  zrange = self->randzrange()
  omega_k = 0.0
  omega_m = par.omega_m
  omega_l = 1.0 - omega_m

  h = par.h

  cm = obj_new('cosmology')

  zrand = cm->genrandz(nrand, zrange[0], zrange[1], $
                       omega_k, omega_m, omega_l, h=h)

  obj_destroy, cm

  return,zrand

END 

; set up the lens input files for randoms
PRO objshear::setuprand, randnum=randnum

  tm = systime(1)

  IF n_elements(randnum) EQ 0 THEN randnum = self->randnum()

  nperfile = self->nperfile()
  
  nrand = n_elements(randnum)
  FOR i=0L, nrand-1 DO BEGIN 

      tmi = systime(1)


      outfile = self->rand_inputfile(randnum[i])
      print,'*********************************************************'
      print,'Output file: ',outFile

      ;; random points on the sky
      print
      rand = self->genrand(nperfile)


      print,'Creating lensum struct'
      lstruct = self->objshear::lensinput_structdef()

      nlens_init = n_elements(rand)
      lstruct = replicate(lstruct, nlens_init)
      
      lstruct.clambda = rand.clambda
      lstruct.ceta = rand.ceta

      csurvey2eq, rand.clambda, rand.ceta, ra, dec
      lstruct.ra  = ra
      lstruct.dec = dec

      delvarx, rand

      print
      print,'Assigning redshifts'

      lstruct.z = self->genrandz(nperfile)
      
      ;; Calculate some stuff and copy int
      self->objshear::calc_cosmo, lstruct

      ;; Make some generic lens cuts
      wlens = self->objshear::lenscuts(lstruct, nkeep)

      ;; remove the unwanted lenses
      lstruct = lstruct[wlens]

      print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'
      
      print
      print,'Writing to file: ',outfile
      write_idlstruct, lstruct, outfile

      print,'Time for this random'
      ptime,systime(1)-tmi

  ENDFOR 
  

  print,'Overall'
  ptime,systime(1)-tm


END 

  






















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; parameter (config) files for the lensing code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;









; No longer need to add randoms
PRO objshear::write_parfile

  p = self->par_struct()

  ;; do lenses and random
  par_file = self->lensfile('par', /createdir)
  lensin_file = self->lensfile('input')
  lensout_file = self->lensfile('output', /createdir)
  if p.dopairs then begin
      pair_file = self->lensfile('pairs', /createdir)
  endif else begin
      pair_file = 'None'
  endelse


  ;; The user can give their own source file name
  ;; to override the defaults
  ms = obj_new('make_scat')

  CASE p.shape_correction_style OF
      3: BEGIN 
          source_file = ms->princeton_source_file(p.source_sample)
          htmrev_file = ms->princeton_htmrev_file(p.source_sample)
      END 
      4: BEGIN 
          source_file = ms->intrinsic_source_file(p.source_sample)
          htmrev_file = ms->intrinsic_htmrev_file(p.source_sample)
      END 
      ELSE: BEGIN 
          source_file = ms->source_file(p.source_sample)
          htmrev_file = ms->htmrev_file(p.source_sample)
      END 
  ENDCASE 

  obj_destroy,ms

  if p.scinv_sample ne -1 then begin
      sc = obj_new('sdss_sigma_crit', p.scinv_sample)
      scinv_file = sc->file('output', 'meanscinv', project='scinv')
      obj_destroy, sc
  endif else begin
      scinv_file='None'
  endelse

  print
      
  print,'Writing to par file: ',par_file
  openw, lun, par_file, /get_lun
  
  printf, lun, 'lens_file       '+lensin_file
  printf, lun, 'source_file     '+source_file
  printf, lun, 'htmrev_file     '+htmrev_file
  printf, lun, 'scinv_file      '+scinv_file
  printf, lun, 'output_file     '+lensout_file
  printf, lun, 'pair_file       '+pair_file
  
  printf, lun, 'pixel_lensing   0'
  printf, lun, 'dopairs         '+string(p.dopairs,format='(i0)')

  printf, lun, 'h               '+string(p.h,format='(g0)')
  printf, lun, 'omega_m         '+string(p.omega_m,format='(g0)')
  printf, lun, 'sigmacrit_style '+string(p.sigmacrit_style,format='(i0)')
  printf, lun, 'shape_correction_style '+$
    string(p.shape_correction_style,format='(i0)')
  
  printf, lun, 'logbin          '+string(p.logbin,format='(i0)')
  printf, lun, 'nbin            '+string(p.nbin,format='(i0)')
  printf, lun, 'binsize         '+string(p.binsize,format='(g0)')
  
  printf, lun, 'rmin            '+string(p.rmin,format='(g0)')
  printf, lun, 'rmax            '+string(p.rmax,format='(g0)')
  
  printf, lun, 'comoving        '+string(p.comoving,format='(i0)')
  printf, lun, 'depth           '+string(p.depth,format='(i0)')
  printf, lun, 'zbuffer         '+string(p.zbuffer,format='(g0)')
  
  printf, lun, 'sample          '+string(p.sample, format='(i0)')
  printf, lun, 'source_sample   '+string(p.source_sample, format='(i0)')
  printf, lun, 'catalog         '+p.catalog
  free_lun, lun


END 





PRO objshear::write_script

  par_dir = self->lensdir('par', /createdir)
  par_file = self->lensfile('par', /nodir)
  
  script_file = 'run_'+self->sample_string()+'.sh'
  script_file = concat_dir(par_dir, script_file)

  address = 'erin.sheldon@gmail.com'

  print
  print,'Writing to file: ',script_file

  openw, lun, script_file, /get
  printf, lun, '#!/bin/sh'
  printf, lun
  printf, lun, '/usr/bin/time -p objshear '+par_file

  mess = 'finished objshear '+par_file
  printf,lun,'dt=`date`'
  printf, lun, $
    'echo "'+mess+' ${dt}" | mail '+address+' -s "'+mess+'"'

  free_lun, lun
  spawn,['chmod','755',script_file],/noshell


return

END 







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; combining the lensum files.  Note, the sub-sample code does this 
; already for the subsamples
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO objshear::combine_lensum, randnum=randnum, zbinsize=zbinsize


  nrand = n_elements(randnum)
  IF nrand EQ 0 THEN BEGIN 

      comb_file = self->lensfile('combined',/createdir)
      lensum = self->lensread('output', hdr=hdr)

      
      print,'combining lensum'
      sh = combine_lensum(lensum, hdr=hdr)
      hdr = idlstruct_hclean(hdr)
      
      print
      print,'Writing to file: ',comb_file
      
      write_idlstruct, sh, comb_file, hdr=hdr
  ENDIF ELSE BEGIN 

      plot_dir = self->plot_dir(/createdir)
      psfile = concat_dir(plot_dir,$
                          'matchzrand_'+self->sample_string()+'_N1.ps')

      ;; don't overwrite
      WHILE fexist(psfile) DO psfile=newname(psfile)

      begplot,name=psfile, /color

      ;; for redshift matching
      lensin = self->lensread('input')
      comb_files = self->randnumfile(randnum, /createdir)

      FOR i=0L, nrand-1 DO BEGIN 

          print,'------------------------------------------------------------------------'
          self->read_random_files, randnum[i], rlensin, rlensum

          rkeep = self->rand_match(lensin, rlensin, zbinsize=zbinsize)

          print,'combining lensum'
          rsh = combine_lensum(rlensum, index=rkeep)
          
          print
          print,'Writing to file: ',comb_files[i]
          
          write_idlstruct, rsh, comb_files[i]
 
          rlensin=0
          rlensum=0
         
      ENDFOR 

      endplot

  ENDELSE 

END 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Tools for selecting lenses and sub sampling them
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; This should get overridden by inherited class
FUNCTION objshear::where_string, subtype, nbin=nbin
  message,'No where string defined. You must override this method'
END 

FUNCTION objshear::where_select, struct, type, nkeep=nkeep, nbin=nbin, labels=labels, average_tags=average_tags, nodisplay=nodisplay

  IF n_params() LT 2 THEN BEGIN 
      message,'-Syntax: keep = obj->where_select(struct, type [, nkeep=, nbin=, labels=, average_tags=])'
  ENDIF 
  where_string = self->where_string(type, nbin=nbin, $
                                    labels=labels, $
                                    average_tags=average_tags, $
									nodisplay=nodisplay)
  keep = self->struct_select(struct, where_string, nkeep)
  return,keep

END 

function objshear::altselect, struct, subtype, nkeep=nkeep, nbin=nbin, $
    average_tags=average_tags

    message,'inherit from objshear and implement this'
end

;; general program for selecting elements of a structure based upon
;; a where string
FUNCTION objshear::struct_select, $
                 struct, where_string, nkeep, $
                 verbose=verbose

  nkeep = 0L
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: keep = o->struct_select(struct, where_string, nkeep, '
      print,'                                 /verbose)'
      print
      print,'where_string is a string sent to the where() function.  It must'
      print,'  refer to "struct" and its tags'
      print,'keep is a pointer (or pointer array) containing indices for '
      print,'  each where string'
      print,'    e.g.  where_string = "struct.z gt 0.05 and struct.x le 23"'
      return,-1
  ENDIF 
  
  nselect = n_elements(where_string)
  IF nselect EQ 0 THEN BEGIN 
      message,'You must enter a valid where string',/inf
      return,-1
  ENDIF 

  IF nselect EQ 1 THEN BEGIN 
      keep = ptr_new()
      nkeep = 0L
  ENDIF ELSE BEGIN 
      keep = ptrarr(nselect)
      nkeep = lonarr(nselect)
  ENDELSE 

  FOR i=0L, nselect-1 DO BEGIN 

      IF keyword_set(verbose) THEN print,where_string[i]

      select_command = 'tkeep = where('+where_string[i]+', tnkeep)'
      IF NOT execute(select_command) THEN BEGIN 
          ptr_free, keep
          nkeep[*] = 0
          message,'Error executing select statement'
      ENDIF 
      
      keep[i] = ptr_new(tkeep, /no_copy)
      nkeep[i] = tnkeep

  ENDFOR 
  
  return,keep

END 





;; match lensum to rand lensum based on the keep array
FUNCTION objshear::rand_match, lenses, rlenses, nozhist=nozhist, input_index=input_index, zbinsize=zbinsize, title=title

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: rmatch = objshear->rand_match(lenses,rlenses,/nozhist,input_index=)'
      print,'  For not /nozhist need redshifts, otherwize zindex'
      print,'  input_index subscripts lenses'
  ENDIF 

  nlens = n_elements(lenses)
  IF n_elements(input_index) EQ 0 THEN use = lindgen(nlens) $
  ELSE use = input_index

  IF n_elements(zbinsize) EQ 0 THEN zbinsize = 0.01

  ;; Are we matching redshift histograms?
  IF not keyword_set(nozhist) THEN BEGIN 
      nrand = n_elements(rlenses)
      zrand = rlenses.z
      IF self->random_sample() EQ 9 OR self->random_sample() EQ 10 THEN BEGIN 
          ;; need to randomize the discrete redshifts
          zwidth = 0.0027
          zrand = zrand + (randomu(seed, nrand) - 0.5)*zwidth
      ENDIF 
      
      rmatch = match_hist(lenses[use].z, zrand, zbinsize, $
                          x1h=x1h, x2h=x2h, new_x2h=new_x2h)
      minz=min(lenses[use].z, max=maxz)
      plothist, lenses[use].z, bin=zbinsize, /norm, $
        xtitle='z', ytitle='P(z)', title=title, $
        aspect=!gratio
      plothist, zrand[rmatch], bin=zbinsize, min=minz,max=maxz,/norm, $
        /overplot, color=c2i('red'), line=3
      legend,['lenses','random'],line=[0,3],$
        color=[!p.color,c2i('red')], /left, box=0, charsize=1
      key = prompt_kbrd('hit a key')
      
  ENDIF ELSE BEGIN 
      match_multi, lenses[use].zindex, rlenses.zindex, rmatch
  ENDELSE 

  rmatch = rmatch[ sort(rmatch) ]
  return,rmatch
  

END 

; This takes a pointer array of "keep" indices and returns the same.
; See rand_match for simpler version.
FUNCTION objshear::rand_match_bins, lenses, rlenses, keep, nozhist=nozhist,$
                 zbinsize=zbinsize

  IF (n_elements(lenses) EQ 0 OR $
      n_elements(rlenses) EQ 0 OR $
      n_elements(keep) EQ 0) THEN BEGIN 

      on_error, 2
      print,'-Syntax: rkeep = o->rand_match_bins(lenses, rlenses, keep, /nozhist, zbinsize=)'
      print
      message,'Halting'
  ENDIF 

  nrand = n_elements(rlenses)

  nbin = n_elements(keep)
  rkeep = ptrarr(nbin)
  FOR i=0L, nbin-1 DO BEGIN 
      IF ptr_valid(keep[i]) THEN BEGIN 
          IF n_elements(*keep[i]) NE 0 THEN BEGIN 

              use = *keep[i]

              title = 'bin = '+strn(i)
              rmatch = self->rand_match(lenses, rlenses, $
				                        nozhist=nozhist, input_index=use, $
                                        zbinsize=zbinsize, title=title)

              IF rmatch[i] NE -1 THEN BEGIN 
                  rkeep[i] = ptr_new(rmatch, /no_copy)
              ENDIF 
          ENDIF 
      ENDIF 
  ENDFOR 

  !p.title = ''
  return, rkeep

END 


FUNCTION objshear::combine_sub_samples, lensum, keep, $
                 hdr=hdr, $
                 struct=struct, average_tags=average_tags


  nbin = n_elements(keep)
;  IF nbin EQ 1 THEN BEGIN 
;      combine_structs = ptr_new()
;  ENDIF ELSE BEGIN 
;      combine_structs = ptrarr(nbin)
;  ENDELSE 

  nlensum = n_elements(lensum)
  nstruct = n_elements(struct)
  ntags = n_elements(average_tags)

  IF (nstruct NE 0 AND ntags NE 0) AND nlensum NE nstruct THEN BEGIN 
      message,'lensum and struct have different number of objects!'
  ENDIF 

  FOR i=0L, nbin-1 DO BEGIN 
      
      IF ptr_valid(keep[i]) THEN BEGIN 
          IF n_elements(*keep[i]) NE 0 THEN BEGIN 
              print,'------------------------------------------------'
              print,'Sub sampling bin: #'+ntostr(i+1)
              print

              combstruct = combine_lensum(lensum, hdr=hdr, index=*keep[i])

              ;; Add some averages 
              IF nstruct NE 0 AND ntags NE 0 THEN BEGIN 
                  meanstruct = self->meanstruct(lensum, struct, average_tags, $
                                                input_index=*keep[i])

                  combstruct = create_struct(combstruct, meanstruct)
              ENDIF 

              IF n_elements(combine_structs) EQ 0 THEN BEGIN 
                  combine_structs = replicate(combstruct, nbin)
              ENDIF 
              combine_structs[i] = combstruct
;              combine_structs[i] = ptr_new(combstruct, /no_copy)
          ENDIF 
      ENDIF 

  ENDFOR 
  return, combine_structs

END 



PRO objshear::write_sub_samples, structs, outfiles

  nbin = n_elements(structs)
  IF n_elements(outfiles) NE nbin THEN BEGIN
      message,'outfiles must be same size as combine_structs'
  ENDIF 

  print
  FOR i=0L, nbin-1 DO BEGIN 

      st = structs[i]

      print,'Writing to file: ',outfiles[i]
      IF tag_exist(st, 'jackknife_id') THEN BEGIN 
          hdr={jackknife_ids: st.jackknife_id}
          write_idlstruct, st, outfiles[i], hdr=hdr
      ENDIF ELSE BEGIN 
          write_idlstruct, st, outfiles[i]
      ENDELSE 
  ENDFOR 
  
END 

;; general program for sub sampling.  Requires general naming scheme
;; for file and selection routines so that we can use the execute()
;; function. 


PRO objshear::sub_sample, subtype, randnum=randnum, $
            nozhist=nozhist, zbinsize=zbinsize

  IF n_elements(subtype) EQ 0 THEN BEGIN 
      print,'-syntax: o->sub_sample, subtype, randnum=, /nozhist, zbinsize='
      message,'You must enter the subtype string'
  ENDIF 

  ;; main catalog
  struct = self->get()
  
  ;; only need those that made it into the lens input catalogs
  lensin = self->lensread('input')
  zindex = lensin.zindex
  struct = struct[zindex]

  print,'running subselect on lens sample'
  if self->select_type(subtype) eq 'alt' then begin
      keep = self->altselect(struct, subtype, nkeep=nkeep, nbin=nbin, $
          average_tags=average_tags)
  endif else begin
      keep = self->where_select(struct, subtype, nkeep=nkeep, nbin=nbin, $
          average_tags=average_tags, /nodisplay)
  endelse

  nrand = n_elements(randnum)  
  if nrand eq 0 then begin 

      ;; lens output file
      lensum = self->lensread('output', hdr=hdr)
      hdr = idlstruct_hclean(hdr)

      ;; will autumatically select all bins
      outfiles = self->lensfile('combined', subtype=subtype, /createdir)
      ;; combine the sub samples
      combine_structs = $
        self->combine_sub_samples(lensum, keep, hdr=hdr, $
                                  struct=struct, $
                                  average_tags=average_tags)

      ;; write the outputs
      self->write_sub_samples, combine_structs, outfiles

  endif else begin 

	  print,'Doing matching to random'
      ;; don't overwrite


	  psdir = self->plot_dir(subtype=subtype,/createdir)
      for ri=0l, nrand-1 do begin 

		  psfile = $
			  self->type()$
			  +'_'+self->sample_string()$
			  +'_'+subtype+'_matchzrand'$
			  +string(randnum[ri],f='(i02)')+'.ps'

		  psfile=path_join(psdir, psfile)
		  begplot,psfile,/color

          ;; automatically selects all bins
          outfiles = self->randnumfile(randnum[ri], subtype=subtype, $
                                       /createdir)

          self->read_random_files, randnum[ri], rlensin, rlensum, hdr, $
            nozhist=nozhist
          nrlensin = n_elements(rlensin) & nrlensum = n_elements(rlensum)
          if nrlensin ne nrlensum then begin 
              message,'rlensin and rlensum have different number of objects!'
          endif 


          ;; Send lensin, may need redshifts
          rkeep = self->rand_match_bins(lensin,rlensin,keep,nozhist=nozhist,$
                                        zbinsize=zbinsize)
          rlensin = 0

          ;; combine the sub samples
          combine_structs = $
            self->combine_sub_samples(rlensum, rkeep, hdr=hdr)
          
          ;; write the outputs
          self->write_sub_samples, combine_structs, outfiles
          

          ;; free memory
          rlensum=0

          ptr_free, rkeep

		  endplot

      endfor 

  endelse 

  ptr_free, keep



END 

;; must be able to subscript both with same index
;; so we must have performed a struct=struct[lensum.zindex]
FUNCTION objshear::meanstruct, lensum, struct, tags, $
                 input_index=input_index, status=status

  IF (n_elements(lensum) EQ 0 OR $
      n_elements(struct) EQ 0 OR $
      n_elements(tags) EQ 0) THEN BEGIN 
      message,'-Syntax: meanstr = oshear->meanstruct(lensum, struct, tags, input_index=, status=)'
  ENDIF 

  ntags = n_elements(tags)

  status = 1
  nlensum = n_elements(lensum)
  IF n_elements(input_index) EQ 0 THEN w = lindgen(nlensum) ELSE w=input_index

  wtot = total( lensum[w].weight )

 
  FOR i=0L, ntags-1 DO BEGIN 

      ;; strip trailing [index] or (index)
      bracket_pos  = stregex(tags[i], '[\[\(]')
      bracket_pos2 = stregex(tags[i], '[]\)]')
      IF bracket_pos NE -1 THEN BEGIN 

          tag = strmid(tags[i], 0, bracket_pos)

          ;; Now extract the element for our name
          element_str = '_'+strmid(tags[i], $
                                   bracket_pos+1, bracket_pos2-bracket_pos-1)
      ENDIF ELSE BEGIN 
          tag = tags[i]
          element_str = ''
      ENDELSE 

      IF tag_exist(struct[0], tag, index=itag) THEN BEGIN 

          mean_command = $
            'mean_tag = total( lensum[w].weight*struct[w].'+tags[i]+' )/wtot'
          sdev_command = $
            'sdev_tag = total( lensum[w].weight*(struct[w].'+tags[i]+'-mean_tag)^2 )/wtot'
          err_command = $
            'err_tag = sqrt( total( lensum[w].weight^2*(struct[w].'+tags[i]+'-mean_tag)^2 )/wtot^2 > 0)'

          IF NOT execute(mean_command) THEN message,'Unable to take mean of tag '+tags[i]
          IF NOT execute(sdev_command) THEN message,'Unable to take sdev of tag '+tags[i]
          IF NOT execute(err_command)  THEN message,'Unable to take err  of tag '+tags[i]

          IF n_elements(meanstruct) EQ 0 THEN BEGIN 
              meanstruct = create_struct('mean_'+tag+element_str, mean_tag, $
                                         'sdev_'+tag+element_str, sdev_tag, $
                                         'err_' +tag+element_str, err_tag)
          ENDIF ELSE BEGIN 
              meanstruct = create_struct(meanstruct, $
                                         'mean_'+tag+element_str, mean_tag, $
                                         'sdev_'+tag+element_str, sdev_tag, $
                                         'err_' +tag+element_str, err_tag)
          ENDELSE 

      ENDIF ELSE BEGIN 
          message,'Tag '+ntostr(tag)+' does not exist in structure',/inf
      ENDELSE 

  ENDFOR 

  IF n_elements(meanstruct) NE 0 THEN BEGIN 
      status = 0
      return, meanstruct
  ENDIF ELSE BEGIN 
      return, -1
  ENDELSE 

END 





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Combine all randoms
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO objshear::combine_allrand, subtype=subtype

  randnum = self->randnum()

  IF n_elements(subtype) EQ 0 THEN BEGIN 

      outfile = self->lensfile('random')
      sharray = self->randnumread(randnum)
      shave = combine_shstructs(sharray)
      
      print
      print,'Writing to file: ',outfile
      write_idlstruct, shave, outfile

  ENDIF ELSE BEGIN 
  
      nbin = self->subtype_nbin(subtype)
      bins = lindgen(nbin)

      FOR i=0L, nbin-1 DO BEGIN 
          
          outfile = self->lensfile('random',subtype=subtype, bin=bins[i])

          sharray = self->randnumread(randnum, $
                                      subtype=subtype, bin=bins[i])
          shave = combine_shstructs(sharray)
          
          print
          print,'Writing to file: ',outfile
          write_idlstruct, shave, outfile
      ENDFOR 


  ENDELSE 



END 















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Correcting profiles
; You must run combine_allrand before this program.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO objshear::corr, subtype=subtype


  IF n_elements(subtype) NE 0 THEN BEGIN 
      nbin = self->subtype_nbin(subtype)
      bin = lindgen(nbin)
  ENDIF 

  sharr    = self->lensread('combined', subtype=subtype)
  rsharr   = self->lensread('random',   subtype=subtype)

  outfiles = self->lensfile('corr',     subtype=subtype, /createdir)

  print,''
  print,'-----------------------------------------------------------------'


  num = n_elements(sharr)
  FOR i=0L, num-1 DO BEGIN 

      sh = sharr[i]
      rsh = rsharr[i]
      outfile = outfiles[i]


      ;; Subtract the random.  
      ;; NOTE: This must be done *before* the clustering correction!
      
      ;; only correct where we have a good measurement.
      w = where(sh.meanr/1000 GT 1, nw)

      ;; kludge of this.  Should really re do the other samples with
      ;; new random matching
      IF self->type() EQ 'maxbcg' AND self->sample() LT 11 THEN BEGIN 
          message,'should get maxbcg random stuff working for other samples'
      ENDIF ELSE BEGIN 
          IF nw NE 0 THEN BEGIN 
              sh.sigma[w] = sh.sigma[w] - rsh.sigma[w]
          ENDIF ELSE BEGIN 
              message,'radius cut failed'
          ENDELSE 
      ENDELSE 
      


      ;; Clustering correction + shear polarizability
      ;; This will preserve allthe info in the original struct
      shc = correct_shear(sh, rsh)

      print,'Writing to file: ',outfile

      write_idlstruct, shc, outfile, hdr=hdr
  ENDFOR 

END 







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Jackknifing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; type, nbin are optional
PRO objshear::jackknife, subtype=subtype

  lensin = $
    self->lensread('input', columns=['zindex','clambda','ceta'])
  lensum = $
    self->lensread('output',$
                   columns=['rsum','sigma','orthosig','wsum','owsum','weight'])
  
  IF n_elements(subtype) EQ 0 THEN BEGIN 
      ;; jackknifing the overall sample
      sh = self->lensread('corr')
      rsh = self->lensread('random')
      outfile = self->lensfile('jackknife', /createdir)

      ;; will inherit the quantities from sh, such 
      ;; as averages
      jackstruct = self->jackknife_sample(lensin, lensum, sh, rsh)

      print
      print,'Writing to jackknife file: ',outfile
      write_idlstruct, jackstruct, outfile

  ENDIF ELSE BEGIN 

      sharr = self->lensread('corr', subtype=subtype)
      rsharr = self->lensread('random', subtype=subtype)
      outfiles = self->lensfile('jackknife', subtype=subtype, /createdir)

      struct = self->get()
      struct = struct[lensin.zindex]

      print
      print,'getting keep array'
      if self->select_type(subtype) eq 'alt' then begin
          keeparray = self->altselect(struct, subtype, nkeep=nkeep, nbin=nbin)
      endif else begin
          keeparray = self->where_select(struct, subtype, nkeep=nkeep, nbin=nbin)
      endelse
      struct = 0



      FOR i=0L, nbin-1 DO BEGIN 

          print,'-------------------------------------------------------------'
          print,'Will write to file: ',outfiles[i]

          keep = *keeparray[i]
          sh = sharr[i]
          rsh = rsharr[i]

          jackstruct = self->jackknife_sample(lensin[keep],lensum[keep],sh,rsh)
          print
          print,'Writing to jackknife file: ',outfiles[i]
          ;;IF fexist(outfiles[i]) THEN message,'Already exists!'
          write_idlstruct, jackstruct, outfiles[i]


      ENDFOR 

      ptr_free, keeparray

  ENDELSE 


END 




; We enter sh because it will have many things calculated which will
; be outside the scope of this program
FUNCTION objshear::jackknife_sample, lensin, lensum, sh_in, rsh

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: jack = ml->jackknife_sample(lensin,lensum,sh,rsh)'
      return,-1
  ENDIF 

  sh = sh_in

  par = self->objshear::par_struct()

  print
  print,'Creating samples'

  IF par.shape_correction_style EQ 3 THEN BEGIN 
      jackknife_file = esheldon_config('princeton_jackknife_file')
  ENDIF 
  jackstruct = create_shearjackknife_samples(lensin, lensum, $
                                             jackknife_ids=jackknife_ids,$
                                             jackknife_file=jackknife_file)
  
  print
  print,'Jackknifing'
  wjackknife, jackstruct.sigsum_sub, jackstruct.wsum_sub, $
    sigma, sigmaerr, covariance
  wjackknife, jackStruct.osigsum_sub, jackStruct.owsum_sub, $
    orthosig, orthosigerr, orthocov

  corr_matrix = cov2corr(covariance, status=status)
  IF status NE 0 THEN message,'Failed to calculate corr matrix for sigma'
  orthocorr_matrix = cov2corr(orthocov,status=status)
  IF status NE 0 THEN message,'Failed to calculate corr matrix for ortho'
  
  
  ;; correct and copy in.  Subtraction must occur first.  See ::corr
  w = where(sh.meanr/1000 GT 1, nw)
  IF nw NE 0 THEN sigma[w] = sigma[w] - rsh.sigma[w]

  sh.sigma = sigma*sh.corr/sh.ssh
  sh.sigmaerr = sigmaerr*sh.corr/sh.ssh

  sh.orthosig = orthosigerr
  sh.orthosigerr = orthosigerr

  ;; This will properly get the errors into the covariance matrix
  covout = corr2cov(corr_matrix, sh.sigmaerr)
  
  jackstruct = create_struct(sh, $
                             'covariance', covout, $
                             'correlation', corr_matrix, $
                             'orthocov', orthocov, $
                             'orthocorrelation', orthocorr_matrix)

  return,jackstruct

END 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Invert to 3D profile
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function objshear::m200interp, rmpc, mass, mass_err, meanz, doplot=doplot

  ;; Units solar masses/Mpc^3
  rhocrit = rhocrit(z=meanz, omega_m=omega_m)

  rhocum = mass/(4d/3d*!dpi*rmpc^3)
  rhocum_err = mass_err/(4d/3d*!dpi*rmpc^3)

  rhorel = rhocum/rhocrit
  rhorel_err = rhocum_err/rhocrit
  
  
  ;; Heavily smooth and make sure not negative
  wbad = where(rhocum lt 0,nbad, comp=wgood)
  if nbad ne 0 then begin 
      rhocum[wbad] = interpol( rhocum[wgood], rmpc[wgood], rmpc[wbad] )
  endif 
    
  ;; Now look for where meanrho crosses 200
  j=0 & nr=n_elements(rhocum)
  while (rhorel[j] gt 200) and (j lt nr) do begin 
      j = j+1
  endwhile 
  
  if j gt 0 and j lt nr then begin 
      ;; First interpolate radius from rho
      r200 = interpol(rmpc, rhorel, 200.0)
      
      ;; Now interpolate mass/err from radius
      m200 = interpol(mass, rmpc, r200)
      m200_err = interpol(mass_err, rmpc, r200)
      
  endif else begin
      r200 = -9999d
      m200 = -9999d
      m200_err = 9999d
  endelse 

  st = { $
         rhocrit: rhocrit, $
         rhocum: rhocum, $
         rhocum_err: rhocum_err, $
         r200: r200, $
         m200: m200, $
         m200_err: m200_err $
       }

  return, st

end 

;; Mass/rhocrit in solar masses
;; Needs to be replaced with a reasonable function rather than a 
;; power law.
function objshear::m200fit, rmpc, mass, mass_err, meanz, doplot=doplot

  par = self->par_struct()
  omega_m = par.omega_m

  ;; Units solar masses/Mpc^3
  rhocrit = rhocrit(z=meanz, omega_m=omega_m)

  rhocum = mass/(4d/3d*!dpi*rmpc^3)
  rhocum_err = mass_err/(4d/3d*!dpi*rmpc^3)

  rhorel = rhocum/rhocrit
  rhorel_err = rhocum_err/rhocrit

  aguess = [1.0d, -2d]
  fitst = fit_r200power(rmpc, rhorel, rhorel_err, aguess,/silent)

  if keyword_set(doplot) then begin 
      pplot, rmpc, rhorel, yerr=rhorel_err, psym=8, $
        /xlog, /ylog, aspect=1, title='Power law r200'
      oplot, rmpc, fitst.rhobyrhocrit_fit, color=c2i('red')
      oplot, [0.001, 100], [200,200]

      oplot, [fitst.r200,fitst.r200],[0.01,1.e20], color=c2i('red')
  endif 

  m200     = 200d*rhocrit*4d/3d*!dpi*fitst.r200^3
  m200_err = 200d*rhocrit*4d/3d*!dpi*2d*fitst.r200^2*fitst.r200_err

  st = { $
         rhocrit: rhocrit, $
         rhocum: rhocum, $
         rhocum_err: rhocum_err, $
         m200: m200, $
         m200_err: m200_err $
       }
  st = create_struct(st, fitst)
  return, st

end 


function objshear::nfwfit, rmpc, mass, mass_err, meanz

  par = self->par_struct()
  omega_m = par.omega_m
  define_nfw_200, z=meanz, omega_m=omega_m
  fit_nfw_mass, rmpc, mass, mass_err, p, perr, $
    mass_yfit = mass_yfit, nfw_yfit=nfw_yfit, bias_yfit=bias_yfit

  ;; bias < 0 try fitting just nfw
  if p[2] lt 0 then begin 
      fit_nfw_mass, rmpc, mass, mass_err, p, perr, $
        mass_yfit = mass_yfit, nfw_yfit=nfw_yfit, bias_yfit=bias_yfit,/positive_bias
  endif 


  r200 = p[0]
  r200_err = perr[0]
  c = p[1]
  c_err = perr[1]
  b = p[2]
  b_err = perr[2]


  rhocrit = rhocrit(z=meanz, omega_m=omega_m)

  m200 = 200d*rhocrit*4d/3d*!dpi*r200^3
  m200_err = 200d*rhocrit*4d/3d*!dpi*2d*r200^2*r200_err

  oplot, [r200], [m200], psym=2, color=c2i('red'), symsize=1.5

  st = { $
         rhocrit: rhocrit, $
         r200: r200, $
         r200_err: r200_err, $
         m200: m200, $
         m200_err: m200_err, $
         c: c, $
         c_err: c_err, $
         b: b, $
         b_err: b_err, $
         mass_yfit: mass_yfit, $
         nfw_yfit: nfw_yfit,$
         linear_yfit: bias_yfit $
       }

  return, st

end 


PRO objshear::invert, subtype=subtype, corrected=corrected, dops=dops

  ;; This will be jackknife when we get that done
  if keyword_set(corrected) then begin 
      dtype = 'corr'
  endif else begin 
      dtype = 'jackknife'
  endelse 

  out_dtype = dtype + '_invert'
  t = self->lensread(dtype, subtype=subtype)
  
  psfile = self->plotfile(out_dtype, subtype=subtype)
  psfile = repstr(psfile, out_dtype, out_dtype+'_qa')
  if keyword_set(dops) then begplot, psfile, /color

  outfiles = self->lensfile(out_dtype, subtype=subtype, /createdir)

  nbin = n_elements(t)
  nrad = n_elements(t[0].meanr)

  nstr = ntostr(nrad-1)
  arrval = 'dblarr('+nstr+')'
  cov_arrval = 'dblarr('+nstr+','+nstr+')'
  newtags = ['ir',$
             'drho',    'drho_err',    'drho_cov', $
             'massin',  'massin_err',  'massin_cov', $
             'massout', 'massout_err', 'massout_cov', $
             'masscomb', 'masscomb_err', 'masscomb_cov', $
             'rhocum', 'rhocum_err', $
             'rhocrit', $
             'r200', 'r200_err', 'm200', 'm200_err', 'c','c_err','b','b_err',$
             'r200_1', 'm200_1', 'm200_1_err', $
             'r200_2', 'r200_2_err', 'm200_2', 'm200_2_err', $
             'vel', 'vel_err', 'vel_cov', $
             'mass_yfit', 'nfw_yfit', 'bias_yfit']
  tagvals = [arrval, $
             arrval, arrval, cov_arrval, $
             arrval, arrval, cov_arrval, $
             arrval, arrval, cov_arrval, $
             arrval, arrval, cov_arrval, $
             arrval, arrval, $
             '-9999d', $
             replicate('-9999d', 8), $
             '-9999d', '-9999d', '-9999d', $
             '-9999d', '-9999d', '-9999d', '-9999d', $
             arrval, arrval, cov_arrval,$
             arrval, arrval, arrval]
  add_tags, t, newtags, tagvals, outst
    
  dp = obj_new('deproject')

  print,'Inverting'
  FOR i=0L, nbin-1 DO BEGIN 

      r = t[i].meanr/1000.0

      dsigma = t[i].sigma

      IF keyword_set(corrected) THEN BEGIN 
          dsigma_cov = diagonal_array(t[i].sigmaerr)
      ENDIF ELSE BEGIN 
          dsigma_cov = t[i].covariance
      ENDELSE 

      ist = dp->invert(r, dsigma, dsigma_cov)

      outst[i].ir = ist.ir

      outst[i].drho = ist.drho
      outst[i].drho_err = ist.drho_err
      outst[i].drho_cov = ist.drho_cov

      outst[i].massin = ist.massin
      outst[i].massin_err = ist.massin_err
      outst[i].massin_cov = ist.massin_cov

      outst[i].massout = ist.massout
      outst[i].massout_err = ist.massout_err
      outst[i].massout_cov = ist.massout_cov

      outst[i].masscomb = ist.masscomb
      outst[i].masscomb_err = ist.masscomb_err
      outst[i].masscomb_cov = ist.masscomb_cov

      outst[i].vel = ist.vel
      outst[i].vel_err = ist.vel_err
      outst[i].vel_cov = ist.vel_cov

      key = prompt_kbrd('hit a key')
      IF key EQ 'q' THEN return

  ENDFOR 

  ;; Run the mass fits separately so we can see
  ;; the plots all in a row

  print,'Fitting NFW+Bias'
  for i=0L, nbin-1 do begin 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Calculate r200,m200
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


      fitst = self->nfwfit(outst[i].ir, outst[i].massout, outst[i].massout_err, t[i].mean_z)

      fitst1 = self->m200interp(outst[i].ir, $
                                outst[i].massout*1.e12, $
                                outst[i].massout_err*1.e12, t[i].mean_z)
      fitst2 = self->m200fit(outst[i].ir, $
                             outst[i].massout*1.e12, $
                             outst[i].massout_err*1.e12, t[i].mean_z)

      
      outst[i].r200 = fitst.r200
      outst[i].r200_err = fitst.r200_err
      outst[i].m200 = fitst.m200
      outst[i].m200_err = fitst.m200_err
      outst[i].c = fitst.c
      outst[i].c_err = fitst.c_err
      outst[i].b = fitst.b
      outst[i].b_err = fitst.b_err
      outst[i].mass_yfit = fitst.mass_yfit
      outst[i].nfw_yfit = fitst.nfw_yfit
      outst[i].bias_yfit = fitst.bias_yfit


      outst[i].r200_1 = fitst1.r200
      outst[i].m200_1 = fitst1.m200
      outst[i].m200_1_err = fitst1.m200_err

      outst[i].r200_2 = fitst2.r200
      outst[i].r200_2_err = fitst2.r200_err
      outst[i].m200_2 = fitst2.m200
      outst[i].m200_2_err = fitst2.m200_err

      outst[i].rhocum = fitst2.rhocum
      outst[i].rhocum_err = fitst2.rhocum_err
      outst[i].rhocrit = fitst2.rhocrit

  endfor 

  simpctable, colorlist=clist
  rho = outst[0].rhocum/outst[0].rhocrit
  rho_err = outst[0].rhocum_err/outst[0].rhocrit

  yt = '<'+!csym.rho+'>/'+!csym.rho+'!Dc!N'
  yt = textoidl('<\rho>/\rho_c')
  pplot, outst[0].ir, rho, yerr=rho_err, $
    psym=8, aspect=1, xtitle=estitle('mpcxtitle3d'), ytitle=yt, $
    /xlog, /ylog
;  oplot, [outst[0].r200], [outst[0].m200

  for i=1L, nbin-1 do begin 
      rho = outst[i].rhocum/outst[i].rhocrit
      rho_err = outst[i].rhocum_err/outst[i].rhocrit
      pplot, outst[i].ir, rho, $
        color=clist[i], /overplot
  endfor 

  if nbin gt 1 then begin 
      ws = self->where_string(subtype, labels=labels)
      legend, labels, /right, line=0, box=0, charsize=0.7, $
        color=clist[0:nbin-1]
  endif 
  oplot, [0.001, 100], [200,200]

  obj_destroy, dp

  if keyword_set(dops) then endplot

  self->write_sub_samples, outst, outfiles

END  












;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Comparing to the intrinsic alignments bounds
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; convert the intrinsic signal to the redshift and 
;; mean sigmacrit of the lens sample, including correction
FUNCTION objshear::convert_intrinsic, lsh

  nbin = n_elements(lsh)
  p = obj_new('percolate',2)
      
  par = self->par_struct()
  
  sh_int = p->lensread('combined')

  t = sh_int
  t.meanr = sh_int.meanr
  t.tmeanr = sh_int.tmeanr
  sh_int_convert = replicate(t, nbin)
  
  ;; scale intrinsic shear to correct redshifts

  scinv = sdss_sigma_crit(par.omega_m,lsh.mean_z)  
  FOR i=0L, nbin-1 DO BEGIN 

      
      ;; Conver to sigma at this redshift and 

;      frac = interpol( lsh[i].frac_clust, lsh[i].meanr, sh_int.meanr ) 
      
      ;; This is the fraction at that redshift times the boost factor
      cMinusOne = lsh[i].corr-1.0
      frac  = interpol( cMinusOne, lsh[i].meanr, sh_int.meanr ) 

      sh_int_convert[i].sigma = sh_int.sigma/scinv[i]*frac
      sh_int_convert[i].sigmaerr = sh_int.sigmaerr/scinv[i]*abs(frac)


      w = where(sh_int_convert[i].sigmaerr GT 0)

      wts = 1.0/sh_int_convert[i].sigmaerr[w]^2
      tsum = total( sh_int_convert[i].sigma[w]*wts, /cumulative )
      wsum = total( wts, /cumulative)

      sh_int_convert[i].tsigma[w] = tsum/wsum
      sh_int_convert[i].tsigmaerr[w] = sqrt(1.0/wsum)


  ENDFOR 
  
  obj_destroy,p

  return,sh_int_convert

END 



PRO objshear::compare_intrinsic, subtype, cumulative=cumulative, $
            dops=dops, color=color

  subtype = 'ngals200_12'

  IF keyword_set(dops) THEN BEGIN 
      plot_dir = self->plot_dir(subtype=subtype,/createdir)
      sstr = self->sample_string()
      file = strlowcase(subtype)+'_compare_intrinsic_'+sstr

      IF keyword_set(cumulative) THEN BEGIN 
          file = file + '_cumulative'
      ENDIF 

      IF keyword_set(color) THEN BEGIN 
          file = file+'_color'
      ENDIF 

      file = file + '.eps'

      file = concat_dir(plot_dir, file)

      begplot, file, /color, /encapsulated

  endif 

  
  ws = self->where_string(subtype, labels=labels)

  sharr = self->lensread('jackknife', subtype=subtype)

  n =  n_elements(sharr)
  nbin = self->subtype_nbin(subtype)

  IF nbin NE n THEN message,'nbin != n'

  IF n_elements(subtype) NE 0 THEN BEGIN 
;      mess = 'bin'+ntostr(1+lindgen(nbin))
      mess = labels
  ENDIF 


  !p.charsize=0.9
  lcharsize=0.6
  thick=2
  !p.thick=thick
  !p.charthick=thick
  !x.thick=thick
  !y.thick=thick
  
  sh_int = self->convert_intrinsic(sharr)

  xtitle = estitle('mpcxtitle2')
  ytitle = estitle('deltaytitle')

  IF keyword_set(color) THEN BEGIN 
      oclr = c2i('red')
  ENDIF ELSE BEGIN 
      oclr = !p.color
  ENDELSE 

  oline = 1

  IF keyword_set(cumulative) THEN ytitle = 'Cumulative '+ytitle

  erase & multiplot, [4,3], /square, $
      mxtitle=xtitle, mytitle=ytitle, $
      mytitoffset=1, $
      xtickf='loglabels', $
      xgap=0.03

  FOR i=0L, nbin-1 DO BEGIN 
      
      IF keyword_set(cumulative) THEN BEGIN 
          rad = sharr[i].tmeanr/1000
          sigma = sharr[i].tsigma
          sigmaerr = sharr[i].tsigmaerr

          irad = sh_int[i].tmeanr/1000
          isigma = sh_int[i].tsigma
          isigmaerr = sh_int[i].tsigmaerr
      ENDIF ELSE BEGIN 
          rad = sharr[i].meanr/1000
          sigma = sharr[i].sigma
          sigmaerr = sharr[i].sigmaerr

          irad = sh_int[i].meanr/1000
          isigma = sh_int[i].sigma
          isigmaerr = sh_int[i].sigmaerr
      ENDELSE 

      yrange = prange(sigma, isigma, sigmaerr, isigmaerr)
     
      pplot, rad, sigma, yerr=sigmaerr, $
        /xlog, xrange=xrange, yrange=yrange, ystyle=3, $
        xticklen=0.04, yticklen=0.04, $
        hat=0

      oplot, [0.001, 100], [0,0]
      pplot, irad, isigma, yerr=isigmaerr, /overplot, color=oclr, hat=0, $
        line=oline


      IF nbin GT 0 THEN BEGIN 
          legend, mess[i], /right, box=0, charsize=lcharsize, margin=0
      ENDIF 

      IF i EQ -1 THEN BEGIN 

          legend, $
              ['lenses', 'intrinsic'], $
              line=[0, oline], $
              color=[!p.color, oclr], $
              /right, /bottom
      ENDIF 

      multiplot, /doyaxis

  ENDFOR 
  multiplot,/default


  print
  print,'---------------------------------------------------------------'
  print,'% cumulative error limites within 100 kpc'
  print

  IF subtype EQ 'ngals200_12' THEN self->ngals200_bins, 12, lowlim, highlim

  FOR i=0L, n-1 DO BEGIN 

      
      ;; interpolate to 100 kpc
      limrad = 100                 ;kpc

      sint = interpol(sh_int[i].tsigma, sh_int[i].tmeanr, limrad)
      sinterr = interpol(sh_int[i].tsigmaerr, sh_int[i].tmeanr, limrad)

      sintlow  = sint - 2.*sinterr
      sinthigh = sint + 2.*sinterr
      

      smeas = interpol(sharr[i].tsigma, sharr[i].tmeanr, limrad)
      smeaserr = interpol(sharr[i].tsigmaerr, sharr[i].tmeanr, limrad)

      IF i NE n-1 THEN cont = '\\\\' ELSE cont = ''

      gam = !csym.gamma+'!DT!N'

      IF n_elements(lowlim) EQ 0 THEN BEGIN 
          format='(%"%d & [%10.3g, %10.3g] & %10.3g $\\pm$ %10.3g '+cont+'")'
          print, format=format, $
            i+1, sintlow, sinthigh, smeas, smeaserr
      ENDIF ELSE BEGIN 

          IF lowlim[i] EQ highlim[i] THEN BEGIN 
              label = 'N_{200} = '+ntostr(lowlim[i])
          ENDIF ELSE BEGIN 
              label = ntostr(lowlim[i])+' \le N_{200} \le '+ntostr(highlim[i])
          ENDELSE 
            
          format='(%"$%s$ & [%10.3g, %10.3g] & %10.3g $\\pm$ %10.3g '+cont+'")'
          print, format=format, $
            label, sintlow, sinthigh, smeas, smeaserr
      ENDELSE 

  ENDFOR 
  print


  if keyword_set(dops) then endplot,/trim

END 












;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plot the positions of lenses and random points
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION objshear::ceta2region, ceta

  stripes = eta2stripenum(ceta)

  regions = bytarr(n_elements(ceta))

  w=where(stripes LE 39, nreg)
  regions[w] = 1

  w=where(stripes GT 39, nreg)
  IF nreg NE 0 THEN regions[w] = 2
  return,regions

END 



FUNCTION objshear::randind, num_orig, fracuse=fracuse

  IF n_elements(fracuse) EQ 0 THEN fracuse = 1.0/10.0

  nkeep = long(num_orig*fracuse)

  keep = long( arrscl(randomu(seed, nkeep), 0, num_orig-1, $
                      arrmin=0.0, arrmax=1.0) )
  
  return,keep

END 


PRO objshear::plot_lensrand, l=l, r=r, dops=dops, dopng=dopng, $
            psym_lens=psym_lens, symsize_lens=symsize_lens, $
            fracuse=fracuse

  IF n_elements(l) EQ 0 OR n_elements(r) EQ 0 THEN BEGIN 

      l = self->lensread('input') 
      
      self->read_random_files, 0, r, rout, /nozhist 
      rout = 0

      nl=n_elements(l)
      nr=n_elements(r)

      IF n_elements(fracuse) NE 0 THEN BEGIN 
          r = r[self->randind(nr,frac=fracuse)]
      ENDIF 
  ENDIF 

  plot_dir = self->plot_dir(/createdir)

  sampstring = self->sample_string()
  pngfile = concat_dir(plot_dir,'lensrand_positions_'+sampstring+'.png')
  psfile = concat_dir(plot_dir,'lensrand_positions_'+sampstring+'.eps')

  IF keyword_set(dops) THEN BEGIN 
      begplot, name=psfile, /color, /landscape
  ENDIF 

  ytitle = textoidl('\eta_c')
  xtitle = textoidl('\lambda_c')
  xrange = [-70, 70]
  botyrange = [-40,50]
  topyrange = [120,170]
  oplotcolor = c2i('red')

  regions = self->ceta2region(l.ceta)
  rregions = self->ceta2region(r.ceta)

  w1 = where(regions EQ 1, nw1)
  w2 = where(regions EQ 2, nw2)
  rw1 = where(rregions EQ 1, nrw1)
  rw2 = where(rregions EQ 2, nrw2)


  IF n_elements(psym_lens) EQ 0 THEN psym_lens=3
  IF nrw1 NE 0 AND nrw2 NE 0 THEN BEGIN 

      plottwo, $
        r[rw2].clambda, r[rw2].ceta, $
        r[rw1].clambda, r[rw1].ceta, $
        toppsym=3, botpsym=3, $
        botyrange=botyrange, topyrange=topyrange, $
        xrange=xrange, xstyle=1, $
        xtitle = xtitle, topytitle=ytitle, botytitle=ytitle, $
        frac1 = 0.25, $
        $
        xoplot = l[w2].clambda, yoplot = l[w2].ceta, $
        oplotsym=3, oplotcolor=oplotcolor, $
        /ynozero, ystyle=1+2
      
      oplot,l[w1].clambda,l[w1].ceta,psym=psym_lens,symsize=symsize_lens,color=oplotcolor
  ENDIF ELSE BEGIN 

      plot,r.clambda,r.ceta,$
        /iso, $
        psym=3,$
        /ynozero, $
        title=title, xtitle=xtitle, ytitle=ytitle, xrange=xrange, xstyle=1
      
;      oplot,l.clambda,l.ceta,psym=8,color=oplotcolor,symsize=0.25
      oplot,l.clambda,l.ceta,color=oplotcolor, psym=psym_lens, symsize=symsize_lens
  ENDELSE 

  legend,$
    ['lens','random'],$
    psym=[8,8],color=[oplotcolor,!p.color],/right,box=0,charsize=1.5

  IF keyword_set(dops) THEN BEGIN 
      endplot, /landfix
  ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
      print
      print,'Writing to file: ',pngfile
      write_png, pngfile, tvrd(/true)
  ENDIF 

END 

PRO objshear::plot_aitoff, l=l, dops=dops, dopng=dopng, sample=sample, $
  colorps=colorps

  if n_elements(l) eq 0 then l=self->lensread('input')

  psfile  = 'lens_aitoff'
  pngfile = psfile

  IF keyword_set(sample) THEN BEGIN 
      plot_dir=self->plot_dir(/createdir)
      psfile  = concat_dir(plot_dir,psfile  + '_'+self->sample_string()+'.eps')
      pngfile = concat_dir(plot_dir,pngfile + '_'+self->sample_string()+'.png')
      IF n_elements(l) EQ 0 THEN l = self->lensread('input')
  ENDIF ELSE BEGIN 
      plot_dir=self->plot_dir(/base,/createdir)
      psfile  = concat_dir(plot_dir,psfile  + '_'+self->catalog()+'.eps')
      pngfile = concat_dir(plot_dir,pngfile + '_'+self->catalog()+'.png')
      IF n_elements(l) EQ 0 THEN l = self->get()
  ENDELSE 

  IF keyword_set(colorps) THEN psfile = repstr(psfile, '.eps', '_color.eps')

  IF keyword_set(dops) THEN BEGIN 
      psfile = psfile
      begplot,name=psfile, /color, /encap, xsize=8.5, ysize=11
      fracuse = 0.02

      !p.charsize = 1
      !p.thick=2
      !x.thick=2
      !y.thick=2
      !p.charthick=2

      IF NOT keyword_set(colorps) THEN BEGIN 
          label_color = c2i('grey50')
      ENDIF ELSE BEGIN 
          label_color = c2i('red')
      ENDELSE 
  ENDIF ELSE BEGIN 
      
      label_color = c2i('red')
      fracuse = 0.05
  ENDELSE 
  grid_color = c2i('grey20')

  simpctable


  plot, [0], xrange=[-200, 200], yrange=[-100, 100], $
        position = aspect(1.0/2.0), $
        xstyle=4, ystyle=4
;  plot, [0], xrange=[-200, 200], yrange=[-100, 100], $
;        /iso;, 
;        xstyle=4, ystyle=4

  ang2pix,l.clambda,l.ceta,ix,iy,pixnum,resolution=16, /survey
  rmd = rem_dup(pixnum)
  display_pixel, pixnum[rmd], /aitoff, resolution=16, /fill, /radec, $
                 /over_plot,color=!p.color

;  myaitoff, l.ra, l.dec, x, y
;  plotrand, x, y, fracuse=fracuse, $
;            xrange = [-200,200], yrange=[-90, 90], psym=3,$
;            xstyle=4, ystyle=4, position=aspect(1.0/!gratio)

  myaitoff_grid, label=2, color=grid_color, lcolor=label_color

  IF keyword_set(dops) THEN BEGIN 
      endplot, /trim_bbox
  ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
      write_png, pngfile, tvrd(/true)
  ENDIF 

END 






































;; Mean critical density from Carlos' deconvolutions
FUNCTION objshear::sigmacrit, photoz=photoz

  par_struct = self->par_struct()

  nzlens = 1000
  minzlens = 0.02
  maxzlens = 1.0
  zlens = arrscl(findgen(nzlens), minzlens, maxzlens)

  scritinv = sdss_sigma_crit(par_struct.omega_m, zlens, $
                             h = par_struct.h, $
                             scinvstruct=scinvstruct, $
                             photoz=photoz)

  

  return,scinvstruct

END 


PRO objshear::sigmacrit_plot, dops=dops, compare=compare


  pz = obj_new('photoz_uchicago')
  pzst = pz->deconv_read()
  psdir = pz->plot_dir(/createdir)
  obj_destroy,pz

  sc = self->sigmacrit_calc()

  IF keyword_set(dops) THEN BEGIN 
      name = concat_dir(psdir,'sigmacrit_inv.eps')
      begplot, name=name, /encap
      oclr = !p.color
  ENDIF ELSE oclr = c2i('green')


  ytitle = textoidl('\Sigma^{-1}_{crit} [10^{-4} pc^2/M_{\sun}]')
  xtitle = textoidl('z_{Lens}')
  sig_inv = sc.sig_inv/1.e-4
  aplot, !gratio, sc.zlens, sig_inv, ytitle=ytitle, xtitle=xtitle, line=3

  maxp = max(sig_inv)/2.0
  dndz = pzst.dndz/max(pzst.dndz)*maxp
  oplot, pzst.z, dndz, psym=10, color=oclr

  legend,[siginv, 'dn/dz source'],line=[3,0], color=[!p.color, oclr],$
         /right, box=0, charsize=1.2

  IF keyword_set(dops) THEN endplot, /trim_bbox

END 








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Generic plotting routines
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION objshear::mplot_value, nbin

  CASE nbin OF
      1: message,"nbin=1 doesn't make sense FOR multiplot"
      2: mplot_value=[1,2]
      3: mplot_value=[3,1]
      4: mplot_value=[2,2]
      5: mplot_value=[3,2]
      6: mplot_value=[3,2]
      8: mplot_value=[4,2]
      9: mplot_value=[3,3]
      12: mplot_value=[4,3]
      16: mplot_value=[4,4]
      ELSE: message,'unknown nbin: '+ntostr(nbin)
  ENDCASE 
  return,mplot_value

END 

FUNCTION objshear::axis_labeled, nbin
  CASE nbin OF
      1: BEGIN 
          xlabeled = 1
          ylabeled = 1
      END 
      2: BEGIN 
          xlabeled = [0,1]
          ylabeled = [1,1]
      END 
      3: BEGIN 
          xlabeled = [1,1,1]
          ylabeled = [1,0,0]
      END 
      6: BEGIN 
          xlabeled = [0,0,0,1,1,1]
          ylabeled = [1,0,0,1,0,0]
      END 
      8: BEGIN 
          xlabeled = [0,0,0,0,1,1,1,1]
          ylabeled = [1,0,0,0,1,0,0,0]
      END 
      12: BEGIN 
          xlabeled = [0,0,0,0,0,0,0,0,1,1,1,1]
          ylabeled = [1,0,0,0,1,0,0,0,1,0,0,0]
      END 
      16: BEGIN 
          xlabeled = [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
          ylabeled = [1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0]
      END 
      ELSE: message,'Unsupported nbin: '+ntostr(nbin)
  ENDCASE 
  st = {xlabeled: xlabeled, $
        ylabeled: ylabeled}
  return,st
END 




FUNCTION objshear::titles, type

  CASE type OF
      ELSE: BEGIN 
          xtitle = estitle('mpcxtitle2')
          ytitle = estitle('deltaytitle')
      END 
  ENDCASE 
  st = {xtitle: xtitle, $
        ytitle: ytitle}
  return,st

END 

FUNCTION objshear::typetags, sh, type

  CASE type OF 
      'mass': BEGIN 
          IF (NOT tag_exist(sh, 'mass', index=tag) OR $
              NOT tag_exist(sh, 'mass_err', index=errtag) ) THEN BEGIN 
              message,'Mass tags not found'
          ENDIF 
      END 
      'clustcorr': BEGIN 
          IF (NOT tag_exist(sh, 'corr', index=tag) OR $
              NOT tag_exist(sh, 'corr_err', index=errtag) ) THEN BEGIN 
              message,'corr/corr_err tags not found'
          ENDIF 
      END 
      ELSE: BEGIN 
          IF (NOT tag_exist(sh, 'sigma', index=tag) OR $
              NOT tag_exist(sh, 'sigmaerr', index=errtag) ) THEN BEGIN 
              message,'sigma/sigmaerr tags not found'
          ENDIF 
      END
  ENDCASE 
  return,{tag:tag, errtag:errtag}
END 


FUNCTION objshear::plotstruct, sh, type, $
                 xlog=xlog, xmin=xmin, xmax=xmax, $
                 ylog=ylog, ymin=ymin, ymax=ymax, $
                 symmetric=symmetric
  
  nbin = n_elements(sh)

  ;; Tag info
  tg=self->typetags(sh, type)
  tit = self->titles(type)

  data = sh.(tg.tag)
  err = sh.(tg.errtag)
  IF type EQ 'clustcorr' THEN BEGIN 
      data = data - 1.0
  ENDIF 

  yrange = prange( data, err, slack=1.5, symmetric=symmetric)
  xrange = [min(sh.meanr), max(sh.meanr)]/1000

  IF keyword_set(xlog) THEN BEGIN 
      xrange[0] = 10.0/1000.0   ; 10 kpc
 

      xrange[1] = xrange[1]*2.0
      xtickformat = 'loglabels'
      xticklen = 0.04
  ENDIF ELSE BEGIN 
      xlog=0
      xtickformat = ''
      xticklen = !x.ticklen
  ENDELSE 
  IF n_elements(xmin) NE 0 THEN xrange[0] = xmin
  IF n_elements(xmax) NE 0 THEN xrange[1] = xmax

  IF keyword_set(ylog) THEN BEGIN 
      yrange[0] = 1.e-2
      yrange[1] = yrange[1]*1.5
      ytickformat = 'loglabels'
      yticklen = 0.04
  ENDIF ELSE BEGIN 
      ylog=0
      ytickformat = ''
      yticklen = !y.ticklen
  ENDELSE 
  IF n_elements(ymin) NE 0 THEN yrange[0] = ymin
  IF n_elements(ymax) NE 0 THEN yrange[1] = ymax

  st = $
    {type: type, $
     tag: tg.tag, $
     errtag: tg.errtag, $
     $
     xlog: xlog, $
     xtickformat: xtickformat, $
     xticklen: xticklen, $
     xrange: xrange, $
     xtitle: tit.xtitle, $
     $
     ylog: ylog, $
     ytickformat: ytickformat, $
     yticklen: yticklen, $
     yrange: yrange, $
     ytitle: tit.ytitle $
    }  

  return,st
END 


PRO objshear::_plot_profile, sh, ps, $
            color=color, $
            psym=psym, linestyle=linestyle, $
            label=label, xlabeled=xlabeled, ylabeled=ylabeled, $
            overplot=overplot, $
            aspect=aspect

  IF n_elements(psym) EQ 0 AND n_elements(linestyle) EQ 0 THEN psym=8

  IF ps.xlog THEN BEGIN 
      xlog = 1
  ENDIF 
  IF ps.ylog THEN BEGIN 
      ylog = 1
  ENDIF 

  data = sh.(ps.tag)
  err = sh.(ps.errtag)

  IF ps.type EQ 'clustcorr' THEN data = data-1.0

  pplot, sh.meanr/1000.0, data, yerror=err, $
    overplot=overplot, $
    aspect=aspect, $
    color=color, $
    xlog=ps.xlog, ylog=ps.ylog, $
    xrange=ps.xrange, yrange=ps.yrange, $
    xticklen=ps.xticklen, yticklen=ps.yticklen, $
    $
    xstyle=3, ystyle=3, hat=0, $
    $
    psym=psym, linestyle=linestyle, symsize=0.5
    

END 

FUNCTION objshear::labels, subtype, bin=bin
  message,'You must override this function in descendant'
END 

PRO objshear::plot_profile, type, subtype=subtype, bin=bin, $
            sample=sample, $
            nsplit=nsplit, $
            xlog=xlog, xmin=xmin, xmax=xmax, $
            ylog=ylog, ymin=ymin, ymax=ymax, symmetric=symmetric, $
            incolor=incolor, $
            charsize=charsize, $
            label_charsize=label_charsize, $
            dops=dops, landscape=landscape, encapsulated=encapsulated

  IF n_elements(type) EQ 0 THEN BEGIN 
      message,'-Syntax: o->plot_profile, type, '+$
        'subtype=, bin=, sample=, /xlog, xmin=, /ylog, ymin=',/inf
      message,'type = "combined", "corr", "jack", etc',/inf
      print
      message,'halting'
  ENDIF 

  IF type EQ 'clustcorr' THEN dtype = 'corr' ELSE dtype=type

  IF keyword_set(dops) THEN BEGIN
      file = self->plotfile(dtype, subtype=subtype, bin=bin, sample=sample, $
                            color=incolor, $
                            encapsulated=encapsulated, $
                            /createdir)

      IF type EQ 'clustcorr' THEN BEGIN 
          file = repstr(file, 'corr.', 'clustcorr.')
      ENDIF 

      begplot, file, $
        encapsulated=encapsulated, color=incolor, landscape=landscape

    thick=2
    !p.thick=thick
    !p.charthick=thick
    !x.thick=thick
    !y.thick=thick
  ENDIF 

  IF n_elements(charsize) NE 0 THEN BEGIN 
      chold = !p.charsize
      !p.charsize = charsize
  ENDIF 

  esheldon_setup



  sharr = self->lensread(dtype, subtype=subtype, bin=bin, sample=sample)

  IF n_elements(subtype) NE 0 THEN labels = self->labels(subtype, bin=bin)

  IF n_elements(nsplit) NE 0 THEN BEGIN
      ;; This means that we should overplot the splits.
      ;; The actual number of plots will be nbin/nsplit
      nbin = n_elements(sharr)/nsplit
      psyms = -[8, 4, 1, 5, 6, 7]
      npsym = n_elements(psyms)
      linestyles=[0,1,3,5]
      nline = n_elements(linestyles)
  ENDIF ELSE BEGIN 
      nsplit = 1
      nbin = n_elements(sharr)
      psym = 8
  ENDELSE 

  simpctable, colorlist=colorlist

  nc = n_elements(colorlist)

  lab = self->axis_labeled(nbin)
  plotstruct = self->plotstruct(sharr, type, $
                                xlog=xlog, xmin=xmin, xmax=xmax, $
                                ylog=ylog, ymin=ymin, ymax=ymax, $
                                symmetric=symmetric)

  IF nbin GT 1 THEN BEGIN 
      mplot_value = self->mplot_value(nbin)
      erase & multiplot, mplot_value, /square, $
          mxtitle=plotstruct.xtitle, $
          mytitle=plotstruct.ytitle, $
          ytickf=plotstruct.ytickformat, xtickf=plotstruct.xtickformat
  ENDIF ELSE BEGIN 
      aspect = 1
  ENDELSE 


  IF nbin GT 8 AND keyword_set(dops) THEN BEGIN 
      !x.thick = 2
      !y.thick = 2
      !p.thick = 2
  ENDIF 

  j = 0
  ish = 0
  FOR i=0L, nbin-1 DO BEGIN 

      FOR j=0, nsplit-1 DO BEGIN 

          IF j GT 0 THEN overplot=1 ELSE overplot=0

          IF keyword_set(incolor) THEN BEGIN 
              color=colorlist[j MOD nc]
              add_arrval, color, lcolors
          ENDIF 
          IF n_elements(psyms) NE 0 THEN BEGIN 
              psym = psyms[j MOD npsym]
              lpsym = psym
              IF nsplit NE 1 THEN lpsym=abs(lpsym)
              add_arrval, lpsym, lpsyms
          ENDIF 
          IF n_elements(linestyles) NE 0 THEN BEGIN 
              linestyle=linestyles[j MOD nline]
              add_arrval, linestyle, llinestyles
          ENDIF 
          self->_plot_profile, sharr[ish], plotstruct, $
            label=label, $
            xlabeled=lab.xlabeled[i], $
            ylabeled=lab.ylabeled[i], $
            overplot=overplot, $
            aspect=aspect, psym=psym, linestyle=linestyle, $
            color=color

          IF n_elements(labels) NE 0 THEN BEGIN 
              add_arrval, labels[ish], llabels
          ENDIF 

          ish = ish+1
      ENDFOR 

      IF n_elements(llabels) NE 0 THEN BEGIN 
          ;;help,llabels,llinestyles,lpsyms,lcolors
          IF n_elements(label_charsize) EQ 0 THEN label_charsize=1

          legend, llabels, line=llinestyles, psym=lpsyms, color=lcolors, $
            /right, box=0, charsize=label_charsize, margin=0
      ENDIF 
      
      delvarx, lcolors, lpsyms, llinestyles, llabels

      IF nbin GT 1 AND i NE nbin-1 THEN multiplot
  ENDFOR 


  IF nbin GT 1 THEN BEGIN 
      multiplot, /reset
  ENDIF 

  IF n_elements(charsize) NE 0 THEN !p.charsize=chold

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encapsulated) THEN trim_bbox=1
      IF keyword_set(landscape) THEN landfix=1
      endplot,trim_bbox=trim_bbox, landfix=landfix
  ENDIF 

END 










FUNCTION objshear::split_ratios, sharr, sharr_split, index, nsplit

  nrad = n_elements(sharr[0].meanr)
  arrval = replicate(-9999.0, nrad)
  struct = replicate( {meanr: arrval, $
                       ratio: arrval, $
                       ratioerr: abs(arrval), $
                       ratio_mean: 0.0, $
                       ratio_meanerr: 0.0}, nsplit)

  radref = sharr[index].meanr
  ref = sharr[index].sigma
  referr = sharr[index].sigmaerr

  FOR i=0L, nsplit-1 DO BEGIN 

      ii = nsplit*index + i
      t = sharr_split[ii].sigma
      terr = sharr_split[ii].sigmaerr

      w = where(ref NE 0.0, nw)

      struct[i].meanr[w] = radref[w]

      ratio = t[w]/ref[w]
      struct[i].ratio[w] = ratio

;      raterr = struct[i].ratio[w]*sqrt( (referr[w]/ref[w]) ^2 + $
;                                        (terr[w]/t[w])^2 )

      raterr = struct[i].ratio[w]*terr[w]/t[w]

      struct[i].ratioerr[w] = raterr

      wmom, ratio, raterr, wmean, wsig, werr
      
      struct[i].ratio_mean = wmean
      struct[i].ratio_meanerr = werr

      print,'ratio_mean['+ntostr(index)+','+ntostr(i)+'] = '+ntostr(wmean)+' '+!plusminus+' '+ntostr(werr)

  ENDFOR 

  return, struct
  
END 

PRO objshear::_plot_ratios, sharr, xrange, psym, linestyle, color, $
  xlabeled, ylabeled


  num = n_elements(sharr)

  w=where(sharr.meanr GT 0)
  yrange=prange(sharr.ratio[w], sharr.ratioerr[w])

  yrange = [-1,3]

  IF xlabeled THEN xtitle=estitle('mpcxtitle2')
  IF ylabeled THEN ytitle='ratio'

  pplot, $
    sharr[0].meanr/1000, sharr[0].ratio, yerr=sharr[0].ratioerr, $
    xrange=xrange, xstyle=3, /xlog, xtitle=xtitle, $
    yrange=yrange, ystyle=3, ytitle=ytitle, $
    psym=psym[0], linestyle=linestyle[0], $
    color=color[0], hat=0

  oplot, [ 1.e-5, 1.e5], [sharr[0].ratio_mean, sharr[0].ratio_mean], $
         color=color[0], line=2
    
  FOR i=1,num-1 DO BEGIN 

      pplot, $
        sharr[i].meanr/1000, sharr[i].ratio, yerr=sharr[i].ratioerr, $
        psym=psym[i], linestyle=linestyle[i], $
        color=color[i], /overplot,hat=0

      oplot, [ 1.e-5, 1.e5], [sharr[i].ratio_mean, sharr[i].ratio_mean], $
             color=color[i], line=2

  ENDFOR 
    
  oplot, [1.e-5, 1.e5], [1,1], thick=2*!p.thick

END 

;; Plot the ratio of splits to the unsplit
PRO objshear::plot_ratio, dtype, subtype, subtype_split, nsplit, $
        bin=bin, sample=sample, $
        incolor=incolor, $
        charsize=charsize, $
        label_charsize=label_charsize, $
        dops=dops, landscape=landscape, encapsulated=encapsulated

  IF n_params() LT 4 THEN BEGIN 
      message,'-Syntax: o->plot_ratio, dtype, subtype, subtype_split, nsplit, '+$
        'subtype=, bin=, /xlog, xmin=, /ylog, ymin=',/inf
      message,'type = "combined", "corr", "jack", etc',/inf
      print
      message,'halting'
  ENDIF 

  IF keyword_set(dops) THEN BEGIN
      file = self->plotfile(dtype, subtype=subtype+'_'+subtype_split, bin=bin, sample=sample, $
                            color=incolor, $
                            encapsulated=encapsulated, $
                            /ratio, $
                            /createdir)

      begplot, file, $
        encapsulated=encapsulated, color=incolor, landscape=landscape

    thick=2
    !p.thick=thick
    !p.charthick=thick
    !x.thick=thick
    !y.thick=thick

  ENDIF 

  IF n_elements(charsize) NE 0 THEN BEGIN 
      chold = !p.charsize
      !p.charsize = charsize
  ENDIF 

  esheldon_setup



  sharr       = self->lensread(dtype, subtype=subtype, sample=sample)
  sharr_split = self->lensread(dtype, subtype=subtype_split, sample=sample)

  xlog = 1
  xrange = [ 0.5*min(sharr[0].meanr), max(sharr[0].meanr)*1.5]/1000

  labels = self->labels(subtype_split, bin=bin)


  ;; This means that we should overplot the splits.
  ;; The actual number of plots will be nbin/nsplit
  nbin = n_elements(sharr)
  psyms = -[8, 4, 1, 5, 6, 7]
  npsym = n_elements(psyms)
  linestyles=[0,1,3,5]
  nline = n_elements(linestyles)

  simpctable, colorlist=colorlist

  nc = n_elements(colorlist)

  IF nbin GT 1 THEN BEGIN 
      mplot_value = self->mplot_value(nbin)
      erase & multiplot, mplot_value, /square
  ENDIF ELSE BEGIN 
      aspect = 1
  ENDELSE 

  lab = self->axis_labeled(nbin)

  IF nbin GT 8 AND keyword_set(dops) THEN BEGIN 
      !x.thick = 2
      !y.thick = 2
      !p.thick = 2
  ENDIF 

  ish = 0
  FOR i=0L, nbin-1 DO BEGIN 

      split_ratios = self->split_ratios(sharr, sharr_split, $
                                        i, nsplit)

      FOR j=0L, nsplit-1 DO BEGIN 

          IF keyword_set(incolor) THEN BEGIN 
              color=colorlist[j MOD nc]
              add_arrval, color, lcolors
          ENDIF ELSE BEGIN 
              add_arrval, !p.color, lcolors
          ENDELSE 
          IF n_elements(psyms) NE 0 THEN BEGIN 
              psym = psyms[j MOD npsym]
              lpsym = psym
              lpsym=abs(lpsym)
              add_arrval, psym, ppsyms
              add_arrval, lpsym, lpsyms
          ENDIF 
          IF n_elements(linestyles) NE 0 THEN BEGIN 
              linestyle=linestyles[j MOD nline]
              add_arrval, linestyle, llinestyles
          ENDIF 

          IF n_elements(labels) NE 0 THEN BEGIN 
              add_arrval, labels[ish], llabels
          ENDIF 
         
          ish = ish+1
      ENDFOR 
      
      self->_plot_ratios, $
        split_ratios, xrange, ppsyms, llinestyles, lcolors, $
        lab.xlabeled[i], lab.ylabeled[i]
      

      IF n_elements(label_charsize) EQ 0 THEN label_charsize=1
      
      legend, llabels, line=llinestyles, psym=lpsyms, color=lcolors, $
              /right, box=0, charsize=label_charsize, margin=0
      
      delvarx, lcolors, lpsyms, llinestyles, llabels

      IF nbin GT 1 AND i NE nbin-1 THEN multiplot
  ENDFOR 


  IF nbin GT 1 THEN BEGIN 
      multiplot, /reset
  ENDIF 

  IF n_elements(charsize) NE 0 THEN !p.charsize=chold

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encapsulated) THEN trim_bbox=1
      IF keyword_set(landscape) THEN landfix=1
      endplot,trim_bbox=trim_bbox, landfix=landfix
  ENDIF 

END 















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Plot the randoms and real for comparison
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO objshear::compare_real_rand, subtype=subtype

  j = self->lensread('combined', subtype=subtype)
  r = self->lensread('random', subtype=subtype)
  
  num = n_elements(j)

  xrange = [0.5*min(j.meanr),max(j.meanr)*1.5]/1000

  FOR i=0L, num-1 DO BEGIN 

      w=where(j[i].meanr/1000 GT 1.0, nw)

      r_yrange = prange(r[i].sigma, r[i].sigmaerr, /slack, /symmetric)

      j_ymax = max(j[i].sigma[w] + j[i].sigmaerr[w])*1.1
      yrange = [r_yrange[0], j_ymax]



      pplot, r[i].meanr/1000.0, r[i].sigma, yerr=r[i].sigmaerr, $
        yrange=yrange, ystyle=3, /xlog, xrange=xrange, xstyle=3, $
        psym=4
      oplot, [1.e-4, 1.e10], [0,0]

      pplot, j[i].meanr/1000, j[i].sigma, yerr=j[i].sigmaerr, $
        color=c2i('red'), psym=-8, /overplot


      key = prompt_kbrd('hit a key')
      IF key EQ 'q' THEN return
  ENDFOR 


END 


; Just plot two of them for comparison
PRO objshear::compare_real_rand_two, subtype=subtype, last=last, dops=dops, $
            xrange=xrange

  IF keyword_set(dops) THEN BEGIN 
      tfile = self->plotfile('combined',subtype=subtype)
      psfile = repstr(tfile, 'combined.ps', 'comparetwo_rand_combined.eps')
      begplot, psfile, /color, /encapsulated
  ENDIF 

  esheldon_setup
  j = self->lensread('combined', subtype=subtype)
  r = self->lensread('random', subtype=subtype)
  
  ws = self->where_string(subtype, labels=labels)

  num = n_elements(j)

  IF n_elements(last) EQ 0 THEN last =  num-1

  j = j[ [0, last] ]
  r = r[ [0, last] ]
  labels = labels[ [0, last] ]

  num = 2

  erase & multiplot, [1,2]

  IF n_elements(xrange) EQ 0 THEN xrange = [0.5*min(j.meanr),max(j.meanr)*1.5]/1000

  ytitle = estitle('deltaytitle')

  FOR i=0L, num-1 DO BEGIN 

      w=where(j[i].meanr/1000 GT 1.0, nw)

      r_yrange = prange(r[i].sigma, r[i].sigmaerr, /slack, /symmetric)

      j_ymax = max(j[i].sigma[w] + j[i].sigmaerr[w])*1.1

;      yrange = [r_yrange[0], j_ymax]
      yrange = [-1, j_ymax]

      yrange=[-1.0,2.5]

      IF i EQ 1 THEN xtitle = estitle('mpcxtitle2') ELSE xtitle=''

      pplot, r[i].meanr/1000.0, r[i].sigma, yerr=r[i].sigmaerr, $
        yrange=yrange, ystyle=3, /xlog, xrange=xrange, xstyle=3, $
        xtitle=xtitle, ytitle=ytitle, $
        psym=4
      oplot, [1.e-4, 1.e10], [0,0]

      pplot, j[i].meanr/1000, j[i].sigma, yerr=j[i].sigmaerr, $
        psym=-8, /overplot

      legend, [labels[i], 'random'], $
        /bottom, /right, box=0, psym=[-8, 4], charsize=1.5

      IF i NE num-1 THEN multiplot

  ENDFOR 

  multiplot, /reset

  IF keyword_set(dops) THEN endplot, /trim_bbox

END 











; Note you can plot all on different panels with plot_profile
PRO objshear::plot_clustcorr_over, subtype, dops=dops, color=color


  IF keyword_set(dops) THEN BEGIN
      file = self->plotfile('corr', subtype=subtype,$
                            color=color, $
                            /encapsulated, $
                            /createdir)

      file = repstr(file, 'corr', 'clustcorr')
      print,file
      begplot, file, /encapsulated, color=color
  ENDIF 



  t = self->lensread('corr', subtype=subtype)
  nt = n_elements(t)

  lines=[0,1,2,3,4,5]
  IF keyword_set(dops) THEN BEGIN 
      thick = [1,1,1,1,1,1,7,7,7,7,7,7,10,10,10,10,10,10]
  ENDIF ELSE BEGIN 
      thick = [1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3]
  ENDELSE 
  nl = n_elements(lines)

  IF keyword_set(color) THEN BEGIN 
      simpctable, colorlist=colors
  ENDIF ELSE BEGIN 
      colors = replicate(!p.color, 100)
  ENDELSE 
  setup_mystuff
  xtitle=estitle('mpcxtitle2')
;  ylog=1
;  xlog=1
;  ytitle='C(R)-1'
;  pplot, t[0].meanr/1000, t[0].corr-1, $
;    xlog=xlog,ylog=ylog, $
;    yrange=[0.01,6], $
;    xrange=[0.01,20], $
;    xstyle=1,ystyle=1,aspect=1,$
;    line=lines[0],xtit=xtitle,ytit=ytitle, thick=thick[0]
  ytitle='C(R)'
  ylog=0
  xlog=1
  pplot, t[0].meanr/1000, t[0].corr, $
    xlog=xlog,ylog=ylog, $
    yrange=[1,4.5], $
    xrange=[0.01,20], $
    xstyle=1,ystyle=3,aspect=1,$
    line=lines[0],xtit=xtitle,ytit=ytitle, thick=thick[0]


  FOR i=0,nt-1 DO BEGIN 

;      pplot, t[i].meanr/1000, t[i].corr, yerr=t[i].corr_err, hat=0, $
;        /overplot,color=colors[i],line=lines[i mod nl], $
;        thick=thick[i], errthick=thick[i]
      pplot, t[i].meanr/1000, t[i].corr, hat=0, $
        /overplot,color=colors[i],line=lines[i mod nl], $
        thick=thick[i], errthick=thick[i]

      add_arrval, lines[i MOD nl], ll
      add_arrval, colors[i], cc
      add_arrval, thick[i], tt
  ENDFOR 


  ws = self->where_string(subtype, labels=labels)

  labels = reverse(labels)
  ll = reverse(ll)
  cc = reverse(cc)
  tt = reverse(tt)
  legend, labels, lines=ll, color=cc, /right,box=0,charsiz=1.5,thick=tt

  IF keyword_set(dops) THEN endplot,/trim_bbox

END 






; Calculate the weighted histogram of source redshifts around
; an input set of centers.  Uses lensing weights but currently 
; ignores redshift errors.
pro objshear::calc_weighted_zdist, lcat, scat, rmax

    

end








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Combine this sample with another one. The way this works, if one has more 
; radial bins than the other, then it is assumed that these can just be
; appended to the sample with fewere bins.
;
; caller is the object type of the caller to this function, so we can 
; make a new object
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function objshear::copy_samples, t1, t2
    nrad1 = n_elements(t1.meanr)
    nrad2 = n_elements(t2.meanr)

    if nrad1 eq nrad2 then message,'Both have same number of radial bins'

    if nrad1 gt nrad2 then begin
        tlong = t1
        tshort = t2
        nshort = nrad2

        sh = tlong
    endif else begin
        tlong = t2
        tshort = t1
        nshort = nrad1

        sh = tlong
    endelse

    zero_struct, sh

    nrad = max([nrad1, nrad2])

    tags = tag_names(sh)
    ntags=n_elements(tags)
    for i=0L, ntags-1 do begin

        ; copy in tags in both that are in all
        if tag_exist(tshort, tags[i], ind=ws) and tag_exist(tlong, tags[i], ind=wl) then begin

            sz=size(tlong.(wl))

            case sz[0] of
                ; scalars just copied from short
                0: sh.(i) = tshort.(ws)
                1: begin
                    ; arrays
                    ; only copy ones with nrad length in tlong. This covers all arrays
                    ; in zshstruct
                    if n_elements(tlong.(wl)) eq nrad then begin
                        array = replicate(tshort.(ws)[0], nrad)
                        array[0:nshort-1] = tshort.(ws)
                        array[nshort:nrad-1] = tlong.(wl)[nshort:nrad-1]

                        sh.(i) = array
                    endif else begin
                        message,'array '+tags[i]+' is not nrad in size',/inf
                    endelse
                end
                2: begin
                    ; matricies. Only do nradxnrad arrays.
                    if sz[1] eq nrad and sz[2] eq nrad then begin
                        ; simplest to copy long then copy in short as a sub-array
                        mat = tlong.(wl)
                        mat[0:nshort-1, 0:nshort-1] = tshort.(ws)
                        sh.(i) = mat
                    endif else begin
                        message,'matrix '+tags[i]+' is not nradxnrad in size',/inf
                    endelse
                end
                else: message,'Arrays with ndim > 2 not supported'
            endcase

        endif
    endfor

    if tag_exist(sh, 'covariance') then begin
        sh.correlation = cov2corr(sh.covariance)
    endif
    if tag_exist(sh, 'orthocov') then begin
        sh.orthocorrelation = cov2corr(sh.orthocov)
    endif


    return, sh
end

pro objshear::combine_samples, sample1, sample2, dtype, subtype=subtype

    if n_params() lt 3 then begin
        print,'-Syntax: o->combine_samples, sample1, sample2, dtype, subtype='
        on_error, 2
        message,'Halting'
    endif

    ; first read the files from this sample
    t1 = self->lensread(dtype, sample=sample1, subtype=subtype, count=count1)
    t2 = self->lensread(dtype, sample=sample2, subtype=subtype, count=count2)

    outfiles = self->lensfile(dtype, subtype=subtype, sample=[sample1,sample2],/createdir)

    if count2 ne count1 then message,'must be same length'


    print
    for i=0L, count1-1 do begin
        sh = self->copy_samples(t1[i], t2[i])
        print,'Writing to file: ',outfiles[i]
        write_idlstruct, sh, outfiles[i]
    endfor

end






function objshear::read_pairs, file, nrows=nrows, status=status

    status=1
    if n_elements(file) eq 0 then begin
        print,'-Syntax: pairs=o->read_pairs(file, nrows=, status=)'
        print,'   send nrows= to read a chunk'
        on_error, 2
        message,'Halting'
    endif

    openr, lun, file, /get_lun, error=status
    if status ne 0 then return, -1

    ;; first is the number of rows in a 64-bit field
    filerows = 0LL
    readu, lun, filerows

    if n_elements(nrows) eq 0 then nrows = filerows

    print
    print,'Reading '+ntostr(nrows)+' rows from file: ',file
    struct = {lindex:0L, sindex:0L, rkpc: 0.0, weight: 0.0}
    struct = replicate(struct, nrows)

    readu, lun, struct

    free_lun, lun

    return, struct
end


pro objshear::analyze_zdist, l, allpz, p

    if n_elements(l) eq 0 then begin
        p=self->par_struct()
        ; the main lens file
        l=self->lensread('input')
        bcg = self->get()
        bcg = bcg[l.zindex]
        w=where(bcg.ngals200 ge 3 and l.z ge 0.1 and l.z le 0.3)
        l = l[w]

        ; the source file
        ms=obj_new('make_scat') 
        source_file = ms->princeton_source_file(p.source_sample)
        s=read_idlstruct(source_file)
        allpz=s.photoz_z
        s=0
        obj_destroy, ms

        ; the pairs file
        p=self->lensread('pairs') 

        help, l, allpz, p
    endif

    ; pair up the z and weight
    zmatch = allpz[p.sindex]

;    zbin=0.01

    ; First the overall distribution
;    plothist, allpz, bin=zbin, peak=1 

;    bs=binner(zmatch, weights=p.weight, bin=zbin)
;    pplot, bs.xcenter, float(bs.hist)/max(bs.hist), psym=10, /over, color=!darkgreen
;    pplot, bs.xcenter, bs.whist/max(bs.whist), psym=10, /over, color=!red
;    legend, 'All',box=0
;    legend, $
;        ['Sources', 'Used Sources', 'Used Sources Weighted'], $
;        line=[0,0,0], color=[!p.color, !darkgreen, !red], $
;        /right, box=0


    ;; Now as a function of lens redshift
    ;nperbin = 6100
    nperbin = 14600
    lbs=binner(l.z, nperbin=nperbin, reverse_indices=rev)
    wn = where(lbs.hist eq nperbin, nbin)

    print,'n(zbin) = ',nbin
    clr=make_rainbow(nbin)
    print,clr

    ;key=prompt_kbrd('hit a key')

    zbin = 0.03
    for i=0L, nbin-1 do begin
        if rev[i] ne rev[i+1] then begin
            wl = rev[ rev[i]:rev[i+1]-1 ]

            tzi = l[wl].zindex

            match_multi, tzi, p.lindex, mp
            print,'npair = ',n_elements(mp)

            ; now plot weighted z dist
            tpz = zmatch[mp]
            tpw = p[mp].weight

            tbs = binner(tpz, weights=tpw, bin=zbin)

            if i eq 0 then begin
                w=where(tbs.whist gt 0.0)
                miny = min(tbs.whist, max=maxy)
                yrange=[miny*0.1, maxy*1.1]
                xrange=[0, max(tbs.xcenter)*1.2]
                pplot, tbs.xcenter, tbs.whist, psym=10, /ylog, yrange=yrange, ystyle=3, xrange=xrange, xstyle=1
            endif
            pplot, tbs.xcenter, tbs.whist, psym=10, color=clr[i], /overplot

        endif
    endfor 

    mess = 'z = '+ntostr(lbs.xmean[0:nbin-2])
    legend, mess, color=clr, line=0, /right,box=0

end

FUNCTION objshear::cleanup

  ptr_free, self.par_structPtr
  return, 1

END 

PRO objshear__define

  struct = { $
             objshear, $
             par_structPtr: ptr_new() $
           }

END 
