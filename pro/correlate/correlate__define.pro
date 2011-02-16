
;+
; NAME:
;  correlate__define
;
; PURPOSE:
;
; CATEGORY:
;
; CALLING SEQUENCE:
;  c=obj_new('correlate', par)
;
; INPUTS:
;  par: parameter structure
;
; METHODS:
;  file methods
;     corrdir_base
;     corrdir_info
;     randinfo
;     rzsub_info
;     sample
;     corrdir_type_exists
;     corrdir_dirtype_exists
;     corrdir_dirtype_hassub
;     corrdir
;     extension
;     add_randnums
;     corrfile
;     list_random
;     corr_read_kcorr
;     corr_read_htm
;     corr_read
;
;     fitfile
;
;  plotfile
;     plotdir_base
;     plotdir
;     plotfile
;
;  config, par, pbs file creation
;     write_pbs
;     write_par
;     write_script
;       output_type
;       corrtype
;
;  mag limits
;     print_max_apparent_imag
;     print_min_apparent_inmgy
;
;  structures
;     primary_struct
;     secondary_struct
;
;     radlumcolor_sumstruct
;     rad_sumstruct
;
;     meanstruct
;     means
;
;  summing outputs
;     _sum_output
;     _sum_match_rows
;     _sum_match_keeparray
;     combine_sumstructs
;     wcombine_sumstructs
;     sums2meanerr
;     calc_sum_derived
;     calc_totaldensity
;     _base_array
;     read_chunk
;     _process_subset_requests
;     sum_radlulmcolor_output
;     sum_radlumcolor_pairs
;     get_pair_index
;     sum_rad_output
;     
;  subsampling
;     rzsub
;       The idea with the rzsub stuff is to do be able to re-use a set of 
;       randoms for many samples.  Do a bunch of randoms and then take the 
;       averages in bins of redshift.  Then to generate the appropriate data 
;       for a sample just add these up weighted by the redshift histogram of 
;       the actual data.
;
;       rzsub_sample - redshift subsampling for rz and rr
;         rzsub_indices - get indices for objects in each histogram bin
;           rzsub_hist - bin according to rz binning for the current sample, 
;           return ind
;       rzsub_combine - combine random samples.  Earch redshift bin gets its 
;           own file.
;
      
;       rzsub_sample_jack
;         rzsub_jack_indices - return pointer array, eac element n_z,n_j
;            (rzsub_indices,rzsub_hist)
;
;       rzsub_fix_jackknifes
;       rzsub_read_jackknife_ids
;       rzsub_read_by_jackknife_id
;       rzsub_jackknife_combine
;
;       get_jackknife_rzsub
;       jackknife_file
;       match_hist_weight
;       rzsub_match
;       rzsub_match_keeparray
;
;     get_input_index
;     sub_sample
;
;  edge
;     model_usedarea
;     model_usedarea_plotfile
;
;  correct
;     _correct
;     correct
;
;  jackknife
;     get_index
;     get_jackknife_indices
;     create_jackknife_sumstructs
;     jackknife_test
;     jackknife_covariance
;
;  invert and model
;     nfwfit
;     invert
;
;  plotting routines, analysis
;     mplot_value
;     plot_usedarea
;     plot_usedarea_multi
;     plot_profile_over
;     radformat
;     plot_lumcolor_vs_rad
;     color_fractions
;     plot_bluefrac
;     find_colorlum_peak
;     plot_redseq
;     color_gaussfit
;     plot_lumfunc
;     run_schechter_fit
;     plot_all_schechter_fits
;     schechter_fit
;     plot_lum_vs_rad
;     plot_kgmr_vs_rad
;     _invert_integrate
;     invert_gettags
;
; EXAMPLES:
;   To sub-sample.  
;   -----------------------------
;     First you must have the rzsub set up, which averages the rd and rr
;     into redshift bins.
;       rzsub_sample, type, primary_randnum=, secondary_randnum=, /jackknife
;       rzsub_combine, type
;       rzsub_jackknife_combine, type
;
;     Now run the sub-sample code.  It will find the randoms from the files 
;     on disk.
;       sub_sample, type, subtype
;       correct, subtype=
;       means, subtype=
;   
;   To jackknife
;   ------------------------------
;     Run rzsub_sample with the /jackknife keyword (see above)
;     sub_sample, type, subtype, /jackknife
;     correct, subtype=, /jackknife
;     
;	To invert:
;	----------------------------
;		::invert, subtype=, corrected=, /dops
;
;
;   
; MODIFICATION HISTORY:
;
;-
function correlate::init, par

  if n_elements(par) eq 0 then begin 
      message,'-Syntax: co = obj_new("correlate", par)',/inf
      return, 0
  endif 
  newpar = self->add2par(par)
  help,newpar,/str, output=helpout
  for i=0,n_elements(helpout)-1 do print,helpout[i]

  self.par_struct = ptr_new(newpar)
  return, 1
end 











; Add the bin information to the parameter structure
function correlate::add2par, par

  ;; add the radial bins
  logRmin = alog10(par.rmin)
  logRmax = alog10(par.rmax)
  logBinsize = ( logRmax - logRmin )/par.nrad

  radbins_min = dblarr(par.nrad)
  radbins_max = dblarr(par.nrad)
  area = dblarr(par.nrad)

  for i=0l, par.nrad-1 do begin 
      radbins_min[i] = 10.0^(logrmin + i*logbinsize)
      radbins_max[i] = 10.0^(logrmin + (i+1)*logbinsize)
      area[i] = !pi*(radbins_max[i]^2 - radbins_min[i]^2)
  endfor 

  ;; add the lum bins
  if par.nlum gt 0 then begin 
      loglbins_min = dblarr(par.nlum)
      loglbins_max = dblarr(par.nlum)
      loglbinsize = (par.loglmax - par.loglmin)/par.nlum

      for i=0l,par.nlum-1 do begin 
          loglbins_min[i] = par.loglmin + i*loglbinsize
          loglbins_max[i] = par.loglmin + (i+1)*loglbinsize
      endfor 
  endif else begin 
      loglbins_min=-1
      loglbins_max=-1
  endelse
  ;; add the gmr bins
  if par.nkgmr gt 0 then begin 
      kgmrbins_min = dblarr(par.nlum)
      kgmrbins_max = dblarr(par.nlum)
      kgmrbinsize = (par.kgmrmax - par.kgmrmin)/par.nkgmr

      for i=0l,par.nlum-1 do begin 
          kgmrbins_min[i] = par.kgmrmin + i*kgmrbinsize
          kgmrbins_max[i] = par.kgmrmin + (i+1)*kgmrbinsize
      endfor 
  endif else begin
      kgmrbins_min = -1
      kgmrbins_max = -1
  endelse

  newpar = create_struct(par, $
                         'radbins_min', radbins_min, $
                         'radbins_max', radbins_max, $
                         'area', area, $
                         'loglbins_min', loglbins_min, $
                         'loglbins_max', loglbins_max, $
                         'kgmrbins_min', kgmrbins_min, $
                         'kgmrbins_max', kgmrbins_max)
  return, newpar

end 
function correlate::par_struct
  return, *self.par_struct
end 



; area of this catalog. Run stand alone
; calculate_area code in sdsspixIDL/app
; note, dr406 now uses the "BOUND" mask
function correlate::area, catalog, radians=radians
  case catalog of
      'dr4plus01': area = 7398.233579d
      'hv1gals': area = 5156.6195
      'rhv1gals': area = 5156.6195
      else: message,'do not have an area for catalog '+ntostr(catalog)
  end 
  if keyword_set(radians) then begin 
      area = area * (!dpi/180d)^2
  endif 
  return, area
end 

function correlate::secondary_density, secondary_randnum=secondary_randnum

  par = self->par_struct()


  if n_elements(secondary_randnum) eq 0 then begin 
      area = self->area(par.DSsample,/radians)
      secondary_file = self->corrfile('secondary','input')
  endif else begin 
      area = self->area(par.RSsample,/radians)
      secondary_file = self->corrfile('secondary_random','input', $
                                      secondary_randnum=secondary_randnum)
  endelse 
  hdr = read_idlheader(secondary_file)
  number = hdr.nrows
  density = number/area

  return, density
end 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; File related methods
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function correlate::corrdir_base
  return,expand_tilde(esheldon_config('correlate_dir'))
end 

; Info each directory type
function correlate::corrdir_info

  ddtypes = ['par','script','pbs','output',$
             'means', $
             'combined','corrected','corrected_invert', $
             'matchdr','matchrd','matchrr',$
             'jackknife','jackknife_invert', $
             'jsample','matchdr_jsample','matchrd_jsample','matchrr_jsample']
  par = self->par_struct()

  info = $
    { $
      kcorr:            {basedir: 'input',  dirtypes: 'input', sub: 0, sample: par.Ksample}, $
      $
      primary:          {basedir: 'input',  dirtypes: 'input', sub: 0, sample: par.DPsample}, $
      primary_random:   {basedir: 'input',  dirtypes: 'input', sub: 0, sample: par.RPsample}, $
      secondary:        {basedir: 'input',  dirtypes: ['input','htm'], sub: 0, sample: par.DSsample}, $
      secondary_random: {basedir: 'input',  dirtypes: ['input','htm'], sub: 0, sample: par.RSsample}, $
      $
      dd:               {basedir: 'output', dirtypes:  ddtypes, sub: 1, sample: par.psample + '_'+par.numerator_output_type+'_'+par.DPsample +'_'+par.DSsample}, $
      dr:               {basedir: 'output', dirtypes: ['par','script','pbs','output'], sub: 0, sample: par.psample+'_'+par.denominator_output_type+'_'+par.DPsample +'_'+par.RSsample}, $
      rd:               {basedir: 'output', dirtypes: ['par','script','pbs','output','rzsub','rzsub_jsample'], sub: 0, sample: par.psample+'_'+par.numerator_output_type+'_'+par.RPsample +'_'+par.DSsample}, $
      rr:               {basedir: 'output', dirtypes: ['par','script','pbs','output','rzsub','rzsub_jsample'], sub: 0, sample: par.psample+'_'+par.denominator_output_type+'_'+par.RPsample +'_'+par.RSsample}  $
    }

  return, info
end 
function correlate::corrdir_info_old

  ddtypes = ['par','script','pbs','output',$
             'means', $
             'combined','corrected','corrected_invert', $
             'matchdr','matchrd','matchrr',$
             'jackknife','jackknife_invert', $
             'jsample','matchdr_jsample','matchrd_jsample','matchrr_jsample']
  par = self->par_struct()

  info = $
    { $
      kcorr:            {basedir: 'input',  dirtypes: 'input', sub: 0, sample: par.Ksample}, $
      $
      primary:          {basedir: 'input',  dirtypes: 'input', sub: 0, sample: par.DPsample}, $
      primary_random:   {basedir: 'input',  dirtypes: 'input', sub: 0, sample: par.RPsample}, $
      secondary:        {basedir: 'input',  dirtypes: ['input','htm'], sub: 0, sample: par.DSsample}, $
      secondary_random: {basedir: 'input',  dirtypes: ['input','htm'], sub: 0, sample: par.RSsample}, $
      $
      dd:               {basedir: 'output', dirtypes:  ddtypes, sub: 1, sample: par.psample + '_'+par.DPsample +'_'+par.DSsample}, $
      dr:               {basedir: 'output', dirtypes: ['par','script','pbs','output'], sub: 0, sample: par.psample+'_'+par.DPsample +'_'+par.RSsample}, $
      rd:               {basedir: 'output', dirtypes: ['par','script','pbs','output','rzsub','rzsub_jsample'], sub: 0, sample: par.psample+'_'+par.RPsample +'_'+par.DSsample}, $
      rr:               {basedir: 'output', dirtypes: ['par','script','pbs','output','rzsub','rzsub_jsample'], sub: 0, sample: par.psample+'_'+par.RPsample +'_'+par.RSsample}  $
    }

  return, info
end 


FUNCTION correlate::randinfo, type, dtype
  st = {prand:0, srand: 0}
  ltype = strlowcase(type)
  ldtype = strlowcase(dtype)
  tst = strlowcase(type)+'-'+strlowcase(dtype)

  CASE ltype OF
      'primary_random': st.prand=1
      'secondary_random': st.srand=1
;      'dd': BEGIN 
;          CASE ldtype OF
;              'matchrd': st.prand = 1
;              'matchdr': st.srand = 1
;              'matchrr': BEGIN 
;                  st.prand = 1
;                  st.srand = 1
;              END 
;              ELSE:
;          ENDCASE 
;      END 
      'rd': st.prand = 1
      'dr': st.srand = 1
      'rr': BEGIN 
          st.prand = 1
          st.srand = 1
      END 
      ELSE: 
  ENDCASE 

  return, st
      
END 

; This rzsample just describes how the redshifts will be binned.
FUNCTION correlate::rzsub_info
  par = self->par_struct()
  CASE par.rzsample OF
      'rz01': BEGIN 
          ;; Note: histogram will give 41 bins for binsize=0.005, but
          ;; will only take the first 40.
          return, {rzmin: 0.1, rzmax: 0.3, rzstep: 0.005, nrz: 40}
      END 
      ELSE: message,'Unsupported rzsample: '+par.rzsample
  ENDCASE 
END 


FUNCTION correlate::sample, type, info=info

  IF n_elements(type) EQ 0 THEN BEGIN 
      on_error,2
      print,'-Syntax: sample = obj->sample(type, info=)'
      print
      message,'halting'
  ENDIF 
  IF self->corrdir_type_exists(type, info=info) THEN BEGIN 
      return,info.sample
  ENDIF ELSE BEGIN 
      message,'No matching type: "'+strlowcase(type)+'"',/inf
      return,''
  ENDELSE 
END 




FUNCTION correlate::corrdir_type_exists, type, info=info

  IF n_elements(type) EQ 0 THEN BEGIN 
      on_error,2
      print,'-Syntax: if co->corrdir_type_exists(type,info=info) then ...'
      print
      message,'Halting'
  ENDIF 
  corrdir_info = self->corrdir_info()
  IF tag_exist(corrdir_info, type, ind=ind) THEN BEGIN 
      info = corrdir_info.(ind)
      return, 1
  ENDIF ELSE BEGIN 
      return, 0
  ENDELSE 
END 

FUNCTION correlate::corrdir_dirtype_exists, type, dtype, info=info

  IF n_elements(type) EQ 0 OR n_elements(dtype) EQ 0 THEN BEGIN 
      on_error,2
      print,'-Syntax: if co->corrdir_dirtype_exists(type,dtype,info=) then ...'
      print
      message,'Halting'
  ENDIF 
  IF self->corrdir_type_exists(type, info=info) THEN BEGIN 
      w=where(info.dirtypes EQ strlowcase(dtype), nw)
      IF nw EQ 0 THEN return,0 ELSE return,1
  ENDIF ELSE BEGIN 
      return,0
  ENDELSE 
END 

FUNCTION correlate::corrdir_dirtype_hassub, type, info=info

  IF n_elements(type) EQ 0 THEN BEGIN 
      on_error,2
      print,'-Syntax: if co->corrdir_dirtype_hassub(type,info=) then ...'
      print
      message,'Halting'
  ENDIF 
  IF self->corrdir_type_exists(type, info=info) THEN BEGIN 
      return,info.sub
  ENDIF ELSE BEGIN 
      return,0
  ENDELSE 
END 


function correlate::corrdir, type, dtype=dtype, subtype=subtype, createdir=createdir


  if n_elements(type) eq 0 then begin 
      on_error,2
      print,'-Syntax: dir=co->corrdir(type, dtype=, subtype=, /createdir)'
      print
      print,'      type = kcorr|primary|primary_random|secondary|secondary_random|dd|dr|rd|rr'
      message,'Halting'
  endif 

  if self->corrdir_type_exists(type, info=info) then begin 
      dir = self->corrdir_base()
      dir = concat_dir(dir, info.basedir)
      dir = concat_dir(dir, strlowcase(type))

      sample = self->sample(type)
      dir = concat_dir(dir, sample)

      if n_elements(dtype) ne 0 then begin 
          if self->corrdir_dirtype_exists(type, dtype, info=info) then begin 

              if n_elements(subtype) ne 0 then begin 
                  if info.sub then begin 
                      dir = concat_dir(dir, 'sub')
                      dir = concat_dir(dir, strlowcase(subtype))
                  endif else begin 
                      message,'type: "'+strlowcase(type)+'" does not support subsamples'
                  endelse 
              endif 

              dir = concat_dir(dir, strlowcase(dtype))
          endif else begin 
              message,'type: "'+strlowcase(type)+'" has no directory: "'+strlowcase(dtype)+'"'
          endelse 
      endif 
  endif else begin 
      message,'unknown type: "'+strlowcase(type)+'"'
  endelse 

  if keyword_set(createdir) then begin 
      if not file_test(dir, /dir) then file_mkdir, dir
  endif 
  return,dir
end 

function correlate::corrdir_old, type, dtype=dtype, subtype=subtype, createdir=createdir
    dir=self->corrdir(type, dtype=dtype, subtype=subtype, createdir=createdir)
    dir = repstr(dir, 'c_r_', '')
    dir = repstr(dir, 'cl_r_', '')
    dir = repstr(dir, 'cl_clr_', '')
    return,dir
end 

function correlate::extension, type, dtype

  if n_elements(type) eq 0 or n_elements(dtype) eq 0 then begin 
      on_error,2
      print,'-Syntax: ext=co->extension(type, dtype)'
      print
      message,'Halting'
  endif 

  ltype = strlowcase(type)
  ldtype = strlowcase(dtype)
  if ldtype eq 'par' then begin 
      ext = '.conf'
  endif else if ldtype eq 'script' then begin 
      ext = '.sh'
  endif else if ldtype eq 'pbs' then begin 
      ext = '.pbs'
  endif else if ltype eq 'kcorr' or ldtype eq 'htm' then begin 
      ext = '.bin'
  endif else begin 
      ext = '.st'
  endelse 
  return, ext
end 


function correlate::add_randnums, oldfiles, randnum
  nf = n_elements(oldfiles)
  nrand = n_elements(randnum)
  for fi=0l, nf-1 do begin 
      for ri=0l,  nrand-1 do begin 
          newfile = oldfiles[fi]+'_'+strn(randnum[ri], length=2, padchar='0')
          add_arrval, newfile, newfiles
      endfor 
  endfor 
  return, newfiles
end 

function correlate::corrfile, type, dtype, subtype=subtype, bin=bin, rzbin=rzbin, primary_randnum=primary_randnum, secondary_randnum=secondary_randnum, createdir=createdir, nodir=nodir


  if n_elements(type) eq 0 or n_elements(dtype) eq 0 then begin 
      on_error,2
      print,'-Syntax: file=co->corrfile(type, dtype, subtype=, bin=, rzbin=, primary_randnum=, secondary_randnum=, /createdir, /nodir)'
      print
      print,'      type = kcorr|primary|primary_random|secondary|secondary_random|dd|dr|rd|rr'
      message,'halting'
  endif 

  par = self->par_struct()

  dir = self->corrdir(type, dtype=dtype, subtype=subtype, createdir=createdir)

  ltype = strlowcase(type)
  ldtype = strlowcase(dtype)

  file = ltype+'_'+self->sample(type)
  
  ;; sub types
  if n_elements(subtype) ne 0 then begin 
      lsubtype = strlowcase(subtype)
      file = file + '_'+lsubtype

      ;; bin numbers for subtype 
      nbin = self->subtype_nbin(subtype)
      if n_elements(bin) ne 0 then begin 
          if bin lt 0 or bin ge nbin then begin 
              message,'bin out of bounds: [0,'+ntostr(nbin-1)+']'
          endif 
          binstr = strn(bin, length=2, padchar='0')
      endif else begin 
          for i=0,nbin-1 do add_arrval, strn(i, length=2, padchar='0'), binstr
      endelse 
          
      file = file + '_'+binstr

  ENDIF 

  if ldtype eq 'rzsub' or ldtype eq 'rzsub_jsample' then begin 

      rzinfo = self->rzsub_info()
      file = file + '_'+par.rzsample

      nrzbin = rzinfo.nrz
      if n_elements(rzbin) ne 0 then begin 
          if rzbin lt 0 or rzbin ge nrzbin then begin 
              message,'rzBin out of bounds: [0,'+ntostr(nrzbin-1)+']'
          endif 
          rzbinstr = strn(rzbin, length=2, padchar='0')
      endif else begin 
          for i=0,nrzbin-1 do add_arrval, strn(i, length=2, padchar='0'), rzbinstr
      endelse 

      file = file +'_'+rzbinstr
  endif 

  file = file + '_'+ldtype

  ;; random numbers for matched stuff
  randinfo = self->randinfo(type, dtype)
  if randinfo.prand then begin 
      if n_elements(primary_randnum) ne 0 then begin 
          file = self->add_randnums(file, primary_randnum)
      endif else begin 
          ;; We now have combined files.
;          on_error, 2
;          message,'You must enter a primary_randnum for "'+ltype+'", "'+ldtype+'"'
      endelse 
  endif 
  if randinfo.srand then begin 
      if n_elements(secondary_randnum) ne 0 then begin 
          file = self->add_randnums(file, secondary_randnum)
      endif else begin 
          ;; We now have combined files.
;          on_error, 2
;          message,'You must enter a secondary_randnum for "'+ltype+'", "'+ldtype+'"'
      endelse 
  endif 

  ;; extension
  ext = self->extension(type,dtype)
  file = file + ext

  ;; add directory?
  if not keyword_set(nodir) then begin 
      file = concat_dir(dir, file)
  endif 

  return, file

end 

;; list the actual random files we have on disk
function correlate::list_random, type, subtype=subtype, bin=bin, rzbin=rzbin, count=count

  if n_elements(type) eq 0 then begin 
      on_error, 2
      print,'-Syntax: rfiles = obj->list_random(type, subtype=, bin=, rzbin=, count=)'
      print,' type = dr|matchrd|matchdr|matchrr|rd_rzsub|rd_rzsub_jsample|rr_rzsub|rr_rzsub_jsample'
  endif 
  case strlowcase(type) of
      'matchrd': begin 
          rfile0 = self->corrfile('dd', 'matchrd', $
                                  subtype=subtype, bin=bin, primary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00.st', '_[0-9][0-9].st')
      end 
      'matchdr': begin 
          rfile0 = self->corrfile('dd', 'matchdr', $
                                  subtype=subtype, bin=bin, secondary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00.st', '_[0-9][0-9].st')
      end 
      'matchrr': begin 
          rfile0 = self->corrfile('dd', 'matchrr', $
                                  subtype=subtype, bin=bin, $
                                  primary_randnum=0, secondary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00_00.st', '_[0-9][0-9]_[0-9][0-9].st')
      end 
      'rd_rzsub': begin 
          IF n_elements(rzbin) EQ 0 THEN message,'You must enter rzbin='
          rfile0 = self->corrfile('rd','rzsub',rzbin=rzbin,primary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00.st', '_[0-9][0-9].st')
      end 
      'rr_rzsub': begin 
          if n_elements(rzbin) eq 0 then message,'You must enter rzbin='
          rfile0 = self->corrfile('rr','rzsub',rzbin=rzbin,$
                                  primary_randnum=0,secondary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00_00.st', '_[0-9][0-9]_[0-9][0-9].st')
      end 
      'rd_rzsub_jsample': begin 
          if n_elements(rzbin) eq 0 then message,'You must enter rzbin='
          rfile0 = self->corrfile('rd','rzsub_jsample',rzbin=rzbin,primary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00.st', '_[0-9][0-9].st')
      end 
      'rr_rzsub_jsample': begin 
          if n_elements(rzbin) eq 0 then message,'you must enter rzbin='
          rfile0 = self->corrfile('rr','rzsub_jsample',rzbin=rzbin,$
                                  primary_randnum=0,secondary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00_00.st', '_[0-9][0-9]_[0-9][0-9].st')
      end 

      'dr': begin 
          rfile0 = self->corrfile('dr', 'output', secondary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00.st', '_[0-9][0-9].st')
      end 
      'rd': begin 
          rfile0 = self->corrfile('rd', 'output', primary_randnum=0)
          rfile_pattern = repstr(rfile0, '_00.st', '_[0-9][0-9].st')
      end 
      'rr': begin 
          rfile0 = self->corrfile('rr', 'output', $
              primary_randnum=0, secondary_randnum=0)
          rfile_pattern = repstr(rfile0, '00_00.st', '[0-9][0-9]_[0-9][0-9].st')
      end 
      else: begin 
          on_error, 2
          message,'type must be dr|rd|rr|matchrd|matchdr|matchrr|rd_rzsub|rr_rzsub'
      end 
  endcase 

  rfiles = findfile(rfile_pattern, count=count)
  return, rfiles
end 

; assumes rr only have same randnums for prim,sec
function correlate::list_randnums, type, subtype=subtype, bin=bin, rzbin=rzbin, count=count

    files = self->list_random(type, subtype=subtype, bin=bin, rzbin=rzbin, count=count)
    ; just extract the last number from file names.
    nums = long( strmid(files, 4, 2,/reverse) )
    return, nums
end


; make reverse indices file; this is generic for all

pro correlate::make_htm_revind, type, secondary_randnum=secondary_randnum

    if type ne 'secondary' and type ne 'secondary_random' then begin
        message,'type must be "secondary" or "secondary_random"'
    endif
    if type eq 'secondary_random' $
            and n_elements(secondary_randnum) eq 0 then begin
        message,'Send secondary_randnum= when type="secondary_random"'
    endif
    files = self->corrfile(type,'input',$
        secondary_randnum=secondary_randnum)
    rev_files = self->corrfile(type,'htm', $
        secondary_randnum=secondary_randnum, /createdir)
    nf = n_elements(files)

    for i=0l, nf-1 do begin 

        print,'-----------------------------------------------'

        file = files[i]
        rev_file = rev_files[i]

        print
        print,'Will write to file: ',rev_file
        print
        print,'Reading from file: ',file

        st = read_idlstruct(file, column='htm_index')
        htm_index = st.htm_index
        st = 0

        print
        print,'Creating reverse indices'

        minid = min(htm_index, max=maxid)
        h = histogram( htm_index-minid, min=0, rev=rev )
        nrev = n_elements(rev)



        print
        print,'Writing to rev file: ',rev_file
        openw, lun, rev_file, /get_lun
        writeu, lun, nrev

        writeu, lun, minid
        writeu, lun, maxid

        writeu, lun, rev

        free_lun, lun

        rev = 0
        htm_index=0
    endfor 

end 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; File reading routines
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function correlate::corr_read_kcorr, file

  openr, lun, file, /get_lun

  nz = 0l
  readu, lun, nz
  z = fltarr(nz)
  readu, lun, z

  ngmr = 0l
  readu, lun, ngmr
  gmr = fltarr(ngmr)
  readu, lun, gmr

  nrmi = 0l
  readu, lun, nrmi
  rmi = fltarr(nrmi)
  readu, lun, rmi

  nband = 0l
  readu, lun, nband
  bands = lonarr(nband)
  readu, lun, bands

  kcorr = fltarr(nz, ngmr, nrmi, nband)

  tmp = 0.0
  for iz=0l, nz-1 do begin 
      for igmr=0l, ngmr-1 do begin 
          for irmi=0l, nrmi-1 do begin 
              for ib=0l, nband-1 do begin 
                  
                  readu, lun, tmp
                  kcorr[iz, igmr, irmi, ib] = tmp

              endfor
          endfor
      endfor
  endfor 

  free_lun, lun

  struct = $
    { $
      nz: nz, $
      zmin: min(z, max=zmax), $
      zmax: zmax, $
      zstep: z[1]-z[0], $
      z: z, $
      $
      ngmr: ngmr, $
      gmrmin: min(gmr, max=gmrmax), $
      gmrmax: gmrmax, $
      gmrstep: gmr[1]-gmr[0], $
      gmr: gmr, $
      $
      nrmi: nrmi, $
      rmimin: min(rmi, max=rmimax), $
      rmimax: rmimax, $
      rmistep: rmi[1]-rmi[0], $
      rmi: rmi, $
      $
      bands: bands, $
      kcorr: kcorr $
    }

  return, struct
END 

function correlate::corr_read_htm, file
  openr, lun, file, /get_lun
  nrev = 0L
  minid = 0L
  maxid = 0L
  readu, lun, nrev
  readu, lun, minid
  readu, lun, maxid
  
  rev = lonarr(nrev)
  readu, lun, rev
  
  struct = {nrev: nrev, $
            min: minid, $
            max: maxid, $
            rev: temporary(rev)}
  return, struct
end 

function correlate::corr_read, type, dtype, subtype=subtype, bin=bin, rzbin=rzbin, primary_randnum=primary_randnum, secondary_randnum=secondary_randnum, columns=columns, silent=silent


  if n_elements(type) eq 0 or n_elements(dtype) eq 0 then begin 
      on_error,2
      print,'-Syntax: file=co->corr_read(type, dtype, subtype=, bin=, primary_randnum=, secondary_randnum=, columns=, hdr=, /silent)'
      print
      message,'Halting'
  endif 

  ltype = strlowcase(type)
  ldtype = strlowcase(dtype)
  file = self->corrfile(type, dtype, subtype=subtype, bin=bin, rzbin=rzbin, primary_randnum=primary_randnum, secondary_randnum=secondary_randnum)
  ext = self->extension(type,dtype)
  case ext of
      '.st': begin 
          if n_elements(file) eq 1 and not keyword_set(silent) then begin 
              print
              print,'Reading file: ',file
          endif 
          struct = read_idlstruct_multi(file, columns=columns, silent=silent)
      end 
      '.bin': begin 
          if ltype eq 'kcorr' and not keyword_set(silent) then begin 
              print
              print,'Reading file: ',file
              struct = self->corr_read_kcorr(file)
          endif else if ltype eq 'htm' and not keyword_set(silent) then begin 
              print
              print,'Reading file: ',file
              struct = self->corr_read_htm(file)
          endif else begin 
              message,'.bin must be of dtype "kcorr" or "htm"'
          endelse 
      end 
      else: message,'Cannot read files with extension: '+ext
  endcase 
  return, struct
end 



;; thse will be used by children of this class
function correlate::primary_struct, num

  struct = { $
             bcg_id: 0L, $
             ra: 0d, $
             dec: 0d, $
             z: 0.0 $
           }
  if n_elements(num) ne 0 then begin 
      struct = replicate(struct, num[0])
  endif 
  return,struct
end 
function correlate::secondary_struct, num, random=random

  if keyword_set(random) then begin 
      struct = { $
                 ra: 0d, $
                 dec: 0d, $
                 htm_index: 0l $
               }
  endif else begin 
      struct = { $
                 ra: 0d, $
                 dec: 0d, $
                 gflux: 0.0, $
                 rflux: 0.0, $
                 iflux: 0.0, $
                 htm_index: 0l $
               }
  endelse 

  if n_elements(num) ne 0 then begin 
      struct = replicate(struct, num[0])
  endif 
  return,struct
end 



function correlate::plotdir_base
    return,expand_tilde( esheldon_config('plot_dir') )
end

function correlate::plotdir, subtype=subtype, base=base, matchzrand=matchzrand, createdir=createdir
  dir = self->plotdir_base()
  dir = concat_dir(dir, 'correlate')

  sample = self->sample('dd')
  dir = concat_dir(dir, sample)

  if n_elements(subtype) ne 0 and not keyword_set(base) then begin 
      dir = concat_dir(dir, strlowcase(subtype))
  endif 

  if keyword_set(matchzrand) then begin 
      dir = concat_dir(dir, 'matchzrand')
  endif 

  if keyword_set(createdir) then begin 
      if not file_test(dir, /dir) then file_mkdir, dir
  endif 

  return,dir  
end 

function correlate::plotfile, dtype, subtype=subtype, bin=bin, color=color, encapsulated=encapsulated, ratio=ratio, nodir=nodir, createdir=createdir


  if n_elements(dtype) eq 0 then begin 
      on_error,2
      print,'-Syntax: psfile = c->plotfile(dtype, subtype=, bin=, /color, /encapsulated, /ratio, /nodir, /createdir)'
      print
      message,'Halting'
  endif  

  ldtype = strlowcase(dtype)

  tfile = self->corrfile('dd', dtype, subtype=subtype, bin=0, /nodir)

  if keyword_set(encapsulated) then ext='.eps' else ext='.ps'
  if keyword_set(color) then cstr = '_color' else cstr=''
  if keyword_set(ratio) then rstr = '_ratio' else rstr=''

  ext = rstr + cstr + ext

  file = repstr(tfile, '.st', ext)

  if n_elements(subtype) ne 0 then begin 
      if n_elements(bin) ne 0 then begin 
          ;; plot for specific bin
          bstr = strn(bin[0], len=2, padchar='0')
          file = repstr(file, '00_'+dtype, bstr+'_'+ldtype)
      endif else begin 
          file = repstr(file, '00_'+dtype, ldtype)
      endelse 
  endif 

  if not keyword_set(nodir) then begin 
      dir = self->plotdir(subtype=subtype, createdir=createdir)  
      file = concat_dir(dir, file)
  endif 

  return, file

end 

; files for fitted parameters
function correlate::fitfile, dtype, subtype=subtype, bin=bin, ratio=ratio, nodir=nodir, createdir=createdir

    psfile = self->plotfile(dtype, subtype=subtype, bin=bin, ratio=ratio, nodir=nodir, createdir=createdir)
    fitfile = repstr(psfile,'.ps','.fit')
    return, fitfile

end







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Creating input configure files and scripts
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; write par, script, and pbs files
pro correlate::write_par_script_pbs, mafalda=mafalda, hours=hours, mem=mem

    sep = mkstr(50, val='-')

    nrand = self->numrand()

    self->write_par, 'dd', mafalda=mafalda
    self->write_script, 'dd', mafalda=mafalda
    self->write_pbs, 'dd', hours=hours, mem=mem, mafalda=mafalda
    print,sep

    for i=0l, nrand-1 do begin 
        self->write_par, 'rd', primary_randnum=i, mafalda=mafalda
        self->write_script, 'rd', primary_randnum=i, mafalda=mafalda

        self->write_pbs, 'rd', primary_randnum=i, hours=hours, mem=mem, mafalda=mafalda
    endfor 
    print,sep

    for i=0l, nrand-1 do begin 
        self->write_par, 'dr', secondary_randnum=i, mafalda=mafalda
        self->write_script, 'dr', secondary_randnum=i, mafalda=mafalda

        self->write_pbs, 'dr', secondary_randnum=i, hours=hours, mem=mem, mafalda=mafalda
    endfor 
    print,sep

    ; We only write when the numbers are equal for now
    for i=0l, nrand-1 do begin 
        self->write_par, 'rr', primary_randnum=i, secondary_randnum=i, mafalda=mafalda
        self->write_script, 'rr', primary_randnum=i, secondary_randnum=i, mafalda=mafalda

        self->write_pbs, 'rr', primary_randnum=i, secondary_randnum=i,$ 
            hours=hours, mem=mem, mafalda=mafalda
    endfor 

end 


pro correlate::write_pbs, type, hours=hours, mem=mem, mafalda=mafalda, $
        primary_randnum=primary_randnum, secondary_randnum=secondary_randnum

    if n_elements(type) eq 0 then begin
        on_error, 2
        print,'-Syntax: c->write_pbs, type, hours=, mem=, primary_randnum=, secondary_randnum='
        print
        print,'for mafalda:'
        print,'    hours:'
        print,'      < 24 for medium queue'
        print,'      168 hours is a week'
        print,'      720 hours is 30 days'
        print,'    mem in megabytes'
        message,'Halting'
    endif 

    parfile = self->corrfile(type, 'par', $
        primary_randnum=primary_randnum, $
        secondary_randnum=secondary_randnum, /createdir)
    parnodir = self->corrfile(type, 'par', $
        primary_randnum=primary_randnum, $
        secondary_randnum=secondary_randnum, $
        /nodir)
    pbs = self->corrfile(type, 'pbs', $
        primary_randnum=primary_randnum, $
        secondary_randnum=secondary_randnum, /createdir)

    output_file = self->corrfile(type, 'output', $
        primary_randnum=primary_randnum, $
        secondary_randnum=secondary_randnum, $
        /createdir)


    ;; where we will copy
    outdir = self->corrdir(type, dtype='output')

    if keyword_set(mafalda) then begin
        parfile  = repstr(parfile, '/mount/early2', '/home')
        parfile  = repstr(parfile, '.conf', '_mafalda.conf')
        parnodir = repstr(parnodir, '.conf', '_mafalda.conf')
        output_file =    repstr(output_file, '/mount/early2', '/scratch')

        scratchdir = repstr(outdir, '/mount/early2', '/scratch')
        datadir = repstr(outdir, '/mount/early2', '/home')
    endif


    ;; The job name
    jobname = strlowcase(type)
    if n_elements(primary_randnum) ne 0 then begin 
        jobname = jobname + '_'+strn(primary_randnum, len=2, padchar='0')
    endif 
    if n_elements(secondary_randnum) ne 0 then begin 
        jobname = jobname + '_'+strn(secondary_randnum, len=2, padchar='0')
    endif 

    print,'Writing to pbs: ',pbs,f='(A-25,A)'
    openw, lun, pbs, /get

    if keyword_set(mafalda) then begin
        if n_elements(hours) eq 0 or n_elements(mem) eq 0 then begin
            message,'For mafalda must specify hours and mem'
        endif

        pbs_machine = repstr(pbs, '/mount/early2','/home')

        hstr = ntostr(hours,f='(I)')
        memstr = ntostr(mem,f='(I)')
        printf, lun, "#!/bin/sh"
        printf, lun, "#PBS -V"               ; Keep environment
        printf, lun, "#PBS -l cput="+hstr+":0:0"   ; at most this many hours
        printf, lun, "#PBS -l mem="+memstr+"mb"    ; at most 1.8Gb
        printf, lun, "#PBS -N "+jobname      ; job name
        printf, lun, "#PBS -m ae"            ; send email on abort and end
        printf, lun
        printf, lun, "# dir in home directory. Only copy into this directory"
        printf, lun, "# after the job has completed"
        printf, lun, "datadir="+datadir 
        printf, lun, "mkdir -p $datadir" ; Make sure the datadir is there.
        printf, lun
        printf, lun, "# the local scratch disk on the node"
        printf, lun, "scratch="+scratchdir
        printf, lun, "mkdir -p $scratch" ; Make sure scratch area is there.
        printf, lun
        printf, lun, "# The output file on scratch disk"
        printf, lun, "output="+output_file
        printf, lun, '/home/esheldon/local/bin/correlate '+parfile+' &> '+pbs_machine+'.out'
        printf, lun
        printf, lun, "# copy the data to the home area and remove the file"
        printf, lun, "cp -f $output $datadir/"
        printf, lun, "if [ $? == 0 ] ;then rm $output; fi"
        printf, lun
    endif else begin
        dirsep, pbs, pbsdir, pbs_bname
        dirsep, parfile, pardir, par_bname

        printf, lun, '#!/bin/bash'
        printf, lun, "#PBS -N "+jobname      ; job name
        printf, lun, "#PBS -M erin.sheldon@gmail.com"
        printf, lun, "#PBS -m ae"            ; send email on abort and end
        printf, lun
        printf, lun, "outdir="+outdir
        printf, lun, "pbsdir="+pbsdir
        printf, lun, "pardir="+pardir
        printf, lun, "mkdir -p $outdir"  ; make sure datadir is there
        printf, lun, "mkdir -p $pbsdir"   ; make sure pbs dir is there
        printf, lun
        printf, lun, '/usr/bin/time -p /home/esheldon/local/bin/correlate $pardir/'+par_bname+' &> $pbsdir/'+pbs_bname+'.out'
    endelse

    free_lun, lun

END 


pro correlate::write_pbs_suball, type, wait=wait

    if n_params() lt 1 then begin
        print,'usage:  m->write_pbs_suball, type, wait='
        print,'  type rd:dr:rr'
        print,'  wait in seconds'
        on_error, 2
        message,'halting'
    endif

    if n_elements(wait) ne 0 then begin
        wtstr = ntostr(wait,f='(I)')
    endif else begin
        wtstr = '60'
    endelse

    nrand = self->numrand()
    randnums = lindgen(nrand)
    case type of
        'rd': begin
            pbs = self->corrfile(type, 'pbs', $
                primary_randnum=randnums, /createdir)
            file = repstr(pbs[0],'_00.pbs', '_suball.sh')
        end
        'dr': begin
            pbs = self->corrfile(type, 'pbs', $
                secondary_randnum=randnums, /createdir)
            file = repstr(pbs[0],'_00.pbs', '_suball.sh')
        end
        'rr': begin
            ; by default it does all pairs, we don't want that
            ; here
            pbs=strarr(nrand)
            for i=0L, nrand-1 do begin
                pbs[i] = self->corrfile(type, 'pbs', $
                    primary_randnum=i, $
                    secondary_randnum=i, /createdir)
            endfor
            file = repstr(pbs[0],'_00_00.pbs', '_suball.sh')
        end
        else: message,'rd, dr, rr'
    endcase

    print,'Writing to file: ',file
    openw,lun,file,/get_lun

    for i=0L, nrand-1 do begin
        printf, lun, 'echo "qsub '+pbs[i]+'"'
        printf, lun, 'qsub '+pbs[i]
        printf, lun, 'sleep '+wtstr
        printf, lun
    endfor

    free_lun, lun

    return

end


pro correlate::write_script, type, $
            primary_randnum=primary_randnum, $
            secondary_randnum=secondary_randnum, $
            mafalda=mafalda


  if n_elements(type) eq 0 then begin 
      on_error, 2
      print,'-Syntax: c->write_script, type, primary_randnum=, secondary_randnum=, /mafalda'
      print
      message,'Halting'
  endif 

  parfile = self->corrfile(type, 'par', $
                           primary_randnum=primary_randnum, $
                           secondary_randnum=secondary_randnum, /createdir)
  parnodir = self->corrfile(type, 'par', $
                          primary_randnum=primary_randnum, $
                          secondary_randnum=secondary_randnum, $
                          /nodir)
  script = self->corrfile(type, 'script', $
                          primary_randnum=primary_randnum, $
                          secondary_randnum=secondary_randnum, /createdir)

  if keyword_set(mafalda) then begin 
      script   = repstr(script, '.sh', '_mafalda.sh')
      parfile  = repstr(parfile, '/mount/early2', '/home')
      parfile  = repstr(parfile, '.conf', '_mafalda.conf')
      parnodir = repstr(parnodir, '/mount/early2', '/home')
      parnodir = repstr(parnodir, '.conf', '_mafalda.conf')
  endif 

  executable = '/usr/bin/time -p correlate'

  print,'Writing to script: ',script,f='(A-25,A)'

  openw, lun, script, /get

  printf, lun, executable+' '+parfile

  mess = 'finished correlate '+parnodir
  printf,lun,'dt=`date`'
  address = 'erin.sheldon@gmail.com'
  printf, lun, $
    'echo "'+mess+' $dt" | mail '+address+' -s "'+mess+'"'

  free_lun, lun

end 

function correlate::output_type, oname
    case strlowcase(oname) of
        'c_r': return, 1 ;; counts binned by radius
        'cl_r': return, 2 ;; counts+lum binned by radius
        'cl_clr': return, 3 ;; counts+lum binned by color-lum-rad
        else: message,'Unknown output typename: '+ntostr(oname)
    endcase
end 
function correlate::corrtype, cname 
  case strlowcase(cname) of
    'real_secondaries': return, 1
    'rand_secondaries': return, 2
    else: message,'Unknown corr typename: '+ntostr(cname)
  endcase
end 

pro correlate::write_par, type, $
            primary_randnum=primary_randnum,$ 
            secondary_randnum=secondary_randnum, $
            mafalda=mafalda


  if n_elements(type) eq 0 then begin 
      on_error, 2
      print,'-Syntax: c->write_par, type, primary_randnum=, secondary_randnum=, /mafalda'
      print,'type = dd|dr|rd|rr'
      print
      message,'Halting'
  endif 

  par = self->par_struct()
  parfile = self->corrfile(type, 'par', $
                           primary_randnum=primary_randnum, $
                           secondary_randnum=secondary_randnum, $
                           /createdir)
  

  ltype = strlowcase(type)
  case ltype of
      'dd': begin 
          primary_file = self->corrfile('primary','input', /createdir)
          secondary_file = self->corrfile('secondary','input', /createdir)
          htmrev_file = self->corrfile('secondary','htm', /createdir)
          kcorr_file = self->corrfile('kcorr','input', /createdir)
          corrtype = self->corrtype('real_secondaries') 
          output_type = self->output_type(par.numerator_output_type)
      end
      'dr': begin 
          primary_file = self->corrfile('primary','input', /createdir)
          secondary_file = self->corrfile('secondary_random','input', $
                                          secondary_randnum=secondary_randnum, $
                                          /createdir)
          htmrev_file = self->corrfile('secondary_random','htm', $
                                       secondary_randnum=secondary_randnum, $
                                       /createdir)
          kcorr_file = self->corrfile('kcorr','input', /createdir)
          corrtype = self->corrtype('rand_secondaries')
          output_type = self->output_type(par.denominator_output_type) 
      end
      'rd': begin 
          primary_file = self->corrfile('primary_random','input', $
                                        primary_randnum=primary_randnum, $
                                        /createdir)
          secondary_file = self->corrfile('secondary','input', /createdir)
          htmrev_file = self->corrfile('secondary','htm', /createdir)
          kcorr_file = self->corrfile('kcorr','input', /createdir)
          corrtype = self->corrtype('real_secondaries') 
          output_type = self->output_type(par.numerator_output_type)
      end
      'rr': begin 
          primary_file = self->corrfile('primary_random','input', $
                                        primary_randnum=primary_randnum, $
                                        /createdir)
          secondary_file = self->corrfile('secondary_random','input', $
                                          secondary_randnum=secondary_randnum, $
                                          /createdir)
          htmrev_file = self->corrfile('secondary_random','htm', $
                                       secondary_randnum=secondary_randnum, $
                                       /createdir)
          kcorr_file = self->corrfile('kcorr','input', /createdir)
          corrtype = self->corrtype('rand_secondaries') 
          output_type=self->output_type(par.denominator_output_type)
      end
  endcase 



  output_file = self->corrfile(type, 'output', $
                               primary_randnum=primary_randnum, $
                               secondary_randnum=secondary_randnum, $
                               /createdir)

  if keyword_set(mafalda) then begin 
      parfile =        repstr(parfile,          '.conf',         '_mafalda.conf')
      primary_file =   repstr(primary_file,     '/mount/early2', '/home')
      secondary_file = repstr(secondary_file,   '/mount/early2', '/home')
      htmrev_file =    repstr(htmrev_file,      '/mount/early2', '/home')
      kcorr_file =     repstr(kcorr_file,       '/mount/early2', '/home')
      output_file =    repstr(output_file,      '/mount/early2', '/scratch')
  endif

  print,'Writing to par: ',parfile,f='(A-25,A)'

;;stop

  openw, lun, parfile, /get_lun
;  lun=-1

  printf, lun, 'version             '+par.version
  printf, lun, 'sample              '+self->sample(type)
  printf, lun, 'primary_file        '+primary_file
  printf, lun, 'secondary_file      '+secondary_file
  printf, lun, 'htmrev_file         '+htmrev_file
  printf, lun, 'kcorr_file          '+kcorr_file
  printf, lun, 'output_file         '+output_file
  printf, lun, 'corrtype            '+ntostr(corrtype)
  printf, lun, 'output_type         '+ntostr(output_type)
  printf, lun, 'h                   '+ntostr(par.h)
  printf, lun, 'omega_m             '+ntostr(par.omega_m)
  printf, lun, 'nrad                '+ntostr(par.nrad)
  printf, lun, 'rmin                '+ntostr(par.rmin)
  printf, lun, 'rmax                '+ntostr(par.rmax)
  printf, lun, 'nlum                '+ntostr(par.nlum)
  printf, lun, 'lumband             '+ntostr(par.lumband)
  printf, lun, 'loglmin             '+ntostr(par.loglmin)
  printf, lun, 'loglmax             '+ntostr(par.loglmax)
  printf, lun, 'nkgmr               '+ntostr(par.nkgmr)
  printf, lun, 'kgmrmin             '+ntostr(par.kgmrmin)
  printf, lun, 'kgmrmax             '+ntostr(par.kgmrmax)
  printf, lun, 'comoving            '+ntostr(par.comoving)
  printf, lun, 'depth               '+ntostr(par.depth)

  free_lun, lun

end 







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Help figuring out mag limits
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; Print the corresponding apparent magnitude given a redshift and
;; an absolute magnitude.  
;;
;; Quick i-band: 
;;   z=0.30,  absmag=-19.0 -> m = 21.0706
;;   z=0.21,  absmag=-18.0 -> m = 21.0807
;;   z=0.135, absmag=-17   -> m = 21.0644

pro correlate::print_max_apparent_mag, band, z, absmag, band_shift, $
        omega_m=omega_m


    if n_params() lt 4 then begin 
        on_error, 2
        print,'-Syntax: mb->print_max_apparent_mag, band, z, absmag, band_shift, omega_m='
        print,'  band_shift should correspond with that in the kcorr struct'
        print
        message,'Halting'
    endif 

    ks = self->corr_read('kcorr','input')
    zint = fix( interpol(lindgen(ks.nz), ks.z, z) )

    tk = reform( ks.kcorr[zint, *, *, band] )

    sz = size(tk,/dim)
    index = lindgen(sz[0]*sz[1])
    x = index mod sz[0]
    y = index/sz[0]

    max_kcorr = max(tk, mind)


    tk = ks.kcorr[zint, *, *, band]
    w=where(tk ne -9999.0)
    min_kcorr = min( tk[w] )

    da = angdist(0.0, z, omega_m=omega_m, DL=dlum)
    distmod = 5.0*( alog10(dlum) + 6 ) - 5.0

    m_max = absmag + distmod + max_kcorr
    ; = 22.5 - 2.5*alog10(nmgy)
    nmgy_min = 10d^( -0.4*(m_max-22.5) )

    lumsolar = sdss_am2lumsolar(absmag, clr=band, band_shift=0.25)
    log_lumsolar = alog10(lumsolar)

    llstr = ntostr(log_lumsolar, f='(F0.2)')

    band_names = ['u','g','r','i','z']
    bname = band_names[band]
    print,'z = ',z
    print,'gmr = ',ks.gmr[x[mind]]
    print,'rmi = ',ks.rmi[y[mind]]
    print,'zint = ',zint
    print,'absmag('+bname+') = ',absmag
    print,'lumsolar('+bname+') = ',lumsolar
    print,'log(lumsolar) = '+llstr
    print,'min '+bname+' kcorr', min_kcorr
    print,'max '+bname+' kcorr', max_kcorr,'  flux: ',10.0^(-0.4*max_kcorr)
    print,'Max apparent '+bname+'-band mag:  ',m_max
    print,'Min apparent '+bname+'-band nmgy: ',nmgy_min

end 


pro correlate::print_min_apparent_inmgy, z, ilum, omega_m=omega_m


  if n_params() lt 2 then begin 
      on_error, 2
      print,'-Syntax: mb->print_min_apparent_inmgy, z, ilum, omega_m='
      print,'ilum in units of 10^10 solar'
      print
      message,'Halting'
  endif 

  sunknmgy = ( self->sunknmgy() )[3]

  ks = self->read_kcorr()
  zint = fix( interpol(lindgen(ks.nz), ks.z, z) )

  tk = reform( ks.kcorr[zint, *, *, 3] )

  sz = size(tk,/dim)
  index = lindgen(sz[0]*sz[1])
  x = index mod sz[0]
  y = index/sz[0]

  max_i_kcorr = max(tk, mind)
  
  min_i_kflux = 10.0^( -0.4*max_i_kcorr )

  da = angdist_lambda(z, omegamat=omega_m, dlum=dlum)

  knmgy_min = ilum*sunknmgy/(dlum^2)
  nmgy_min = knmgy_min*min_i_kflux

  print,'z = ',z
  print,'gmr = ',ks.gmr[x[mind]]
  print,'rmi = ',ks.rmi[y[mind]]
  print,'zint = ',zint
  print,'ilum (10^10) = ',ilum
  print,'min i kflux', min_i_kflux,'  max i kcorr: ',max_i_kcorr
  print,'Min apparent i-band nmgy: ',nmgy_min
end 




































;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Sum rows of the output file
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function correlate::radlumcolor_sumstruct, corrected=corrected

  par = self->par_struct()

  nrad = par.nrad
  nlum = par.nlum
  ngmr = par.nkgmr

  radarray = dblarr(nrad)
  radarray_cumlum = dblarr(nrad,nlum)

  if par.numerator_output_type eq 'cl_clr' then begin
      biglon64 = lon64arr(nrad, nlum, ngmr)
      bigdbl = dblarr(nrad,nlum,ngmr)
  endif else begin 
      biglon64 = lon64arr(nrad)
      bigdbl = dblarr(nrad)
  endelse
  
  if not keyword_set(corrected) then begin 
      rcounts = lon64arr(nrad)
	  rcounts_cumlum = lon64arr(nrad,nlum)
      counts = biglon64
  endif else begin 
      rcounts = dblarr(nrad)
      counts = bigdbl
  endelse 


  ilum = bigdbl 

  if not keyword_set(corrected) then begin 
      sumstruct = $
        { $
          npoints: 0LL, $
          nrad: nrad, $
          nlum: nlum, $
          ngmr: ngmr, $
          totpairs: 0LL, $
          jackknife_id: 0L, $
          $
          radbins_min: par.radbins_min, $
          radbins_max: par.radbins_max, $
          area: par.area, $
          loglbins_min: par.loglbins_min, $
          loglbins_max: par.loglbins_max, $
          kgmrbins_min: par.kgmrbins_min, $
          kgmrbins_max: par.kgmrbins_max, $
          $
          rsum:  radarray, $
          rsum2: radarray, $
          $
          kgflux:  radarray, $      ; sum of flux over all primaries
          kgflux2: radarray, $      ; sum of flux^2 over all primaries, for color errors
          krflux:  radarray, $
          krflux2: radarray, $
          kiflux:  radarray, $
          kiflux2: radarray, $
          $
          radcounts:   rcounts, $ ; totaled over color,lum
		  $ ; totaled over color, cumulative sum over lum from low to high
		  radcounts_cumlum: rcounts_cumlum,    $ 
          ilumcounts:  ilum,    $ ; sum if i-band luminosity over all primaries
          ilumcounts2: ilum,    $ ; sum of square of above, for errors
          counts:      counts,  $ ; Total counts
          $                
          $
          $                     ;   -- Derived quantities: means/errs over primaries
          $
          r:         radarray, $    ; mean radius
          r_err:     radarray, $
          kgmr:      radarray, $    ; mean total color
          kgmr_err:  radarray, $ 
          krmi:      radarray, $
          krmi_err:  radarray, $
          $
          radilum:         radarray, $ ; mean luminosity in radial bins
          radilum_err:     radarray, $
          radilumdens:     radarray, $
          radilumdens_err: radarray, $
          radnumdens:      radarray, $ ; mean density in radial bins, per Mpc^2
          radnumdens_err:  radarray, $
		  $ ; mean density in radial bins, but cumulative across luminosity
          radnumdens_cumlum: radarray_cumlum, $ 
          radnumdens_cumlum_err:  radarray_cumlum, $
          $
          ilum:         bigdbl, $ ; mean luminosity in all bins
          ilum_err:     bigdbl, $
          ilumdens:     bigdbl, $
          ilumdens_err: bigdbl, $
          numdens:      bigdbl, $ ; mean density in all bins, per Mpc^2
          numdens_err:  bigdbl  $
        }
      return, sumstruct
  endif else begin 
      corrected = $
        { $
          npoints: 0LL, $
          nrad: nrad, $
          nlum: nlum, $
          ngmr: ngmr, $
          totpairs: 0LL, $
          jackknife_id: 0L, $
          $
          radbins_min: par.radbins_min, $
          radbins_max: par.radbins_max, $
          area: par.area, $
          loglbins_min: par.loglbins_min, $
          loglbins_max: par.loglbins_max, $
          kgmrbins_min: par.kgmrbins_min, $
          kgmrbins_max: par.kgmrbins_max, $
          $
          $
          $  ;   -- Derived quantities: means/errs over primaries
          $
          r:        radarray, $     ; mean radius
          r_err:    radarray, $
          kgmr:     radarray, $     ; mean total color
          kgmr_err: radarray, $
          krmi:     radarray, $
          krmi_err: radarray, $
          $; mean luminosity density in rad bins, per Mpc^2
          radilumdens:     radarray, $ 
          radilumdens_err: radarray, $
		  $; mean density in radial bins, per Mpc^2
          radnumdens:      radarray, $ 
          radnumdens_err:  radarray, $
          radnumdens_cumlum: radarray_cumlum, $ 
          radnumdens_cumlum_err:  radarray_cumlum, $
          $
          ilumdens:     bigdbl, $ ; mean lum density in all bins, per Mpc^2
          ilumdens_err: bigdbl, $
          numdens:      bigdbl, $ ; mean density in all bins, per Mpc^2
          numdens_err:  bigdbl, $
          $
		  $ ; mean luminosity density in rad bins, per Mpc^2
          radilumdens_field:     radarray, $ 
          radilumdens_field_err: radarray, $
		  $ ; mean density in radial bins, per Mpc^2
          radnumdens_field:      radarray, $ 
          radnumdens_field_err:  radarray, $
		  $ ; now with cumulative luminosities
          radnumdens_cumlum_field:      radarray_cumlum, $ 
          radnumdens_cumlum_field_err:  radarray_cumlum, $
          $
		  $ ; mean lum density in all bins, per Mpc^2
          ilumdens_field:     bigdbl, $ 
          ilumdens_field_err: bigdbl, $
          numdens_field:      bigdbl, $ ; mean density in all bins, per Mpc^2
          numdens_field_err:  bigdbl  $

        }
      return, corrected
  endelse 
end 


function correlate::rad_sumstruct, corrected=corrected

  par = self->par_struct()

  nrad = par.nrad
  nlum = par.nlum
  ngmr = par.nkgmr

  radarray = dblarr(nrad)
  l64rad = lon64arr(nrad)

  if not keyword_set(corrected) then begin 
      rcounts = lon64arr(nrad)
      counts = l64rad
  endif else begin 
      rcounts = dblarr(nrad)
      counts = radarray
  endelse 

  if not keyword_set(corrected) then begin 
      sumstruct = $
        { $
          npoints: 0LL, $
          nrad: nrad, $
          nlum: nlum, $
          ngmr: ngmr, $
          totpairs: 0LL, $
          jackknife_id: 0L, $
          $
          radbins_min: par.radbins_min, $
          radbins_max: par.radbins_max, $
          area: par.area, $
          loglbins_min: par.loglbins_min, $
          loglbins_max: par.loglbins_max, $
          kgmrbins_min: par.kgmrbins_min, $
          kgmrbins_max: par.kgmrbins_max, $
          $
          rsum:  radarray, $
          rsum2: radarray, $
          $
          totaldensity:       0d, $ ; total(secondary density/DA^2)
          usedarea:     radarray, $ ; counts/totaldensity
          usedarea_err: radarray, $
          radcounts:    rcounts,  $ ; These two totals over color,lum
          counts:       rcounts,  $ ; Total counts
          $
          $                     ;   -- Derived quantities: means/errs over primaries
          $
          r:         radarray, $    ; mean radius
          r_err:     radarray, $
          $
          radnumdens:      radarray, $ ; mean density in radial bins, per Mpc^2
          radnumdens_err:  radarray  $
        }
      return, sumstruct
  endif else begin 
      corrected = $
        { $
          npoints: 0LL, $
          nrad: nrad, $
          nlum: nlum, $
          ngmr: ngmr, $
          totpairs: 0LL, $
          jackknife_id: 0L, $
          $
          radbins_min: par.radbins_min, $
          radbins_max: par.radbins_max, $
          area: par.area, $
          loglbins_min: par.loglbins_min, $
          loglbins_max: par.loglbins_max, $
          kgmrbins_min: par.kgmrbins_min, $
          kgmrbins_max: par.kgmrbins_max, $
          $                     ;   -- Derived quantities: means/errs over primaries
          $
          totaldensity:       0d, $ ; total(secondary density/DA^2)
          usedarea:     radarray, $ ; counts/totaldensity
          usedarea_err: radarray, $
          $
          r:        radarray, $     ; mean radius
          r_err:    radarray, $
          $
          radnumdens:      radarray, $ ; mean density in radial bins, per Mpc^2
          radnumdens_err:  radarray  $
        }
      return, corrected
  endelse 
end 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; For summing over input primaries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro correlate::_sum_output, instruct, outstruct, index=index

    if n_elements(instruct) eq 0 or n_elements(outstruct) eq 0 then begin 
        on_error, 2
        print,'-Syntax: c->_sum_output, instruct, outstruct, index='
        print
        message,'Halting'
    endif 

    ninstruct = n_elements(instruct)

    nind = n_elements(index)
    if nind eq 0 then begin 
        nind = ninstruct
        index = lindgen(nind)
    endif

	par = self->par_struct()
    if par.numerator_output_type eq 'cl_clr' then begin
        cube=1
    endif else begin
        cube=0
    endelse

    lumthere = $
        tag_exist(outstruct, 'ilumcounts') or tag_exist(outstruct, 'radlum')

    outstruct.npoints = outstruct.npoints + nind

    ;; It's faster to do these with total
    if nind gt 1 then begin 
        outstruct.totpairs = $
            outstruct.totpairs + total_int( instruct[index].totpairs )
  
        outstruct.rsum  = $
            outstruct.rsum + total( instruct[index].rsum, 2, /double )
        outstruct.rsum2 = $
            outstruct.rsum2 + total( instruct[index].rsum^2, 2, /double )
      
        if lumthere then begin 
            outstruct.kgflux = $
                outstruct.kgflux + total( instruct[index].kgflux, 2, /double )
            outstruct.krflux = $
                outstruct.krflux + total( instruct[index].krflux, 2, /double )
            outstruct.kiflux = $
                outstruct.kiflux + total( instruct[index].kiflux, 2, /double )
          
            outstruct.kgflux2 = $
                outstruct.kgflux2 + total( instruct[index].kgflux^2,2,/double )
            outstruct.krflux2 = $
                outstruct.krflux2 + total( instruct[index].krflux^2,2,/double )
            outstruct.kiflux2 = $
                outstruct.kiflux2 + total( instruct[index].kiflux^2,2,/double )
        endif 
    endif else begin 
        outstruct.totpairs = outstruct.totpairs + instruct[index].totpairs

        outstruct.rsum  = outstruct.rsum + instruct[index].rsum
        outstruct.rsum2 = outstruct.rsum2 + instruct[index].rsum^2

        if lumthere then begin 
            outstruct.kgflux = outstruct.kgflux + instruct[index].kgflux
            outstruct.krflux = outstruct.krflux + instruct[index].krflux
            outstruct.kiflux = outstruct.kiflux + instruct[index].kiflux

            outstruct.kgflux2 = outstruct.kgflux2 + instruct[index].kgflux^2
            outstruct.krflux2 = outstruct.krflux2 + instruct[index].krflux^2
            outstruct.kiflux2 = outstruct.kiflux2 + instruct[index].kiflux^2
        endif 
    endelse 

    ; its faster to do this in a loop! presumably it has to do 
    ; with the strides..
    if cube then begin
        nrad = n_elements(outstruct.r)
        radcounts = lon64arr(nrad)
		nlum = n_elements(instruct[0].counts[0,*,0] )
		radcounts_cumlum = lon64arr(nrad, nlum)
    endif
    for ii=0l, nind-1 do begin 
        i = index[ii]


        if not cube then begin
            outstruct.radcounts = $
                outstruct.radcounts + instruct[index].radcounts
            if lumthere then begin
                outstruct.ilumcounts  = $
                    outstruct.ilumcounts  + instruct[i].radlum
                outstruct.ilumcounts2 = $
                    outstruct.ilumcounts2 + instruct[i].radlum^2
            endif
        endif else begin

            outstruct.counts = outstruct.counts + instruct[i].counts

            ;; needed to do this here because of the weighted summing
            ;; for randoms; small numbers get killed in the rounding
            for ir=0l,nrad-1 do begin 
                radcounts[ir] = total( instruct[i].counts[ir,*,*], /int )
            endfor 
            outstruct.radcounts = outstruct.radcounts + radcounts

            if lumthere then begin 
				; Now just summing over color, and doing a cumulative sum
				; over luminosity
				for ir=0L, nrad-1 do begin
					; pick a subarray for this radius
					subarray = reform(instruct[i].counts[ir,*,*])
					; total over color. This is now 2-d
					subarray = total(subarray, 2, /int)
					; now cumulatively sum over luminosity, but such that the 
					; value is the sum over all indices >= i
					for ilum=0L, nlum-1 do begin
						radcounts_cumlum[ir, ilum] = $
							total( subarray[ilum:nlum-1], /int)
					endfor
				endfor
				outstruct.radcounts_cumlum = $
					outstruct.radcounts_cumlum + radcounts_cumlum

				outstruct.ilumcounts  = $
					outstruct.ilumcounts  + instruct[i].ilum
				outstruct.ilumcounts2 = $
					outstruct.ilumcounts2 + instruct[i].ilum^2
            endif 

        endelse


    endfor 

end 




; indices is index for this chunk
; rows is the requested rows
pro correlate::_sum_match_rows, instruct, indices, rows, outstruct
  if n_params() lt 4 then begin 
      print,'-Syntax: obj->_sum_match_rows, instruct, indices, rows, outstruct'
      print
      message,'halting'
  endif 
  match, indices, rows, mind, mrows
  if mrows[0] ne -1 then begin 
      self->_sum_output, instruct, outstruct, index=mind
  endif 
end 

pro correlate::_sum_match_keeparray, instruct, indices, keeparray, outstructs
  if n_params() lt 4 then begin 
      print,'-Syntax: obj->_sum_match_keeparray, instruct, indices, keeparray, outstructs'
      print
      message,'halting'
  endif 

  nbin = n_elements(keeparray)
  if n_elements(outstructs) ne nbin then message,'outstructs must be same size as keeparray'
  for i=0l, nbin-1 do begin 

      if ptr_valid(keeparray[i]) then begin 
          match, indices, *keeparray[i], mind, mrows
          if mrows[0] ne -1 then begin 
              t = outstructs[i]
              self->_sum_output, instruct, t, index=mind
              outstructs[i] = t
          endif 
      endif 

  endfor 

end 







;; for summing over sumstructs
function correlate::combine_sumstructs, sumstructs

	nsum = n_elements(sumstructs)
	if nsum eq 0 then begin 
		on_error, 2
		print,'-Syntax: st = obj->combine_sumstructs(sumstructs)'
		print
		message,'Halting'
	endif 

	if nsum eq 1 then return, sumstructs

	par = self->par_struct()
	if par.numerator_output_type eq 'cl_clr' then begin
		cube=1
	endif else begin
		cube=0
	endelse

	lumthere = $
		tag_exist(sumstructs, 'ilumcounts') or tag_exist(sumstructs, 'radlum')

	if lumthere then begin 
		outstruct= self->radlumcolor_sumstruct()
	endif else begin 
		outstruct= self->rad_sumstruct()
	endelse 


	if cube then begin
		nrad = n_elements(sumstructs[0].r)
		radcounts = lon64arr(nrad)
		nlum = n_elements(sumstructs[0].counts[0,*,0] )
		radcounts_cumlum = lon64arr(nrad, nlum)
	endif

	for index=0l, nsum-1 do begin 

		outstruct.npoints = outstruct.npoints + sumstructs[index].npoints

		outstruct.totpairs = outstruct.totpairs + sumstructs[index].totpairs

		outstruct.rsum  = outstruct.rsum + sumstructs[index].rsum
		outstruct.rsum2 = outstruct.rsum2 + sumstructs[index].rsum2

		outstruct.counts = outstruct.counts + sumstructs[index].counts

		if cube then begin
			;; needed to do this here because of the weighted summing
			;; for randoms; small numbers get killed in the rounding
			for ir=0l,nrad-1 do begin 
				radcounts[ir] = total(sumstructs[index].counts[ir,*,*],/int)
			endfor 
			outstruct.radcounts = outstruct.radcounts + radcounts

			if lumthere then begin
				; Now just summing over color, and doing a cumulative sum
				; over luminosity
				for ir=0L, nrad-1 do begin
					; pick a subarray for this radius
					subarray = reform(sumstructs[index].counts[ir,*,*])
					; total over color. This is now 2-d
					subarray = total(subarray, 2, /int)
					; now cumulatively sum over luminosity, but such that the 
					; value is the sum over all indices >= i
					for ilum=0L, nlum-1 do begin
						radcounts_cumlum[ir, ilum] = $
							total( subarray[ilum:nlum-1], /int)
					endfor
				endfor
				outstruct.radcounts_cumlum = $
					outstruct.radcounts_cumlum + radcounts_cumlum
			endif

		endif else begin
			outstruct.radcounts = $
				outstruct.radcounts + sumstructs[index].radcounts
		endelse
      
		if lumthere then begin 
			outstruct.kgflux = outstruct.kgflux + sumstructs[index].kgflux
			outstruct.krflux = outstruct.krflux + sumstructs[index].krflux
			outstruct.kiflux = outstruct.kiflux + sumstructs[index].kiflux

			outstruct.kgflux2 = outstruct.kgflux2 + sumstructs[index].kgflux2
			outstruct.krflux2 = outstruct.krflux2 + sumstructs[index].krflux2
			outstruct.kiflux2 = outstruct.kiflux2 + sumstructs[index].kiflux2

			outstruct.ilumcounts  = $
				outstruct.ilumcounts  + sumstructs[index].ilumcounts
			outstruct.ilumcounts2 = $
				outstruct.ilumcounts2 + sumstructs[index].ilumcounts2
		endif else begin 
			outstruct.totaldensity = $
				outstruct.totaldensity + sumstructs[index].totaldensity
		endelse 
	endfor 

	self->calc_sum_derived, outstruct

  return, outstruct
end 


;; for summing over sumstructs with weights
function correlate::wcombine_sumstructs, sumstructs, weights

  nsum = n_elements(sumstructs)
  nweight = n_elements(weights)
  if nsum eq 0 or nweight eq 0 then begin 
      on_error, 2
      print,'-Syntax: st = obj->wcombine_sumstructs(sumstructs, weights)'
      print
      message,'Halting'
  endif 

  if nsum eq 1 then return, sumstructs

  if nweight ne nsum then message,'Weights must be same size as sumstructs'

  ;; make sure floating point
  wt = double(weights)

  ;; make sure max is unity so the weight*counts gives the percentage of
  ;; the counts to keep; just easier to think about
  perc_keep = wt/max(wt)

  par = self->par_struct()
  if par.numerator_output_type eq 'cl_clr' then begin
      cube=1
  endif else begin
      cube=0
  endelse

  lumthere = $
      tag_exist(sumstructs, 'ilumcounts') or tag_exist(sumstructs, 'radlum')

  if lumthere then begin 
      outstruct= self->radlumcolor_sumstruct()
  endif else begin 
      outstruct= self->rad_sumstruct()
  endelse 


  if cube then begin
      nrad = n_elements(sumstructs[0].r)
      radcounts = lon64arr(nrad)
	  nlum = n_elements(sumstructs[0].counts[0,*,0] )
	  radcounts_cumlum = lon64arr(nrad, nlum)
  endif


  for index=0l, nsum-1 do begin 

      outstruct.npoints = outstruct.npoints + $
          round(  sumstructs[index].npoints*perc_keep[index], /l64)

      outstruct.totpairs = outstruct.totpairs + $
          round( sumstructs[index].totpairs*perc_keep[index], /l64)
  
      outstruct.rsum  = outstruct.rsum + sumstructs[index].rsum*perc_keep[index]
      outstruct.rsum2 = outstruct.rsum2 + sumstructs[index].rsum2*perc_keep[index]

      outstruct.counts = outstruct.counts + $
          round( sumstructs[index].counts*perc_keep[index], /l64)
      
      if cube then begin
          ;; need to do this here for precision issues in weighting
          for ir=0l,nrad-1 do begin 
              radcounts[ir] = total_int( sumstructs[index].counts[ir,*,*] )
          endfor 
          outstruct.radcounts = outstruct.radcounts + $
              round( radcounts*perc_keep[index], /l64)

			if lumthere then begin
				; Now just summing over color, and doing a cumulative sum
				; over luminosity
				for ir=0L, nrad-1 do begin
					; pick a subarray for this radius
					subarray = reform(sumstructs[index].counts[ir,*,*])
					; total over color. This is now 2-d
					subarray = total(subarray, 2, /int)
					; now cumulatively sum over luminosity, but such that the 
					; value is the sum over all indices >= i
					for ilum=0L, nlum-1 do begin
						radcounts_cumlum[ir, ilum] = $
							total( subarray[ilum:nlum-1], /int)
					endfor
				endfor
				outstruct.radcounts_cumlum = $
					outstruct.radcounts_cumlum + $
					round( radcounts_cumlum*perc_keep[index], /l64 )
			endif

       endif else begin
          outstruct.radcounts = outstruct.radcounts + $
              round( sumstructs[index].radcounts*perc_keep[index], /l64)
       endelse


      if lumthere then begin 
          outstruct.kgflux = outstruct.kgflux + sumstructs[index].kgflux*perc_keep[index]
          outstruct.krflux = outstruct.krflux + sumstructs[index].krflux*perc_keep[index]
          outstruct.kiflux = outstruct.kiflux + sumstructs[index].kiflux*perc_keep[index]
          
          outstruct.kgflux2 = $
              outstruct.kgflux2 + sumstructs[index].kgflux2*perc_keep[index]
          outstruct.krflux2 = $
              outstruct.krflux2 + sumstructs[index].krflux2*perc_keep[index]
          outstruct.kiflux2 = $
              outstruct.kiflux2 + sumstructs[index].kiflux2*perc_keep[index]
          
          outstruct.ilumcounts  = $
              outstruct.ilumcounts  + sumstructs[index].ilumcounts*perc_keep[index]
          outstruct.ilumcounts2 = $
              outstruct.ilumcounts2 + sumstructs[index].ilumcounts2*perc_keep[index]
      endif else begin 
          outstruct.totaldensity = $
              outstruct.totaldensity + sumstructs[index].totaldensity*perc_keep[index]
      endelse 
  endfor 



  self->calc_sum_derived, outstruct

  return, outstruct
end 




pro correlate::sums2meanerr, xsum, x2sum, nx, xmean, xvar, xerr

  if n_params() lt 5 then begin 
      on_error, 2
      print,'-Syntax: obj->sums2meanerr, xsum, x2sum, nx,  xmean, var, xerr'
      print
      message,'Halting'
  endif 
  xmean = xsum/nx
  xvar = (x2sum - 2.0*xmean*xsum + nx*xmean^2)/(nx-1.0)
  xerr = sqrt(xvar/nx)

end 

pro correlate::calc_adderr, x1, x1err, x2, x2err, res, err
  res = x1+x2
  err = sqrt( x1err^2 + x2err^2 )
end 

function correlate::fluxerr2colorerr, f1mean, f1err, f2mean, f2err

  if n_params() lt 4 then begin 
      on_error,2
      print,'-Syntax: cerr = obj->fluxerr2colorerr(f1mean,f1err,f2mean,f2err)'
      print
      message,'Halting'
  endif 
  ;; sig_(r-i) = 2.5/ln(10) sqrt( sig_r^2/f_r^2 + sig_i^2/f_i^2 )
  cerr = 2.5/alog(10.)*sqrt( (f1err/f1mean)^2 + (f2err/f2mean)^2 )
  return, cerr
end 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some derived quantities from the summed data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro correlate::calc_sum_derived, st

  if n_elements(st) eq 0 then begin 
      on_error, 2
      print,'-Syntax: co->calc_sum_derived, st'
      print
      message,'Halting'
  endif 

  if tag_exist(st, 'ilumcounts') then begin 
      radlumcolor=1 
	  nlum = st.nlum
	  lumthere=1
  endif else begin 
      radlumcolor=0
	  lumthere=0
  endelse 


  nst = n_elements(st)
  for j=0l, nst-1 do begin 

      for i=0l, st[j].nrad-1 do begin 

          ;; mean lum in all bins. over primaries, so npoints here.
          if radlumcolor then begin 
              self->sums2meanerr, $
                st[j].ilumcounts[i,*,*], st[j].ilumcounts2[i,*,*], st[j].npoints, tm, tv, terr
              st[j].ilum[i,*,*] = tm
              st[j].ilum_err[i,*,*] = terr
              
              st[j].radilum[i] = total( st[j].ilum[i,*,*] )
              st[j].radilum_err[i] = sqrt( total(st[j].ilum_err[i,*,*]^2 ) )

              ;; mean lum in radial bins.  We don't have enough precision to
              ;; do it this way:
              ;;self->sums2meanerr, $
              ;;  total( st[j].ilumcounts[i,*,*] ), total( st[j].ilumcounts2[i,*,*] ), $
              ;;  st[j].npoints, tm, tv, terr
              ;;st[j].radilum[i] = tm
              ;;st[j].radilum_err[i] = terr
              
              ;; Luminosity density in all bins and radial bins
              st[j].ilumdens[i,*,*] = st[j].ilum[i,*,*]/st[j].area[i]
              st[j].ilumdens_err[i,*,*] = st[j].ilum_err[i,*,*]/st[j].area[i]
              
              st[j].radilumdens[i] = st[j].radilum[i]/st[j].area[i]
              st[j].radilumdens_err[i] = st[j].radilum_err[i]/st[j].area[i]

              ;; Number density over all rad,lum,color bins. 
              dcounts = double( st[j].counts[i,*,*] )
              st[j].numdens[i,*,*]     = dcounts/st[j].npoints/st[j].area[i]
              st[j].numdens_err[i,*,*] = sqrt(dcounts)/st[j].npoints/st[j].area[i]

          endif 

          ;; This is not mean per primary!  This notation works over simple radial
          ;; bins too
;          st[j].radcounts[i] = total_int( st[j].counts[i,*,*] )

          if not radlumcolor then begin 
              if st[j].totaldensity gt 0 then begin 
                  st[j].usedarea[i] = st[j].radcounts[i]/st[j].totaldensity
                  st[j].usedarea_err[i] = sqrt(st[j].radcounts[i])/st[j].totaldensity
              endif 
          endif 

          ;; now mean over primaries
          dradcounts = double( st[j].radcounts[i] )
          st[j].radnumdens[i] = dradcounts/st[j].npoints/st[j].area[i]
          st[j].radnumdens_err[i] = sqrt(dradcounts)/st[j].npoints/st[j].area[i]

		  if lumthere then begin
			  ;; now mean over primaries, cumulative across luminosity
			  dradcounts_cumlum = double( st[j].radcounts_cumlum[i,*] )
			  st[j].radnumdens_cumlum[i,*] = $
				  dradcounts_cumlum/st[j].npoints/st[j].area[i]
			  st[j].radnumdens_cumlum_err[i,*] = $
				  sqrt(dradcounts_cumlum)/st[j].npoints/st[j].area[i]
		  endif



          ;; Mean radius *per secondary* not per primary.
          self->sums2meanerr, $
            st[j].rsum[i], st[j].rsum2[i], st[j].radcounts[i], tm, tv, terr
          st[j].r[i] = tm
          st[j].r_err[i] = terr

          


      endfor 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Mean total color and error
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      if radlumcolor then begin 
          self->sums2meanerr, st[j].kgflux, st[j].kgflux2, st[j].npoints, kgfluxmean, var, kgfluxerr
          self->sums2meanerr, st[j].krflux, st[j].krflux2, st[j].npoints, krfluxmean, var, krfluxerr
          self->sums2meanerr, st[j].kiflux, st[j].kiflux2, st[j].npoints, kifluxmean, var, kifluxerr
          
          st[j].kgmr = -2.5*alog10( st[j].kgflux/st[j].krflux )
          st[j].kgmr_err = self->fluxerr2colorerr(kgfluxmean, kgfluxerr, krfluxmean, krfluxerr )
          
          st[j].krmi = -2.5*alog10( st[j].krflux/st[j].kiflux )
          st[j].krmi_err = self->fluxerr2colorerr(krfluxmean, krfluxerr, kifluxmean, kifluxerr )
      endif 

  endfor 

end 


function correlate::calc_totaldensity, z, secondary_randnum=secondary_randnum

  if n_elements(z) eq 0 then begin 
      on_error, 2
      print,'-Syntax: td = c->calc_totaldensity(z, secondary_randnum=)'
      print
      message,'Halting'
  endif 

  par = self->par_struct()

  ;; number of secondaries per steradian
  sdens = self->secondary_density(secondary_randnum=secondary_randnum)

  da = angdist_lambda( z, omegam=par.omega_m, h=par.h ) ; mpc
  totaldensity = sdens*total(1.0/da^2, /double)
  return, totaldensity

end 

function correlate::_base_array, type
  if n_elements(type) eq 0 then type='float'
  par=self->par_struct()
  if par.nlum gt 0 then begin
    basearr=fltarr(par.nrad,par.nlum,par.nkgmr)
  endif else begin
    basearr=fltarr(par.nrad)
  endelse
  if type eq 'long' then return,long(basearr) else return,basearr
end
function correlate::radlumcolor_base_struct
  par = self->par_struct()
  ilumarr = self->_base_array()
  carr = self->_base_array('long')
  st = { $
         index: 0L, $
         totpairs: 0L, $
         rsum: fltarr(par.nrad), $
         kgflux: fltarr(par.nrad), $
         krflux: fltarr(par.nrad), $,
         kiflux: fltarr(par.nrad), $
         counts: carr, $
         ilum: ilumarr $
       }
  return, st
end 
function correlate::read_chunk, lun, base_struct, num


  if n_params() lt 3 then begin 
      on_error, 2
      print,'-Syntax: st = c->read_chunk(lun, base_struct, num)'
      print
      message,'Halting'
  endif 
  out = replicate(base_struct, num)
  readu, lun, out

  return, out

end 



;; This is written assuming that, most of the time, we are interested
;; in most of the file.  
pro correlate::_process_subset_requests, nrows, rowsin, keeparray, rows, st, rows_input, keep_input

  rows_input=0
  keep_input=0
  if n_elements(keeparray) ne 0 then begin 
      nbin = n_elements(keeparray)
      for i=0l, nbin-1 do begin

          if ptr_valid(keeparray[i]) then begin 

              rmd = rem_dup(*keeparray[i])
              if n_elements(rmd) ne n_elements(*keeparray[i]) then begin 
                  message,'keeparray['+ntostr(i)+'] is not unique'
              endif 
              kmin = min(*keeparray[i], max=kmax)
              if kmin lt 0 or kmin ge nrows then begin 
                  message,'keeparray['+ntostr(i)+'] is out of range [0,'+ntostr(nrows-1)+']'
              endif 
          endif 

      endfor 

      st = replicate(st, nbin)
      keep_input = 1
  endif else if n_elements(rowsin) ne 0 then begin 
      w=where(rowsin lt 0 or rowsin ge nrows, nw)
      if nw ne 0 then message,'input rows out of range [0,'+ntostr(nrows-1)+']'

      ;; will also sort
      rows = rowsin[ rem_dup(rowsin) ]
      rows_input = 1
  endif 


end 

function correlate::sum_radlumcolor_output, file, rows=rowsin, keeparray=keeparray


  if n_elements(file) eq 0 then begin 
      on_error, 2
      print,'-Syntax: sumstruct = c->sum_radlumcolor_output(file, rows=, keeparray=)'
      print
      message,'Halting'
  endif 

  t0 = systime(1)

  ;; Get the output structure
  st = self->radlumcolor_sumstruct()


  ;; open the file
  openr, lun, file, /get_lun

  ;; Read the header
  hdr = read_idlheader(lun)
  nrows = hdr.nrows

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Deal with requested subsets
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  self->_process_subset_requests, nrows, rowsin, keeparray, rows, st, rows_input, keep_input

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Split into chunks
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  chunk = 10000L

  chunkdiv = nrows/chunk
  chunkmod = nrows mod chunk

  if chunkdiv ne 0 then begin 
      nchunk = chunkdiv
      nleft = chunkmod
  endif else begin 
      nchunk = 0
      nleft = nrows
  endelse 

  print
  print,'Summing radlumcolor rows from file '+file  
  print
  print,'Chunksize = ',chunk
  print,'Nchunk    = ',nchunk
  print,'Remainder = ',nleft


  base_struct = self->radlumcolor_base_struct()

  ;; Read and sum the chunks
  beg = 0L
  nchstr = '/'+ntostr(nchunk)
  indices = lindgen(nrows)


  for i=0l, nchunk-1 do begin 
      if i eq 0 then begin 
          print,'  Reading Chunk: '+ntostr(i+1)+nchstr, format='($,a)'
      endif else begin 
          print,' '+ntostr(i+1)+nchstr, format='($,a)'
      endelse 

      tmp = self->read_chunk(lun, base_struct, chunk)

      print,'*',format='($,a)'

      if keep_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+chunk-1]
          self->_sum_match_keeparray, tmp, matchid, keeparray, st
      endif else if rows_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+chunk-1]
          self->_sum_match_rows, tmp, matchid, rows, st
      endif else begin 
          self->_sum_output, tmp, st
      endelse 

      tmp = 0
      beg = beg+chunk
  endfor 

  if nleft ne 0 then begin 
      print,'  Reading remainder',format='($,a)'

      tmp = self->read_chunk(lun, base_struct, nleft)
      print,'*',format='($,a)'

      if keep_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+nleft-1]
          self->_sum_match_keeparray, tmp, matchid, keeparray, st
      endif else if rows_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+nleft-1]
          self->_sum_match_rows, tmp, matchid, rows, st
      endif else begin 
          self->_sum_output, tmp, st
      endelse 
      print

      tmp = 0
  endif 

  free_lun, lun

  ;; some derived quantities
  self->calc_sum_derived, st

  ptime,systime(1)-t0

  return, st

end 



function correlate::sum_radlumcolor_pairs, file, rows=rowsin, keeparray=keeparray


  if n_elements(file) eq 0 then begin 
      on_error, 2
      print,'-Syntax: sumstruct = c->sum_radlumcolor_pairs(file, rows=, keeparray=)'
      print
      message,'Halting'
  endif 

  t0 = systime(1)

  ;; Get the output structure
  st = self->radlumcolor_sumstruct()


  ;; open the file
  openr, lun, file, /get_lun

  ;; Read the header
  hdr = read_idlheader(lun)
  nrows = hdr.nrows

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Deal with requested subsets
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  self->_process_subset_requests, nrows, rowsin, keeparray, rows, st, rows_input, keep_input

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Split into chunks
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  chunk = 10000L

  chunkdiv = nrows/chunk
  chunkmod = nrows MOD chunk

  if chunkdiv ne 0 then begin 
      nchunk = chunkdiv
      nleft = chunkmod
  endif else begin 
      nchunk = 0
      nleft = nrows
  endelse 

  print
  print,'Summing radlumcolor rows from file '+file  
  print
  print,'Chunksize = ',chunk
  print,'Nchunk    = ',nchunk
  print,'Remainder = ',nleft


  base_struct = self->radlumcolor_base_struct()

  ;; Read and sum the chunks
  beg = 0L
  nchstr = '/'+ntostr(nchunk)
  indices = lindgen(nrows)


  for i=0l, nchunk-1 do begin 
      if i eq 0 then begin 
          print,'  Reading Chunk: '+ntostr(i+1)+nchstr, format='($,a)'
      endif else begin 
          print,' '+ntostr(i+1)+nchstr, format='($,a)'
      endelse 

      tmp = self->read_chunk(lun, base_struct, chunk)

      print,'*',format='($,a)'

      if keep_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+chunk-1]
          self->_sum_match_keeparray, tmp, matchid, keeparray, st
      endif else if rows_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+chunk-1]
          self->_sum_match_rows, tmp, matchid, rows, st
      endif else begin 
          self->_sum_output, tmp, st
      endelse 

      tmp = 0
      beg = beg+chunk
  endfor 

  if nleft ne 0 then begin 
      print,'  Reading remainder',format='($,a)'

      tmp = self->read_chunk(lun, base_struct, nleft)
      print,'*',format='($,a)'

      if keep_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+nleft-1]
          self->_sum_match_keeparray, tmp, matchid, keeparray, st
      endif else if rows_input then begin 
          ;; matchid matches up with the chunk we have read
          matchid = indices[beg:beg+nleft-1]
          self->_sum_match_rows, tmp, matchid, rows, st
      endif else begin 
          self->_sum_output, tmp, st
      endelse 
      print

      tmp = 0
  endif 

  free_lun, lun

  ;; some derived quantities
  self->calc_sum_derived, st

  ptime,systime(1)-t0

  return, st

end 



function correlate::get_pair_index, file

  openr, lun, file, /get_lun

  nrows = 0l
  readu, lun, nrows

  print,'nrows = ',nrows

  nrows=nrows-1
  st = {index:0l, npairs:0l}
  st = replicate(st, nrows)

  readu, lun, st
  free_lun, lun
  return, st

end 






function correlate::sum_rad_output, file, rows=rowsin, keeparray=keeparray

  if n_elements(file) eq 0 then begin 
      on_error, 2
      print,'-Syntax: sumstruct = c->sum_rad_output(file, rows=, keeparray=)'
      print
      message,'Halting'
  endif 

  st = self->rad_sumstruct()

  hdr = read_idlheader(file)
  nrows = hdr.nrows

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Deal with requested subsets
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  self->_process_subset_requests, nrows, rowsin, keeparray, rows, st, rows_input, keep_input  

  print
  print,'Summing rad output from file: ',file
  struct = read_idlstruct(file, status=status)
  if status ne 0 then message,'Failed to read file'
  if keep_input then begin 
      indices = lindgen(nrows)
      self->_sum_match_keeparray, struct, indices, keeparray, st
  endif else if rows_input then begin 
      struct = struct[rows]
      self->_sum_output, struct, st
  endif else begin 
      self->_sum_output, struct, st
  endelse 

  ;; some derived quantities
  self->calc_sum_derived, st

  return, st
end













;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Perform subsampling
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; This will take a set of redshifts, bin according to the rz binning for this
;; sample, and return the indices in each bin.  

function correlate::rzsub_hist, rz, reverse_indices=reverse_indices
  if n_elements(rz) eq 0 then begin 
      on_error, 2
      print,'hist = corr->rzsub_hist(z, reverse_indices=)'
      print
      message,'Halting'
  endif 

  rzinfo = self->rzsub_info()
  h = histogram(rz, min=rzinfo.rzmin, max=rzinfo.rzmax, binsize=rzinfo.rzstep, $
                reverse_indices=reverse_indices)

  nh = n_elements(h)
  if nh lt rzinfo.nrz then message,'Found '+ntostr(nh)+' bins, expecting at least '+ntostr(rzinfo.nrz)

  ;; There is always one last bin we don't want
  return, h[0:rzinfo.nrz-1]
end 

; get indices for objects in each histogram bin
function correlate::rzsub_indices, rz, hist=hist

  if n_elements(rz) eq 0 then begin 
      on_error, 2
      print,'keeparray = corr->rzsub_indices(z, hist=)'
      print
      message,'Halting'
  endif 

  hist = self->rzsub_hist(rz, rev=rev)

  nbin = n_elements(hist)
  keep = ptrarr(nbin)

  for i=0l, nbin-1 do begin 

      if rev[i] ne rev[i+1] then begin 
          w = rev[ rev[i]:rev[i+1]-1 ]
          hist[i] = n_elements(w)
          keep[i] = ptr_new(w, /no_copy)
      endif 

  endfor 

  return, keep
end 

;; returns pointer array, each element corresponds to a n_z,n_j 
pro correlate::rzsub_jack_indices, rz, clam, ceta, $
             allbins, z_ids, jackknife_ids, jackknife_file=jackknife_file
             

  rkeep = self->rzsub_indices(rz)

  nr = n_elements(rkeep)

  for i=0l, nr-1 do begin 

      print
      print,'Getting jackknife for zbin = '+ntostr(i+1)+'/'+ntostr(nr)
      ind = *rkeep[i]
      tjk = self->get_jackknife_indices(ind, clam[ind], ceta[ind], $
                                       jackknife_file=jackknife_file, $
                                       jackknife_ids=tjid)

      nj = n_elements(tjk)
      for j=0l, nj-1 do begin 
          add_arrval, ptr_new(*tjk[j]), allbins
      endfor 
      add_arrval, tjid, jackknife_ids
      add_arrval, replicate(i, nj), z_ids

      ptr_free, tjk
  endfor 

  ptr_free, rkeep

end 

;; This is just for doing the redshift subsampling for rz and rr
pro correlate::rzsub_sample_old, type, $
             primary_randnum=primary_randnum, $
             secondary_randnum=secondary_randnum

  case strlowcase(type) of
      'rd': begin 

          nrand = n_elements(primary_randnum) 
          if nrand eq 0 then message,'You must enter a primary_randnum for "rd"'

          for i=0l, nrand-1 do begin 

              print,'------------------------------------------------------------'
              
              rfile = self->corrfile('rd', 'output', $
                                     primary_randnum=primary_randnum[i])
              
              outfiles = $
                self->corrfile('rd', 'rzsub', primary_randnum=primary_randnum[i], $
                               /createdir)
              
              rcorrin = self->corr_read('primary_random','input',primary_randnum=primary_randnum[i])

              print
              print,'Sub sampling into rzsub redshift bins'
              rkeep = self->rzsub_indices(rcorrin.z)

              ;; combine the sub samples
              sumstructs = self->sum_radlumcolor_output(rfile, keeparray=rkeep)
              
              ;; write the outputs
              self->objshear::write_sub_samples, sumstructs, outfiles

              ptr_free, rkeep

          endfor 

      end 
      'rr': begin 

          nprand = n_elements(primary_randnum)
          nsrand = n_elements(secondary_randnum)

          ;; These must pair up
          if nprand eq 0 then message,'You must enter a primary_randnum for "rr"'
          if nsrand eq 0 then message,'You must enter a secondary_randnum for "rr"'

          if nprand ne nsrand then message,'primary_randnum and secondary_randnum must be same size'


          for i=0l, nprand-1 do begin 

              print,'------------------------------------------------------------'
              print,'prand = '+ntostr(primary_randnum[i])+' srand = '+ntostr(secondary_randnum[i])
              rfile = self->corrfile('rr', 'output', $
                                     primary_randnum=primary_randnum[i], $
                                     secondary_randnum=secondary_randnum[i])

              outfiles = $
                self->corrfile('rr', 'rzsub', $
                               primary_randnum=primary_randnum[i], $
                               secondary_randnum=secondary_randnum[i], $
                               /createdir)

              rcorrin = self->corr_read('primary_random','input',$
                                        primary_randnum=primary_randnum[i])

              print
              print,'Sub sampling into rzsub redshift bins'
              rkeep = self->rzsub_indices(rcorrin.z)

              ;; combine the sub samples
              sumstructs = self->sum_rad_output(rfile, keeparray=rkeep)

              ;; Get the total density
              nbin = n_elements(sumstructs)
              for j=0l, nbin-1 do begin 
                  ind = *rkeep[j]
                  sumstructs[j].totaldensity = $
                    self->calc_totaldensity(rcorrin[ind].z, $
                                            secondary_randnum=secondary_randnum[i])
                  sumstructs[j].usedarea = $
                    sumstructs[j].radcounts/sumstructs[j].totaldensity
                  sumstructs[j].usedarea_err = $
                    sqrt(sumstructs[j].radcounts)/sumstructs[j].totaldensity
              endfor 

              
              ;; write the outputs
              self->objshear::write_sub_samples, sumstructs, outfiles

              ptr_free, rkeep
          endfor 



      end 
      else: message,'rzsub_sample works only with "rd" and "rr"'
  endcase 

end 

;; Add a header item which is the jackknife regions used
pro correlate::rzsub_fix_jackknifes

  info = self->rzsub_info()
  
  for zi=0l, info.nrz-1 do begin 
      
      print
      print,'zi = '+ntostr(zi+1)+'/'+ntostr(info.nrz)
      print,'--------------------------------------------------------------'

      files = self->list_random('rd_rzsub_jsample',rzbin=zi)
      nf = n_elements(files)
      for fi=0l, nf-1 do begin 
          t = read_idlstruct(files[fi])
          hdr = {jackknife_id: t.jackknife_id}
          print,'Writing to file: ',files[fi]
          write_idlstruct, t, files[fi], hdr=hdr
          t = 0
;          key = prompt_kbrd('hit a key')
      endfor 

      files = self->list_random('rr_rzsub_jsample',rzbin=zi)
      nf = n_elements(files)
      for fi=0l, nf-1 do begin 
          t = read_idlstruct(files[fi])
          hdr = {jackknife_id: t.jackknife_id}
          print,'Writing to file: ',files[fi]
          write_idlstruct, t, files[fi], hdr=hdr
          t = 0
;          key = prompt_kbrd('hit a key')
      endfor 


  endfor 

end 

;; get all the unique ids
function correlate::rzsub_read_jackknife_ids
    info = self->rzsub_info()  

    for zi=0l, info.nrz-1 do begin 
        files = self->list_random('rr_rzsub_jsample',rzbin=zi)
        nf = n_elements(files)
        for fi=0l, nf-1 do begin 
            t = read_idlstruct(files[fi], /silent)
            if tag_exist(t, 'jackknife_id') then begin
                tids = t.jackknife_id
            endif else begin       
                h = read_idlheader(files[fi], /silent)
                if tag_exist(h, 'jackknife_id') then begin
                    tids = h.jackknife_id
                endif else begin
                    message,'No jackknife ids found'
                endelse
            endelse
            add_arrval, tids, jids
        endfor 
    endfor 

    rmd = rem_dup(jids)
    jids = jids[rmd]
    return, jids
end 

;; read the rzsub redshift bins, but only a given jackknife id
function correlate::rzsub_read_by_jackknife_id, type, jid

  if n_elements(type) eq 0 or n_elements(jid) ne 1 then begin 
      on_error, 2
      print,'-Syntax: st = c->rzsub_read_by_jackknife_id(type, jid)'
      print,' type = "rd" or "rr"'
      message,'Halting'
  endif 

  if type ne 'rd' and type ne 'rr' then message,'type must be "rd" or "rr"'

  info = self->rzsub_info()

  for zi=0l, info.nrz-1 do begin 
      rz = self->corr_read(type, 'rzsub_jsample', rzbin=zi, /silent)
      
      w=where(rz.jackknife_id eq jid, nw)

      if nw eq 0 then begin 
          message,'jid '+ntostr(jid)+' did not match',/inf
      endif else begin 
          add_arrval, rz[w], outstruct
      endelse 

  endfor 
  if n_elements(outstruct) eq 0 then begin 
      message,'jackknife id '+ntostr(jid)+' did not match any files'
  endif 
  return, outstruct

end 

pro correlate::rzsub_combine, type

  case type of
      'rd': list_type = 'rd_rzsub'
      'rr': list_type = 'rr_rzsub'
      else: begin 
          on_error, 2
          message,'type must be "rd" or "rr"'
      end 
  endcase 

  info = self->rzsub_info()

  outfiles = self->corrfile(type, 'rzsub')

  ;; Each redshift bin will have it's own file still
  for zi=0l, info.nrz-1 do begin 
      
      print
      print,'zi = '+ntostr(zi+1)+'/'+ntostr(info.nrz)
      print,'--------------------------------------------------------------'

      files = self->list_random(list_type,rzbin=zi)
      nf = n_elements(files)

      ;; loop over the files
      for fi=0l, nf-1 do begin 
          print,'Reading file: ',files[fi]
          t = read_idlstruct(files[fi])

          ;; Only need to combine if this is not the first time
          ;; we have done this jackknife sample
          if fi eq 0 then begin 
              outst = t
          endif else begin 
              concat_dstructs, outst, t, tt

			  outst = self->combine_sumstructs(tt)
          endelse 
          t = 0

      endfor 

      ;; write out the combined file
      print
      print,'writing to file: ',outfiles[zi]
;      if fexist(outfiles[zi]) then message,'file exists'
      write_idlstruct, outst, outfiles[zi]


  endfor 

end 

pro correlate::rzsub_jackknife_combine, type

  case type of
      'rd': list_type = 'rd_rzsub_jsample'
      'rr': list_type = 'rr_rzsub_jsample'
      else: begin 
          on_error, 2
          message,'type must be "rd" or "rr"'
      end 
  endcase 

  jids = self->rzsub_read_jackknife_ids()
  nj = n_elements(jids)

  info = self->rzsub_info()

  outfiles = self->corrfile(type, 'rzsub_jsample')

  ;; each redshift bin will have it's own file still
  for zi=0l, info.nrz-1 do begin 
      
      print
      print,'zi = '+ntostr(zi+1)+'/'+ntostr(info.nrz)
      print,'--------------------------------------------------------------'

      files = self->list_random(list_type,rzbin=zi)
      nf = n_elements(files)

      ;; Loop over the files
      for fi=0l, nf-1 do begin 
          print,'Reading file: ',files[fi]
          t = read_idlstruct(files[fi])

          ;; Loop over the jackknife regions
          for ij=0,nj-1 do begin 

              ;; does this jackknife region exist in this file?
              w = where(t.jackknife_id eq jids[ij],nw)
              if nw ne 0 then begin 

                  w = w[0]

                  ;; first match, create output struct
                  if n_elements(outst) eq 0 then begin 
                      tmp = t[0]
                      zero_struct, tmp
                      outst = replicate(tmp, n_elements(jids))
                  endif 

                  ;; Only need to combine if this is not the first time
                  ;; we have done this jackknife sample
                  if outst[ij].npoints eq 0 then begin 
                      outst[ij] = t[w]
                  endif else begin 
                      tt = [outst[ij], t[w]]
                      outst[ij] = self->combine_sumstructs(tt)
                  endelse 

              endif 

          endfor 

          t = 0

      endfor 

      ;; Keep track of the jackknife id
      outst.jackknife_id = jids

      ;; Write out the combined file
      print
      print,'Writing to file: ',outfiles[zi]
;      IF fexist(outfiles[zi]) THEN message,'File exists'
      write_idlstruct, outst, outfiles[zi]

      ;; Need to clear this for the if npoints eq 0 statement above
      delvarx, outst

  endfor 

end 

function correlate::get_jackknife_rzsub, type, jackknife_id

  ;; This is the list of files for each rzbin.  We will just get
  ;; the particular jackknife regions because of memory problems.
  files = self->corrfile(type, 'rzsub_jsample')

  rzinfo = self->rzsub_info()
  case strlowcase(type) of
      'rd': begin 
          rzsubs = replicate(self->radlumcolor_sumstruct(), rzinfo.nrz)
          rztype = 'rd_rzsub'
      end 
      'rr': begin 
          rzsubs = replicate(self->rad_sumstruct(), rzinfo.nrz)
          rztype = 'rr_rzsub'
      end 
      else: message,'type must be rd or rr'
  endcase 

  for rzbin=0,rzinfo.nrz-1 do begin 
      rfiles = self->list_random(rztype, rzbin=rzbin)
      tmp = read_idlstruct_multi(rfiles)
      rzsubs[rzbin] = self->combine_sumstructs(tmp)
  endfor 
  return, rzsubs
end 



function correlate::jackknife_file
  tf = esheldon_config('jackknife_file')
  dirsep, tf, dir, file
  file = concat_dir(dir, 'jackknife_samples_dr4_res256_nstripe4.st') 
  return, file
end 

;; this is just for doing the redshift subsampling for rd and rr
pro correlate::rzsub_sample, type, $
             primary_randnum=primary_randnum, $
             secondary_randnum=secondary_randnum, $
             jackknife=jackknife

  if n_elements(type) eq 0 then begin 
      on_error, 2
      print,'-Syntax: c->rzsub_sample_jack, type, primary_randnum=, secondary_randnum='
      print
      message,'Halting'
  endif 

  case strlowcase(type) of
      'rd': begin 

          nrand = n_elements(primary_randnum) 
          if nrand eq 0 then message,'You must enter a primary_randnum for "rd"'

          for i=0l, nrand-1 do begin 

              print,'------------------------------------------------------------'

              ;; The input file to the correlation code to get redshifts
              rcorrin = self->corr_read('primary_random','input',primary_randnum=primary_randnum[i])
              
              ;; Output file from correlation code
              rfile = self->corrfile('rd', 'output', $
                                     primary_randnum=primary_randnum[i])
              
              ;; Output files.  If jackknife, we need to do each zsub
              ;; separately and ask for the jackknife samples for each zsub

              if not keyword_set(jackknife) then begin 
                  outfiles = $
                    self->corrfile('rd', 'rzsub', $
                                   primary_randnum=primary_randnum[i], /createdir)

                  ;; Get indices for each redshift bin                  
                  print
                  print,'Sub sampling into rzsub redshift bins'
                  rkeep = self->rzsub_indices(rcorrin.z)
                  
                  ;; combine the sub samples
                  sumstructs = self->sum_radlumcolor_output(rfile, keeparray=rkeep)
                  
                  ;; write the outputs
                  self->objshear::write_sub_samples, sumstructs, outfiles

                  ptr_free, rkeep
              endif else begin 
                  
                  eq2csurvey, rcorrin.ra, rcorrin.dec, clam, ceta
                  num = n_elements(rcorrin)

                  outfiles = $
                    self->corrfile('rd', 'rzsub_jsample', $
                                   primary_randnum=primary_randnum[i], $
                                   /createdir)

                  print
                  print,'Creating rzsub+jackknife indices'
                  self->rzsub_jack_indices, rcorrin.z, clam, ceta, $
                    allbins, z_ids, j_ids

                  ;; combine the sub samples
                  sumstructs = $
                    self->sum_radlumcolor_output(rfile, keeparray=allbins)

                  ;; Now output by redshift bin.
                  uz_ids = z_ids[ rem_dup(z_ids) ]
                  nuz = n_elements(uz_ids)
                  for iuz=0l, nuz-1 do begin 
                      zid = uz_ids[iuz]
                      w = where(z_ids eq zid, nw)

                      tstructs = sumstructs[w]
                      jackstructs = self->create_jackknife_sumstructs(tstructs)
                      jackstructs.jackknife_id = j_ids[w]
                      
                      print
                      print,'Writing to file: ',outfiles[iuz]
                      write_idlstruct, jackstructs, outfiles[iuz]
                  endfor 

                  ptr_free, allbins

              endelse 
          endfor 

      end 
      'rr': begin 

          nprand = n_elements(primary_randnum)
          nsrand = n_elements(secondary_randnum)

          ;; these must pair up
          if nprand eq 0 then message,'You must enter a primary_randnum for "rr"'
          if nsrand eq 0 then message,'You must enter a secondary_randnum for "rr"'

          if nprand ne nsrand then message,'primary_randnum and secondary_randnum must be same size'


          for i=0l, nprand-1 do begin 

              print,'------------------------------------------------------------'
              print,'prand = '+ntostr(primary_randnum[i])+' srand = '+ntostr(secondary_randnum[i])
              rfile = self->corrfile('rr', 'output', $
                                     primary_randnum=primary_randnum[i], $
                                     secondary_randnum=secondary_randnum[i])

              rcorrin = self->corr_read('primary_random','input',$
                                        primary_randnum=primary_randnum[i])

              if not keyword_set(jackknife) then begin 
                  outfiles = $
                    self->corrfile('rr', 'rzsub', $
                                   primary_randnum=primary_randnum[i], $
                                   secondary_randnum=secondary_randnum[i], $
                                   /createdir)

                  
                  print
                  print,'Sub sampling into rzsub redshift bins'
                  rkeep = self->rzsub_indices(rcorrin.z)
                  
                  ;; combine the sub samples
                  sumstructs = self->sum_rad_output(rfile, keeparray=rkeep)
                  
                  ;; Get the total density
                  nbin = n_elements(sumstructs)
                  for j=0l, nbin-1 do begin 
                      if ptr_valid(rkeep[j]) then begin 
                          ind = *rkeep[j]
                          sumstructs[j].totaldensity = $
                            self->calc_totaldensity(rcorrin[ind].z, $
                                                    secondary_randnum=secondary_randnum[i])
                          sumstructs[j].usedarea = $
                            sumstructs[j].radcounts/sumstructs[j].totaldensity
                          sumstructs[j].usedarea_err = $
                            sqrt(sumstructs[j].radcounts)/sumstructs[j].totaldensity
                      endif 
                  endfor 

              
                  ;; write the outputs
                  self->objshear::write_sub_samples, sumstructs, outfiles
                  
                  ptr_free, rkeep
              endif else begin 

                  eq2csurvey, rcorrin.ra, rcorrin.dec, clam, ceta
                  num = n_elements(rcorrin)

                  outfiles = $
                    self->corrfile('rr', 'rzsub_jsample', $
                                   primary_randnum=primary_randnum[i], $
                                   secondary_randnum=secondary_randnum[i],$
                                   /createdir)

                  print
                  print,'Creating rzsub+jackknife indices'
                  self->rzsub_jack_indices, rcorrin.z, clam, ceta, $
                    allbins, z_ids, j_ids

                  ;; combine the sub samples
                  sumstructs = self->sum_rad_output(rfile, keeparray=allbins)

                  print
                  print,'Calculating total density'
                  ;; Get the total density
                  num = n_elements(sumstructs)
                  for j=0l, num-1 do begin 
                      if ptr_valid(allbins[j]) then begin 
                          ind = *allbins[j]
                          sumstructs[j].totaldensity = $
                            self->calc_totaldensity(rcorrin[ind].z, $
                                                    secondary_randnum=secondary_randnum[i])
                          sumstructs[j].usedarea = $
                            sumstructs[j].radcounts/sumstructs[j].totaldensity
                          sumstructs[j].usedarea_err = $
                            sqrt(sumstructs[j].radcounts)/sumstructs[j].totaldensity
                      endif 
                  endfor 


                  ;; now output by redshift bin.
                  uz_ids = z_ids[ rem_dup(z_ids) ]
                  nuz = n_elements(uz_ids)
                  for iuz=0l, nuz-1 do begin 
                      zid = uz_ids[iuz]
                      w = where(z_ids eq zid, nw)
                      
                      tstructs = sumstructs[w]
                      jackstructs = self->create_jackknife_sumstructs(tstructs)

                      jackstructs.jackknife_id = j_ids[w]

                      print
                      print,'Writing to file: ',outfiles[iuz]
                      write_idlstruct, jackstructs, outfiles[iuz]

                  endfor 

                  ptr_free, allbins


              endelse 
          endfor 



      end 
      else: message,'rzsub sampling works only with "rd" and "rr"'
  endcase 

end 





;; Returns the weight needed such that two histograms have the same shape
;; h1,h2 are histograms for sample 1 and sample 2.
;; For us, sample will be the redshifts we want to match, sample2 will be
;; randoms
;; This weight is the percentage of h2 objects to keep from each bin.

function correlate::match_hist_weight, h1, h2

  n1 = n_elements(h1)
  n2 = n_elements(h2)

  if n1 eq 0 or n2 eq 0 then begin 
      on_error, 2
      print,'-Syntax: w = corr->match_hist_weight(h1, h2)'
      print
      message,'Halting'
  endif 

  if n1 ne n2 then message,'Two histograms must be same size'

  wh = where(h1 gt 0 and h2 gt 0, ngood)

  if ngood eq 0 then begin 
      message,'No good values found in histograms'
  endif 

  ;; The percentage of x2 to keep in each bin
  perc_keep = fltarr(n1)
  ratio = fltarr(n1)

  ;; normalize by the maximum in x1 histogram
  ratio[wh] = double( h1[wh] )/double( h2[wh] )
  maxratio = max(ratio)

  ;; The percentage to keep
  perc_keep[wh] = ratio[wh]/maxratio

  return, perc_keep

end 


function correlate::rzsub_match, z, rzsubs, doplot=doplot
  if n_elements(z) eq 0 or n_elements(rzsubs) eq 0 then begin 
      on_error, 2
      print,'-Syntax: st = c->rzsub_match(z, rzsubs, /doplot)'
      print
      message,'Halting'
  endif 

  ;; The number of points used is just the redshift histogram
  rzh = rzsubs.npoints

  ;; histogram zs
  zh = self->rzsub_hist(z)

  ;; Weight to get proper correction factors for input redshift hist
  weight = self->match_hist_weight(zh, rzh)  

  if keyword_set(doplot) then begin 
      plot, zh/total(zh), psym=10
      oplot, rzh*weight/total(rzh*weight), color=c2i('green'), psym=10, line=2
      key = prompt_kbrd('hit a key')
  endif 

  ;; Now combine with weight
  struct = self->wcombine_sumstructs(rzsubs, weight)
  return, struct

end 
function correlate::rzsub_match_keeparray, z, keep, rzsubs

  if n_params() lt 3 then begin 
      on_error, 2
      print,'-Syntax: sums = corr->rzsub_match(z, keep, rzsubs)'
  endif 
  nbin = n_elements(keep)

  if keyword_set(radlumcolor) then begin 
      out = self->radlumcolor_sumstruct()
  endif else begin 
      out = self->rad_sumstruct()
  endelse 
  out = replicate(out, nbin)

  ;; The number of points used is just the redshift histogram
  rzh = rzsubs.npoints

  for i=0l, nbin-1 do begin 

      if ptr_valid(keep[i]) then begin 
          ind = *keep[i]

          tout = self->rzsub_match(z[ind], rzsubs, /doplot)
          
          if i eq 0 then begin 
              out = replicate(tout, nbin)
          endif 
          out[i] = tout
      endif 
  endfor 

  return, out
end 



















function correlate::get_input_index, inputstruct
  if tag_exist(inputstruct, 'bcg_id', ind=ind) then begin 
      return, ind
  endif else if tag_exist(inputstruct, 'index', ind=ind) then begin 
      return, ind
  endif else begin 
      message,'no index found'
  endelse 
end 




;; This is for sub-sampling relative to the primaries.
;; There is another sub-sampler for sub-sampling the random
;; redshifts into redshift bins, see above.
;;
;; For /jackknife, the jackknife samples are created.  These
;; are the total-jsample, not the samples themselves, and can
;; be used directly in the covariance calculation.
PRO correlate::sub_sample, type, subtype, $
             zbinsize=zbinsize, jackknife=jackknife

  IF n_elements(type) EQ 0 OR n_elements(subtype) EQ 0 THEN BEGIN 
      on_error, 2
      print,'-Syntax: c->sub_sample, type, subtype, /jackknife'
      print
      message,'You must enter the type and subtype strings.  Halting'
  ENDIF 

  ;; main catalog
  struct = self->get()

  ;; lens input file. Just to make sure input and output are the same size
  corrin = self->corr_read('primary','input')
  
  ;; only need those that made it into the lens input catalogs
  indind = self->get_input_index(corrin)
  index = corrin.(indind)
  struct = struct[index]
  
  keep = self->objshear::where_select(struct, subtype, nkeep=nkeep, nbin=nbin)
  struct = 0

  CASE strlowcase(type) OF
      'dd': BEGIN 
          ;; will autumatically select all bins
          file = self->corrfile('dd', 'output')

          IF NOT keyword_set(jackknife) THEN BEGIN 
              outfiles = self->corrfile(type, 'combined', subtype=subtype, /createdir)
			  ;ptr_free, keep
			  ;help,outfiles
			  ;stop
          
              ;; combine the sub samples
              sumstructs = self->sum_radlumcolor_output(file, keeparray=keep)
              
              ;; write the outputs
              self->objshear::write_sub_samples, sumstructs, outfiles
          ENDIF ELSE BEGIN 
              outfiles = self->corrfile(type, 'jsample', $
                                        subtype=subtype, /createdir)

              eq2csurvey, corrin.ra, corrin.dec, clam, ceta

              FOR bin=0L, nbin-1 DO BEGIN 
                  print
                  print,'---------------------------------------------------'
                  IF ptr_valid(keep[bin]) THEN BEGIN 

                      ind = *keep[bin]
                      jk = self->get_jackknife_indices(ind, clam[ind], ceta[ind], $
                                                       jackknife_ids=jids)
                      njack = n_elements(jk)

                      print,'Extracting '+ntostr(njack)+' jackknife regions'
                      
                      ;; Sum structs in each jackknife region
                      dd_jsamp = self->sum_radlumcolor_output(file, keeparray=jk)
                      ptr_free, jk

                      ;; Now create summed minus one 
                      print,'Creating summed jackknife samples'
                      dd_jsamp = self->create_jackknife_sumstructs(dd_jsamp)
                      ;; Keep track of jackknife ids for later
                      dd_jsamp.jackknife_id = jids

                      print
                      print,'Writing to file: ',outfiles[bin]
                      write_idlstruct, dd_jsamp, outfiles[bin]

                 ENDIF 

              ENDFOR 
          ENDELSE 

      END 

      'dr': BEGIN 

          ;; Any primary index arrays will be same for each of the random
          ;; secondaries.

          drfiles = self->list_random('dr', count=nfiles)

          IF NOT keyword_set(jackknife) THEN BEGIN 
              
              starray = replicate( self->rad_sumstruct(), nfiles, nbin )
              
              
              FOR i=0L, nfiles-1 DO BEGIN 
                  
                  print,'------------------------------------------------------------'
                  file = drfiles[i]
                  
                  secondary_randnum = long( stregex(file,'[0-9][0-9].st',/extract) )
                  
                  print,'Secondary randnum = ',secondary_randnum
                  ;; combine the sub samples
                  sumstructs = self->sum_rad_output(file, keeparray=keep)
                  
                  ;; Get the total density
                  FOR j=0L, nbin-1 DO BEGIN 
                      IF ptr_valid(keep[j]) THEN BEGIN 
                          ind = *keep[j]
                          sumstructs[j].totaldensity = $
                            self->calc_totaldensity(corrin[ind].z, $
                                                    secondary_randnum=secondary_randnum)
                          sumstructs[j].usedarea = $
                            sumstructs[j].radcounts/sumstructs[j].totaldensity
                          sumstructs[j].usedarea_err = $
                            sqrt(sumstructs[j].radcounts)/sumstructs[j].totaldensity
                      ENDIF 
                  ENDFOR 
                  
                  starray[i,*] = sumstructs
                  
              ENDFOR 
              
              sumstructs = replicate( self->rad_sumstruct(), nbin )
              FOR i=0L, nbin-1 DO BEGIN 
                  sumstructs[i] = self->combine_sumstructs(starray[*,i])
              ENDFOR 
              
              ;; write the outputs
              outfiles = self->corrfile('dd', 'matchdr', subtype=subtype, /createdir)
              self->objshear::write_sub_samples, sumstructs, outfiles
              

          ENDIF ELSE BEGIN 
              
              outfiles = self->corrfile('dd', 'matchdr_jsample', $
                                        subtype=subtype, /createdir)
              
              eq2csurvey, corrin.ra, corrin.dec, clam, ceta
              FOR bin=0L, nbin-1 DO BEGIN 
                  
                  print
                  print,'-----------------------------------------------------------'
                  
                  IF ptr_valid(keep[bin]) THEN BEGIN 
                      
                      ind = *keep[bin]

                      jk = self->get_jackknife_indices(ind, clam[ind], ceta[ind], $
                                                       jackknife_ids=jackids)
                      njack = n_elements(jk)

                      print,'Extracting '+ntostr(njack)+' jackknife regions from files'
                      
                      ;; Loop over files and keep the running sum
                      FOR fi=0L, nfiles-1 DO BEGIN 
                          file = drfiles[fi]
                          
                          secondary_randnum = long( stregex(file,'[0-9][0-9].st',/extract) )

                          ;; Sum structs in each jackknife region
                          tdr_jsamp = self->sum_rad_output(file, keeparray=jk)

                          FOR ji=0L, njack-1 DO BEGIN 
                              IF ptr_valid(jk[ji]) THEN BEGIN 

                                  jind = *jk[ji]

                                  tdr_jsamp[ji].totaldensity = $
                                    self->calc_totaldensity(corrin[jind].z, $
                                                            secondary_randnum=secondary_randnum)
                                  tdr_jsamp[ji].usedarea = $
                                    tdr_jsamp[ji].radcounts/tdr_jsamp[ji].totaldensity
                                  tdr_jsamp[ji].usedarea_err = $
                                    sqrt(tdr_jsamp[ji].radcounts)/tdr_jsamp[ji].totaldensity
                              ENDIF 
                          ENDFOR 
                          
                          ;; Now create summed minus one 
                          print,'Creating summed jackknife samples'
                          tdr_jsamp = self->create_jackknife_sumstructs(tdr_jsamp)
                          
                          ;; Keep the running sum
                          IF fi EQ 0 THEN BEGIN 
                              dr_jsamp = tdr_jsamp
                          ENDIF ELSE BEGIN 
                              FOR ji=0L, njack-1 DO BEGIN 
                                  tt = [dr_jsamp[ji], tdr_jsamp[ji]]
                                  dr_jsamp[ji] = self->combine_sumstructs(tt)
                              ENDFOR 
                          ENDELSE 
                          
                      ENDFOR 
                      ptr_free, jk
                      
                      ;; Keep track of jackknife ids for later
                      dr_jsamp.jackknife_id = jackids
                      print
                      print,'Writing to file: ',outfiles[bin]
                      write_idlstruct, dr_jsamp, outfiles[bin]
                  ENDIF 
                  
              ENDFOR 
              
              
          ENDELSE 
      END 


      ;; These require histogram matching
      'rd': BEGIN 

          ;; Start a plot
          pdir = self->plotdir(subtype=subtype,/matchzrand,/createdir)
          psfile = 'rd_'+self->sample('dd')+'_'+subtype+'_matchzrand'
          IF keyword_set(jackknife) THEN psfile = psfile + '_jack'
          psfile = psfile + '_N1.ps'
          psfile = concat_dir(pdir, psfile)

          ;; don't overwrite
          WHILE fexist(psfile) DO psfile=newname(psfile)
          
          begplot,name=psfile, /color
          !p.multi = [0,0,2]


          IF NOT keyword_set(jackknife) THEN BEGIN 
              rzsubs = self->corr_read('rd','rzsub')
              
              ;; Now combine over rzbins with the appropriate weight
              print
              print,'Matching redshift histograms'
              sumstructs = self->rzsub_match_keeparray(corrin.z, keep, rzsubs)
              
              ;; write the outputs
              outfiles = self->corrfile('dd', 'matchrd', subtype=subtype, /createdir)
              
              self->objshear::write_sub_samples, sumstructs, outfiles
          ENDIF ELSE BEGIN 

              outfiles = self->corrfile('dd', 'matchrd_jsample', $
                                        subtype=subtype, /createdir)


              ;; Get jackknife id for each primary
              eq2csurvey, corrin.ra, corrin.dec, clam, ceta
              jids = csurvey2jack(clam, ceta, file=self->jackknife_file())


              ;; These are already combined to have a jackknife region
              ;; *removed*
              print,'Reading jackknife/redshift samples'
              rzsubs = self->corr_read('rd', 'rzsub_jsample')

              ;; Unique jackknife ids
              ujids = rzsubs[ rem_dup(rzsubs.jackknife_id) ].jackknife_id
              njack = n_elements(ujids)



              ;; Loop over subsample bins
              FOR bin=0L, nbin-1 DO BEGIN 
                  print,'bin = '+ntostr(bin+1)+'/'+ntostr(nbin)
                  IF ptr_valid(keep[bin]) THEN BEGIN 

                      ind = *keep[bin]

                      ;; Loop over jackknife regions.
                      back = backspace(5)
                      maxstr = strn(max(ujids)+1, len=2, padchar='0')
                      print,'Jackknife id = 00/'+maxstr, format='($,a)'
                      FOR ji=0L, njack-1 DO BEGIN 

                          ;; We need to remove the jackknife region to
                          ;; get proper redshift distribution
                          wjack = where( jids[ind] NE ujids[ji] )
                          wjack = ind[wjack]

                          print,back,format='($,a)'
                          print,strn(ujids[ji]+1,len=2,padchar='0')+'/'+maxstr, $
                            format='($,a)'

                          ;; This equal because they already have the jackknife
                          ;; removed
                          w = where(rzsubs.jackknife_id EQ ujids[ji])
                          IF ji EQ 0 THEN doplot=1 ELSE doplot=0
                          tout = self->rzsub_match(corrin[wjack].z, rzsubs[w], doplot=doplot)

                          add_arrval, tout, outstructs
                      ENDFOR 
                      print
                  ENDIF 

                  ;; Write out this subsample
                  print,'Writing file: ',outfiles[bin]
                  write_idlstruct, outstructs, outfiles[bin]

                  delvarx, outstructs
              ENDFOR 

          ENDELSE 

          !p.multi=0
          endplot

      END 

      'rr': BEGIN 

          ;; Start a plot
          pdir = self->plotdir(subtype=subtype,/matchzrand,/createdir)
          psfile = 'rr_'+self->sample('dd')+'_'+subtype+'_matchzrand'
          IF keyword_set(jackknife) THEN psfile = psfile + '_jack'
          psfile = psfile + '_N1.ps'
          psfile = concat_dir(pdir, psfile)

          ;; don't overwrite
          WHILE fexist(psfile) DO psfile=newname(psfile)
          
          begplot,name=psfile, /color
          !p.multi = [0,0,2]

          IF NOT keyword_set(jackknife) THEN BEGIN 

              ;; Get all the randoms in each rzbin, total them up.
              print,'Combining rzsubs'
              rzsubs = self->corr_read('rr','rzsub')              

              ;; Now combine over rzbins with the appropriate weight
              print
              print,'Matching redshift histograms'
              sumstructs = self->rzsub_match_keeparray(corrin.z, keep, rzsubs)

          
              ;; write the outputs
              outfiles = self->corrfile('dd', 'matchrr', subtype=subtype, /createdir)
              
              self->objshear::write_sub_samples, sumstructs, outfiles

          ENDIF ELSE BEGIN 

              ;; This is exactly the same code as the rd sample except
              ;; replaced with 'rr'
              outfiles = self->corrfile('dd', 'matchrr_jsample', $
                                        subtype=subtype, /createdir)

              ;; Get jackknife id for each primary
              eq2csurvey, corrin.ra, corrin.dec, clam, ceta
              jids = csurvey2jack(clam, ceta, file=self->jackknife_file())

              print,'Reading jackknife/redshift samples'
              rzsubs = self->corr_read('rr', 'rzsub_jsample')
              ujids = rzsubs[ rem_dup(rzsubs.jackknife_id) ].jackknife_id
              njack = n_elements(ujids)

              ;; Loop over subsample bins
              FOR bin=0L, nbin-1 DO BEGIN 
                  print,'bin = '+ntostr(bin+1)+'/'+ntostr(nbin)
                  IF ptr_valid(keep[bin]) THEN BEGIN 

                      ind = *keep[bin]

                      ;; Loop over jackknife regions.
                      back = backspace(5)
                      maxstr = strn(max(ujids)+1, len=2, padchar='0')
                      print,'Jackknife id = 00/'+maxstr, format='($,a)'
                      FOR ji=0L, njack-1 DO BEGIN 

                          ;; We need to remove the jackknife region to
                          ;; get proper redshift distribution
                          wjack = where( jids[ind] NE ujids[ji] )
                          wjack = ind[wjack]

                          print,back,format='($,a)'
                          print,strn(ujids[ji]+1,len=2,padchar='0')+'/'+maxstr, $
                            format='($,a)'

                          ;; This equal because they already have the jackknife
                          ;; removed
                          w = where(rzsubs.jackknife_id EQ ujids[ji])
                          IF ji EQ 0 THEN doplot=1 ELSE doplot=0
                          tout = self->rzsub_match(corrin[wjack].z, rzsubs[w], doplot=doplot)

                          add_arrval, tout, outstructs
                      ENDFOR 
                      print
                  ENDIF 

                  ;; Write out this subsample
                  print,'Writing file: ',outfiles[bin]
                  write_idlstruct, outstructs, outfiles[bin]

                  delvarx, outstructs
              ENDFOR 

          ENDELSE 

          !p.multi=0
          endplot

      END 


      ELSE: message,'type not supported: '+ntostr(type)
  ENDCASE 

  ;; Free sub-sample pointers memory
  ptr_free, keep

END 



function correlate::meanstruct, struct, tags, index=index

  nstruct = n_elements(struct)
  ntags = n_elements(tags)
  if nstruct eq 0 or ntags eq 0 then begin 
      on_error, 2
      print,'-Syntax: meanst = c->meanstruct(struct,tags,index=)'
      message,'Halting'
  endif 

  if n_elements(index) eq 0 then index=lindgen(n_elements(struct))

  for i=0l, ntags-1 do begin 

      ;; strip trailing [index] or (index)
      bracket_pos  = stregex(tags[i], '[\[\(]')
      bracket_pos2 = stregex(tags[i], '[]\)]')
      if bracket_pos ne -1 then begin 

          tag = strmid(tags[i], 0, bracket_pos)

          ;; Now extract the element for our name
          element_str = '_'+strmid(tags[i], $
                                   bracket_pos+1, bracket_pos2-bracket_pos-1)
      endif else begin 
          tag = tags[i]
          element_str = ''
      endelse 

      IF tag_exist(struct[0], tag, index=itag) THEN BEGIN 

          command = 'mom = moment( struct[index].'+tags[i]+', sdev=sdev_tag, /nan)'
          IF NOT execute(command) THEN message,'Unable to take moments of tag '+tags[i]

          command = 'min_tag = min( struct[index].'+tags[i]+', max=max_tag, /nan)'
          IF NOT execute(command) THEN message,'Unable to take min/max of tag '+tags[i]

          command = 'median_tag = median( struct[index].'+tags[i]+')'
          IF NOT execute(command) THEN message,'Unable to take median of tag '+tags[i]


          mean_tag = mom[0]
          err_tag = sdev_tag/sqrt(nstruct)


          mean_name = tag+element_str+'_mean'
          sdev_name = tag+element_str+'_sdev'
          err_name  = tag+element_str+'_err'
          median_name = tag+element_str+'_median'
          min_name  = tag+element_str+'_min'
          max_name  = tag+element_str+'_max'

          IF n_elements(meanstruct) EQ 0 THEN BEGIN 
              meanstruct = create_struct(mean_name, mean_tag, $
                                         sdev_name, sdev_tag, $
                                         err_name, err_tag, $
                                         median_name, median_tag, $
                                         min_name, min_tag, $
                                         max_name, max_tag)

          ENDIF ELSE BEGIN 
              meanstruct = create_struct(meanstruct, $
                                         mean_name, mean_tag, $
                                         sdev_name, sdev_tag, $
                                         err_name, err_tag, $
                                         median_name, median_tag, $
                                         min_name, min_tag, $
                                         max_name, max_tag)
          ENDELSE 

      ENDIF ELSE BEGIN 
          message,'Tag '+ntostr(tag)+' does not exist in structure',/inf
      ENDELSE 

  ENDFOR 

  return, meanstruct

end 

pro correlate::means, subtype=subtype, tags=tags

  struct = self->get()
  
  if n_elements(tags) eq 0 then tags=self->average_tags(subtype=subtype)

  if n_elements(subtype) ne 0 then begin 
      keep = self->objshear::where_select(struct, subtype, nkeep=nkeep, $
                                          nbin=nbin, labels=labels)
  endif else begin 
      keep = ptr_new( lindgen(nstruct) )
  endelse 

  outfiles = self->corrfile('dd','means', subtype=subtype, /createdir)
  n = n_elements(outfiles)

  for i=0L, n-1 do begin 
      
      tmeanstruct = self->meanstruct(struct, tags, index=*keep[i])
      add_arrval, tmeanstruct, meanstruct

  endfor 

  ptr_free, keep
  self->objshear::write_sub_samples, meanstruct, outfiles

end 























FUNCTION correlate::model_usedarea, struct, max_model_rad, max_use_rad, doplot=doplot

  IF n_params() LT 3 THEN BEGIN 
      on_error, 2
      print,'-Syntax: struct = c->model_usedarea(structarray, max_model_rad, max_use_rad)'
      print,'Returns a copy of input struct with modeled usedarea for small radii'
      print
      message,'Halting'
  ENDIF 

  degree = 5
  nst = n_elements(struct)
  newstruct = struct


  xtit=textoidl('R [h^{-1} Mpc]')
  FOR i=0L, nst-1 DO BEGIN 

      ;; Get radii
      r = struct[i].r
      wfit = where(r LT max_model_rad, nwfit)
      IF nwfit EQ 0 THEN message,'No good model radii for i='+ntostr(i)

      wuse = where(r LT max_use_rad, nwuse)
      IF nwuse EQ 0 THEN message,'No good use radii for i='+ntostr(i)

      fa = struct[i].usedarea/struct[i].area
      faerr = struct[i].usedarea_err/struct[i].area

      ;; Fit polynumial
      res = FitAreaPoly(r[wfit], fa[wfit], faerr[wfit], degree=degree, yfit=yfit)

      ;; Fix the usedarea
      newstruct[i].usedarea[wuse] = yfit[wuse]*struct[i].area[wuse]

      ;; Plot the result
      IF keyword_set(doplot) THEN BEGIN 
          rp = arrscl( findgen(1000), min(r[wuse]), max(r[wuse]) )
          yfit = AreaPoly(rp, res)

          pplot, r, fa, yerr=faerr, psym=8, /xlog, aspect=1.0, $
              ytitle='Fractional Area', xtitle=xtit
          oplot, rp, yfit, thick=2, color=c2i('darkGreen')
          key = prompt_kbrd('hit a key')
      ENDIF 
  ENDFOR 

  return, newstruct

END 
FUNCTION correlate::model_usedarea_plotfile, subtype=subtype
  pdir = self->plotdir(subtype=subtype,/createdir)
  psfile = 'model_usedarea_dr_'+self->sample('dd')
  IF n_elements(subtype) NE 0 THEN psfile = psfile +'_'+subtype
  psfile = psfile + '.ps'
  psfile = concat_dir(pdir, psfile)
  return, psfile
END 
FUNCTION correlate::usedarea_plotfile, subtype=subtype, bin=bin
  pdir = self->plotdir(subtype=subtype,/createdir)
  psfile = 'usedarea_dr_'+self->sample('dd')
  ext='.ps'
  IF n_elements(subtype) NE 0 THEN begin
      psfile = psfile +'_'+subtype
      nb=n_elements(bin)
      if nb eq 1 then begin
          psfile = psfile + '_'+ntostr(bin)
          ext='.eps'
      endif else begin
          ext='.ps'
      endelse
  endif
  psfile = psfile + ext
  psfile = path_join(pdir, psfile)
  return, psfile
END 




FUNCTION correlate::_correct, dd, rd, dr, rr

	if n_params() lt 4 then begin 
		on_error, 2
		print,'-Syntax: corrstruct = c->_correct(dd, rd, dr, rr)'
		print
		message,'Halting'
	endif 

	st = self->radlumcolor_sumstruct(/corrected)

	copy_struct, dd, st

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Get the mean color and error
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	self->sums2meanerr, dd.kgflux, dd.kgflux2, dd.npoints, dd_kgf, v, dd_kgferr
	self->sums2meanerr, dd.krflux, dd.krflux2, dd.npoints, dd_krf, v, dd_krferr
	self->sums2meanerr, dd.kiflux, dd.kiflux2, dd.npoints, dd_kif, v, dd_kiferr

	self->sums2meanerr, rd.kgflux, rd.kgflux2, rd.npoints, rd_kgf, v, rd_kgferr
	self->sums2meanerr, rd.krflux, rd.krflux2, rd.npoints, rd_krf, v, rd_krferr
	self->sums2meanerr, rd.kiflux, rd.kiflux2, rd.npoints, rd_kif, v, rd_kiferr

	;; subtract counts and keep track of errors
	self->calc_adderr, dd_kgf, dd_kgferr, -rd_kgf, rd_kgferr, kgf, kgferr
	self->calc_adderr, dd_krf, dd_krferr, -rd_krf, rd_krferr, krf, krferr
	self->calc_adderr, dd_kif, dd_kiferr, -rd_kif, rd_kiferr, kif, kiferr  

	st.kgmr = -2.5*alog10(kgf/krf)
	st.kgmr_err = self->fluxerr2colorerr(kgf, kgferr, krf, krferr )
	st.krmi = -2.5*alog10(krf/kif)
	st.krmi_err = self->fluxerr2colorerr(krf, krferr, kif, kiferr )

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; Densities
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	;; Radial 
	;; lum

	st.radilumdens_field = rd.radilum/rr.usedarea
	st.radilumdens_field_err = rd.radilum_err/rr.usedarea
	st.radilumdens     = dd.radilum/dr.usedarea - st.radilumdens_field
	st.radilumdens_err = $
		sqrt( (dd.radilum_err/dr.usedarea)^2 + st.radilumdens_field_err^2 )


	;; number
	ddcounts = double(dd.radcounts)/dd.npoints
	ddcounts_err = sqrt(double(dd.radcounts))/dd.npoints

	rdcounts = double(rd.radcounts)/rd.npoints
	rdcounts_err = sqrt(double(rd.radcounts))/rd.npoints

	st.radnumdens_field = rdcounts/rr.usedarea
	st.radnumdens_field_err = rdcounts_err/rr.usedarea
	st.radnumdens      = ddcounts/dr.usedarea - st.radnumdens_field
	st.radnumdens_err  = $
		sqrt( (ddcounts_err/dr.usedarea)^2 + st.radnumdens_field_err^2 )

	;; now number but with cumulative sums over luminosity
	ddcounts_cl = double(dd.radcounts_cumlum)/dd.npoints
	ddcounts_cl_err = sqrt(double(dd.radcounts_cumlum))/dd.npoints

	rdcounts_cl = double(rd.radcounts_cumlum)/rd.npoints
	rdcounts_cl_err = sqrt(double(rd.radcounts_cumlum))/rd.npoints

	for ilum=0L, dd[0].nlum-1 do begin
		st.radnumdens_cumlum_field[*,ilum] = rdcounts_cl[*,ilum]/rr.usedarea
		st.radnumdens_cumlum_field_err[*,ilum] = $
			rdcounts_cl_err[*,ilum]/rr.usedarea
		st.radnumdens_cumlum[*,ilum] = $
			ddcounts_cl[*,ilum]/dr.usedarea - st.radnumdens_cumlum_field[*,ilum]
		st.radnumdens_cumlum_err[*,ilum]  = $
			sqrt( (ddcounts_cl_err[*,ilum]/dr.usedarea)^2 + $
				  st.radnumdens_cumlum_field_err[*,ilum]^2 )
	endfor



	;; means in all bins
	for i=0l, dd.nrad-1 do begin 

		st.ilumdens_field[i,*,*] = rd.ilum[i,*,*]/rr.usedarea[i]
		st.ilumdens_field_err[i,*,*] = rd.ilum_err[i,*,*]/rr.usedarea[i]
		st.ilumdens[i,*,*] = $
			dd.ilum[i,*,*]/dr.usedarea[i] - st.ilumdens_field[i,*,*]
		st.ilumdens_err[i,*,*] = $
			sqrt( (dd.ilum_err[i,*,*]/dr.usedarea[i])^2 + $
			st.ilumdens_field_err[i,*,*]^2 )

		ddcounts = double(dd.counts[i,*,*])/dd.npoints
		ddcounts_err = sqrt(double(dd.counts[i,*,*]))/dd.npoints

		rdcounts = double(rd.counts[i,*,*])/rd.npoints
		rdcounts_err = sqrt(double(rd.counts[i,*,*]))/rd.npoints

		st.numdens_field[i,*,*] = rdcounts/rr.usedarea[i]
		st.numdens_field_err[i,*,*] = rdcounts_err/rr.usedarea[i]
		st.numdens[i,*,*]  = ddcounts/dr.usedarea[i] - st.numdens_field[i,*,*]
		st.numdens_err[i,*,*]  = $
			sqrt( (ddcounts_err/dr.usedarea[i])^2 + $
			st.numdens_field_err[i,*,*]^2 )


	endfor 

  return, st
END 



;; For now just dd-rd
PRO correlate::correct, subtype=subtype, $
             max_model_rad=max_model_rad, max_use_rad=max_use_rad, $
             jackknife=jackknife

  dd = self->corr_read('dd','combined',subtype=subtype)
  dr = self->corr_read('dd','matchdr', subtype=subtype)
  rd = self->corr_read('dd','matchrd', subtype=subtype)
  rr = self->corr_read('dd','matchrr', subtype=subtype)  

  print
  print,'Modeling small scale UsedArea'
  psfile = self->model_usedarea_plotfile(subtype=subtype)

  begplot, name=psfile, /color
  IF n_elements(max_model_rad) EQ 0 THEN max_model_rad = 4.0 ; Mpc
  IF n_elements(max_use_rad) EQ 0 THEN max_use_rad = 3.0
  dr = self->model_usedarea(dr, max_model_rad, max_use_rad, /doplot)
  endplot

  IF NOT keyword_set(jackknife) THEN BEGIN 
      outfiles = self->corrfile('dd','corrected',subtype=subtype,/createdir)
  ENDIF ELSE BEGIN 
      outfiles = self->corrfile('dd','jackknife',subtype=subtype,/createdir)
  ENDELSE 

  nbin = n_elements(dd)
  FOR bin=0L, nbin-1 DO BEGIN 

      print,'--------------------------------------------------------------------'

      ddcorr = self->_correct(dd[bin], rd[bin], dr[bin], rr[bin])

      IF keyword_set(jackknife) THEN BEGIN 
          jdd = self->corr_read('dd', 'jsample', subtype=subtype, bin=bin)
          jdr = self->corr_read('dd', 'matchdr_jsample', subtype=subtype, bin=bin)
          jrd = self->corr_read('dd', 'matchrd_jsample', subtype=subtype, bin=bin)
          jrr = self->corr_read('dd', 'matchrr_jsample', subtype=subtype, bin=bin)

          jdr = self->model_usedarea(jdr, max_model_rad, max_use_rad)          

          ncov = self->jackknife_covariance(jdd, jdr, jrd, jrr, err=nerr)
          lcov = self->jackknife_covariance(jdd, jdr, jrd, jrr, err=lerr, $
                                            /luminosity)

          ddcorr.radnumdens_err = nerr
          ddcorr.radilumdens_err = lerr

          ddcorr = create_struct(ddcorr, $
                                 'radnumdens_cov', ncov, $
                                 'radilumdens_cov', lcov)

      ENDIF 

      print
      print,'Writing to file: ',outfiles[bin]
      write_idlstruct, ddcorr, outfiles[bin]

  ENDFOR 

END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Jackknifing
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION correlate::get_index, struct
  ;; I stupidly made the index called bcg_id the first time around
  IF tag_exist(struct, 'bcg_id') THEN BEGIN 
      index = struct.bcg_id
  ENDIF ELSE IF tag_exist(index) THEN BEGIN 
      index = struct.index
  ENDIF ELSE BEGIN 
      message,'No index found'
  ENDELSE 
  return, index
END 


;; Return the indices of objects in each jackknife region
;; optionaly return the jackknife ids for each region found
FUNCTION correlate::get_jackknife_indices, inputid, clambda, ceta, $
                  jackknife_file=jackknife_file, jackknife_ids=jackknife_ids

  IF n_elements(jackknife_file) EQ 0 THEN jackknife_file = self->jackknife_file()
  jid = csurvey2jack(clambda, ceta, file=jackknife_file)

  minjack = min(jid, max=maxjack)
  IF minjack LT 0 THEN message,'Some not found in a jackknife region'

  h = histogram(jid-minjack, min=0, rev=rev)

  ;; How many unique jackknife regions do we have?
  wh = where(h NE 0, Nsub)

  jack_keeparray = ptrarr(Nsub)
  jackknife_ids = lonarr(Nsub)

  ;; Now aggregate the indices by jackknife region
  FOR sub=0L, Nsub-1 DO BEGIN 

      iSub = wh[sub]
      wSub = rev[ rev[iSub]:rev[iSub+1]-1 ]
      ind = inputid[ wSub ]

      jack_keeparray[sub] = ptr_new(ind, /no_copy)
      jackknife_ids[sub] = jid[ wSub[0] ]
  ENDFOR 

  return, jack_keeparray
END 


;; Take a set of jackknife sumstructs (either radlumcolor or rad) and
;; produce the sumstructs for the total samples (with each jackknife region
;; removed in turn)
FUNCTION correlate::create_jackknife_sumstructs, struct

  IF n_elements(struct) EQ 0 THEN BEGIN 
      on_error, 2
      print,'-Syntax: covariance=c->jackknife_covariance(struct)'
      print
      message,'Halting'
  ENDIF 

  tt=systime(1)

  Njack = n_elements(struct)
  Nrad = n_elements(struct[0].radcounts)
  
  allind = lindgen(Njack)
  FOR i=0L, Njack-1 DO BEGIN 

      ind = allind
      remove, i, ind

      ;; Get sums with this jackknife region removed. This works for both
      ;; radlumcolor and rad
      tstruct = self->combine_sumstructs(struct[ind])

      IF i EQ 0 THEN BEGIN 
          outstruct = replicate(tstruct, Njack)
      ENDIF 
      
      outstruct[i] = tstruct

  ENDFOR 

  return, outstruct

END 

;; Test the jackknife samples for a subtype
PRO correlate::jackknife_test, corrmat
  
  subtype = 'ilum200_16'
  bin = 0
  jdd = self->corr_read('dd','jsample',subtype=subtype, bin=bin)
  dd = self->corr_read('dd','combined',subtype=subtype, bin=bin)


  jdr = self->corr_read('dd','matchdr_jsample',subtype=subtype, bin=bin)
  dr = self->corr_read('dd','matchdr',subtype=subtype, bin=bin)

  jrd = self->corr_read('dd','matchrd_jsample',subtype=subtype, bin=bin)
  rd = self->corr_read('dd','matchrd',subtype=subtype, bin=bin)

  jrr = self->corr_read('dd','matchrr_jsample',subtype=subtype, bin=bin)
  rr = self->corr_read('dd','matchrr',subtype=subtype, bin=bin)

  
  simpctable, colorlist=clist
  nclr = n_elements(clist)

  ;; dd
  pplot, dd.r, dd.radnumdens, yerr=dd.radnumdens_err, psym=8, /xlog,/ylog, $
    title = 'dd', charsize=2

  FOR i=0L, n_elements(jdd)-1 DO BEGIN 
      pplot, jdd[i].r, jdd[i].radnumdens, color=clist[i MOD nclr], /overplot
  ENDFOR 

  key = prompt_kbrd('hit a key')

  ;; dr
  pplot, dr.r, dr.radnumdens, yerr=dr.radnumdens_err, psym=8, /xlog, $
    title = 'dr', charsize=2, /ynozero

  FOR i=0L, n_elements(jdr)-1 DO BEGIN 
      pplot, jdr[i].r, jdr[i].radnumdens, color=clist[i MOD nclr], /overplot
  ENDFOR 

  key = prompt_kbrd('hit a key')


  ;; rd
  pplot, rd.r, rd.radnumdens, yerr=rd.radnumdens_err, psym=8, /xlog, $
    title = 'rd', charsize=2, /ynozero

  FOR i=0L, n_elements(jrd)-1 DO BEGIN 
      pplot, jrd[i].r, jrd[i].radnumdens, color=clist[i MOD nclr], /overplot
  ENDFOR 

  key = prompt_kbrd('hit a key')

  ;; rr
  pplot, rr.r, rr.radnumdens, yerr=rr.radnumdens_err, psym=8, /xlog, $
    title = 'rr', charsize=2, /ynozero

  FOR i=0L, n_elements(jrr)-1 DO BEGIN 
      pplot, jrr[i].r, jrr[i].radnumdens, color=clist[i MOD nclr], /overplot
  ENDFOR 

  key = prompt_kbrd('hit a key')

  ;; corr
  corr = self->corr_read('dd','corrected', subtype=subtype, bin=bin)
  cov = self->jackknife_covariance_new(jdd,jdr,jrd,jrr,meanvals=meanvals,err=err)

  pplot, corr.r, corr.radnumdens, yerr=corr.radnumdens_err, psym=8, /xlog, /ylog,$
    title = 'rr', charsize=2
  pplot, corr.r, corr.radnumdens, yerr=err, color=c2i('green'), /overplot


  colprint,corr.radnumdens,meanvals,corr.radnumdens_err,err

  key = prompt_kbrd('hit a key')
  calc_correlation_matrix, cov, corrmat

  image_contour, corrmat, $
    levels=lindgen(10)/10.0,proc='tvim2', color=c2i('red'), $
    c_labe=replicate(1,18), /scale

  
  key = prompt_kbrd('hit a key')

  pplot, corr.r, err/corr.radnumdens, /xlog
  pplot, corr.r, corr.radnumdens_err/corr.radnumdens, $
    /xlog, color=c2i('green'), /overplot
return


  sz = size(corrmat)
  index=lindgen(sz[0]*sz[1])

  x=index MOD sz[0]
  y=index/sz[1]
  r = sqrt(x^2 + y^2)

  FOR i=0L, sz[1]-1 DO BEGIN 
      IF i EQ 0 THEN BEGIN 
          plot, corrmat[i,*], psym=8
          oplot, corrmat[i,*]
      ENDIF ELSE BEGIN 
          oplot, corrmat[i,*], color=clist[i], psym=8
          oplot, corrmat[i,*], color=clist[i]
      ENDELSE 
  ENDFOR 
END 

;; All structs must match up by jackknife id
;; If modeling usedare, must be done before calling this function
FUNCTION correlate::jackknife_covariance, dd, dr, rd, rr, meanvals=meanvals, err=err, luminosity=luminosity, sigma_clip=sigma_clip

  np = n_elements(dd) + n_elements(dr)+n_elements(rd)+n_elements(rr)
  IF np LT 4 THEN BEGIN 
      on_error, 2
      print,'-Syntax: covariance=c->jackknife_covariance(dd, dr, err=, meanvals=meanvals, /luminosity, /sigma_clip)'
      print
      message,'Halting'
  ENDIF 

  tt=systime(1)

  NJack = n_elements(dd)
  Nrad = n_elements(dd[0].radcounts)
  
  print,'NJack = ',NJack
  print,'Nrad = ',Nrad

  ;; Outputs: mean, covariance and diagonal terms
  err = dblarr(Nrad)
  covariance = dblarr(Nrad, Nrad)

  sampval = dblarr(Nrad,NJack)

  FOR i=0L, NJack-1 DO BEGIN 

      IF NOT keyword_set(luminosity) THEN BEGIN 
          dd_radcounts = double(dd[i].radcounts)
          rd_radcounts = double(rd[i].radcounts)
          sampval[*,i] = $
            dd_radcounts[*]/dd[i].npoints/dr[i].usedarea[*] - $
            rd_radcounts[*]/rd[i].npoints/rr[i].usedarea[*]
      ENDIF ELSE BEGIN 
          sampval[*,i] = $
            dd[i].radilum[*]/dr[i].usedarea[*] - rd[i].radilum[*]/rr[i].usedarea[*]
      ENDELSE 

  ENDFOR 

  ;; Now jackknife covariance between variables
  nsigma1 = 4.0
  nsigma2 = 3.5
  niter=3
  FOR jRad=0L, Nrad-1 DO BEGIN 

      IF keyword_set(sigma_clip) THEN BEGIN 
          sigma_clip, sampval[jRad,*], smean, ssigma, nsig=nsigma1, nIter=3, /silent
          jw = where(sampval[jRad,*] GT smean-nsigma2*ssigma AND $
                     sampval[jRad,*] LT smean+nsigma2*ssigma)
      ENDIF ELSE BEGIN 
          jw = lindgen(NJack)
      ENDELSE 
      jmean = mean_check(sampval[jRad,jw])

      add_arrval, jmean, meanvals

      FOR iRad=jRad, Nrad-1 DO BEGIN 

          IF keyword_set(sigma_clip) THEN BEGIN 
              sigma_clip, sampval[iRad,*], smean, ssigma, nsig=nsigma1, nIter=niter, /silent
              iw = where(sampval[iRad,*] GT smean-nsigma2*ssigma AND $
                         sampval[iRad,*] LT smean+nsigma2*ssigma)
          ENDIF ELSE BEGIN 
              iw = lindgen(NJack)
          ENDELSE 
          imean = mean_check(sampval[iRad,iw])

          tmp = total( (sampval[iRad,iw]-imean)*(sampval[jRad,jw]-jmean) )
         
          covariance[jRad, iRad] = tmp*(NJack - 1.)/NJack

          IF jRad NE iRad THEN covariance[iRad, jRad] = covariance[jRad, iRad]
          IF jRad EQ iRad THEN err[iRad] = sqrt( covariance[iRad, iRad] )

      ENDFOR 
  ENDFOR 

  return, covariance

END 








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Invert to 3D profile
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Need to adopt this to work with number density or lum somehow
function correlate::nfwfit, rmpc, mass, mass_err, meanz

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
         bias_yfit: bias_yfit $
       }

  return, st

end 


pro correlate::invert, subtype=subtype, corrected=corrected, dops=dops

  if keyword_set(corrected) then begin 
      dtype = 'corrected'
  endif else begin 
      dtype = 'jackknife'
  endelse 

  out_dtype = dtype + '_invert'
  t = self->corr_read('dd', dtype, subtype=subtype)

  struct = self->get(nrows=nstruct)
  if n_elements(subtype) ne 0 then begin 
      keep = self->objshear::where_select(struct, subtype, nkeep=nkeep, nbin=nbin)
      ws = self->where_string(subtype, labels=labels)
  endif else begin 
      keep = ptr_new( lindgen(nstruct) )
  endelse 
  
  psfile = self->plotfile(out_dtype, subtype=subtype)
  psfile = repstr(psfile, out_dtype, out_dtype+'_qa')
  if keyword_set(dops) then begplot, psfile, /color

  outfiles = self->corrfile('dd', out_dtype, subtype=subtype, /createdir)

  nbin = n_elements(t)
  nrad = n_elements(t[0].r)

  nstr = ntostr(nrad-1)
  arrval = 'dblarr('+nstr+')'
  cov_arrval = 'dblarr('+nstr+','+nstr+')'
  newtags = ['ir',$
             'iradnumdens', 'iradnumdens_err', 'iradnumdens_cov', 'iradcounts', 'iradcounts_err', 'iradcounts_cov', $
             'iradilumdens', 'iradilumdens_err', 'iradilumdens_cov', 'iradilum', 'iradilum_err', 'iradilum_cov', $
             'bcg_ilum_mean','bcg_ilum_err','label']
  tagvals = [arrval, $
             arrval, arrval, cov_arrval, arrval, arrval, cov_arrval, $
             arrval, arrval, cov_arrval, arrval, arrval, cov_arrval, $
             '0d','0d','""']
  
  add_tags, t, newtags, tagvals, outst
    
  dp = obj_new('deproject')
  
  for i=0L, nbin-1 do begin 

      r = t[i].r
      proj_numdens = t[i].radnumdens
      proj_lumdens = t[i].radilumdens
      if keyword_set(corrected) then begin 
          proj_numdens_cov = diagonal_array(t[i].radnumdens_err)
          proj_lumdens_cov = diagonal_array(t[i].radilumdens_err)
      endif else begin 
          proj_numdens_cov = t[i].radnumdens_cov
          proj_lumdens_cov = t[i].radilumdens_cov
      endelse 


      nst = dp->invert(r, proj_numdens, proj_numdens_cov, /corr2d)
      lst = dp->invert(r, proj_lumdens, proj_lumdens_cov, /corr2d)

      outst[i].ir = nst.ir

      outst[i].iradnumdens = nst.drho
      outst[i].iradnumdens_err = nst.drho_err
      outst[i].iradnumdens_cov = nst.drho_cov

      ; Note, adding 1 for BCG.  This is not consistent and should be
	  ; changed
      outst[i].iradcounts = nst.massout
      outst[i].iradcounts_err = nst.massout_err
      outst[i].iradcounts_cov = nst.massout_cov

      outst[i].iradilumdens = lst.drho
      outst[i].iradilumdens_err = lst.drho_err
      outst[i].iradilumdens_cov = lst.drho_cov


      outst[i].iradilum = lst.massout
      outst[i].iradilum_err = lst.massout_err
      outst[i].iradilum_cov = lst.massout_cov

      ;; Adding mean lum of central BCGs
      w = *keep[i]
      nw = n_elements(w)

      mom = moment(struct[w].bcg_ilum,/double)
      meanlum = mom[0]
      meanlum_err = sqrt( mom[1]/nw )

      outst[i].bcg_ilum_mean = meanlum
      outst[i].bcg_ilum_err = meanlum_err

      if n_elements(labels) ne 0 then outst[i].label = labels[i]

      key = prompt_kbrd('hit a key')
      if key eq 'q' then begin 
          ptr_free, keep
          return
      endif 
  ENDFOR 

  obj_destroy, dp
  ptr_free, keep

  if keyword_set(dops) then endplot

  self->objshear::write_sub_samples, outst, outfiles

end 







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Some plotting routines.  Need to be generalized
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION correlate::mplot_value, nbin

  CASE nbin OF
      1: message,"nbin=1 doesn't make sense FOR multiplot"
      2: mplot_value=[1,2]
      3: mplot_value=[3,1]
      6: mplot_value=[3,2]
      8: mplot_value=[4,2]
      9: mplot_value=[3,3]
      12: mplot_value=[4,3]
      14: mplot_value=[4,4]
      16: mplot_value=[4,4]
	  18: mplot_value = [5,4]
      ELSE: message,'unknown nbin: '+ntostr(nbin)
  ENDCASE 
  return,mplot_value

END 

PRO correlate::plot_usedarea, dtype, subtype=subtype, bin=bin, dops=dops, _extra=_extra

    IF n_elements(dtype) EQ 0 THEN BEGIN 
        on_error, 2
        print,'-Syntax: co->plot_usedarea, dtype, subtype=, bin=, /dops, _extra='
        print
        message,'Halting'
    ENDIF 

    IF dtype NE 'matchdr' AND dtype NE 'matchrr' THEN BEGIN 
        message,'dtype must be matchdr or matchrr'
    ENDIF 

    IF n_elements(subtype) NE 0 AND n_elements(bin) EQ 0 THEN BEGIN 
        message,'for subtypes, send a bin number or use ::plot_usedarea_multi'
    ENDIF 

    if keyword_set(dops) then begin
        psfile=self->usedarea_plotfile(subtype=subtype, bin=bin)
        begplot, psfile, /encapsulated
    endif
    stsum = self->corr_read('dd',dtype,subtype=subtype,bin=bin)
    nst=n_elements(stsum)

	xtit=textoidl('R [h^{-1} Mpc]')
    pplot, stsum[0].r, stsum[0].usedarea/stsum[0].area, $
        yerr=stsum[0].usedarea_err/stsum[0].area, $
        psym=8, $
        /xlog,/ynozero,yrange=yrange,thick=2,xrange=xrange, xstyle=3,ystyle=3, $
        xtitle=xtit, ytitle='Fractional Area', $
        xtickf='loglabels', aspect=1, charsize=2, _extra=_extra

    max_model_rad = 4.0 ; Mpc
    max_use_rad = 3.0
    ;; new usedarea will be the model
    dr = self->model_usedarea(stsum[0], max_model_rad, max_use_rad)

    w=where(dr.r le max_use_rad)
    pplot, dr.r[w], dr.usedarea[w]/dr.area[w], /overplot

    if keyword_set(dops) then endplot, /trim_bbox

END 





PRO correlate::plot_usedarea_multi, dtype, subtype

  IF n_elements(dtype) EQ 0 OR n_elements(subtype) EQ 0 THEN BEGIN 
      on_error, 2
      print,'-Syntax: co->plot_usedarea, dtype, subtype, bin='
      print
      message,'Halting'
  ENDIF 

  IF dtype NE 'matchdr' AND dtype NE 'matchrr' THEN BEGIN 
      message,'dtype must be matchdr or matchrr'
  ENDIF 

  nbin = self->subtype_nbin(subtype)

  mpv = self->mplot_value(nbin)

  ws=self->where_string(subtype, labels=labels)

  st = self->corr_read('dd',dtype,subtype=subtype,bin=bin)

  !p.multi = [0,mpv] & !p.charsize=2
  simpctable, color=color
  FOR bin=0L, nbin-1 DO BEGIN 

      IF bin EQ 0 THEN BEGIN 
          ratio = st.usedarea/st.area
          r = st.r
          w=where(r GT 0 AND ratio GT 0)
          ;;yrange = [min(ratio[w]) > 0.8, max(ratio[w]) < 1.2]
          xrange = [0.5*min(r[w]), max(r[w])*1.5]
          yrange=[0.9,1.05]
      ENDIF 
      
      sti = st[bin]

	  xtit=textoidl('R [h^{-1} Mpc]')
      mplot, sti.r, sti.usedarea/sti.area, $
        /xlog,/ynozero,yrange=yrange,thick=2,xrange=xrange, xstyle=1,ystyle=1, $
        mxtitle=xtit, mytitle='Used Area/Area', $
        xtickf='loglabels', xticklen=0.04
      oploterror, sti.r, sti.usedarea/sti.area, sti.usedarea_err/sti.area

      oplot, [1.e-5, 1.e5], [1,1]
      legend, labels[bin],/right,box=0, charsize=1

  ENDFOR 

  !p.multi=0 & !p.charsize=1

END 






PRO correlate::plot_profile, dtype, subtype, bin, tag, reverse_labels=reverse_labels, noerr=noerr, dops=dops, dopng=dopng,ytitle=ytitle, yrange=yrange, xrange=xrange, aspect=aspect, lcharsize=lcharsize, _extra=_extra

    if n_params() lt 4 then begin
        on_error, 2
        print,'-Syntax: m->plot_profile, dtype, subtype, bin, tag, /reverse_labels, /noerr, /dops, /dopng, ytitle=, yrange=, xrange=, aspect=, lcharsize=, _extra='
        print
		print,'dtype = "dd"'
		print,'subtype = "ngals200_12"'
        message,'Halting'
    endif

    IF keyword_set(dops) or keyword_set(dopng) THEN BEGIN 
        message,'ps png not yet supported'
    ENDIF 

    ddc=self->corr_read('dd', dtype, subtype=subtype, bin=bin)

    IF NOT tag_exist(ddc,tag,index=ti) THEN message,'No such tag: '+ntostr(tag)
    IF NOT tag_exist(ddc,tag+'_err',index=tei) THEN message,'No such error tag: '+ntostr(tag+'_err')

    xlog=1
    ylog=1
    IF tag EQ 'kgmr' AND tag EQ 'krmi' THEN ylog=0

    r = ddc.r
    t = ddc.(ti)
    w=where( r GT 0 AND t GT 0 )
    if n_elements(xrange) eq 0 then xrange = [0.5*min(r[w]), 1.5*max(r[w])]
    if n_elements(yrange) eq 0 then yrange = [0.5*min(t[w]), 1.5*max(t[w])]
    if n_elements(aspect) eq 0 then aspect=1

    if n_elements(ytitle) eq 0 then ytitle=tag

	xtit=textoidl('R [h^{-1} Mpc]')
    if not keyword_set(noerr) then yerr=ddc[0].(tei)
    pplot, ddc.r, ddc.(ti), yerr=yerr, $
        /xlog, /ylog, xtickf='loglabels', ytickf='loglabels', $
        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
        xtitle=xtit, ytitle=ytitle, aspect=aspect, $
        _extra=_extra

    ws = self->where_string(subtype, labels=labels)

    IF keyword_set(reverse_labels) THEN BEGIN 
        labels=reverse(labels)
        lcolors=reverse(lcolors)
    ENDIF 

    if n_elements(lcharsize) eq 0 then lcharsize=1
    legend, labels[bin], /right, box=0, charsize=lcharsize, lin=0


END 





PRO correlate::plot_profile_loop, dtype, subtype, tag, color=color, reverse_labels=reverse_labels, noerr=noerr, dops=dops, dopng=dopng,ytitle=ytitle, yrange=yrange, xrange=xrange, aspect=aspect, lcharsize=lcharsize, _extra=_extra

    ws = self->where_string(subtype, labels=labels)
    nbin=n_elements(ws)
    for i=0L, nbin-1 do begin
        self->plot_profile, dtype, subtype, i, tag, reverse_labels=reverse_labels, noerr=noerr, dops=dops, dopng=dopng,ytitle=ytitle, yrange=yrange, xrange=xrange, aspect=aspect, lcharsize=lcharsize, _extra=_extra
        key = prompt_kbrd('hit a key')
    endfor
end


; Plot a single cluster bin but all luminosity thresholds
pro correlate::plot_radnumdens_cumlum, dtype, subtype, bin, band, band_shift, absmag=absmag, reverse_labels=reverse_labels, noerr=noerr, dops=dops, dopng=dopng,ytitle=ytitle, yrange=yrange, xrange=xrange, aspect=aspect, lcharsize=lcharsize, _extra=_extra

    if n_params() lt 3 then begin
        on_error, 2
        print,'-Syntax: m->plot_profile, dtype, subtype, bin, tag, /noerr, /dops, /dopng, ytitle=, yrange=, xrange=, aspect=, lcharsize=, _extra='
        print
		print,'dtype = "dd"'
		print,'subtype = "ngals200_12"'
        message,'Halting'
    endif

    if keyword_set(dops) or keyword_set(dopng) then begin 
        message,'ps png not yet supported'
    endif 

    ddc=self->corr_read('dd', dtype, subtype=subtype, bin=bin)

	tag = 'radnumdens_cumlum'
    if not tag_exist(ddc,tag, index=ti) then begin
		message,'No such tag: '+ntostr(tag)
	endif
	errtag = tag+'_err'
    if not tag_exist(ddc,errtag,index=tei) then begin
		message,'No such error tag: '+ntostr(errtag)
	endif

    xlog=1
    ylog=1
    IF tag EQ 'kgmr' AND tag EQ 'krmi' THEN ylog=0

    r = ddc.r
    t = ddc.(ti)
    w=where( r GT 0 AND t[*,0] GT 0 )
    if n_elements(xrange) eq 0 then $
		xrange = [0.5*min(r[w]), 50*max(r[w])]
    if n_elements(yrange) eq 0 then $
		yrange = [0.001*min(t[w,0]), 1.5*max(t[w,0])]
    if n_elements(aspect) eq 0 then aspect=1

	xtit=textoidl('R [h^{-1} Mpc]')
	ytitle=textoidl('Projected N(R)')

	y = ddc.(ti)
    if not keyword_set(noerr) then yerr=ddc[0].(tei)
    pplot, ddc.r, y[*,0], yerr=yerr[*,0], $
        /xlog, /ylog, xtickf='loglabels', ytickf='loglabels', $
        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
        xtitle=xtit, ytitle=ytitle, aspect=aspect, $
        _extra=_extra, /nohat

	colors=make_rainbow(ddc.nlum)
	for ilum=0L, ddc.nlum-1 do begin
		pplot, ddc.r, y[*,ilum], yerr=yerr[*,ilum], color=colors[ilum], $
			/overplot, /nohat
	endfor




	limits = ddc.loglbins_min
    if keyword_set(absmag) then begin
		if n_elements(band) eq 0 or n_elements(band_shift) eq 0 then begin
			message,'You must enter a band and band_shift for absmag'
		endif
        solarmags = k_solar_magnitudes(band_shift=band_shift, /silent)
        limits = -2.5*limits + solarmags[band]
		labels = ntostr(limits, f='(f6.2)')
		labels[0] = 'M < '+labels[0]
    endif else begin
		labels = ntostr(limits, f='(f6.2)')
		labels[0] = 'L > '+labels[0]
	endelse
    legend, labels, /right, box=0, charsize=lcharsize, lin=0, color=colors

    ws = self->where_string(subtype, labels=blabels)
    legend, blabels[bin], /right, /bottom, box=0, charsize=lcharsize

end 





PRO correlate::plot_profile_over, dtype, subtype, tag, color=color, reverse_labels=reverse_labels, noerr=noerr, dops=dops, dopng=dopng,ytitle=ytitle, yrange=yrange, xrange=xrange, aspect=aspect, lcharsize=lcharsize, _extra=_extra

    if n_params() lt 3 then begin
        on_error, 2
        print,'-Syntax: m->plot_profile_over, dtype, subtype, tag, /color, /reverse_labels, /noerr, /dops, /dopng, ytitle=, yrange=, xrange=, aspect=, lcharsize=, _extra='
        print
        message,'Halting'
    endif

    ;; !!! NEED TO UPDATE DIRECTORY NAMEING SCHEME AND THIS TOO
    par = self->par_struct()
    dir = self->plotdir(subtype=subtype)
    plotfile = par.DPsample+'_'+par.DSsample+'_'+subtype+'_'+tag
    plotfile = concat_dir(dir, plotfile)
    IF keyword_set(dops) THEN BEGIN 
        if keyword_set(color) then plotfile=plotfile+'_color'
        plotfile = plotfile+'.eps'
        begplot, plotfile, /encap, color=color, xsize=8.5
        !p.charsize=2
    ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
        plotfile = plotfile+'.png'
    ENDIF 

    ddc=self->corr_read('dd', dtype, subtype=subtype)
    nbin = n_elements(ddc)

    IF NOT tag_exist(ddc,tag,index=ti) THEN message,'No such tag: '+ntostr(tag)
    IF NOT tag_exist(ddc,tag+'_err',index=tei) THEN message,'No such error tag: '+ntostr(tag+'_err')

    xlog=1
    ylog=1
    IF tag EQ 'kgmr' AND tag EQ 'krmi' THEN ylog=0

    r = ddc.r
    t = ddc.(ti)
    w=where( r GT 0 AND t GT 0 )
    if n_elements(xrange) eq 0 then xrange = [0.5*min(r[w]), 1.5*max(r[w])]
    if n_elements(yrange) eq 0 then yrange = [0.5*min(t[w]), 1.5*max(t[w])]
    if n_elements(aspect) eq 0 then aspect=1

    if keyword_set(color) or !d.name eq 'X' then begin
        colors=make_rainbow(nbin)
    endif else begin
        loadct, 0
        colors = arrscl(findgen(nbin), 0, 200)
    endelse

    if n_elements(ytitle) eq 0 then ytitle=tag
    if n_elements(ymin) eq 0 then ymin = 1.e-5

    if not keyword_set(noerr) then begin
        yerr=ddc[0].(tei)
    endif

    xtitle = textoidl('R [h^{-1} Mpc]')
    pplot, ddc[0].r, ddc[0].(ti), yerr=yerr, $
        /xlog, /ylog, xtickf='loglabels', ytickf='loglabels', $
        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
        xtitle=xtitle, ytitle=ytitle, aspect=aspect, $
        _extra=_extra

    FOR i=0L, nbin-1 DO BEGIN 
        if not keyword_set(noerr) then yerr=ddc[i].(tei)
        pplot, ddc[i].r, ddc[i].(ti), yerr=yerr, /overplot, $
            color=colors[i],/ylog

        add_arrval, colors[i], lcolors
    ENDFOR 

    ws = self->where_string(subtype, labels=labels)

    IF keyword_set(reverse_labels) THEN BEGIN 
        labels=reverse(labels)
        lcolors=reverse(lcolors)
    ENDIF 

    if n_elements(lcharsize) eq 0 then lcharsize=1
    legend, labels, /bottom, box=0, charsize=lcharsize, lin=0, color=lcolors


    IF keyword_set(dops) THEN BEGIN 
        endplot, /trim_bbox
    ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
        write_png, plotfile, tvrd(/true)
    ENDIF 


END 



function correlate::radformat, rad
    if rad lt 0.1 then begin
        format = '(F5.3)'
    endif else if rad lt 1 then begin
        format = '(F4.2)'
    endif else if rad lt 10 then begin
        format = '(F3.1)'
    endif else if rad lt 100 then begin
        format = '(F4.1)' 
    endif else begin
        format = '(I)'
    endelse
    return, format
end

PRO correlate::plot_lumcolor_vs_rad, dtype, subtype, bin, band_shift, $
        dops=dops

    np = n_elements(dtype) + n_elements(subtype) + n_elements(bin) + $
            n_elements(band_shift)
    if np lt 4 then begin
        on_error, 2
        print,'-Syntax: contour_lumcolor_vs_rad, dtype, subtype, bin, band_shift, /dops'
        print
        message,'Halting'
    endif

    if keyword_set(dops) then begin
        file = self->plotfile(dtype, sub=subtype, bin=bin, $
            /createdir,/encapsulated)
        file = repstr(file,'.eps','_colorlum_vs_rad.eps')
        begplot,file,/encapsulated,xsize=8,ysize=7
    endif

    ddcorr=self->corr_read('dd',dtype,subtype=subtype, bin=bin)

    xrange = [min(ddcorr.loglbins_min), max(ddcorr.loglbins_max)]
    yrange = [min(ddcorr.kgmrbins_min), max(ddcorr.kgmrbins_max)]

    bshift = ntostr(band_shift,format='(f4.2)')

    mytitle = textoidl('^{'+bshift+'}(g-r)')
    mxtitle = textoidl('log(L_{'+bshift+'i}/L_{\odot}h^2)')
    if !d.name eq 'PS' then begin
        !p.charsize=0.9
        !p.charthick=2
        !x.thick=2
        !y.thick=2
        !p.thick=2  
    endif else begin
        !p.charsize=1
    endelse
  
    levels = 2+2*findgen(50)
    ;levels = [1,2,3,4,5,6,7,8,9]
    ;levels = [levels, 10*(1+lindgen(10))]
    ;levels = lindgen(100)
    ;levels = 1+lindgen(10)
    ;;labels = ntostr( long(levels) )
    colors=c2i('red')

    loglum = (ddcorr.loglbins_min+ddcorr.loglbins_max)/2.0
    gmr = (ddcorr.kgmrbins_min+ddcorr.kgmrbins_max)/2.0

    nx = ddcorr.nlum
    ny = ddcorr.ngmr
    xmin = (ddcorr.loglbins_min[0]+ddcorr.loglbins_max[0])/2
    xmax = (ddcorr.loglbins_min[nx-1]+ddcorr.loglbins_max[nx-1])/2
    ymin = (ddcorr.kgmrbins_min[0]+ddcorr.kgmrbins_max[0])/2
    ymax = (ddcorr.kgmrbins_min[ny-1]+ddcorr.kgmrbins_max[ny-1])/2

    titsize=2
    ncols = 5
    nrows = 4
	mpval = self->mplot_value( ddcorr.nrad )
    erase & multiplot, mpval, /square, $
        mxTitle=mxtitle, mxTitSize=titsize, mxTitOffset=0.8, $
        mytitle=mytitle, myTitSize=titsize

    for i=0,ddcorr.nrad-1 do begin 

        im = reform(ddcorr.numdens[i,*,*]) > 0
        err = reform(ddcorr.numdens_err[i,*,*])
        w=where(im gt 0)
        ;im[w] = im[w]/err[w]

        im = im/max(im)
        ;alpha = 0.1
        implot, im, type=type, image_scaled=image_scaled,$
            alpha=alpha, xrange=[xmin,xmax], yrange=[ymin,ymax], $
            xticklen=0.04, yticklen=0.04


        ; Accumulate across the luminosity axis
        cvals = (ddcorr.kgmrbins_min + ddcorr.kgmrbins_max)/2
        cdist = total(im, 1)
        histogram_points, cvals, cdist, xp, yp,  type='vertical'
        minlum = 11.0
        xp = xmax - xp/max(xp)*(xmax-minlum)
        oplot, xp, yp

        ; Across the g-r axis
        lvals = (ddcorr.loglbins_min+ddcorr.loglbins_max)/2
        ldist = total(im,2)
        ldist = ldist/max(ldist)

        minshow = 0.001
        w=where(ldist gt minshow,ncomp=ncomp,comp=comp)
        ldist[w] = alog10(ldist[w])
        histogram_points, lvals[w], ldist[w], xp, yp

        maxgmr = 0.5
        yp = yp - min(yp)
        yp = yp/max(yp)*maxgmr
        oplot, xp, yp


        ;rad = (ddcorr.radbins_min[i]+ddcorr.radbins_max[i])/2
        rad = rnd(ddcorr.r[i], 3)
        format = self->radformat(rad)
        rad = ntostr(rad,format=format)
        if rad eq '0.100' then rad = '0.10'
        legend, $
            'r = '+rad+' Mpc', charsize=0.75,/right,box=0,margin=0

        if i lt (ddcorr.nrad-1) then begin
            multiplot
        endif

    endfor 

    multiplot,/default

    if keyword_set(dops) then endplot,/trim_bbox
end 

function correlate::color_fractions, dd, redcut

    cvals = (dd.kgmrbins_min + dd.kgmrbins_max)/2
        
    st = { $
            allcount:0.0, allerr:0.0, $
            bluecount:0.0, blueerr:0.0,$
            redcount:0.0, rederr:0.0, $
            bluefrac:0.0, bluefracerr: 0.0 $
         }
    st=replicate(st, dd.nrad)

    for i=0,dd.nrad-1 do begin 

        numdens = reform(dd.numdens[i,*,*])
        err = reform(dd.numdens_err[i,*,*])
 
        cdist = total(numdens, 1)
        cdisterr = sqrt( total(err^2,1) )

        wblue = where(cvals lt redcut, comp=wred)

        st[i].allcount = total(cdist)
        st[i].allerr = sqrt( total(cdisterr^2) )

        st[i].bluecount = total(cdist[wblue])
        st[i].blueerr = sqrt( total(cdisterr[wblue]^2) )

        st[i].redcount = total(cdist[wred])
        st[i].rederr = sqrt( total(cdisterr[wred]^2) )

        st[i].bluefrac = st[i].bluecount/st[i].allcount
        st[i].bluefracerr = $
            st[i].bluefrac*sqrt( (st[i].blueerr/st[i].bluecount)^2 + $
                                 (st[i].rederr/st[i].redcount)^2 )
    endfor

    return, st

end


; Plot the mean galaxy luminosity versus radius
; Calculate mean lum in one of two ways: 
;  1) ratio the lumdens to the numdens.  
;  2) Calculate the mean numdens histogram vs. lum.  Should have less variance.

pro correlate::calc_meanlum_ratio, dd, ratio, ratio_err
    ratio=dd.radilumdens/dd.radnumdens*1.e10
    ratio_err = abs(ratio)*sqrt( (dd.radilumdens_err/dd.radilumdens)^2 + $
                                 (dd.radnumdens_err/dd.radnumdens)^2 )
end

function correlate::get_lumfunc, dd, radbin
    if n_params() lt 2 then begin
        print,'-Syntax: st = c->get_lumfunc(cstruct, radbin)'
        on_error, 2
        message,'Halting'
    endif
    lmin = 10d^dd.loglbins_min
    lmax = 10d^dd.loglbins_max
    lum = (lmin + lmax)/2d
    loglum = alog10(lum)

    numdens = reform(dd.numdens[radbin,*,*])
    err = reform(dd.numdens_err[radbin,*,*])
    lumfunc = total(numdens, 2)
    lumfunc_err = sqrt( total(err^2,2) )

    st = {                           $
            lum:         lum,        $
            logl:        loglum,     $
            lumfunc:     lumfunc,    $
            lumfunc_err: lumfunc_err $
        }
    return, st
end
pro correlate::calc_meanlum_hist, dd, meanlum, meanlum_err

    lmin = 10d^dd.loglbins_min
    lmax = 10d^dd.loglbins_max
    lmean = (lmin + lmax)/2d

    nrad = n_elements(dd.r)
    meanlum = dblarr(nrad)
    meanlum_err = dblarr(nrad)
    for i=0L, nrad-1 do begin
        lst = self->get_lumfunc(dd, i)
        ; Get the mean luminosity from the luminosity function
        ; If more than 10% of this histogram is negative, we will just
        ; set to -1
        lf = lst.lumfunc
        lferr = lst.lumfunc_err

        ;; only count bins not consistent with zero
        wbad=where(lf+lferr lt 0, nbad, comp=wgood, ncomp=ngood)
        if nbad ne 0 then lbad = total(lf[wbad]) else lbad=0.0
        if ngood ne 0 then lgood = total(lf[wgood]) else lgood=0.0
        
;        if abs(lbad)/lgood gt 0.02 then begin
        if float(nbad)/float(ngood) gt 0.05 then begin
            meanlum[i] = -1
        endif else begin
            meanlum[i] = total( lf*lst.lum )/total(lf)
            pplot, [meanlum[i],meanlum[i]], [0.0, 1.e10],/overplot
        endelse
        pplot, lmean, lf, yerr=lferr, psym=10, /xlog
        pplot, [1,1.e15], [0,0], /overplot

        key=prompt_kbrd('hit a key')
    endfor
end
pro correlate::calc_meanlum, dd, meanlum, meanlum_err, ratio=ratio
    if keyword_set(ratio) then begin
        self->calc_meanlum_ratio, dd, meanlum, meanlum_err
    endif else begin
        self->calc_meanlum_hist, dd, meanlum, meanlum_err
    endelse
end
pro correlate::plot_meanlum, dtype, subtype, error=error, ratio=ratio, $
        dops=dops

    dd=self->corr_read('dd', dtype, sub=subtype)

    n=n_elements(dd)

    colors=make_rainbow(n)
    yrange=[0,2]*1.e10

    str=mkstr(150,val='-')
    print, dd[0].r, f='(22F6.2)'
    print,str
    for i=0, n-1 do begin

        self->calc_meanlum, dd[i], meanlum, mlerr, ratio=ratio
        if keyword_set(error) then meanlum_err=mlerr
        if i eq 0 then overplot=0 else overplot=1
        pplot, dd[i].r, meanlum, yerr=meanlum_err, yrange=yrange, $
            /xlog, ystyle=3, overplot=overplot
        pplot, dd[i].r, meanlum, yerr=meanlum_err, color=colors[i], /overplot

        print, meanlum/1.e10, f='(22F6.2)'
        if keyword_set(error) then begin
            print, meanlum_err/1.e10, f='(22F6.2)'
            print, meanlum/meanlum_err, f='(22F6.2)'
        endif
        print,str
    endfor

end


pro correlate::plot_bluefrac, dtype, subtype, dops=dops, fill=fill

    np = n_elements(dtype) + n_elements(subtype)
    if np lt 2 then begin
        on_error, 2
        print,'-Syntax: plot_bluefrac, dtype, subtype, /dops, /fill'
        print
        message,'Halting'
    endif

    if keyword_set(dops) then begin
        file = self->plotfile(dtype, sub=subtype, bin=bin, $
            /createdir,/encapsulated)
        file = repstr(file,'.eps','_bluefrac.eps')
        if keyword_set(fill) then file=repstr(file,'bluefrac','bluefrac_fill')
        begplot,file,/color,/encapsulated,xsize=8,ysize=7
    endif

    dd=self->corr_read('dd',dtype,subtype=subtype)
    ws=self->where_string(subtype, labels=labels)
    nbin=n_elements(dd)

    redcut = 1.2
    maxerr = 0.1

    xtitle = textoidl('R [h^{-1} Mpc]')
    ytitle = 'Blue Fraction'

    clist = make_rainbow(nbin)

    ind=findgen(nbin)
    ind = reverse( findgen(nbin) )
    yrange = [0.0, 0.8]
    xrange = [min(dd[0].r), max(dd[0].r)]
    for ii=0L, nbin-1 do begin
        i = ind[ii]

        cf = self->color_fractions(dd[i], redcut)

        w=where(cf.bluefracerr lt maxerr)
        if ii eq 0 then begin

            pplot, dd[0].r[w], cf[w].bluefrac, aspect=1, $
                yrange=yrange, ystyle=1, $
                /xlog, xrange=xrange, xstyle=3, $
                xtitle=xtitle, ytitle=ytitle, $
                overplot=overplot
        endif 

        if not keyword_set(fill) then begin
            pplot, dd[0].r[w], cf[w].bluefrac, yerr=cf[w].bluefracerr, color=clist[i], /overplot  
            ;pplot, dd[0].r[w], cf[w].bluefrac, color=clist[i], /overplot  
        endif else begin
            errfill, dd[0].r[w], cf[w].bluefrac, cf[w].bluefracerr, $
                color=clist[i]
        endelse

    endfor

    legend,labels,color=clist,line=0, $
        /right,/bottom,box=0,charsize=1
    legend, 'g-r < '+ntostr(redcut, format='(F0.1)'),$
        /left,box=0

    if keyword_set(dops) then endplot,/trim_bbox
end


; Smooth the color-mag image and find the peak
function correlate::find_colorlum_peak, dd

    nc = n_elements(dd.kgmrbins_min)
    nl = n_elements(dd.loglbins_min)
    cvals = (dd.kgmrbins_min + dd.kgmrbins_max)/2
    lvals = (dd.loglbins_min + dd.loglbins_max)/2
        
    lrange = [min(lvals),max(lvals)]
    crange = [min(cvals),max(cvals)]
    st = { $
            radius: 0.0, clrmax: 0.0, lmax: 0.0 $
         }
    st=replicate(st, dd.nrad)

    for i=0,dd.nrad-1 do begin 

        numdens = reform(dd.numdens[i,*,*])
        err = reform(dd.numdens_err[i,*,*])
 

        nsmooth = smooth(numdens,2)

        maxn = max(nsmooth,wmax)
        
        wc = wmax/nc
        wl = wmax mod nc 

        implot, nsmooth > 0, type='asinh', xrange=lrange, yrange=crange

        print, lvals[wl], cvals[wc]
        pplot, [lvals[wl]], [cvals[wc]], psym=7, /overplot, color=c2i('red')
        key = prompt_kbrd('hit a key')
        if key eq 'q' then return, -1

    endfor

    return, st

end

; Find and plot the location of the red sequence as a function of
; radius
pro correlate::plot_redseq, dtype, subtype, dops=dops, fill=fill

    np = n_elements(dtype) + n_elements(subtype)
    if np lt 2 then begin
        on_error, 2
        print,'-Syntax: plot_redseq, dtype, subtype, /dops, /fill'
        print
        message,'Halting'
    endif

    if keyword_set(dops) then begin
        file = self->plotfile(dtype, sub=subtype, bin=bin, $
            /createdir,/encapsulated)
        file = repstr(file,'.eps','_redseq.eps')
        if keyword_set(fill) then file=repstr(file,'redseq','redseq_fill')
        begplot,file,/color,/encapsulated,xsize=8,ysize=7
    endif

    dd=self->corr_read('dd',dtype,subtype=subtype)
    ws=self->where_string(subtype, labels=labels)
    nbin=n_elements(dd)

    redcut = 1.2
    maxerr = 0.1

    xtitle = textoidl('R [h^{-1} Mpc]')
    ytitle = 'Red Peak'

    clist = make_rainbow(nbin)

    ind=findgen(nbin)
    ind = reverse( findgen(nbin) )
    yrange = [0.0, 0.8]
    xrange = [min(dd[0].r), max(dd[0].r)]
    for ii=0L, nbin-1 do begin
        i = ind[ii]

        rs = self->find_colorlum_peak(dd[i])
        for i=0L, nrad-1 do begin

        endfor

        stop
        w=where(cf.bluefracerr lt maxerr)
        if ii eq 0 then begin

            pplot, dd[0].r[w], cf[w].bluefrac, aspect=1, $
                yrange=yrange, ystyle=1, $
                /xlog, xrange=xrange, xstyle=3, $
                xtitle=xtitle, ytitle=ytitle, $
                overplot=overplot
        endif 

        if not keyword_set(fill) then begin
            pplot, dd[0].r[w], cf[w].bluefrac, yerr=cf[w].bluefracerr, color=clist[i], /overplot  
            ;pplot, dd[0].r[w], cf[w].bluefrac, color=clist[i], /overplot  
        endif else begin
            errfill, dd[0].r[w], cf[w].bluefrac, cf[w].bluefracerr, $
                color=clist[i]
        endelse

    endfor

    legend,labels,color=clist,line=0, $
        /right,/bottom,box=0,charsize=1
    legend, 'g-r < '+ntostr(redcut, format='(F0.1)'),$
        /left,box=0

    if keyword_set(dops) then endplot,/trim_bbox
end




PRO correlate::color_gaussfit, dtype, subtype, bin, dops=dops

    np = n_elements(dtype) + n_elements(subtype) + n_elements(bin)
    if np lt 3 then begin
        on_error, 2
        print,'-Syntax: color_gaussfit, dtype, subtype, bin, /dops'
        print
        message,'Halting'
    endif

    if keyword_set(dops) then begin
        file = self->plotfile(dtype, sub=subtype, bin=bin, $
            /createdir,/encapsulated)
        file = repstr(file,'.eps','_color_gaussfit.eps')
        begplot,file,/encapsulated,xsize=8,ysize=7
    endif

    ddcorr=self->corr_read('dd',dtype,subtype=subtype, bin=bin)

    xrange = [min(ddcorr.loglbins_min), max(ddcorr.loglbins_max)]
    yrange = [min(ddcorr.kgmrbins_min), max(ddcorr.kgmrbins_max)]

    !p.charsize=1
  
    gmr = (ddcorr.kgmrbins_min+ddcorr.kgmrbins_max)/2.0

    nx = ddcorr.nlum
    ny = ddcorr.ngmr
    xmin = (ddcorr.kgmrbins_min[0]+ddcorr.kgmrbins_max[0])/2
    xmax = (ddcorr.kgmrbins_min[ny-1]+ddcorr.kgmrbins_max[ny-1])/2


    cvals = (ddcorr.kgmrbins_min + ddcorr.kgmrbins_max)/2

    redcut = 1.2
    wred = where(cvals ge 1.2)

    ncols = 5
    nrows = 4
    erase & multiplot, [5, 4], /square

    simpctable, colorlist=clr
    nclr = n_elements(clr)
    clr = clr[1:nclr-1]
    for i=0,ddcorr.nrad-1 do begin 


        numdens = reform(ddcorr.numdens[i,*,*])
        err = reform(ddcorr.numdens_err[i,*,*])
 
        ; Accumulate across the luminosity axis
        cdist = total(numdens, 1)
        ;err = 1.0/total(1.0/err^2,1)
        err = total(err^2,1)
        err = sqrt(err)
        norm = qgauss(cdist, cvals, 20)
        cdist = cdist/norm
        err = err/norm
        pplot, cvals, cdist, yerr=err, /nohat, $
            yrange=[0, 3.8], ystyle=1, psym=10  

        w=where(cdist gt 0,nw)
        if nw gt 9 then begin
            aguess = [0.5, 1.32, 0.15, $
                        0.5, 0.9, 0.2,$
                        0.1, 0.6, 1.0 ]

            pvals=fitngauss(cvals[w], cdist[w], $
                err=err[w], yfit=yfit, aguess=aguess)

            xx = arrscl(findgen(1000), min(cvals), max(cvals))
            yfit = ngauss(xx, pvals)
            pplot, xx, yfit, /overplot, color=c2i('darkgreen'), thick=2

            ng = n_elements(aguess)/3
            for j=0L, ng-1 do begin
                jj = j*3
                gauss = pvals[jj]*exp( -(xx-pvals[jj+1])^2/2.0/pvals[jj+2]^2 )
                pplot, xx, gauss, /overplot, color=clr[j]
            endfor
        endif

        ;rad = (ddcorr.radbins_min[i]+ddcorr.radbins_max[i])/2
        rad = rnd(ddcorr.r[i], 3)
        format = self->radformat(rad)
        rad = ntostr(rad,format=format)
        if rad eq '0.100' then rad = '0.10'
        legend, $
            'r = '+rad+' Mpc', charsize=0.7,/left,box=0,margin=0

        if i lt (ddcorr.nrad-1) then begin
            multiplot
        endif

    endfor 
;    multiplot,/reset

    add_title, 'xtitle', 'Rest Frame g-r', $
        offset=0.7, charsize=1.5
    add_title, 'ytitle', 'Number Density [# h!U2!N Mpc!U-2!N]', $
        offset=-1.5, charsize=1.5

    multiplot,/reset

    if keyword_set(dops) then endplot,/trim_bbox
end 


pro correlate::plot_lumfunc, dtype, subtype, bin, band_shift, absmag=absmag, linear=linear, lum=lum, color=color, dops=dops, dopng=dopng

    np = n_elements(dtype) + n_elements(subtype) + n_elements(bin) + $
        n_elements(band_shift)
    if np lt 4 then begin
        on_error, 2
        print,'-Syntax: plot_lumfunc, dtype, subtype, bin, band_shift, /linear, /dops, /dopng'
        print
        message,'Halting'
    endif

    bshift = ntostr(band_shift,format='(f4.2)')

    psfile = self->plotfile(dtype, sub=subtype, bin=bin, $
        /createdir,/encapsulated)
    if keyword_set(linear) then begin
        psfile = repstr(psfile,'.eps','_lumfunc_linear.eps')
    endif else begin
        psfile = repstr(psfile,'.eps','_lumfunc.eps')
    endelse
    if keyword_set(lum) then psfile=repstr(psfile,'.eps','_lumweight.eps')

    if keyword_set(color) then psfile=repstr(psfile, '.eps', '_color.eps')
 
    if keyword_set(dops) then begin
       begplot,psfile,color=color,/encapsulated,xsize=8,ysize=7
    endif
    if keyword_set(dopng) then begin
        file = self->plotfile(dtype, sub=subtype, bin=bin, $
            /createdir,/encapsulated)
        file = repstr(file,'.eps','_lumfunc.png')
    endif

    if keyword_set(color) then oclr = !red else oclr=!grey50

    !p.charsize=1
    if !d.name eq 'PS' then begin
        fitcolor=!blue
        !p.charsize=0.9
        lchar=0.75
        lmarg=-0.5

        !p.thick=2
        !x.thick=2
        !y.thick=2
        !p.charthick=2
    endif else begin
        fitcolor=c2i('green')
        ;!p.charsize=1
        ;!p.charsize=2
        lchar=1
        lmarg=0
    endelse
    ddcorr=self->corr_read('dd',dtype,subtype=subtype, bin=bin, /silent)

    nx = ddcorr.nlum
    lvals = (ddcorr.loglbins_min + ddcorr.loglbins_max)/2

    ; in absmag
    if keyword_set(absmag) then begin
        solarmags = k_solar_magnitudes(band_shift=band_shift, /silent)
        absmag = -2.5*lvals + solarmags[3]
        ;xrange = reverse([min(absmag), max(absmag)])
        xrange = [-18.9, -24.8]
        xtitle = 'M!S!Di!R!U0.25!N - 5log(h)'

        xvals = absmag
        xstyle=1
    endif else begin
        xrange = [min(ddcorr.loglbins_min), max(ddcorr.loglbins_max)]
        ;xtitle = 'log( L/L'+sunsymbol()+' h!U2!N)'
        ;xtitle = 'log( L!D'+bshift+'i!N/L'+sunsymbol()+' h!U2!N)'

        xtitle = textoidl('log(L_{'+bshift+'i}/L_{\odot}h^2)')

        xrange = [min(lvals),max(lvals)]
        xvals = lvals
        xstyle=1
    endelse
    
    ncols = 5
    nrows = 4
    if not keyword_set(linear) then begin
        yrange = [1.e-3, 3.8]
        ylog=1
        ystyle=1
        ytickformat='loglabels'
    endif else begin
        ylog=0 
        ystyle=3
        yrange=[-0.1, 2]
    endelse
    if keyword_set(lum) then begin
        ytitle = textoidl('Normalized L \times \phi')
    endif else begin
        ytitle = textoidl('Normalized \phi')
    endelse

    titsize=2
    erase & multiplot, [ncols, nrows], /square, $
        ytickformat=ytickformat, $
        mxTitle=xtitle, mxTitSize=titsize, mxTitOffset=0.8, $
        myTitle=ytitle, myTitSize=titsize


    ws=self->where_string(subtype, labels=labels, /nodisplay)
    xs = !csym.chi+'!U2!N/'+!csym.nu
    simpctable, colorlist=clr
    nclr = n_elements(clr)
    clr = clr[1:nclr-1]
    for i=0,ddcorr.nrad-1 do begin 

        if keyword_set(lum) then begin
            dens = reform(ddcorr.ilumdens[i,*,*])
            err = reform(ddcorr.ilumdens_err[i,*,*])
        endif else begin
            dens = reform(ddcorr.numdens[i,*,*])
            err = reform(ddcorr.numdens_err[i,*,*])
        endelse
 
        ; Accumulate across the luminosity axis
        ldist = total(dens, 2)
        err = sqrt( total(err^2,2) )
        norm = qgauss(ldist, lvals, 20)
        ldist = ldist/norm
        err = err/norm
 
        pplot, xvals, ldist, yerr=err, $
            xrange=xrange, xstyle=xstyle, yticklen=0.04, $
            ylog=ylog, yrange=yrange, ystyle=ystyle, psym=10, $
            hat=0

        if keyword_set(linear) then begin
            pplot, xrange, [0,0],/overplot, color=oclr
        endif

        rad = rnd(ddcorr.r[i], 3)
        format = self->radformat(rad)
        rad = ntostr(rad,format=format)
        if rad eq '0.100' then rad = '0.10'
        legend, $
            'r = '+rad+' Mpc', charsize=lchar,/right,/top,box=0,margin=0

        multiplot
    endfor 
    multiplot, /default

    if keyword_set(dops) then endplot,/trim_bbox
    if keyword_set(dopng) then write_png, file, tvrd(/true)

    dirsep, psfile, tdir, tfile 
    print,'\begin{figure}[p]'
    print,'\centering'
    print,'\includegraphics[scale=0.9]{plots/'+tfile+'}'
    print,'\caption{$'+labels[bin]+'$}'
    print,'\end{figure}'


    !p.multi=0

end 



;
; Schechter fitting
;

function correlate::schechter_file, dtype, subtype, ftype, bin=bin, all=all, cumulative=cumulative, reduced=reduced, addstring=addstring

    post='_schechter'

    if keyword_set(cumulative) then post=post+'_cum' 
    if keyword_set(reduced) then post=post+'_reduced'
    
    if n_elements(bin) eq 0 then post=post+'_all'

    if n_elements(addstring) ne 0 then post=post+addstring

    file = self->plotfile(dtype, sub=subtype, bin=bin, $
        /createdir,/encapsulated)
    file = repstr(file,'.eps',post)
    case ftype of
        'ps': file=file+'.ps'
        'eps': file=file+'.eps'
        'png': file=file+'.png'
        'fits': file=file+'.fits'
    endcase

    return, file
end


; run the schechter fit on all bins
function correlate::run_schechter_fit, dtype, subtype, band_shift, cumulative=cumulative, labels=labels, dops=dops, allps=allps

    if keyword_set(allps) and keyword_set(dops) then message,'do not send /allps and /dops'
    if keyword_set(allps) then begin
        file=self->schechter_file(dtype,subtype,'ps',cumulative=cumulative)
        begplot,file,/color
    endif


    ddcorr=self->corr_read('dd',dtype,subtype=subtype, bin=bin, /silent)
    ws = self->where_string(subtype, labels=labels)
    r = ddcorr[0].r
    nbin=n_elements(ddcorr)
    ddcorr=0

    
    for bin=0L, nbin-1 do begin
        fitst = self->schechter_fit(dtype, subtype, bin, band_shift, $
                    cumulative=cumulative, dops=dops)

        if bin eq 0 then begin
            allfitst = replicate(fitst, nbin)
        endif
        allfitst[bin] = fitst
    endfor

    if keyword_set(allps) then endplot

    fitfile=self->schechter_file(dtype, subtype, 'fits',cumulative=cumulative)
    print,'Writing to file: ',fitfile
    mwrfits, allfitst, fitfile, /create
    return, allfitst

end

; this is much simpler than plot_radbin_schechter_fits
pro correlate::plot_radbin_schechter_fits, dtype, subtype, radbin, meanz, $
        cumulative=cumulative, dops=dops

    nd=n_elements(dtype) & ns=n_elements(subtype) & nr=n_elements(radbin)
    nz=n_elements(meanz)
    if (nd+ns+nr+nz) lt 4 then begin
        print,'-Syntax: c->plot_radbin_schechter_fits, dtype, subtype, radbin, meanz, /cumulative, /dops'
        on_error, 2
        message, 'Halting'
    endif

    lumfile=self->schechter_file(dtype, subtype, 'eps', /reduced, $
        cumulative=cumulative, addstring='_lumfunc')
    bestfitfile=self->schechter_file(dtype, subtype, 'eps', /reduced, $
        cumulative=cumulative, addstring='_bestfit')

    if keyword_set(dops) then begin
        begplot,lumfile,/color,/encap
    endif

    !p.charsize=1.7
    ddcorr=self->corr_read('dd',dtype,subtype=subtype, /silent)

    ws=self->where_string(subtype, labels=labels)
    fitfile=self->schechter_file(dtype,subtype,'fits',cumulative=cumulative)
    print,'Reading fitfile: ',fitfile
    t=mrdfits(fitfile,1)
    n=n_elements(t)
    nrad=n_elements(t.mstar)

    colors = make_rainbow(n)
    rcolors = reverse(colors)

    ; plot the luminosity functions for this radial bin
    lvals = (ddcorr[0].loglbins_min + ddcorr[0].loglbins_max)/2
    lvals = 10^lvals
    yrange=[7.e-4, 0.7]

    ;xtitle = 'log( L/L'+sunsymbol()+' h!U2!N)'
    xtitle = 'L [ h!U-2!N L'+sunsymbol()+' ]'
    ytitle = '# [ h!U2!N Mpc!U-2!N ]'

    for i=n-1,0,-1 do begin

        if keyword_set(cumulative) then begin
            cumstruct = self->cumulative_densities(ddcorr[i])
            numdens = reform(cumstruct.numdens[radbin,*,*])
            err = reform(cumstruct.numdens_err[radbin,*,*])
        endif else begin
            numdens = reform(ddcorr[i].numdens[radbin,*,*])
            err = reform(ddcorr[i].numdens_err[radbin,*,*])
        endelse
        ; Accumulate across the luminosity axis
        ldist = total(numdens, 2)
        err = sqrt( total(err^2,2) )

        psym=0
        if i eq (n-1) then begin
            pplot, lvals, ldist, psym=psym, $
                /xlog, xstyle=3, xtickf='loglabels', xrange=[3.e9, 7.e11], $
                xtitle = xtitle, xticklen=0.04, $
                ytitle = ytitle, $
                yrange=yrange, ystyle=3, /ylog, $
                aspect=1, ytickf='loglabels', yticklen=0.04
        endif 
        pplot, lvals, ldist, psym=psym, /overplot, $
            color=rcolors[i]
    endfor

    if n_elements(labels) ne 0 then begin
        legend, reverse(labels), /right, charsize=1, color=rcolors, line=0
    endif


    key=prompt_kbrd('hit a key')
    if key eq 'q' then return

    if keyword_set(dops) then endplot, /trim_bbox

    if keyword_set(dops) then begin
        begplot,bestfitfile,/color, /encap
    endif

    colors = make_rainbow(n)
    mstar = t.mstar[radbin]
    mstar_err=t.mstar_err[radbin]
    alpha = t.alpha[radbin]
    alpha_err=t.alpha_err[radbin]

    lstar = sdss_am2lumsolar(mstar, clr=3, amerr=mstar_err, lumerr=lstar_err, $
                             band_shift=meanz)
    colprint, mstar, mstar_err, lstar, lstar_err

;    yrange = [-1.0,-0.5]
;    xrange = [-20.2, -21.5]
;    pplot, [-20], [-0.7], $
;        aspect=1, $
;        xtitle='M!D*!N - 5 log(h)', $
;        ytitle=!csym.alpha, $
;        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3

;    for i=0L, n-1 do begin
;        pplot, [mstar[i]], xerr=mstar_err[i], [alpha[i]], yerr=alpha_err[i], $
;            psym=8, $
;            xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
;            color=colors[i], /overplot
;    endfor
    
;    if n_elements(labels) ne 0 then begin
;        legend, reverse(labels), /right, charsize=1, color=reverse(colors), $
;            psym=8
;    endif

;    if keyword_set(dops) then endplot, /trim_bbox


    !p.charsize=2
    yrange = [-1.0,-0.5]
    xrange = [9.e9, 1.8e10]/1.e10
    xtitle='L!D*!N [ 10!U10!N h!U-2!N L'+sunsymbol()+' ]'
    ytitle=!csym.alpha
    pplot, [-20], [-0.7], $
        aspect=1, $
        xtitle=xtitle, ytitle=ytitle, $
        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3

    for i=0L, n-1 do begin

        pplot, $
            [lstar[i]]/1.e10, xerr=lstar_err[i]/1.e10, $
            [alpha[i]], yerr=alpha_err[i], $
            psym=8, $
            xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
            color=colors[i], /overplot
    endfor
    
    if n_elements(labels) ne 0 then begin
        legend, reverse(labels), /right, charsize=1, color=reverse(colors), $
            psym=8
    endif


    if keyword_set(dops) then endplot, /trim_bbox

    colprint, $
        ntostr(t.meanlum[radbin]*1.e10)+' '+!plusminus+' '+ ntostr(t.meanlum_err[radbin]*1.e10)+'  '+ ntostr(lstar)+' '+!plusminus+' '+ ntostr(lstar_err)

    pplot, lindgen(n), lstar/1.e10, yerr=lstar_err/1.e10, $
        yrange=[0,2]
    pplot, lindgen(n), t.meanlum[radbin], yerr=t.meanlum_err[radbin], $
        /overplot, color=!blue

end 

pro correlate::plot_all_schechter_fits, dtype, subtype, $
        cumulative=cumulative, dops=dops

    if keyword_set(dops) then begin
        file=self->schechter_file(dtype,subtype,'ps', /reduced, $
            cumulative=cumulative)
        begplot,file,/color
    endif

    ws=self->where_string(subtype, labels=labels)
    fitfile=self->schechter_file(dtype,subtype,'fits',cumulative=cumulative)
    allfitst=mrdfits(fitfile,1)

    nbin = n_elements(allfitst)
    clist = make_rainbow(nbin)

    mstarytitle = 'M!S!D*!N!R!U0.25!N - 5 log(h)'
    xtitle = textoidl('R [h^{-1} Mpc]')
    r = allfitst[0].r
    pplot, [0], xrange=[min(r),max(r)], xstyle=3, /xlog, $
        yrange=[-18,-22], ystyle=3, $
        xtitle=xtitle, ytitle=mstarytitle, $
        aspect=1
    for bin=0L, nbin-1 do begin
        w=where(allfitst[bin].chi2per lt 3)
        pplot, allfitst[bin].r[w], allfitst[bin].mstar[w], $
            yerr=allfitst[bin].mstar_err[w], $
            color=clist[bin], /overplot
    endfor

    if n_elements(labels) ne 0 then begin
        legend, labels, /right, /bottom, box=0, charsize=1, color=clist, line=0
    endif

    key = prompt_kbrd('hit a key')
    
    alphaytitle = !csym.alpha
    r = allfitst[0].r
    pplot, [0], xrange=[min(r),max(r)], xstyle=3, /xlog, $
        yrange=[-1.5,1.5], ystyle=3, $
        xtitle=xtitle, ytitle=alphaytitle, $
        aspect=1
    for bin=0L, nbin-1 do begin
        w=where(allfitst[bin].chi2per lt 3)
        pplot, allfitst[bin].r[w], allfitst[bin].alpha[w], $
            yerr=allfitst[bin].alpha_err[w], $
            color=clist[bin], /overplot
    endfor

    if n_elements(labels) ne 0 then begin
        legend, labels, /right, box=0, charsize=1, color=clist, line=0
    endif


    key = prompt_kbrd('hit a key')
    ; The joint fit as a function of radius for each.  This gets kind
    ; of confusing
    xrange = [-19, -22]
    yrange = [-1.5,1.5]
    pplot, [0], xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, /nodata, $
        xtitle=mstarytitle, ytitle=alphaytitle, aspect=1
    for bin=0L, nbin-1 do begin

        w=where(allfitst[bin].chi2per lt 5)
        pplot, allfitst[bin].mstar[w], allfitst[bin].alpha[w], $
            color=clist[bin], /overplot

    endfor
    if n_elements(labels) ne 0 then begin
        legend, labels, /right, box=0, charsize=1, color=clist, line=0
    endif

    if keyword_set(dops) then endplot

end

function correlate::cumulative_densities, st
    numdens = st.numdens
    numdens_err = st.numdens_err
    ilumdens = st.ilumdens
    ilumdens_err = st.ilumdens_err

    ; for each color/lum bin accumulate counts over radius and devide by
    ; cumulative area
    cumarea = !dpi*(st.radbins_max^2-st.radbins_min[0]^2)
    for i=0L, st.nlum-1 do begin
        for j=0L, st.ngmr-1 do begin
            numdens[*,i,j] = $
                total(st.numdens[*,i,j]*st.area,/cumulative)/cumarea
            err = st.numdens_err[*,i,j]*st.area
            numdens_err[*,i,j] = sqrt(total(err^2,/cumulative))/cumarea

            ilumdens[*,i,j] = $
                total(st.ilumdens[*,i,j]*st.area,/cumulative)/cumarea
            err = st.ilumdens_err[*,i,j]*st.area
            ilumdens_err[*,i,j] = sqrt(total(err^2,/cumulative))/cumarea
        endfor
    endfor

    radnumdens = total(st.radnumdens*st.area, /cumulative)/cumarea 
    err = st.radnumdens_err*st.area
    radnumdens_err = sqrt(total(err^2,/cumulative))/cumarea

    radilumdens = total(st.radilumdens*st.area, /cumulative)/cumarea 
    err = st.radilumdens_err*st.area
    radilumdens_err = sqrt(total(err^2,/cumulative))/cumarea

    out={numdens:numdens, numdens_err:numdens_err, $
         ilumdens:ilumdens, ilumdens_err:ilumdens_err, $
         radnumdens:radnumdens, radnumdens_err:radnumdens_err, $
         radilumdens:radilumdens, radilumdens_err:radilumdens_err}
    return, out

end

function correlate::schechter_fit, dtype, subtype, bin, band_shift, cumulative=cumulative, dops=dops, dopng=dopng

    np = n_elements(dtype) + n_elements(subtype) + n_elements(bin) + $
        n_elements(band_shift)
    if np lt 4 then begin
        on_error, 2
        print,'-Syntax: schechter_fit, dtype, subtype, bin, band_shift, /dops, /dopng'
        print
        message,'Halting'
    endif

    if keyword_set(cumulative) then begin
        post='_schechter_cum' 
    endif else begin
        post='_schechter'
    endelse
    if keyword_set(dops) then begin
        file=self->schechter_file(dtype, subtype, 'eps', bin=bin, $
            cumulative=cumulative)
        begplot,file,/color,/encapsulated,xsize=8,ysize=7
    endif
    if keyword_set(dopng) then begin
        file=self->schechter_file(dtype, subtype, 'png', bin=bin, $
            cumulative=cumulative)
    endif

    if !d.name eq 'PS' then begin
        fitcolor=!blue
        !p.charsize=0.6
        lchar=0.6
        lmarg=-0.5
    endif else begin
        fitcolor=c2i('green')
        !p.charsize=1
        lchar=1
        lmarg=0
    endelse
    ddcorr=self->corr_read('dd',dtype,subtype=subtype, bin=bin, /silent)

    nx = ddcorr.nlum
    lvals = (ddcorr.loglbins_min + ddcorr.loglbins_max)/2


    arrval = dblarr(ddcorr.nrad)
    arrval2 = dblarr(ddcorr.nrad,3,3)
    fitstruct = $
        {r: ddcorr[0].r, $
         phistar: arrval, phistar_err: arrval, $
         mstar: arrval, mstar_err: arrval, $
         alpha: arrval, alpha_err: arrval, $
         covar: arrval2, $
         chi2: arrval, dof: arrval, chi2per: arrval, $
         meanlum: arrval, meanlum_err:arrval $
        }

    ; in absmag
    solarmags = k_solar_magnitudes(band_shift=band_shift, /silent)
    absmag = -2.5*lvals + solarmags[3]

    ;xrange = [min(ddcorr.loglbins_min), max(ddcorr.loglbins_max)]
    xrange = reverse([min(absmag), max(absmag)])
    
    xrange = [-18.9, -24.8]
    ncols = 5
    nrows = 4
    erase & multiplot, [5, 4], /square, $
        mxTitle='M!S!Di!N!R!U0.25!N - 5log h', mxTitSize=1.5, $
        myTitle='Normalized Number Density', myTitSize=1.5, $
        ytickf='loglabels'

    xs = !csym.chi+'!U2!N/'+!csym.nu
    simpctable, colorlist=clr
    nclr = n_elements(clr)
    clr = clr[1:nclr-1]
    for i=0,ddcorr.nrad-1 do begin 

        if keyword_set(cumulative) then begin
            cumstruct = self->cumulative_densities(ddcorr)
            numdens = reform(cumstruct.numdens[i,*,*])
            err = reform(cumstruct.numdens_err[i,*,*])

            mlum=cumstruct.radilumdens[i]/cumstruct.radnumdens[i]
            mlum_err=mlum*sqrt( $
        (cumstruct.radilumdens_err[i]/cumstruct.radilumdens[i])^2 + $
        (cumstruct.radnumdens_err[i]/cumstruct.radnumdens[i])^2 )

        endif else begin
            numdens = reform(ddcorr.numdens[i,*,*])
            err = reform(ddcorr.numdens_err[i,*,*])

            mlum=ddcorr.radilumdens[i]/ddcorr.radnumdens[i]
            mlum_err=mlum*sqrt( $
        (ddcorr.radilumdens_err[i]/ddcorr.radilumdens[i])^2 + $
        (ddcorr.radnumdens_err[i]/ddcorr.radnumdens[i])^2 )
        endelse
 
        ; Accumulate across the luminosity axis
        ldist = total(numdens, 2)
        err = sqrt( total(err^2,2) )
        norm = qgauss(ldist, lvals, 20)
        ldist = ldist/norm
        err = err/norm
        ;pplot, lvals, ldist, yerr=err, /nohat, $
        ;    /ylog, yrange=[1.e-2, 3.8], ystyle=1, psym=10  
        pplot, absmag, ldist, yerr=err, /nohat, $
            xrange=xrange, xstyle=1, yticklen=0.04, $
            /ylog, yrange=[1.e-3, 3.8], ystyle=1, psym=10  


        ;w=where(ldist gt 0,nw)
        w=where(ldist/err gt 1,nw)
        aguess = [1.0d, -20.0d, -1.0d]

        fit = fitschechter(absmag[w], ldist[w], err[w], params=pvals, $
                                aguess=aguess, sigma=sigma)

        xx = arrscl(findgen(1000), min(absmag), max(absmag))
        yfit = mpfit_schechter_absmag(xx, pvals)
        pplot, xx, yfit, /overplot, color=fitcolor, thick=2

        rad = rnd(ddcorr.r[i], 3)
        format = self->radformat(rad)
        rad = ntostr(rad,format=format)
        if rad eq '0.100' then rad = '0.10'
        legend, $
            'r = '+rad+' Mpc', charsize=lchar,/left,/bottom,box=0,margin=0

        mmess = 'M!D*!N '+ntostr(fit.mstar, format='(F0.2)')+!csym.plusminus+ntostr(fit.mstar_err,format='(F0.2)')
        amess = !csym.alpha+' '+ntostr(fit.alpha, format='(F0.2)')+!csym.plusminus+ntostr(fit.alpha_err,format='(F0.2)')

        chimess = xs+' '+ntostr(fit.chi2per,format='(F0.1)')
        legend,[mmess,amess,chimess],/right,box=0,charsize=lchar,margin=lmarg


        if i lt (ddcorr.nrad-1) then begin
            multiplot
        endif

        fitstruct.phistar[i] = fit.phistar
        fitstruct.phistar_err[i] = fit.phistar_err
        fitstruct.mstar[i] = fit.mstar
        fitstruct.mstar_err[i] = fit.mstar_err
        fitstruct.alpha[i] = fit.alpha
        fitstruct.alpha_err[i] = fit.alpha_err
        fitstruct.covar[i,*,*] = fit.covar

        fitstruct.chi2[i] = fit.chi2
        fitstruct.dof[i] = fit.dof
        fitstruct.chi2per[i] = fit.chi2per

        fitstruct.meanlum[i] = mlum
        fitstruct.meanlum_err[i] = mlum_err
    endfor 

    multiplot,/reset

    if keyword_set(dops) then endplot,/trim_bbox
    if keyword_set(dopng) then write_png, file, tvrd(/true)

    return, fitstruct
end 











pro correlate::plot_lum_vs_rad, ddcorr

 !p.multi = [0,5,4] 
 !p.charsize=2

;  erase & multiplot, [5,4]
  
  yrange=[1.e-4,100]
  
  xtitle = 'log(L/L'+sunsymbol()+')'
  
  loglvals = (ddcorr.loglbins_min+ddcorr.loglbins_max)/2.0
  xrange = [0.99*min(ddcorr.loglbins_min), max(ddcorr.loglbins_max)*1.01]
;  xrange = [min(loglvals),max(loglvals)]

  for i=0,ddcorr.nrad-1 do begin 
      
      numdens = reform( ddcorr.numdens[i,*,*] )
      numdens_err = reform( ddcorr.numdens_err[i,*,*])
      
      ldens = total(numdens, 2, /double)
      ldens_err= sqrt( total(numdens_err^2, 2, /double) )
      
;     yrange = [ 1.e-4, 1.5*max(ldens) ]
      
      mplot, loglvals, ldens, $
        /ylog, psym=10, mxtitle=xtitle, mytitle='h!U2!N Mpc!U-2!N', $
        yrange = yrange, ystyle=3, $
        xrange=xrange, xstyle=3, $
        ytickf='loglabels', yticklen=0.04
      
      oploterror, loglvals, ldens, ldens_err, hat=0, psym=10
      
      legend, 'r = '+ntostr(ddcorr.r[i],4,/round)+' h!U-1!NMpc', $
        /right, box=0, charsize=1


;      IF i NE ddcorr.nrad-1 THEN multiplot
      
  endfor 

;  multiplot, /reset

  !p.multi=0
  !p.charsize=1

end 



pro correlate::plot_kgmr_vs_rad, ddcorr

 !p.multi = [0,5,4] 
 !p.charsize=2


 xmold = !x.margin
 ymold = !y.margin

 !x.margin = [5,1.5]
 !y.margin = [3,1.5]


;  erase & multiplot, [5,4]
  
;  yrange=[1.e-4,100]
  
  xtitle = 'g-r'
  
  cvals = (ddcorr.kgmrbins_min+ddcorr.kgmrbins_max)/2.0
  xrange = [0.99*min(ddcorr.kgmrbins_min), max(ddcorr.kgmrbins_max)*1.01]
;  xrange = [min(loglvals),max(loglvals)]

  for i=0,ddcorr.nrad-1 do begin 
      
      numdens = reform( ddcorr.numdens[i,*,*] )
      numdens_err = reform( ddcorr.numdens_err[i,*,*])
      
      cdens = total(numdens, 1, /double)
      cdens_err= sqrt( total(numdens_err^2, 1, /double) )
      
;     yrange = [ 1.e-4, 1.5*max(ldens) ]
      
;      mplot, cvals, cdens, $
;        psym=10, mxtitle=xtitle, mytitle='Mpc!U-2!N', $
;        yrange = yrange, ystyle=3, $
;        xrange=xrange, xstyle=3      
;      oploterror, cvals, cdens, cdens_err, hat=0, psym=10

      pplot, cvals, cdens, yerr=cdens_err, hat=0, $
        psym=10, xtitle=xtitle, ytitle='h!U2!N Mpc!U-2!N', $
        yrange = yrange, ystyle=3, $
        xrange=xrange, xstyle=3, aspect=1, $
        title = 'r = '+ntostr(ddcorr.r[i],4,/round)+' h!U-1!NMpc'


;      IF i NE ddcorr.nrad-1 THEN multiplot
      
  endfor 

;  multiplot, /reset

  !p.multi=0
  !p.charsize=1

  !x.margin = xmold
  !y.margin = ymold

end 





; Obsolete
FUNCTION correlate::_invert_integrate, r, density, rmax

  ;; Cumulative integral: repeat for various rmax
  npoints = 100

  nrad = n_elements(rmax)
  
  integral = dblarr(nrad)

  rmin = min(r)
  FOR i=0L, nrad-1 DO BEGIN 

      w = where(r LE rmax[i], nw)

      IF nw EQ 0 THEN BEGIN 
          message,'rmax value outside range of input r'
      ENDIF 


      gauleg, rmin, rmax[i], npoints, RR, WW

      iDensity = interpol(density, r, RR)
      integrand = 4d*!dpi*iDensity*RR^2

      integral[i] = total( integrand*WW )

  ENDFOR 
  
  return, integral
  

END 

pro correlate::invert_gettags, tagname_input, struct, tag, covtag, corrected=corrected
  tagname = strlowcase(tagname_input)
  tagnames = strlowcase( tag_names(struct) )

  covtagname = keyword_set(corrected) ? tagname+'_err' : tagname+'_cov'

  tag = where(tagnames eq tagname, ntag)
  if ntag eq 0 then message,'Tag not found: '+tagname
  covtag = where(tagnames eq covtagname, ncov)
  if ncov eq 0 then message,'Covariance tag not found: '+covtagname
end 


pro correlate__define
  struct = {$
             correlate, $
             par_struct: ptr_new() $
           }
end 
