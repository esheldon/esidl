;+
; NAME:
;  MAXBCG__DEFINE
;
;
; PURPOSE:
;  Class file containing routines for generating the galaxy catalogs for 
;  ben's maxbcg algorithm. Also includes tools for reading the files.
;
; CATEGORY:
;  SDSS Cluster finding.
;
; CALLING SEQUENCE:
;  mb = obj_new('maxbcg', catalog)
;
; For a new catalog, add it to the ::catname and ::get methods.
; Then run:
;    ::convert
;
; MODIFICATION HISTORY:
;   Created: Early May, 2005 from separate programs.  
;   Author: Erin Sheldon, UofChicago
;
;-


; note this is different from the lensing sample
function maxbcg::init, catalog

	if n_elements(catalog) eq 0 then message,'You must initialize the catalog'

	if not in(self->catalogs(), catalog) then begin
		message,'Unknown catalog '+strn(catalog)
	endif

	self.catalog = catalog
	print,'catalog =                      ',catalog

	return,1

end 

function maxbcg::type
  return,'maxbcg'
end 
function maxbcg::catalog
  return,self.catalog
end 
function maxbcg::catalogs
  return,['public','dr3','dr4plus','dr406','m2n','n2m','hv1','hv1halo',$
	  'gmbcg10']
end 


function maxbcg::maxbcg_maskfile
	case self->catalog() of
		'dr406': f=sdssidl_config('pixel_mask_bound')
		'gmbcg10': f=sdssidl_config('pixel_mask_bound')
		else: message,'bad catalog: '+self->catalog()
	endcase
  return, f
end 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ben's maxbcg output catalog
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function maxbcg::catdir

    dir = expand_path(esheldon_config('maxbcg_dir'))
    if strmatch(self.catalog, 'hv1*') then dir = path_join(dir,'sim')
    return,dir

end 

function maxbcg::catname, $
               fits=fits, $     ; the original fits file
               idit=idit, $
               neighbor=neighbor, $
               gals=gals, $
               spec=spec, $
               centerclass=centerclass

    catalog = self->catalog()

    case catalog of 
        'public': name = 'maxbcg_public_catalog'
        'dr406':  name = 'catalog_full_bcgs'
        'm2n':    name = 'noref/maxbcg_match_to_noref_3_12_06'
        'n2m':    name = 'noref/noref_match_to_maxbcg_3_12_06'
        'hv1': begin
            if keyword_set(gals) then begin
                name = 'hv_v1.00ja_bright_galaxies'
            endif else begin
                name = 'hv_v1.00ja_bright_bcgs_cts_z'
            endelse
            ext = '.fit'
        end
        'hv1halo': begin
            name = 'hv_v1.00ja_halos'
            ext='.fit'
        end
		; gmbcg; ngals =10 dr7
		'gmbcg10': name='gmbcg10'
        ELSE: message,'Unknown input sample: '+catalog
    endcase 

    if catalog eq 'public' then begin
        ext='.fit'
	endif else if catalog eq 'gmbcg10' then begin
		if keyword_set(fits) then begin
			ext='.fits'
		endif else begin
			ext='.st'
		endelse
    endif else if keyword_set(fits) or keyword_set(centerclass) then begin 
        ext = '.fit'
    endif else if keyword_set(spec) then begin 
        ext = '_spec.st'
    endif else if n_elements(ext) eq 0 then begin 
        ext = '.st'
    endif
    name = name+ext

    if keyword_set(neighbor) then begin 
        name = repstr(name, 'bcgs', 'bcggals')
    endif else if keyword_set(centerclass) then begin
        name = repstr(name, 'bcgs', 'bcgs_centerclass')
    endif


    dir=self->catdir()
    name = concat_dir(dir,name)
    return,name

end 

function maxbcg::get, gals=gals, neighbor=neighbor, centerclass=centerclass, $
        spec=spec, fits=fits, $
        columns=columns, rows=rows, nrows=nrows, _extra=_extra

    file = self->catname(gals=gals, neighbor=neighbor, centerclass=centerclass,$
        spec=spec, fits=fits)
    print,'Reading file: '+file

    cat = self->catalog()

    if stregex(file, '.*\.fits?$',/bool) then begin
        bcg = mrdfits(file, 1, columns=columns, rows=rows, _extra=_extra)
    endif else if stregex(file, '.*\.st$',/bool) then begin
        bcg = read_idlstruct(file, columns=columns, rows=rows, _extra=_extra)
    endif else begin
        message,'Do not know how to read files of this type: '+ntostr(file)
    endelse

    nrows = n_elements(bcg)

    return,bcg

end 



function maxbcg::zcuts

  ;; use with strict gt and lt
	catname=self->catalog()
	case catname of
		'dr406': return,[0.04, 0.3]
		'gmbcg10': return, [0.1,0.35]
		else: message,'bad catalog name: '+catname
	endcase

end 

function maxbcg::zselect, z, nkeep

	zcuts = self->zcuts()
	print,'Cutting redshifts at [',zcuts[0],', ',zcuts[1],']',$
		f='(a,f0.2,a,f0.2,a)'
	keep  = where(z GT zcuts[0] AND z LT zcuts[1], nkeep)

	return, keep

end 

function maxbcg::add_csurvey, orig

	n = n_elements(orig)
	newstruct = create_struct(orig[0], 'clambda', 0d, 'ceta', 0d)
	newstruct = replicate(newstruct, n)

	copy_struct, orig, newstruct

	eq2csurvey, newstruct.ra, newstruct.dec, clambda, ceta
	newstruct.clambda = clambda
	newstruct.ceta = ceta
	return, newstruct

end 

function maxbcg::add_tags, struct

  ;; Wish IDL had dictionaries.
  newtags = $
    ['bcg_vagc_index', $
     'bcg_iabsmag','bcg_iabsmag_ivar', $
     'bcg_ilum',   'bcg_ilum_ivar', $
     'iabsmag200', 'iabsmag200_ivar', $
     'ilum200',    'ilum200_ivar', $
     'bcg_kgmr', 'bcg_kgmr_ivar', $ ;; not assuming lrg
     'bcg_krmi', 'bcg_krmi_ivar', $
     $
     'bcg_spec_iabsmag', 'bcg_spec_iabsmag_ivar', $
     'bcg_spec_ilum',    'bcg_spec_ilum_ivar', $
     'spec_iabsmag200',  'spec_iabsmag200_ivar', $
     'spec_ilum200',     'spec_ilum200_ivar']

  ;; Values default to -9999.0, ivar to 0.0
  newtagdef = replicate('0.0', 20)
  val_id = lindgen(10)*2
  newtagdef[val_id] = '-9999.0'

  ;; long int for vagc_index
  newtagdef = ['-9999L', newtagdef]
  add_tags, struct, newtags, newtagdef, newstruct

  return, newstruct

end 

pro maxbcg::add_kcorrect, bcg, neigh

  ;; Calculate k-corrected mags and total luminosities
  print,'-----------------------------------------------------------------'
  print
  print,'k-correcting, luminosities, etc'
  extra_struct=$
    self->calculate_extra(bcg, $
                          bcg.photoz_cts, $
                          neigh, $
                          bcg[neigh.bcg_id].photoz_cts)

  print
  print,'Combining extra info'
  bcg.bcg_iabsmag = extra_struct.bcg_absmag[3]
  bcg.bcg_iabsmag_ivar = extra_struct.bcg_absmag_ivar[3]
  bcg.bcg_ilum = extra_struct.bcg_lumsolar[3]
  bcg.bcg_ilum_ivar = extra_struct.bcg_lumsolar_ivar[3]

  bcg.iabsmag200 = extra_struct.absmag200[3]
  bcg.iabsmag200_ivar = extra_struct.absmag200_ivar[3]
  bcg.ilum200 = extra_struct.lumsolar200[3]
  bcg.ilum200_ivar = extra_struct.lumsolar200_ivar[3]
  delvarx, extra_struct

  print,'-----------------------------------------------------------------'
  print
  print,'Doing spectra extra info'
  w=where(bcg.bcg_spec_z GT 0.0)
  match_multi, bcg[w].bcg_id, neigh.bcg_id, mneigh
  extra_struct=$
    self->calculate_extra(bcg[w], $
                          bcg[w].bcg_spec_z, $
                          neigh[mneigh], $
                          bcg[neigh[mneigh].bcg_id].bcg_spec_z)

  bcg[w].bcg_spec_iabsmag = extra_struct.bcg_absmag[3]
  bcg[w].bcg_spec_iabsmag_ivar = extra_struct.bcg_absmag_ivar[3]
  bcg[w].bcg_spec_ilum = extra_struct.bcg_lumsolar[3]
  bcg[w].bcg_spec_ilum_ivar = extra_struct.bcg_lumsolar_ivar[3]

  bcg[w].spec_iabsmag200 = extra_struct.absmag200[3]
  bcg[w].spec_iabsmag200_ivar = extra_struct.absmag200_ivar[3]
  bcg[w].spec_ilum200 = extra_struct.lumsolar200[3]
  bcg[w].spec_ilum200_ivar = extra_struct.lumsolar200_ivar[3]
  delvarx, extra_struct

end 

pro maxbcg::add_centerclass, bcg

    print,'Adding tag: kde_class'
    add_tags, bcg, 'kde_class', '-9999', newbcg
    delvarx, bcg
    bcg = temporary(newbcg)

    kde = self->get(/centerclass)
    w=where(kde.kde_class lt 0, nw)
    if nw ne 0 then kde[w].kde_class = -9999


    print,'Matching to center classification file'
    match, bcg.mem_match_id, kde.mem_match_id, mbcg, me, /sort
    if mbcg[0] eq -1 then message,'Failed to match centerclass file'
    bcg[mbcg].kde_class = kde[me].kde_class

end


;; Redshift cuts, rename tags, calculate new stuff, write as .st file
function maxbcg::gmbcg10_keeptags
	keeptags = $
		['objid',$
		 'mem_match_id', $
		 'ra', $
		 'dec', $
		 'model_counts', $
		 'model_counts_err', $
		 'photoz', $
		 'photoz_err', $
		 'spz', $
		 'ngals', $
		 'gm_scaledr', $
		 'gm_scaled_ngals', $
		 'gm_ngals_weighted']
	return, keeptags
end
function maxbcg::gmbcg10_extract_tags, str
	keeptags=self->gmbcg10_keeptags()
	newstr = extract_tags(str, keeptags)
	if n_tags(newstr) ne n_elements(keeptags) then begin
		match, keeptags, tag_names(str), mkeep, mstr
		remove, mkeep, keeptags
		print,'unmatched tags: ',keeptags
		on_error,2
		message,'halting'
	endif

	return, newstr
end

function maxbcg::gmbcg10_rename_tags, str
	oldnames = $
		['gm_scaledr', 'gm_scaled_ngals', 'gm_ngals_weighted']
	newnames = $
		['scale', 'sngals', 'wngals']
	newstr = rename_tags(str, oldnames, newnames)
	return, newstr
end
pro maxbcg::convert_gmbcg10, overwrite=overwrite

	fitsfile = self->catname(/fits)
	stfile   = self->catname()

	print,'fits file: ',fitsfile
	print,'new st file: ',stfile

	orig = mrdfits(fitsfile, 1)
	norig = n_elements(orig)

	keep = self->zselect(orig.photoz, nkeep)
	str = orig[keep]
	print,'kept ',nkeep,'/',norig,f='(a,i0,a,i0)'

	str = self->gmbcg10_extract_tags(str)
	str = self->gmbcg10_rename_tags(str)

	if not tag_exist(str[0], 'clambda') then begin 
		print
		print,'Adding clambda/ceta'

		str = self->add_csurvey(temporary(str))
	endif 

	if not keyword_set(overwrite) and fexist(stfile) then begin 
		key = ''
		read,'File '+stfile+' exists. Overwrite? (y/n) ', key
		if strlowcase(key) ne 'y' then return
	endif 

	print,'Writing file: ',stfile
	write_idlstruct, str, stfile

end 




;; Redshift cuts, rename tags, calculate new stuff, write as .st file
pro maxbcg::convert, overwrite=overwrite


  fitsfile = self->catname(/fits)
  stfile   = self->catname()
  nfitsfile  = self->catname(/fits, /neigh)
  nstfile    = self->catname(/neigh)

  print,'fits file: ',fitsfile
  print,'new st file: ',stfile
  print,'neighbor file: ',nfitsfile
  print,'new neighbor file: ',nstfile

  orig = mrdfits(fitsfile, 1)
  neigh = mrdfits(nfitsfile, 1)

  norig = n_elements(orig)

  print
  print,'Making redshift cuts'
  keep = self->zselect(orig.photoz_cts, nkeep)
  orig = orig[keep]


  ;; Color of BCG NOT assuming lrg shape
  ;print
  ;print,'Getting kgmr,krmi NOT assuming LRG'
  ;ft = self->fake_tsobj(orig)
  ;kc = sdss_kcorrect(orig.photoz_cts, $
  ;                   tsobj=ft, flux='model', absmag=absmag, amivar=amivar, $
  ;                   band_shift=self->band_shift(), omega0=0.27, omegal0=0.73)




  print,'Keeping '+ntostr(nkeep)+'/'+ntostr(norig)


  if not tag_exist(orig[0], 'clambda') then begin 
      print
      print,'Adding clambda/ceta'

      orig = self->add_csurvey(temporary(orig))
  endif 

  if self->catalog() ne 'dr406' then begin 
      message,'Check to see if this program works for catalog: '+self->catalog()
  endif 

  oldtags = $
    ['id', $
     'ngals_r200','tngals_r200',$
     'i_lum', 'r_lum', $
     'mean_kgr', 'mean_kri', $
     'med_kgr', 'med_kri', $
     'kgmr', 'krmi']
  newtags = $
    ['bcg_id', $
     'ngals200',  'tngals200',  $
     'ilum200','rlum200', $
     'mean_kgmr200','mean_krmi200',$
     'med_kgmr200','med_krmi200',$
     'bcg_kgmr','bcg_krmi']
  
  print
  newstruct = rename_tags(temporary(orig), oldtags, newtags)
  
  ;; remove run,rerun,camcol,field since ben overwrite id with 
  ;; Some of these will get added back in, just for clarity of code
  rmtags = ['run','rerun','camcol','field',$
            'bcg_rlum', 'bcg_ilum', $
            'rlum200', 'ilum200', $
            'bcg_kgmr', 'bcg_krmi', $
            'med_kgmr200', 'med_krmi200', $
            'mean_kgmr200', 'mean_krmi200',$
            'spectro_r_lum','spectro_i_lum', $
            'spectro_bcg_rlum', 'spectro_bcg_ilum']
  print
  print,'Removing tags: ',rmtags
  newstruct = remove_tags(temporary(newstruct), rmtags)
  
  print
  print,'Adding tags'
  newstruct = self->add_tags(temporary(newstruct))

  ;; Match to spectra
  print
  print,'Matching to lss spectro sample'
  match_struct = self->match2lss(newstruct)
  w = match_struct.mbcg
  newstruct.bcg_spec_z = -9999.0
  newstruct[w].bcg_spec_z = match_struct.z
  newstruct[w].bcg_vagc_index = match_struct.mlss

  ;; Re-index to a unique bcg_id
  print
  print,'Re-indexing'
  maxbcg_reindex, newstruct, neigh

  ;; Add center classification "kde_class"
  self->add_centerclass, newstruct

  ;; k-corrected stuff.  Assumes lrg for all
  self->add_kcorrect, newstruct, neigh

  ;; Color of BCG NOT assuming lrg shape
  print
  print,'Getting kgmr,krmi NOT assuming LRG'
  ft = self->fake_tsobj(newstruct)
  ; some problem with common blocks, need to call a separate one
  kc = self->sdss_kcorrect(newstruct.photoz_cts, $
                     tsobj=ft, flux='model', absmag=absmag, amivar=amivar, $
                     band_shift=self->band_shift(), omega0=0.27, omegal0=0.73)

  newstruct.bcg_kgmr = reform( absmag[1,*] - absmag[2,*] )
  newstruct.bcg_kgmr_ivar = reform( 1.0/(1.0/amivar[1,*] + 1.0/amivar[2,*] ) )
  newstruct.bcg_krmi = reform( absmag[2,*] - absmag[3,*] )
  newstruct.bcg_krmi_ivar = reform( 1.0/(1.0/amivar[2,*] + 1.0/amivar[3,*] ) )

  delvarx, ft, kc, absmag, amivar
  print
  print,'Re-ordering tags'
  front_tags = $
    ['bcg_id', $
     'photoid', $
     'stripe', $
     'ra','dec','clambda','ceta',$
     $
     'z','photoz_cts','photoz_z', 'photoz_zerr', $
     $
     'ngals','ngals200',$
     $
     'bcg_iabsmag',       'bcg_iabsmag_ivar', $
     'bcg_ilum',          'bcg_ilum_ivar', $
     'bcg_kgmr', 'bcg_kgmr_ivar', $
     'bcg_krmi', 'bcg_krmi_ivar', $
     'bcg_vagc_index', $
     'bcg_spec_z', $
     'bcg_spec_iabsmag',  'bcg_spec_iabsmag_ivar', $
     'bcg_spec_ilum',     'bcg_spec_ilum_ivar', $
     $
     'iabsmag200',        'iabsmag200_ivar', $
     'ilum200',           'ilum200_ivar', $
     'spec_iabsmag200',   'spec_iabsmag200_ivar', $
     'spec_ilum200',      'spec_ilum200_ivar', $
     $
     'cmodel_counts', $
     'imag', $
     'gmr', 'gmr_err', 'rmi', 'rmi_err', $
     'objc_prob_gal', $
     'corrselect_flags', $
     'maskflags','edge','lrg_flags']

  newstruct = reorder_tags( temporary(newstruct), front_tags )


  if not keyword_set(overwrite) and fexist(stfile) then begin 
      key = ''
      read,'File '+stfile+' exists. Overwrite? (y/n) ', key
      if strlowcase(key) ne 'y' then return

  endif 

  print,'Writing file: ',stfile
  write_idlstruct, newstruct, stfile



  if not keyword_set(overwrite) and fexist(nstfile) then begin 
      key = ''
      read,'File '+nstfile+' exists. Overwrite? (y/n) ', key
      if strlowcase(key) ne 'y' then return
  endif 

  print,'Writing neighbor file: ',nstfile
  write_idlstruct, neigh, nstfile


end 

pro maxbcg::stuff

  bcg = self->get()
  
  tmpdir = self->catdir()
  table = 'maxbcg'
  self->postgres::struct2table, bcg, table, $
    primary_key='bcg_id', $
    conn='user=postgres', $
    tmpdir=tmpdir, $
    status=status

  if status ne 0 then begin 
      message,'Error stuffing'
  endif 

  index_cols = ['photoid', 'stripe', $
                'z','photoz_cts', 'photoz_z', $
                'ngals', 'ngals200', $
                'bcg_vagc_index', 'bcg_spec_z', $
                'bcg_iabsmag', 'bcg_ilum', $
                'iabsmag200', 'ilum200','ra','dec']

  self->postgres::create_index, 'maxbcg', index_cols, conn='user=postgres'

  query = 'ANALYZE '+table
  print,query
  self->query, query, conn='user=postgres', status=status
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Error analyzing'
  endif 

  query = 'GRANT SELECT ON '+table+' TO sdss'
  print,query
  self->query, query, conn='user=postgres', status=status
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Error granting'
  endif 



end 



pro maxbcg::stuff_neighbors

  neigh = self->get(/neigh)
  
  w=where(neigh.bcg_id GE 0)
  neigh=neigh[w]

  neigh = reorder_tags(temporary(neigh), $
                       ['bcg_id', 'photoid', 'stripe', 'bcg_stripe'])


  tmpdir = self->catdir()
  table = 'maxbcg_neigh'
  self->postgres::struct2table, neigh, table, $
    conn='user=postgres', $
    tmpdir=tmpdir, $
    status=status

  if status ne 0 then begin 
      message,'Error stuffing'
  endif 

  index_cols = ['bcg_id','photoid', 'stripe', 'bcg_stripe']

  self->postgres::create_index, 'maxbcg_neigh', index_cols, $
    conn='user=postgres'

  query = 'ANALYZE maxbcg_neigh'
  print,query
  self->query, query, conn='user=postgres', status=status
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Error analyzing'
  endif 

  query = 'GRANT SELECT ON maxbcg_neigh TO sdss'
  print,query
  self->query, query, conn='user=postgres', status=status
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Error granting'
  endif 



end 



;; Obsolete; info now in main catalog
;; This makes a catalog of just the spectroscopic matched 
;; objects and makes measurement based on the spec z
pro maxbcg::remeasure_spec

  newfile = self->catname(/spec)
  print
  print,'Will write to file: ',newfile

  bcg = self->get()
  neigh = self->get(/neigh)

  ;; First extract redshifts for neighbors. This is
  ;; by far the simplest way

  zbcg = bcg.bcg_spec_z
  zneigh = bcg[neigh.bcg_id].bcg_spec_z

  w = where(bcg.bcg_spec_z GT 0.0, nw)

  bcg = bcg[w]
  zbcg = zbcg[w]

  ;; Now get neighbors that match.
  match_multi, bcg.bcg_id, neigh.bcg_id, mneigh
  neigh = neigh[mneigh]
  zneigh = zneigh[mneigh]


  wbad = where(zneigh LE 0, nbad)
  if nbad ne 0 then message,'Some bad neighbor redshifts'

  extra_struct=self->calculate_extra(bcg, zbcg, neigh, zneigh)

  copy_struct, extra_struct, bcg

  if fexist(newfile) then begin 
      key = ''
      read,'File '+newfile+' exists. Overwrite? (y/n) ', key
      if strlowcase(key) ne 'y' then return
  endif 
  print
  print,'Writing to file: ',newfile
  write_idlstruct, bcg, newfile


end 




; Designed to create a list of good objects in all the basic
; richness estimators.  Will set all values to -9999
pro maxbcg::badvals_fix, st

  badval = -9999

  wbad = $
    where((st.ilum200 ne st.ilum200) OR (st.bcg_ilum ne st.bcg_ilum) OR $
          (st.rlum200 ne st.rlum200) OR (st.bcg_rlum ne st.bcg_rlum), $
          nbad, comp=wgood, ncomp=ngood)

  if nbad ne 0 then begin 
      st[wbad].ngals200 = badval
      st[wbad].ilum200  = badval
      st[wbad].rlum200  = badval

      st[wbad].bcg_ilum = badval
      st[wbad].bcg_rlum = badval

      ;; This too, just for consistency
      st[wbad].tngals = badval
      st[wbad].tngals200 = badval


      st[wbad].mean_kgmr200 = badval
      st[wbad].mean_krmi200 = badval

      st[wbad].med_kgmr200 = badval
      st[wbad].med_krmi200 = badval

      st[wbad].num2_imag = badval
      st[wbad].num3_imag = badval
      st[wbad].num4_imag = badval
  endif 

  ;; bug in colors for ngals200 = 1
  w2 = where(st[wgood].ngals200 eq 1, n2)
  if n2 ne 0 then begin 
      st[wgood[w2]].mean_kgmr200 = st[wgood[w2]].bcg_kgmr
      st[wgood[w2]].mean_krmi200 = st[wgood[w2]].bcg_krmi

      st[wgood[w2]].med_kgmr200 = st[wgood[w2]].bcg_kgmr
      st[wgood[w2]].med_krmi200 = st[wgood[w2]].bcg_krmi
  endif 

end 






pro maxbcg::reindex, bcg, neigh, $
          keepneigh=keepneigh, verbose=verbose

  nneigh = n_elements(neigh)
  keepneigh = lonarr(nneigh)

  ;; Ben indexed by stripe, then id within stripe
  ustripes = bcg[ rem_dup(bcg.stripe) ].stripe
  nstripe = n_elements(ustripes)
  
  for ist = 0L, nstripe-1 do begin 

      stripe = ustripes[ist]

      print,'//////////////////////////////////////////////////////////////'
      print,'Doing stripe: ',stripe


      ;; Must mach bcg stripe and bcg id
      wneighst = where( neigh.bcg_stripe eq stripe, nneighst )
      wbcgst   = where( bcg.stripe eq stripe, nbcgst)

      ;; histogram the neighbor bcgid
      print
      print,'Histogramming neighbor bcg_id'
      h = histogram(neigh[wneighst].bcg_id, $
                    min=0, max=max(bcg.bcg_id), $
                    rev=rev)

      print
      print,'Matching to bcgs'
      print

      if keyword_set(verbose) then begin 
          print,$
            'bcg_id','nmatch','Ngals200',$
            format='(A10,A10,A10)'
      endif 

      for iibcg=0L, nbcgst-1 do begin 

          ibcg = wbcgst[iibcg]

          bcg_id = bcg[ibcg].bcg_id

          bcg[ibcg].bcg_id = ibcg
     
          ;; get neighbors
          if rev[bcg_id] ne rev[bcg_id+1] then begin 


              wmatch = rev[ rev[bcg_id]:rev[bcg_id+1]-1 ]
              wmatch = wneighst[wmatch]
              nmatch = n_elements(wmatch)
                  
              ;; re-index
              neigh[wmatch].bcg_id = ibcg
              keepneigh[wmatch] = 1

              if keyword_set(verbose) then begin 

                  print,$
                    strn(bcg_id),$
                    strn(nmatch),$
                    strn(bcg[ibcg].ngals200),$
                    format='(A10,A10,A10)'
                  
                  if nmatch ne bcg[ibcg].ngals200 then begin 
                      print,' ----  Ngals different'
                  endif 
                               
                  key = prompt_kbrd()
                  if key eq 'q' then return
              endif 
              
          endif else begin 
;              print,'NO NEIGHBORS FOUND.  ID: ',ibcg,$
;                    '  Ngals200: ',bcg[ibcg].ngals200
          endelse 

      endfor ;; over bcgs in this stripe

  endfor ;; over stripes

  keep = where(keepneigh ne 0, nkeep, comp=comp,  ncomp=ncomp)
  print
  print,'Neighbors that matched: '+ntostr(nkeep)+'/'+ntostr(nneigh)
  if ncomp ne 0 then begin 
      neigh[comp].bcg_id = -9999
  endif 
;  neigh = neigh[keep]


end 


; 
; Try to pick out some bad objects
;

; This really just gets high bcg luminosity objects
function maxbcg::getbadbcg, bcg=bcg, badz=badz

    if n_elements(bcg) eq 0 then begin
        bcg = self->get()
    endif

    if n_elements(badz) eq 0 then badz=-9999999.0


    logl=10.0 + alog10(bcg.bcg_ilum)
    logl_lim = 11.0
    w=where(bcg.bcg_spec_z ne badz $
            and bcg.photoz_cts gt 0.1 and bcg.photoz_cts lt 0.3 $
            and bcg.ngals200 ge 3 $
            and logl gt logl_lim)

    return, w

end

pro maxbcg::make_badbcg_table, bcg=bcg
    table = 'badbcg'
    w=self->getbadbcg(bcg=bcg, badz=badz)

    tbcg1 = bcg[w]
    self->badbcg_imsize, tbcg1, pixscale, wpix, hpix

    newtags = $
        ['index','imsize','pixscale',$
         'erinmark','johnmark','timmark','benmark','eduardomark','mattmark',$
         'pubmark']
    newtypes = ['0L','0L','0.0',$
                '0','0','0','0','0','0','0']
    tbcg = struct_addtags(tbcg1, newtags, newtypes)

    tbcg.imsize = wpix
    tbcg.pixscale = pixscale
    tbcg = reorder_tags(tbcg, 'index')

    tbcg = tbcg[reverse(sort(tbcg.bcg_ilum))]
    id = lindgen(n_elements(tbcg))
    tbcg.index = id

    self->struct2table, tbcg, table, primary_key='index'

    icols = ['bcg_id', 'bcg_ilum', 'ilum200', 'ngals200', $
            'photoz_cts', 'bcg_spec_z',$
            'erinmark','johnmark','timmark','benmark','eduardomark','pubmark']
    self->create_index, table, icols

end

pro maxbcg::badbcg_imsize, badbcg, pixscale, wpix, hpix
    ; pixscale: half resolution of SDSS images arcsec per pixel
    pixscale = 0.396

    ; physical scale of image in kpc
    width = 1000
    D = angdist(0.0, badbcg.photoz_cts, omega_m=0.25)*1000 ; kpc/h
    ; radians
    wrad = width/D
    ; arcsec
    warc = wrad*180.0/!pi*3600.0
    ; pixels
    wpix = warc/pixscale

    ; now use half the resolution
    pixfac=1.5
    pixscale = pixscale*pixfac
    wpix=wpix/pixfac
    hpix=wpix

    wpix=long(wpix)
    hpix=long(hpix)
end
pro maxbcg::make_badbcg_html_php, file, nper=nper, badspec=badspec, bcg=bcg

    if keyword_set(badspec) then begin
        badz=-9999.0
    endif

    if n_elements(nper) eq 0 then nper=50

    w=self->getbadbcg(bcg=bcg, badz=badz)

    s=reverse(sort(bcg[w].bcg_ilum))
    w=w[s]
    nw = n_elements(w)

    tags = ['bcg_id','photoz_cts', 'bcg_spec_z', 'bcg_ilum', 'ilum200']
    extracol = extract_tags(bcg[w], tags)
    extracol = rename_tags(extracol, ['bcg_ilum', 'ilum200'], ['logbcglum', 'logl200'])
    extracol.logbcglum = 10.0+alog10(extracol.logbcglum)
    extracol.logl200 = 10.0+alog10(extracol.logl200)

    ; pixscale: half resolution of SDSS images arcsec per pixel
    pixscale = 0.396

    ; physical scale of image in kpc
    width = 1000
    ;D = angdist(0.0, bcg[w].bcg_spec_z, omega_m=0.25)*1000 ; kpc/h
    D = angdist(0.0, bcg[w].photoz_cts, omega_m=0.25)*1000 ; kpc/h
    ; radians
    wrad = width/D
    ; arcsec
    warc = wrad*180.0/!pi*3600.0
    ; pixels
    wpix = warc/pixscale

    ; now use half the resolution
    pixfac=1.5
    pixscale = pixscale*pixfac
    wpix=wpix/pixfac
    hpix=wpix

    wpix=long(wpix)
    hpix=long(hpix)

    openw, lun, file, /get_lun
    printf, lun, '<html><body>'
    printf, lun, '<table border=1>'
    hline = '  <tr><th>number</th><th>bcg id</th><th>photoz</th><th>specz</th><th>log(lbcg)</th><th>log(l200)</th><th>Images</th></tr>'
    printf, lun, hline
    for i=0L, nw-1 do begin
        line = '  <tr>'

        url = 'bcg.php?'
        ii = w[i]
        url = url + 'id='+ntostr(bcg[ii].bcg_id)
        url = url + '&ra='+ntostr(bcg[ii].ra)
        url = url + '&dec='+ntostr(bcg[ii].dec)

        url = url + '&width='+ntostr(wpix[i])
        url = url + '&height='+ntostr(hpix[i])

        url = url + '&scale='+ntostr(pixscale)

        url = url + '&photoz='+ntostr(bcg[ii].photoz_cts)
        url = url + '&specz='+ntostr(bcg[ii].bcg_spec_z)
        url = url + '&loglbcg='+ntostr(10.0+alog10(bcg[ii].bcg_ilum))
        url = url + '&logl200='+ntostr(10.0+alog10(bcg[ii].ilum200))

        url = 'bcg.php?index='+ntostr(i)

        line = line+'<td>'+ntostr(i)+'</td>'
        line = line+'<td>'+ntostr(bcg[ii].bcg_id)+'</td>'
        line = line+'<td>'+ntostr(bcg[ii].photoz_cts)+'</td>'
        line = line+'<td>'+ntostr(bcg[ii].bcg_spec_z)+'</td>'
        line = line+'<td>'+ntostr(10.0+alog10(bcg[ii].bcg_ilum))+'</td>'
        line = line+'<td>'+ntostr(10.0+alog10(bcg[ii].ilum200))+'</td>'
        
        line = line+'<td><a href="'+url+'" target="_blank">images</a></td>'
        line = line + '</tr>'
        printf, lun, line

        if ( ((i+1) mod 10) eq 0) then printf, lun, hline
    endfor
    printf, lun, '</table>'
    free_lun, lun


end


pro maxbcg::make_badbcg_html, filefront, nper=nper, badspec=badspec, bcg=bcg

    if keyword_set(badspec) then begin
        specdiff=0.05 
        badz=-9999.0
    endif else begin
        specdiff=0.0
        badz=9999999.0
    endelse

    if n_elements(nper) eq 0 then nper=50

    w=self->getbadbcg(bcg=bcg, badz=badz, specdiff=specdiff)

    s=reverse(sort(bcg[w].bcg_ilum))
    w=w[s]
    nw = n_elements(w)

    tags = ['bcg_id','photoz_cts', 'bcg_spec_z', 'bcg_ilum', 'ilum200']
    extracol = extract_tags(bcg[w], tags)
    extracol = rename_tags(extracol, ['bcg_ilum', 'ilum200'], ['logbcglum', 'logl200'])
    extracol.logbcglum = 10.0+alog10(extracol.logbcglum)
    extracol.logl200 = 10.0+alog10(extracol.logl200)

    ; pixscale: half resolution of SDSS images arcsec per pixel
    pixscale = 0.396

    ; physical scale of image in kpc
    width = 1000
    ;D = angdist(0.0, bcg[w].bcg_spec_z, omega_m=0.25)*1000 ; kpc/h
    D = angdist(0.0, bcg[w].photoz_cts, omega_m=0.25)*1000 ; kpc/h
    ; radians
    wrad = width/D
    ; arcsec
    warc = wrad*180.0/!pi*3600.0
    ; pixels
    wpix = warc/pixscale

    ; now use half the resolution
    pixscale = pixscale*2
    wpix=wpix/2
    hpix=wpix

    ndiv = nw/nper
    nmod = nw mod nper
    if nmod ne 0 then ndiv = ndiv + 1

    num = lonarr(ndiv)
    num[*] = nper
    num[ndiv-1] = nmod
    beg = 0L
    for i=0L, ndiv-1 do begin
        file = filefront+ntostr(i, f='(i02)')+'.html'
        print,'Writing to file: ',file
        ww = w[beg:beg+num[i]-1]    
        casimagelinks, bcg[ww].ra, bcg[ww].dec, file=file, $
            scale=pixscale, width=wpix[ww], height=hpix[ww], $
            extra_columns=extracol[beg:beg+num[i]-1], /embed, headmod=1
        beg = beg+num[i]
    endfor

end







; analysis of the high bcg luminosity sample
pro maxbcg::plot_highbcglum_deltaz, ngals_min, bcg=bcg


    case ngals_min of
        3: nper2=1000
        10: nper2=200
        else: message,'Unsupported ngals min: '+ntostr(ngals_min)
    endcase

    dir = esheldon_config('plot_dir')
    dir = path_join(dir, 'maxbcg/badbcg')
    psfile = path_join(dir, 'deltaz-ngals200-'+ntostr(ngals_min,format='(i02)')+'.eps')
    begplot, psfile, /encapsulated, /color


    w=self->getbadbcg(bcg=bcg)
    s=reverse(sort(bcg[w].bcg_ilum))
    w=w[s]

    ww=where(bcg[w].bcg_spec_z gt 0 and bcg[w].ngals200 ge ngals_min)

    ww2=where(bcg.photoz_cts gt 0.1 $
                and bcg.photoz_cts lt 0.3 $
                and bcg.bcg_spec_z gt 0 $
                and bcg.ngals200 ge ngals_min)

    bs=binner(bcg[w[ww]].photoz_cts,$
                bcg[w[ww]].photoz_cts - bcg[w[ww]].bcg_spec_z, nper=100)
    bs2=binner(bcg[ww2].photoz_cts, bcg[ww2].photoz_cts - bcg[ww2].bcg_spec_z, $
                nper=nper2)

    hcolor=!blue
    pplot, bs2.xmean, bs2.ymean, yerr=bs2.yerr, psym=4, yrange=[-0.025, 0.05],$
        xtitle='photoz', ytitle='photoz-specz', aspect=1
    pplot, bs2.xmean, bs2.ymean+bs2.ysdev, /over
    pplot, bs2.xmean, bs2.ymean-bs2.ysdev, /over
    pplot, bs.xmean, bs.ymean, yerr=bs.yerr, psym=8,/over,color=hcolor

    legend, ['all', 'high bcglum'], psym=[4,8], color=[!p.color, hcolor]
    legend, 'N200 >= '+ntostr(ngals_min), /right

    colprint, bs2.xmean, bs2.ysdev
    print,mean(bs2.ysdev)

    endplot,/trim_bbox

end


pro maxbcg::plot_highbcglum_marked, user, topn=topn

    dir = esheldon_config('plot_dir')
    dir = path_join(dir, 'maxbcg/badbcg')
    psfile = path_join(dir, 'markfrac_'+user+'.eps')
    begplot, psfile, /encapsulated, /color

    if n_elements(topn) eq 0 then topn=500
    table = 'badbcg'
    markcol = user+'mark'
    topnstr = ntostr(topn)

    ; Order by index -> order reverse by bcg luminosity
    query = 'select * from '+table+' where index < '+topnstr+' order by index'
	pg=obj_new('postgres')
    st = pg->query(query)

    help, st

    if not tag_exist(st, markcol, index=tag) then begin
        message,'user '+user+' has no mark column'
    endif

    ;wdust = where((st.markcol and 2^0) ne 0, ndust)
    ;wstar = where((st.markcol and 2^1) ne 0, nstar)
    ;wspur = where((st.markcol and 2^2) ne 0, nspur)
    ;wfollow = where((st.markcol and 2^3) ne 0, nfollow)
    ldust = (st.(tag) and 2^0) ne 0
    lstar = (st.(tag) and 2^1) ne 0
    lspur = (st.(tag) and 2^2) ne 0
    lfollow = (st.(tag) and 2^3) ne 0

    ii = lindgen(topn)+1
    cumfrac_dust = total(ldust, /cum)/lindgen(topn)
    cumfrac_star = total(lstar, /cum)/lindgen(topn)
    cumfrac_spur = total(lspur, /cum)/lindgen(topn)
    cumfrac_follow = total(lfollow, /cum)/lindgen(topn)

    cumfrac_dust_err = sqrt(total(ldust, /cum))/lindgen(topn)
    cumfrac_star_err = sqrt(total(lstar, /cum))/lindgen(topn)
    cumfrac_spur_err = sqrt(total(lspur, /cum))/lindgen(topn)
    cumfrac_follow_err = sqrt(total(lfollow, /cum))/lindgen(topn)


    !p.multi=[0,2,2]
    !p.charsize=1
    !p.thick=1
    ytitle = 'Fraction With Dust (> Lum)'
    xtitle='BCG Lum'
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_dust, yerr=cumfrac_dust_err, $
        psym=8, yrange=[0,1], ytitle=ytitle, xtitle=xtitle, aspect=1, hat=0
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_dust, color=!blue, /over

    ytitle = 'Fraction of Stars (> Lum)'
    xtitle='BCG Lum'
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_star, yerr=cumfrac_star_err, $
        psym=8, yrange=[0,1], ytitle=ytitle, xtitle=xtitle, aspect=1, hat=0
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_star, color=!blue, /over
    ;axis, yaxis=1, title=ytitle, /noerase

    ytitle = 'Fraction Spurious (> Lum)'
    xtitle='BCG Lum'
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_spur, yerr=cumfrac_spur_err, $
        psym=8, yrange=[0,1], ytitle=ytitle, xtitle=xtitle, aspect=1, hat=0
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_spur, color=!blue, /over

    ytitle = 'Fraction for Follow-up (> Lum)'
    xtitle='BCG Lum'
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_follow, yerr=cumfrac_follow_err, $
        psym=8, yrange=[0,1], ytitle=ytitle, xtitle=xtitle, aspect=1, hat=0
    pplot, 10.0+alog10(st.bcg_ilum), cumfrac_follow, color=!blue, /over

    endplot,/trim_bbox

    !p.multi=0
end










;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; My re-measurement stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function maxbcg::defval
  return, -9999.0
end 
function maxbcg::band_shift
  return, 0.25
end 

function maxbcg::sunabsmag, band_shift=band_shift

  common maxbcg_sunabsmag_block, band_shift_save, sunabsmag_save

  if n_elements(band_shift) eq 0 then band_shift = self->band_shift()

  ;; First call
  if n_elements(band_shift_save) eq 0 then begin 
      band_shift_save = band_shift
      command='sunabsmag_save = k_solar_magnitudes(band_shift=band_shift, /silent)'
	  if not execute(command) then begin
		  message,'Could not execute command: '+command
	  endif
  endif 

  if band_shift eq band_shift_save then begin 
      sun = sunabsmag_save
  endif else begin 
      command='sun = k_solar_magnitudes(band_shift=band_shift, /silent)'
	  if not execute(command) then begin
		  message,'Could not execute command: '+command
	  endif

      band_shift_save = band_shift
      sunabsmag_save = sun
  endelse

  return, sun

end 

function maxbcg::sunknmgy, band_shift=band_shift
  sunabsmag = self->sunabsmag(band_shift=band_shift)
  sunknmgy = 10.0^( -0.4*sunabsmag )*1.e9
  return,sunknmgy
end 






pro maxbcg::fake_modelcounts, struct, mags, magerr

  ;; k_sdssfix checks for -9999.0 in mags and magerr and for those objects
  ;; sets ivar = 0
  def = self->defval()
  num = n_elements(struct)
  mags = replicate(def, 5, num)
  magerr = replicate(def, 5, num)

  nmgy = replicate(def,5, num)


  imag = struct.imag
  rmag = imag + struct.rmi
  gmag = rmag + struct.gmr

  wgood = where(gmag GT 10 AND gmag LT 24 AND (gmag eq gmag), ngood)
  if ngood ne 0 then begin 
      mags[1,wgood] = gmag[wgood]
      magerr[1,wgood] = struct[wgood].gmr_err
  endif 
  wgood = where(rmag GT 10 AND rmag LT 24 AND (rmag eq rmag), ngood)
  if ngood ne 0 then begin 
      mags[2,wgood] = rmag[wgood]
      magerr[2,wgood] = struct[wgood].gmr_err
  endif 
  wgood = where(imag GT 10 AND imag LT 24 AND (imag eq imag), ngood)
  if ngood ne 0 then begin 
      mags[3,wgood] = imag[wgood]
      magerr[3,wgood] = struct[wgood].rmi_err
  endif 


end 

function maxbcg::fake_tsobj, struct

  arrval = replicate(-9999.0, 5)
  tsobj = {counts_model: arrval, counts_modelerr: arrval, $
           reddening: fltarr(5)}

  num = n_elements(struct)
  tsobj = replicate(tsobj, num)

  self->fake_modelcounts, struct, mags, magerr
  tsobj.counts_model = temporary(mags)
  tsobj.counts_modelerr = temporary(magerr)
  return, tsobj
end 

pro maxbcg::lup2nmgy, struct, nmgy, ivar

  self->fake_modelcounts, struct, mags, magerr
  
  ;; convert to maggies. Will take the default errors
  ;; above and turn into ivar=0.  Also converts to AB system
  k_sdssfix, mags, magerr, maggies, ivar

  wbad = where(maggies ne maggies OR ivar ne ivar,nbad, $
               comp=wgood, ncomp=ngood)

  if nbad ne 0 then begin 
      maggies[wbad] = self->defval()
      ivar[wbad] = 0.0
  endif 

  nmgy = maggies
  nmgy[*] = self->defval()
  if ngood ne 0 then begin 
      nmgy[wgood] = maggies[wgood]*1.e9
      ivar[wgood] = ivar[wgood]*1.e-18
  endif 

end 






;; We assume good values in all these conversions
function maxbcg::lumsolar_ivar2amivar, lumsolar, lumsolar_ivar

  on_error, 2
  if n_elements(lumsolar) eq 0 OR n_elements(lumsolar_ivar) eq 0 then begin 
      print,'-syntax: mb->lumsolar_ivar2amivar(lumsolar, lumsolar_ivar)'
      print
      message,'Halting'
  endif 
  amivar = lumsolar^2*lumsolar_ivar*(0.4*alog(10.))^2
  return, amivar
end 
function maxbcg::lumsolar2absmag, lumsolar, filter=filter, band_shift=band_shift, amivar=amivar, lumsolar_ivar=lumsolar_ivar

  on_error, 2
  if n_elements(lumsolar) eq 0 then begin 
      print,'-syntax: mb->lumsolar2absmag(lumsolar, filter=, band_shift=, amivar=, lumsolar_ivar=)'
      print
      message,'Halting'
  endif 

  sun = self->sunabsmag(band_shift=band_shift)
  ;; assumes units of 10^10 solar
  if n_elements(filter) eq 0 then begin 
      ;; assume 5 filters
      absmag = lumsolar
      replicate_inplace, absmag, self->defval()
      for filter=0,4 do begin 
          tlum = reform( lumsolar[filter,*] )
          logLumSolar = 10.0 + alog10(tlum)
          absmag[filter,*] = -2.5*logLumSolar + sun[filter]
      endfor 
  endif else begin 
      logLumSolar = 10.0 + alog10(lumsolar)
      absmag = -2.5*logLumSolar + sun[filter]
  endelse 

  if arg_present(amivar) then begin 
      if n_elements(lumsolar_ivar) eq 0 then begin 
          message,'To calculate amivar you must also enter lumsolar_ivar',/inf
      endif 
      amivar = self->lumsolar_ivar2amivar(lumsolar, lumsolar_ivar)
  endif 

  return,absmag
end 

function maxbcg::amivar2lumsolar_ivar, amivar, lumsolar
  on_error, 2
  if n_elements(amivar) eq 0 OR n_elements(lumsolar) eq 0 then begin 
      print,'-syntax: mb->amivar2lumsolar_ivar(amivar, lumsolar)'
      print
      message,'Halting'
  endif 

  lumsolar_ivar = amivar
  replicate_inplace, lumsolar_ivar, 0.0
  w=where(lumsolar GT 0, nw)
  if nw ne 0 then begin 
      lumsolar_ivar[w] = amivar[w]/lumsolar[w]^2/(0.4*alog(10.))^2
  endif 
  return, lumsolar_ivar
end 



function maxbcg::absmag2lumsolar, absmag, filter=filter, band_shift=band_shift, lumsolar_ivar=lumsolar_ivar, amivar=amivar, verbose=verbose
  ;; assumes units of 10^10 solar

  on_error, 2
  if n_elements(absmag) eq 0 then begin 
      print,'-syntax: mb->absmag2lumsolar(absmag, filter=, band_shift=, lumsolar_ivar=, amivar=)'
      print
      message,'Halting'
  endif 

  sun = self->sunabsmag(band_shift=band_shift)
  if n_elements(filter) eq 0 then begin 
      if keyword_set(verbose) then message,'Doing all filters',/inf
      lumsolar = absmag
      replicate_inplace, lumsolar, self->defval()
      for filter=0,4 do begin 
          tabsmag = reform(absmag[filter,*])
          logLumSolar = (tabsmag - sun[filter])/(-2.5)
          lumsolar[filter,*] = 10.0^(logLumSolar - 10.0)
      endfor 
  endif else begin 
      if keyword_set(verbose) then message,'Doing single filter',/inf
      logLumSolar = (absmag - sun[filter])/(-2.5)
      lumsolar = 10.0^(logLumSolar - 10.0)
  endelse 

  if arg_present(lumsolar_ivar) then begin 
      if n_elements(amivar) eq 0 then begin 
          message,'To calculate lumsolar_ivar you must also enter amivar',/inf
          message,'Ignoring',/inf
      endif else begin 
          lumsolar_ivar = self->amivar2lumsolar_ivar(amivar, lumsolar)
      endelse 
  endif 


  return, lumsolar
end 


function maxbcg::kcorrect, struct, z

  print
  print,'Making fake tsobj'
  tsobj = self->fake_tsobj(struct)

  ;; note, defined so that kmag = mag - kcorrect
  ;; Thus, kmaggies = maggies/kmaggies_correct
  print,'Doing lrg kcorr'
  tm = systime(1)

  band_shift = self->band_shift()

  command='kcorrect_lrg = '+$
    'sdss_kcorrect(z, tsobj=tsobj, '+$
                  'absmag=absmag_lrg, amivar=amivar_lrg, '+$
                  'band_shift=band_shift, omega0=0.27, omegal0=0.73, '+$
                  '/lrg)'

  if not execute(command) then begin
	  message,'failed to execute command: '+command
  endif
  ptime,systime(1)-tm

  arrval = replicate(self->defval(), 5)
  kstruct = {absmag: arrval, $
             absmag_ivar: fltarr(5), $
             lumsolar: arrval, $
             lumsolar_ivar: fltarr(5), $
             kcorrect: arrval}
  kstruct = replicate(kstruct, n_elements(struct))
  kstruct.absmag = absmag_lrg
  kstruct.absmag_ivar = amivar_lrg
  kstruct.kcorrect = kcorrect_lrg

  kstruct.lumsolar = self->absmag2lumsolar(absmag_lrg, $
                                           band_shift=band_shift, $
                                           amivar=amivar_lrg, $
                                           lumsolar_ivar=lumsolar_ivar)
  ;; will only do lumsolar_ivar for the filters we gave, 
  ;; g,r,i   
  for filter=1,3 do begin 
      kstruct.lumsolar_ivar[filter] = reform(lumsolar_ivar[filter,*])
  endfor 

  return, kstruct

end 






;; These must have been run through reindex 
function maxbcg::calculate_extra, bcg, zbcg, neigh, zneigh, $
               kstruct_bcg=kstruct_bcg, kstruct_neigh=kstruct_neigh

  nbcg = n_elements(bcg)
  nneigh = n_elements(neigh)
  nzb = n_elements(zbcg)
  nzn = n_elements(zneigh)
  if nbcg eq 0 OR nneigh eq 0 OR nzb eq 0 OR nzn eq 0 then begin 
      print,'-Syntax: est=mb->calculate_extra(bcg, zbcg, neigh, zneigh)'
      print
      message,'Halting'
  endif 

  if n_elements(kstruct_bcg) eq 0 then begin 
      print
      print,'k-correcting bcgs'
      kstruct_bcg   = self->kcorrect(bcg, zbcg)
  endif 

  if n_elements(kstruct_neigh) eq 0 then begin 
      print
      print,'k-correcting neighbors'
      kstruct_neigh = self->kcorrect(neigh, zneigh)
  endif 

  arrval = replicate(self->defval(), 5)

  estruct = {$
              bcg_absmag: arrval, $
              bcg_absmag_ivar: fltarr(5), $
              bcg_lumsolar: arrval, $
              bcg_lumsolar_ivar: fltarr(5), $
              $
              absmag200: arrval, $
              absmag200_ivar: fltarr(5), $
              lumsolar200: arrval,$
              lumsolar200_ivar: fltarr(5) $
            }

  estruct = replicate(estruct, nbcg)

  estruct.bcg_absmag = kstruct_bcg.absmag
  estruct.bcg_absmag_ivar = kstruct_bcg.absmag_ivar

  estruct.bcg_lumsolar = kstruct_bcg.lumsolar
  estruct.bcg_lumsolar_ivar = kstruct_bcg.lumsolar_ivar
             

  ;; For ngals200 == 1, just copy
  w1 = where(bcg.ngals200 eq 1,n1, comp=w, ncomp=nw)

  if n1 ne 0 then begin 

      estruct[w1].absmag200 = estruct[w1].bcg_absmag
      estruct[w1].absmag200_ivar = estruct[w1].bcg_absmag_ivar

      estruct[w1].lumsolar200 = estruct[w1].bcg_lumsolar
      estruct[w1].lumsolar200_ivar = estruct[w1].bcg_lumsolar_ivar

  endif 
  
  band_shift = self->band_shift()

  print
  print,'Getting cumulative values'
  tm = systime(1)
  maxid = max(bcg.bcg_id)
  hist = histogram(neigh.bcg_id, min=0, max=maxid, rev=rev)

  for ii=0L, nw-1 do begin 

      ;; w is complement of ngals200 = 1
      i = w[ii]

      bcg_id = bcg[i].bcg_id

      if rev[bcg_id] ne rev[bcg_id+1] then begin 
          wm = rev[ rev[bcg_id]:rev[bcg_id+1]-1 ]

          wbad = where(zneigh[wm] ne zbcg[i], nbad)
          if nbad ne 0 then message,'Bad redshifts'

          for filter=0,4 do begin               
              lum = [reform(kstruct_neigh[wm].lumsolar[filter]), $
                     kstruct_bcg[i].lumsolar[filter]]

              ivar = [reform(kstruct_neigh[wm].lumsolar_ivar[filter]), $
                      kstruct_bcg[i].lumsolar_ivar[filter]]


              wivar = where(ivar ne 0, ngood)
              if ngood GT 0 then begin 
                  
                  totlum = total(lum[wivar])
                  totivar = 1.0/total( 1.0/ivar[wivar] )

                  estruct[i].lumsolar200[filter] = totlum
                  estruct[i].lumsolar200_ivar[filter] = totivar

                  estruct[i].absmag200[filter] = $
                    self->lumsolar2absmag(totlum, $
                                          filter=filter, $
                                          band_shift=band_shift, $
                                          amivar=totamivar, $
                                          lumsolar_ivar=totivar)
                  estruct[i].absmag200_ivar[filter] = totamivar
              
              endif 

          endfor 


      endif 

  endfor 

  ptime, systime(1)-tm

  return, estruct
end 










;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Convert Ben's and my random catalogs to ascii versions for Idit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro maxbcg::make_idit, random=random

  ;; Convert to ascii with limited tags for idit

  if not keyword_set(random) then begin 
      bcg = self->get()
      nbcg = n_elements(bcg)

      outstruct = {ra:0d, dec:0d, $
                   ngals:0, $
                   z:0.0, bcg_spec_z:0.0, $
                   id:0L, $
                   maskflags:0}

      outstruct = replicate(outstruct, nbcg)
      copy_struct, bcg, outstruct

      name = self->catname(/idit)
      outfile = repstr(name, '.fit', '_idit.st')
      print
      print,'Writing to outfile: ',outfile
      write_idlstruct, outstruct, outfile, /ascii
      
  endif else begin 

      numrand = long( sdssidl_config('maxbcg_numrand') )
      newstruct = {ra:0d, dec:0d, maskflags:0}

      for i=0L, numrand-1 do begin 

          file = self->catname(randnum=i, /idit)
          print,'Reading file: ',file

          newfile = repstr(files[i], '.fit', '.st')
          newfile = repstr(newfile, 'random', 'random-idit-eq')

      
          outstruct = replicate(newstruct, n_elements(instruct))
          csurvey2eq, instruct.clambda, instruct.ceta, ra, dec
          
          outstruct.ra = ra
          outstruct.dec = dec
          outstruct.maskflags = instruct.maskflags
          
          print,'Writing file: ',newfile
          write_idlstruct, outstruct, newfile, /ascii

      endfor 

  endelse 

end 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The input catalog for Ben's algorithm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function maxbcg::input_stripes, nstripes=nstripes

  reg1stripes = 9+lindgen(31)
  reg2stripes = [76, 82, 86]

  stripes = [reg1stripes, reg2stripes]


  nStripes = n_elements(stripes)
  return,stripes

end 

function maxbcg::input_dir
  dir = concat_dir(esheldon_config('lensinput_dir'),'maxbcg_input')
  catname=self->catalog()
  dir = path_join(esheldon_config('lensinput_dir'),['maxbcg_input',catname])
  return,dir
end 

function maxbcg::ryan_file, stripe

  if n_elements(stripe) eq 0 then begin 
      print,'-Syntax: file = ml->ryan_file(stripe)'
      return,''
  endif 

  dir = self->input_dir()
  stripeString = stripe2string(stripe)
  file = concat_dir(dir, 'stripe'+stripeString+'_cluster.fit')
  
  return,file
end 

function maxbcg::ryan_get, stripe
  file = self->ryan_file(stripe)
  print
  print,'Reading file: ',file
  struct = mrdfits(file,1)
  return,struct
end 

function maxbcg::input_file, stripe, database=database

  if n_elements(stripe) eq 0 then begin 
      print,'-Syntax: file = ml->input_file(stripe, /database)'
      return,''
  endif 

  dir = self->input_dir()
  stripeString = stripe2string(stripe)
  if keyword_set(database) then begin 
      inFile = 'maxbcg_input_stripe'+stripeString+'.pgsql'
  endif else begin 
      inFile = 'maxbcg_input_stripe'+stripeString+'.st'
  endelse 

  inFile = concat_dir(dir, inFile)
  return,inFile
end 



function maxbcg::getrunphotoz, run, rerun, nyu=nyu, status=status

  if n_params() LT 2 then begin 
      message,'-Syntax: mb->getrunphotoz(run, rerun, /nyu, status=)'
  endif 

  if keyword_set(nyu) then begin 
      query = $
        'SELECT '+$
          'photoid, photoz_z, '+$
          'photoz_zerr4 as photoz_zerr, photoz_zwarning '+$
        'FROM '+$
          'zphot '+$
        'WHERE '+$
          'run = '+ntostr(run)+' AND rerun = '+ntostr(rerun)
  endif else begin 
      minid = ntostr( sdss_photoid(run, rerun, 1, 1, 0) )
      maxid = ntostr( sdss_photoid(run, rerun, 6, 5000, 5000) )
      
      query = $
        'SELECT '+$
          'photoid, photoz_z, '+$
          'photoz_zerr4 as photoz_zerr, photoz_zwarning '+$
        'FROM '+$
          'zphot '+$
        'WHERE '+$
          'photoid BETWEEN '+minid+' AND '+maxid
  endelse 
      
  print,query
  pg=obj_new('postgres')
  zst = pg->query(query,status=status)
  
  return,zst

end 

function maxbcg::getrunz, run, rerun, status=status

  if n_params() LT 2 then begin 
      message,'-Syntax: mb->getrunz(run, rerun, status=)'
  endif 

  minid = ntostr( sdss_photoid(run, rerun, 1, 1, 0) )
  maxid = ntostr( sdss_photoid(run, rerun, 6, 5000, 5000) )
      
  query = $
    'SELECT '+$
    '  match_photoid as photoid, z, '+$
    'FROM '+$
    '  specgal '+$
    'WHERE '+$
    '  match_photoid BETWEEN '+minid+' AND '+maxid
      
  print,query
  pg=obj_new('postgres')
  zst = pg->query(query,status=status)
  
  return,zst

end 


function maxbcg::input_structdef, ryan_struct, num

  struct = $
    create_struct(ryan_struct[0],        $
                  'photoid',        0LL, $
                  'stripe',          0b, $
                  'ra',              0d, $
                  'dec',             0d, $
                  'gflux',      -9999.0, $
                  'rflux',      -9999.0, $
                  'iflux',      -9999.0, $
                  'photoz',     -9999.0, $
                  'photoz_ivar',    0.0, $
                  'maskflags',        0)

  frontags = $
    ['photoid', $
     'run','rerun','camcol','field','id', $
     'stripe', $
     'ra', 'dec', $
     'clambda','ceta',$
     'photoz', 'photoz_ivar', $
     'cmodel_counts', $
     'counts_model', $
     'gflux', 'rflux', 'iflux', $
     'modelflux', $
     'gmr','dgmr','rmi','drmi', $
     'objc_prob_gal', $
     'flags', 'flags2', 'objc_flags','objc_flags2', $
     'corrselect_flags', $
     'maskflags', $
     'lrg_flags', 'value_flags']
     
  newstruct = reorder_tags(struct, frontags) 

  if n_elements(num) ne 0 then begin 
      newstruct = replicate(newstruct, num)
  endif 

  return,newstruct
end 


function maxbcg::lup2nmgy, lups, band

  if n_elements(lups) eq 0 OR n_elements(band) eq 0 then begin 
      print,'-Syntax: nmgy = mb->lup2nmgy(lups, band)'
  endif 

  bvalues=[1.4D-10, 0.9D-10, 1.2D-10, 1.8D-10, 7.4D-10]

  nmgy = $
    2.D*bvalues[band]*sinh(-alog(bvalues[band])-0.4D*alog(10.D)*lups)*1.e9
  
  return,nmgy

end 

pro maxbcg::input_stuff, nyu=nyu

  ;; We do this on a stripe-by-stripe basis
  stripes = self->input_stripes()
  dir = self->input_dir()

  nst = n_elements(stripes)
  for i=0L, nst-1 do begin 

      print,'----------------------------------------------------------------'

      stripe = stripes[i]

      stripeString = stripe2string(stripe)

      ;; These are Ryan's files made from the primary area
      inStruct = self->ryan_get(stripe)
      nInStruct = n_elements(inStruct)

      help,inStruct


      ;; spatial masks
      maskFlags = intarr(nInStruct)

      print,'Applying basic mask'
      apply_pixel_mask, inStruct.clambda, inStruct.ceta, $
        basicMasked, basicUnmasked, /basic
      help, basicUnmasked

      print,'Applying simple mask'
      apply_pixel_mask, inStruct.clambda, inStruct.ceta, $
        simpleMasked, simpleUnmasked, /simple
      help,simpleUnmasked

      print,'Applying bound mask'
      apply_pixel_mask, inStruct.clambda, inStruct.ceta, $
        boundMasked, boundUnmasked, /bound
      help,boundUnmasked

      if basicMasked[0] ne -1 then begin 
          maskFlags[basicMasked] = $
            maskFlags[basicMasked] + !FLAGS_MASKED_BASIC
      endif 
      if simpleMasked[0] ne -1 then begin 
          maskFlags[simpleMasked] = $
            maskFlags[simpleMasked] + !FLAGS_MASKED_SIMPLE
      endif 
      if boundMasked[0] ne -1 then begin 
          maskFlags[boundMasked] = $
            maskFlags[boundMasked] + !FLAGS_MASKED_BOUND
      endif 

      ;; There is a cmodel_counts[3] < 21.3 cut.  We must apply
      ;; the colors cuts, mask cuts, and sg-sep cuts
      keep = where(instruct.gmr LT 2.0 AND $
                   instruct.rmi LT 1.0 AND $
                   (maskFlags AND !FLAGS_MASKED_BOUND) eq 0 AND $
                   inStruct.objc_prob_gal GT 0.8, ngood)

      help,keep

      inStruct = inStruct[keep]
      maskFlags = maskFlags[keep]
      nInStruct = n_elements(inStruct)


      ;; output struct
      outStruct = self->input_structdef(inStruct[0], nInStruct)

      copy_struct, inStruct, outStruct


      ;; Calculate the g and r fluxes from i and the colors
      imag = outStruct.cmodel_counts[3]
      rmag = imag + outStruct.rmi
      gmag = rmag + outStruct.gmr

      outStruct.gflux = self->lup2nmgy(gmag, 1)
      outStruct.rflux = self->lup2nmgy(rmag, 2)
      outStruct.iflux = self->lup2nmgy(imag, 3)

      outStruct.photoid = sdss_photoid(inStruct)
      outStruct.stripe = stripe      
      outStruct.maskFlags = maskFlags

      ;; Add ra/dec for 
      csurvey2eq, outStruct.clambda, outStruct.ceta, ra, dec
      outStruct.ra = ra
      outStruct.dec = dec

      print
      print,'Matching to photoz'
      rmd_run = rem_dup(inStruct.run)
      runs = inStruct[rmd_run].run
      reruns = inStruct[rmd_run].rerun
      nruns = n_elements(runs)
      for irun=0L, nruns-1 do begin 

          run = runs[irun]
          rerun = reruns[irun]

          zst = self->getrunphotoz(run, rerun, nyu=nyu, status=status)
          if status ne 0 then message,'error'

          w = where(zst.photoz_zwarning eq 0, nw)
          if nw ne 0 then begin 
              match, outStruct.photoid, zst[w].photoid, mout, mzst, /sort

              if mout[0] ne -1 then begin 
                  nmatch = n_elements(mout)
                  outStruct[mout].photoz = zst[w[mzst]].photoz_z
                  outStruct[mout].photoz_ivar = $
                    1.0/zst[w[mzst]].photoz_zerr^2
              endif else nmatch = 0
              print,'Matched '+ntostr(nmatch)
          endif 
      endfor 

      self->postgres::struct2table, $
        outStruct, 'maxbcg_input', primary_key='photoid', $
        connect_info='user=postgres', status=stuff_status, $
        tmpdir=self->input_dir()
      if stuff_status ne 0 then begin 
          message,'Failed to stuff file'
      endif 

  endfor 

  query = 'GRANT SELECT on maxbcg_input TO sdss'
  print,query
  self->postgres::query, query, conn='user=postgres', status=status
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Could not GRANT ON maxbcg_input'
  endif 

  self->input_create_index

end 

pro maxbcg::input_create_index

  table = 'maxbcg_input'

  ;; multicolumn index
  query = $
    'CREATE INDEX maxbcg_input_rrcfi_index ON '+table+$
    ' (run,rerun,camcol,field,id)'
  print,query
  self->postgres::query, query, status=status, conn='user=postgres'
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Could not create rrcfi index'
  endif 

  index_cols = ['stripe','photoz','photoz_ivar']
  self->postgres::create_index, table, index_cols, conn='user=postgres'

  query = 'ANALYZE '+table
  print,query
  self->postgres::query,query, conn='user=postgres', status=status
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Could not ANALYZE '+table
  endif 




end 

;; Create the meta table
pro maxbcg::input_meta_stuff

  table = 'maxbcg_input_meta'

  stripes = self->input_stripes()

  query = 'SELECT count(*) FROM maxbcg_input'
  print,query
  st = self->postgres::query(query, conn='user=postgres')

  input_struct = $
    {nrows: st.count, $
     stripes: stripes}

  self->struct2table, input_struct, table, conn='user=postgres', status=status
  if status ne 0 then begin 
      message,'Failed to create table: '+table
  endif 


end 











function maxbcg::corrselect_bit, name
  CASE strlowcase(name) OF
      'good_u': return, 0
      'good_g': return, 1
      'good_r': return, 2
      'good_i': return, 3
      'good_z': return, 4
      ELSE: message,'Band corrselect name: '+ntostr(name)
  endcase 
end 
function maxbcg::corrselect_flag, name
  bit = self->corrselect_bit(name)
  return, 2^bit
end 

function maxbcg::input_maxmag
  return, 22.5
end 

pro maxbcg::input_stuff_old, nyu=nyu

  ;; We do this on a stripe-by-stripe basis
  stripes = self->input_stripes()
  dir = self->input_dir()

  maxmag = self->input_maxmag()
  print
  print,'maxmag (g,r,i) = ',maxmag
  nst = n_elements(stripes)
  for i=0L, nst-1 do begin 

      print,'----------------------------------------------------------------'

      stripe = stripes[i]

      stripeString = stripe2string(stripe)

      ;; These are Ryan's files made from the primary area
      inStruct = self->ryan_get(stripe)
      nInStruct = n_elements(inStruct)

      help,inStruct


      ;; spatial masks
      maskFlags = intarr(nInStruct)

      print,'Applying basic mask'
      apply_pixel_mask, inStruct.clambda, inStruct.ceta, $
        basicMasked, basicUnmasked, /basic
      help, basicUnmasked

      print,'Applying simple mask'
      apply_pixel_mask, inStruct.clambda, inStruct.ceta, $
        simpleMasked, simpleUnmasked, /simple
      help,simpleUnmasked

      print,'Applying bound mask'
      apply_pixel_mask, inStruct.clambda, inStruct.ceta, $
        boundMasked, boundUnmasked, /bound
      help,boundUnmasked

      if basicMasked[0] ne -1 then begin 
          maskFlags[basicMasked] = $
            maskFlags[basicMasked] + !FLAGS_MASKED_BASIC
      endif 
      if simpleMasked[0] ne -1 then begin 
          maskFlags[simpleMasked] = $
            maskFlags[simpleMasked] + !FLAGS_MASKED_SIMPLE
      endif 
      if boundMasked[0] ne -1 then begin 
          maskFlags[boundMasked] = $
            maskFlags[boundMasked] + !FLAGS_MASKED_BOUND
      endif 




      ;; Cut bad measurements.
      gflag = self->corrselect_flag('good_g')
      rflag = self->corrselect_flag('good_r')
      iflag = self->corrselect_flag('good_i')

      ;; some things that we don't want pass corrselect flags since
      ;; I always keep targets.  So apply the mag cut.
      keep = where(instruct.cmodel_counts[1] LT maxmag AND $
                   instruct.cmodel_counts[2] LT maxmag AND $
                   instruct.cmodel_counts[3] LT maxmag AND $
                   (inStruct.corrselect_flags AND gflag) ne 0 AND $
                   (inStruct.corrselect_flags AND rflag) ne 0 AND $
                   (inStruct.corrselect_flags AND iflag) ne 0 AND $
                   (maskFlags AND !FLAGS_MASKED_BOUND)   eq 0 AND $
                   inStruct.objc_prob_gal GT 0.8, ngood)

      help,keep

      inStruct = inStruct[keep]
      maskFlags = maskFlags[keep]
      nInStruct = n_elements(inStruct)



      outStruct = self->input_structdef(inStruct[0])
      outStruct = replicate(outStruct, nInStruct)
      
      copy_struct, inStruct, outStruct

      outStruct.photoid = sdss_photoid(inStruct)
      outStruct.stripe = stripe      
      outStruct.maskFlags = maskFlags



      print
      print,'Matching to photoz'
      rmd_run = rem_dup(inStruct.run)
      runs = inStruct[rmd_run].run
      reruns = inStruct[rmd_run].rerun
      nruns = n_elements(runs)
      for irun=0L, nruns-1 do begin 

          run = runs[irun]
          rerun = reruns[irun]

          zst = self->getrunphotoz(run, rerun, nyu=nyu, status=status)
          if status ne 0 then message,'error'

          w = where(zst.photoz_zwarning eq 0, nw)
          if nw ne 0 then begin 
              match, outStruct.photoid, zst[w].photoid, mout, mzst, /sort

              if mout[0] ne -1 then begin 
                  nmatch = n_elements(mout)
                  outStruct[mout].photoz_z = zst[w[mzst]].photoz_z
                  outStruct[mout].photoz_zerr = zst[w[mzst]].photoz_zerr
              endif else nmatch = 0
              print,'Matched '+ntostr(nmatch)
          endif 
      endfor 

      self->postgres::struct2table, $
        outStruct, 'maxbcg_input', primary_key='photoid', $
        connect_info='user=postgres', status=stuff_status, $
        tmpdir=self->input_dir()
      if stuff_status ne 0 then begin 
          message,'Failed to stuff file'
      endif 

  endfor 

  query = 'GRANT SELECT on maxbcg_input TO sdss'
  print,query
  self->postgres::query,query, conn='user=postgres', status=status
  if status ne self->postgres::status_val('no_result') then begin 
      message,'Could not GRANT ON maxbcg_input'
  endif 

  self->input_create_index

end 








;; Ben's input files
pro maxbcg::input_create

  stripes = self->input_stripes()
  dir = self->input_dir()

  pg=obj_new('postgres')
  nst = n_elements(stripes)
  for i=0L, nst-1 do begin 

      stripe = stripes[i]
      outfile = self->input_file(stripe)
      fitfile = repstr(outfile, '.st','.fit')
      sstr = ntostr(stripe)
      print,'-------------------------------------------------------------'

      query = 'select * from maxbcg_input where stripe = '+sstr

      ;; Need to index stripe before running this
      ;; everything but petrocounts in the current database
      query = $
        'SELECT '+$
        '  photoid, run, rerun, camcol, field, id, '+$
        '  stripe, '+$
        '  clambda, ceta, '+$
        '  cmodel_counts, '+$
        '  objc_prob_gal, '+$
        '  gmr, dgmr as gmr_err, rmi, drmi as rmi_err, '+$
        '  photoz_z, photoz_zerr, '+$
        '  corrselect_flags, lrg_flags, value_flags, maskflags '+$
        'FROM '+$
        '  maxbcg_input '+$
        'WHERE '+$
        '  stripe = '+sstr

      print,query
      outstruct = pg->query(query, status=status)


      if status ne 0 then message,'Query failed'

      w = where((outstruct.maskflags AND !FLAGS_MASKED_BASIC) eq 0, $
                n_basic)
      w = where((outstruct.maskflags AND !FLAGS_MASKED_SIMPLE) eq 0, $
                n_simple)
      w = where((outstruct.maskflags AND !FLAGS_MASKED_BOUND) eq 0, $
                n_bound)
      w = where((outstruct.maskflags AND !FLAGS_MASKED_COMBINED) eq 0, $
                n_combined)

      hdr = {stripe: stripe, n_basic: n_basic, n_simple:n_simple, $
             n_bound: n_bound, n_combined: n_combined}
      print_struct, hdr
      print
      print,'Writing output file: ',outFile
      write_idlstruct, outStruct, outFile, hdr=hdr
      print,'Writing fits file: ',fitFile
      mwrfits, outStruct, fitFile, /create

  endfor 

end 



pro maxbcg::get_input_maskcut, cat, keep, $
                              basic=basic, simple=simple, $
                              bound=bound, combined=combined

  if keyword_set(basic) then begin 
      keep=where( (cat.maskFlags AND !FLAGS_MASKED_BASIC) eq 0, nw)
  endif else if keyword_set(simple) then begin 
      keep=where( (cat.maskFlags AND !FLAGS_MASKED_SIMPLE) eq 0, nw)      
  endif else if keyword_set(bound) then begin 
      keep=where( (cat.maskFlags AND !FLAGS_MASKED_BOUND) eq 0, nw)
  endif else if keyword_set(combined) then begin 
      keep=where( (cat.maskFlags AND !FLAGS_MASKED_COMBINED) eq 0, nw)
  endif  

end 

pro maxbcg::get_input_nobj, hdr, nObj, $
                           basic=basic, simple=simple, $
                           bound=bound, combined=combined
  
  if keyword_set(basic) then begin 
      nObj = hdr.n_basic
  endif else if keyword_set(simple) then begin 
      nObj = hdr.n_simple
  endif else if keyword_set(bound) then begin 
      nObj = hdr.n_bound
  endif else if keyword_set(combined) then begin 
      nObj = hdr.n_combine
  endif else begin
      nObj = hdr.nrows
  endelse 

end 

pro maxbcg::get_input_files, files, stripes

  stripePosition = 19
  stripeLength = 2

  ;; Standard input directory
  inDir = concat_dir(esheldon_config('lensinput_dir'), 'maxbcg_input')
  
  pattern = concat_dir(indir, 'maxbcg_input_stripe[0-9][0-9].st')

  ;; get the files that match
  if !version.release GE 5.5 then begin 
      files = file_search(pattern,count=nf)
  endif else begin 
      files = findfile(pattern,count=nf)
  endelse 

  for i=0L, nf-1 do begin 
      dirsep, files[i], td, tf
      stripe = long( strmid(tf, stripePosition, stripeLength) )
      
      add_arrval, stripe, stripes
  endfor 

end 

; convert to fits for Ben
pro maxbcg::input_convert

  self->get_input_files, files, stripes

  fitfiles = repstr(files, '.st', '.fit')


  idlstruct2fits, files, fitfiles

end 
function maxbcg::input_get, $
               stripes=stripes, $
               basic=basic, simple=simple, $
               bound=bound, combined=combined, $
               columns=columns, $
               silent=silent, $
               count=count, $
               status=status, $
               help=help

  status = 1

  if keyword_set(help) then begin 
      print,'-Syntax: '
      print,'  cat = obj->input_get(stripes=, '
      print,'                       /basic, /simple, /bound, /combined, '
      print,'                       columns=, '
      print,'                       /silent, '
      print,'                       count=, '
      print,'                       status=, '
      print,'                       /help)'
      print
      print,'Default is to read all objects from all stripes'
      return, -1
  endif 

  count = 0

  nStripes = n_elements(stripes)
  if nStripes eq 0 then begin 
      stripes = self->input_stripes(nStripes=nStripes)
  endif 

  self->get_input_files, files, fileStripes

  ;; Now match to available stripes
  match, stripes, fileStripes, mget, mfile, /sort

  if n_elements(mget) ne nStripes then begin 
      if mget[0] eq -1 then begin 
          print
          print,'These stripes were not found: '
          print,stripes
          return, -1
      endif else begin 
          remove, mget, stripes
          print
          print,'These stripes were not found: '
          print,stripes
          return, -1
      endelse 
  endif else begin 
      files = files[mfile]
      fileStripes = fileStripes[mfile]
  endelse 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read the headers to find out how many objects there are
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ntotal = 0L
  numlist = lonarr(nStripes)
  for i=0L, nStripes-1 do begin 
      hdr = read_idlheader(files[i])

      self->get_input_nobj, hdr, nObj, $
        basic=basic, simple=simple, $
        bound=bound, combined=combined

      numlist[i] = nObj
      nTotal = ntotal + numlist[i]
  endfor
 
  if nTotal eq 0 then begin
      print
      print,'No objects found!'
      count=0
      return, -1
  endif 
  count = nTotal

  if not keyword_set(silent) then begin
      print
      print,'Total number of objects: '+ntostr(nTotal)
  endif 

  delvarx, cat
  beg = 0L
  for i=0L, nStripes-1 do begin 


      if numlist[i] ne 0 then begin 

          file = files[i]

          if not keyword_set(silent) then begin 
              print
              print,'Reading file: ',file
          endif 

          t = read_idlstruct(file, silent=silent, columns=columns)
          nObj = n_elements(t)

          if nObj GT numlist[i] then begin 
              self->get_input_maskcut, t, keep, $
                basic=basic, simple=simple, $
                bound=bound, combined=combined
              
              if n_elements(keep) ne numlist[i] then begin 
                  print
                  print,'Mismatch between maskcut and numlist'
                  return, -1
              endif 
              t = t[keep]
          endif 

          if n_elements(cat) eq 0 then begin 
              cat=replicate(t[0], nTotal)
              cat[beg:beg+numlist[i]-1] = t
              beg = beg+numlist[i]
          endif else begin 
              cat[beg:beg+numlist[i]-1] = t
              beg = beg+numlist[i]
          endelse 

          ;; Clear these variables each time
          delvarx, t 
          
      endif else begin  
          print,'File is empty: '+ntostr(file)
      endelse 


  endfor 
  
  status = 0
  return,cat

end 











;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generate randoms covering the same space as the input catalog
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function maxbcg::numrand
  return,24
end 
function maxbcg::randnum
  numrand = self->numrand()
  return,lindgen(numrand)
end 
function maxbcg::nperfile
  return,1000000
end 


;; Generate random lenses.  The points are generated from the 
;; BOUND mask for dr4plus.  The edges are checked versus the
;; basic mask for the source sample we are using; that step comes
;; during the setuplens stage for the randoms.


function maxbcg::genrand, nperfile=nperfile

  self->genrand, rlam, reta, nperfile=nperfile

  ;; Output struct
  print
  print,'Copying to random struct'
  ss = create_struct('clambda', 0d, $
                     'ceta', 0d)
  outstruct = replicate(ss, nperfile)
      
  outstruct.clambda = temporary(rlam)
  outstruct.ceta    = temporary(reta)

  return,outstruct
end 

pro maxbcg::genrand, rlam, reta, nperfile=nperfile, rra=rra, rdec=rdec

    if n_elements(nperfile) eq 0 then nperfile = self->nperfile()

    cat = self->catalog()
    if strmatch(cat, 'hv1*') then cat='hv1'
    case cat of 
        'dr406': begin
            stripes = self->input_stripes()

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; Bound mask has missing fields/pieces plus bright stars
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

            maskfile = self->maxbcg_maskfile()
            print,'Generating from maskfile: ',maskfile
            sdss_genrand, stripes, nperfile, rlam, reta, $
                phot_maskfile=maskfile

        end
        'hv1': begin
            ; For the sims this is just a quadrant
            ;   0 < ra < 90
            ;   0 < dec < 90

            minra = 0d
            maxra = 90d
            mindec = 0d
            maxdec = 90d

            ; generate uniformly in sin(dec)
            sinfdec = sin( minra*!d2r )
            sinldec = sin( maxra*!d2r )

            rsindec = arrscl( randomu(seed, nperfile, /double), $
                sinfdec, sinldec, arrmin=0d, arrmax=1d )
            rdec = asin(rsindec)*!r2d
            rra = arrscl( randomu(seed, nperfile, /double), $
                minra, maxra, arrmin=0d, arrmax=1d )
            eq2csurvey, rra, rdec, rlam, reta
        end
		'gmbcg10': begin
			; eventually we will use the dr7 area, but for now just use
			; the old area

            stripes = self->input_stripes()

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; Bound mask has missing fields/pieces plus bright stars
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

            maskfile = self->maxbcg_maskfile()
            print,'Generating from maskfile: ',maskfile
            sdss_genrand, stripes, nperfile, rlam, reta, $
                phot_maskfile=maskfile


		end
        else: message,'unknown catalog: '+self->catalog()
    endcase

end 

function maxbcg::random_dir
  input_dir = self->input_dir()
  random_dir = path_join(input_dir, 'random')
  ;random_dir = concat_dir(input_dir, 'random')
  return, random_dir
end 
function maxbcg::random_file, randnum

  dir = self->random_dir()
  file = 'maxbcg_random'+strn(randnum, len=2, padchar='0')+'.st'
  file = concat_dir(dir, file)
  return,file
end 
function maxbcg::random_read, randnum
  files = self->random_file(randnum)
  print,'Reading file: ',files
  struct = read_idlstruct_multi(files)
  return, struct
end 

pro maxbcg::write_randoms

	dir = self->random_dir()
	if not file_test(dir) then file_mkdir, dir

	randnums = self->randnum()
	nrand = n_elements(randnums)

	for i=0L, nrand-1 do begin 

		print,'----------------------------------------------------------'
		rnum = randnums[i]
		file = self->random_file(rnum)
		st = self->genrand()

		print
		print,'Writing to file: ',file
		write_idlstruct, st, file
	endfor 

end 






























pro maxbcg::zbin, bcg, zbin, s, low, high, num

  s = sort(bcg.z)
  z = bcg[s].z

  nbcg = n_elements(bcg)

  uz = rem_dup(z)
  nuz = n_elements(uz)
  
  low = lonarr(nuz)
  high = lonarr(nuz)
  num = lonarr(nuz)
  zbin = fltarr(nuz)
  
  zold = z[0]
  znew = z[0]

  j = 0L
  for i=0L, nuz-1 do begin 

      low[i] = j
      zbin[i] = znew
      WHILE znew eq zold do begin 
          high[i] = j

          j = j+1
          if (j GT nbcg-1) then break

          znew = z[j]

      endwhile 
      zold = znew

  endfor 

  num = lonarr(nuz)
  for i=0L, nuz-1 do num[i] = high[i]-low[i]+1


end 







pro maxbcg::plothist_ngals200, bcg=bcg, dops=dops

  if n_elements(bcg) eq 0 then bcg = self->get()


  plotdir = '~/plots/maxbcg/'
  psfile = plotdir + $
    'maxbcg_'+self->catalog()+'_ngals200_func'
  psfile = psfile + '.eps'

  if keyword_set(dops) then begin 
      begplot, name=psfile, /encapsulated
  endif 




  yrange = [0.9, 2.e5]
  w = where(bcg.ngals200 GE 3)

  xrange = [0.5*min(bcg[w].ngals200), max(bcg[w].ngals200)*1.5]

  xtitle = 'N!D200!N'
  ytitle = 'N'

  plothist, bcg[w].ngals200, $
    /xlog, xrange=xrange, xstyle=3, $
    /ylog, yrange=yrange, ystyle=3, $
    xtitle=xtitle, ytitle=ytitle, $
    xticklen=0.03, $
    ytickf='loglabels', aspect=!gratio

  if keyword_set(dops) then endplot, /trim_bbox

end 



pro maxbcg::plothist_ilum200_cuts, st, ylog=ylog

  loglum = 10.0 + alog10(st.ilum200)

  bin = 0.1

  if n_elements(ylog) eq 0 then ylog = 1
  plothist, loglum, bin=bin, $
    ylog=ylog,yrange=[0.9,5.e5],ysty=3, $
    ytitle='Number', xtitle='log( ilum200 )'

  simpctable, colorlist=colors

  ncut = 10
  for ngalscut = 1,10 do begin 

      w=where(st.ngals200 GT ngalscut)

      plothist, loglum[w], bin=bin, /overplot, color=colors[ngalscut]

      add_arrval, strn(ngalscut, len=2), cuts

  endfor 

  cuts = 'ngals200 > '+cuts
  legend, cuts, /right, colors=colors[lindgen(ncut)+1], line=0, box=0

end 


;; plot the ilum200 histogram in ngals 200 bins

pro maxbcg::plothist_ngals200_ilum200, st, set, $
          fill=fill, arrfill=arrfill, specz=specz, dops=dops, $
          xrange=xrange, yrange=yrange, nolegend=nolegend

;  on_error, 2
  if n_params() LT 2 then begin 
      print,'-Syntax: mb->plothist_ngals200_ilum200, struct, set, /dops'
  endif 

  if keyword_set(dops) then begin 
      setstr = ntostr(set)
      dir = self->lensdir('plot', /base)

      if keyword_set(specz) then addstr='_specz' else addstr=''
      file = $
        'maxbcg_'+self->catalog()+$
        '_ngals200_ilum200'+addstr+'_hist_set'+setstr+'.ps'
      file = concat_dir(dir, file)
      begplot, file, /color, /landscape
  endif 


  lcharsize=1
  if set eq 1 then begin 
      min_ngals = 1
      max_ngals = 9
      lumbin = 0.1
      if n_elements(yrange) eq 0 then yrange = [0.9, 4.e5]
      bin = 1
  endif else if set eq 2 then begin 
      min_ngals = 10
      max_ngals = 40
      lumbin = 0.05  
      if n_elements(yrange) eq 0 then yrange = [0.9, 1.e3]
      bin = 1
  endif else if set eq 3 then begin 
      min_ngals = 41
      max_ngals = 99
      lumbin = 0.03
      if n_elements(yrange) eq 0 then yrange = [0.9, 100]
      bin = 10
  endif else if set eq 4 then begin 
      min_ngals = 100
      max_ngals = max(st.ngals200)
      lumbin = 0.05
      if n_elements(yrange) eq 0 then yrange = [0.9, 10]
      bin = max_ngals-min_ngals
  endif 

  sttags = tag_names(st)
  if keyword_set(specz) then begin 
      ztag = ( where(sttags eq 'BCG_SPEC_Z') )[0]
      ilumtag = ( where(sttags eq 'SPEC_ILUM200') )[0]
      ivartag = ( where(sttags eq 'SPEC_ILUM200_IVAR') )[0]
  endif else begin 
      ztag = ( where(sttags eq 'PHOTOZ_CTS') )[0]
      ilumtag = ( where(sttags eq 'ILUM200') )[0]
      ivartag = ( where(sttags eq 'ILUM200_IVAR') )[0]
  endelse 


  wz = where(st.(ztag) GT 0.1 AND st.(ztag) LT 0.3, nz)

;  wz = lindgen(n_elements(st))
  wgood = where(st[wz].ngals200 GE min_ngals AND $
                st[wz].ngals200 LE max_ngals AND $
                st[wz].(ivartag) GT 0.0 AND $
                st[wz].(ilumtag) GT 0.0)
  wgood = wz[wgood]

  loglum = alog10(st[wgood].(ilumtag)) + 10.0
  minlum = min(loglum, max=maxlum)


  h = histogram(st[wgood].ngals200, $
                min=min_ngals, max=max_ngals, bin=bin, rev=rev)
  wh = where(h ne 0, nh)
  
  simpctable, colorlist=clist

  wcl=where(clist ne !p.color)
  clist = clist[wcl]
  nclist = n_elements(clist)

  pj = 0
  for ii=0L, nh-1 do begin 
      i = wh[ii]

      w = rev[ rev[i]:rev[i+1]-1 ]
      w = wgood[w]
      nw = n_elements(w)


      if bin eq 1 then begin 
          tn = st[w[0]].ngals200
          tnstr = ntostr(tn)
      endif else begin 
          minn = min(st[w].ngals200, max=maxn)
          tnstr = '['+ntostr(minn)+','+ntostr(maxn)+']'
      endelse 

      loglum = 10.0 + alog10(st[w].(ilumtag))

      if nw GE 2 then begin 
          tmedlum = median(st[w].(ilumtag))

          ivar = st[w].(ivartag)

          ivartot = total(ivar)
          tmeanlum = $
            total(st[w].(ilumtag)*ivar)/ivartot
          print,tmedlum,tmeanlum,sqrt(1.0/ivartot)
          
          
          if pj eq 0 then begin 
              clr = !p.color
          endif else begin 
              clr=clist[pj MOD nclist]
          endelse 
          
          add_arrval, clr, pcolors
          add_arrval, tmedlum, medlum
          add_arrval, tnstr, nstr
;          fill = 1
          
          if n_elements(arrfill) ne 0 then begin 
              fill = arrfill[ii]
          endif 

          if n_elements(xrange) eq 0 then begin 
              xrng=[minlum, 1.05*maxlum]
          endif else begin 
              xrng=xrange
          endelse 

          ylog = 1
          if ylog ne 0 then ytickf='loglabels'

          if not keyword_set(fill) then begin 
              if pj MOD 2 eq 0 then thick=2*!p.thick else thick=!p.thick
          endif 

          if pj eq 0 then begin 

              plothist, loglum, bin=lumbin, $
                        ylog=ylog, yrange=yrange, ystyle=3, ytickf=ytickf,$
                        xrange = xrng, xstyle=3, $
                        xtitle='log!D10!N( ilum200 )', ytitle='N', $
                        fill=fill, fcolor=!p.color, thick=thick
          endif else begin 
              plothist, loglum, bin=lumbin, /overplot, color=clr, $
                        fill=fill, fcolor=clr, thick=thick
          endelse 
          
          pj = pj+1
      endif 

  endfor 

  if not keyword_set(nolegend) then begin 
      mess = 'ngals200 = '+nstr ;+$
                                ;' med(lum200) = '+ntostr(medlum,4)
      legend, mess, line=0, colors=pcolors, /right, box=0, charsize=lcharsize
  endif 

  if keyword_set(dops) then endplot, /landfix

end 



pro maxbcg::plothist_ngals200_ilum200_compare_specz, st, set, dops=dops

  nst = n_elements(st) 
  nset = n_elements(set)
  if nst eq 0 OR nset eq 0 then begin 
      print,'-Syntax: mb->plothist_ngals200_ilum200_compare_specz, st, set'
      print
      message,'Halting'
  endif 

  if keyword_set(dops) then begin 
      setstr = ntostr(set)
      dir = self->lensdir('plot', /base)

      file = $
        'maxbcg_'+self->catalog()+$
        '_ngals200_ilum200_hist_compare_specz_set'+setstr+'.ps'
      file = concat_dir(dir, file)
      begplot, file, /color
  endif 



  !p.multi=[0,0,2]

  xrange = [9,14]

  w = where(st.bcg_spec_z GT 0.1 AND st.bcg_spec_z LT 0.3)
  self->plothist_ngals200_ilum200, st[w], set, xrange=xrange
  legend,'maxBCG z',/left,box=0

  self->plothist_ngals200_ilum200, st[w], set,/spec, xrange=xrange,/nolegend
  legend,'spec z',/left,box=0

  if keyword_set(dops) then endplot
  
  !p.multi=0

end 


pro maxbcg::plot_ngals200_ilum200_compare_specz, st, dops=dops

  nst = n_elements(st) 
  if nst eq 0 then begin 
      print,'-Syntax: mb->plot_ngals200_ilum200_compare_specz, st, /dops'
      print
      message,'Halting'
  endif 

  if keyword_set(dops) then begin 
      dir = self->lensdir('plot', /base)

      file = $
        'maxbcg_'+self->catalog()+$
        '_ilum200_vs_ngals200_compare_specz.ps'
      file = concat_dir(dir, file)
      begplot, file, /color
  endif 


  w=where(st.ngals200 GE 1 AND st.ilum200 GT 0)
  self->plot_tag_vs_tag, st[w], 'ngals200', 'ilum200', $
    /xinteger, $
    ytitle = 'i-band Lum [10!U10!N L'+sunsymbol()+']', $
    xtitle = 'ngals!D200!N'

  w=where(st.ngals200 GE 1 AND st.bcg_spec_z GT 0.1 AND st.bcg_spec_z LT 0.3)
  self->plot_tag_vs_tag, st[w], 'ngals200', 'spec_ilum200', /xinteger, $
    /overplot,color=!darkGreen

  legend, ['maxBCG z', 'spec z'], line=[0,0], color=[!p.color, !darkGreen],$
    /left, box=0

  if keyword_set(dops) then endplot

end 













;; plot the color histogram in ngals 200 bins
;; cstring should be gr or ri
pro maxbcg::plothist_ngals200_color, st, set, cstring, dops=dops

  on_error, 2
  if n_params() LT 3 then begin 
      print,'-Syntax: mb->plothist_ngals200_color, struct, set, cstring, /dops'
  endif 

  if cstring eq 'gr' then begin 
      xtitle = 'median( g-r )'
      if not tag_exist(st, 'med_kgmr200', index=ctag) then begin 
          message,ctag+' does not exist'
      endif 
  endif else if cstring eq 'ri' then begin 
      xtitle = 'median( r-i )'
      if not tag_exist(st, 'med_krmi200', index=ctag) then begin 
          message,ctag+' does not exist'
      endif 
  endif else begin 
      message,"cstring should be 'gr' or 'ri'"
  endelse 


  if keyword_set(dops) then begin 
      setstr = ntostr(set)
      dir = self->lensdir('plot', /base)
      file = $
        'maxbcg_'+self->catalog()+'_ngals200_'+cstring+'_hist_set'+setstr+'.ps'
      file = concat_dir(dir, file)
      begplot, file, /color, /landscape
  endif 



  lcharsize=1
  if set eq 1 then begin 
      min_ngals = 1
      max_ngals = 4
      yrange = [0.9, 1.e6]
      bin = 1
      cbin = 0.1
      fill = 0
  endif else if set eq 2 then begin 
      min_ngals = 5
      max_ngals = 9
      yrange = [0.9, 4.e5]
      bin = 1
      cbin = 0.1
      fill = 0
  endif else if set eq 3 then begin 
      min_ngals = 10
      max_ngals = 40
      yrange = [0.9, 2.e3]
      bin = 1
      cbin = 0.02
      fill = 1
  endif else if set eq 4 then begin 
      min_ngals = 41
      max_ngals = 90
      yrange = [0.9, 200]
      bin = 10
      cbin = 0.02
          fill = 1
  endif else if set eq 5 then begin 
      min_ngals = 91
      max_ngals = max(st.ngals200)
      yrange = [0.9, 15]
      bin = max_ngals-min_ngals
      cbin = 0.012
      fill = 1
  endif 

  wgood = where(st.photoz_cts GE 0.1 AND $
                st.ngals200 GE min_ngals AND $
                st.ngals200 LE max_ngals)
  fluxcolor = st[wgood].(ctag)
  minfluxcolor = min(fluxcolor, max=maxfluxcolor)
  xrange = [0.95*minfluxcolor, 1.1*maxfluxcolor]

  h = histogram(st[wgood].ngals200, $
                min=min_ngals, max=max_ngals, bin=bin, rev=rev)
  wh = where(h ne 0, nh)
  
  simpctable, colorlist=clist

  wcl=where(clist ne !p.color)
  clist = clist[wcl]
  nclist = n_elements(clist)

  ytitle = 'N'
  pj = 0
  for ii=0L, nh-1 do begin 
      i = wh[ii]

      w = rev[ rev[i]:rev[i+1]-1 ]
      w = wgood[w]
      nw = n_elements(w)


      if bin eq 1 then begin 
          tn = st[w[0]].ngals200
          tnstr = ntostr(tn)
          if tn eq 1 then begin 
              fillold = fill
              fill = 1
          endif 
      endif else begin 
          minn = min(st[w].ngals200, max=maxn)
          tnstr = '['+ntostr(minn)+','+ntostr(maxn)+']'
      endelse 

      fluxcolor = st[w].(ctag)

      if nw GE 2 then begin 
          tmedfluxcolor = median(fluxcolor)
          
          if pj eq 0 then begin 
              clr = !p.color
          endif else begin 
              clr=clist[pj MOD nclist]
          endelse 
          
          add_arrval, clr, pcolors
          add_arrval, tmedfluxcolor, medfluxcolor
          add_arrval, tnstr, nstr

;          if bin eq 1 then begin 
;              if tn eq 4 then fill=0 
;          endif  

          if pj eq 0 then begin 
              plothist, fluxcolor, bin=cbin, $
                        /ylog, yrange=yrange, ystyle=3, ytickf='loglabels',$
                        xrange=xrange, xstyle=3, $
                        xtitle=xtitle, ytitle=ytitle, $
                        fill=fill, fcolor=!p.color
          endif else begin 
              plothist, fluxcolor, bin=cbin, /overplot, color=clr, $
                        fill=fill, fcolor=clr
          endelse 
          
          pj = pj+1
      endif 

      if bin eq 1 then begin 
          if tn eq 1 then fill = fillold
      endif 

  endfor 

  mess = 'ngals200 = '+nstr
  legend, mess, line=0, colors=pcolors, /right, box=0, charsize=lcharsize

  if keyword_set(dops) then endplot, /landfix

end 


;; plot the color histogram in ngals 200 bins
;; cstring should be gr or ri
pro maxbcg::plothist_ngals200_totcolor, st, set, cstring, dops=dops

  on_error, 2
  if n_params() LT 3 then begin 
      print,'-Syntax: mb->plothist_ngals200_totcolor, struct, set, cstring, /dops'
  endif 

  if cstring eq 'gr' then begin 
      xtitle = 'total( g-r )'
      b1 = 1
      b2 = 2
  endif else if cstring eq 'ri' then begin 
      xtitle = 'total( r-i )'
      b1 = 2
      b2 = 3
  endif else begin 
      message,"cstring should be 'gr' or 'ri'"
  endelse 


  if keyword_set(dops) then begin 
      setstr = ntostr(set)
      dir = self->lensdir('plot', /base)
      file = $
        'maxbcg_'+self->catalog()+'_ngals200_total_'+cstring+'_hist_set'+setstr+'.ps'
      file = concat_dir(dir, file)
      begplot, file, /color, /landscape
  endif 



  lcharsize=1
  if set eq 1 then begin 
      min_ngals = 1
      max_ngals = 4
      yrange = [0.9, 1.e6]
      bin = 1
      cbin = 0.1
      fill = 0
  endif else if set eq 2 then begin 
      min_ngals = 5
      max_ngals = 9
      yrange = [0.9, 4.e5]
      bin = 1
      cbin = 0.1
      fill = 0
  endif else if set eq 3 then begin 
      min_ngals = 10
      max_ngals = 40
      yrange = [0.9, 2.e3]
      bin = 1
      cbin = 0.1
      fill = 1
  endif else if set eq 4 then begin 
      min_ngals = 41
      max_ngals = 90
      yrange = [0.9, 200]
      bin = 10
      cbin = 0.02
          fill = 1
  endif else if set eq 5 then begin 
      min_ngals = 91
      max_ngals = max(st.ngals200)
      yrange = [0.9, 15]
      bin = max_ngals-min_ngals
      cbin = 0.1
      fill = 1
  endif 

  wgood = where(st.photoz_cts GE 0.1 AND $
                st.ngals200 GE min_ngals AND $
                st.ngals200 LE max_ngals AND $
                st.lumsolar200_ivar[b1] GT 0 AND $
                st.lumsolar200_ivar[b2] GT 0)

  band_shift = self->band_shift()
  absmag200_1 = self->lumsolar2absmag(st.lumsolar200[b1],filter=b1, $
                                      band_shift=band_shift)
  absmag200_2 = self->lumsolar2absmag(st.lumsolar200[b2],filter=b2, $
                                      band_shift=band_shift)

  fluxcolor = absmag200_1[wgood] - absmag200_2[wgood]
  minfluxcolor = min(fluxcolor, max=maxfluxcolor)
  xrange = [0.95*minfluxcolor, 1.1*maxfluxcolor]

  h = histogram(st[wgood].ngals200, $
                min=min_ngals, max=max_ngals, bin=bin, rev=rev)
  wh = where(h ne 0, nh)
  
  simpctable, colorlist=clist

  wcl=where(clist ne !p.color)
  clist = clist[wcl]
  nclist = n_elements(clist)

  ytitle = 'N'
  pj = 0
  for ii=0L, nh-1 do begin 
      i = wh[ii]

      w = rev[ rev[i]:rev[i+1]-1 ]
      w = wgood[w]
      nw = n_elements(w)


      if bin eq 1 then begin 
          tn = st[w[0]].ngals200
          tnstr = ntostr(tn)
          if tn eq 1 then begin 
              fillold = fill
              fill = 1
          endif 
      endif else begin 
          minn = min(st[w].ngals200, max=maxn)
          tnstr = '['+ntostr(minn)+','+ntostr(maxn)+']'
      endelse 

      fluxcolor = absmag200_1[w] - absmag200_2[w]      

      if nw GE 2 then begin 
          
          if pj eq 0 then begin 
              clr = !p.color
          endif else begin 
              clr=clist[pj MOD nclist]
          endelse 
          
          add_arrval, clr, pcolors
          add_arrval, tnstr, nstr

;          if bin eq 1 then begin 
;              if tn eq 4 then fill=0 
;          endif  

          if pj eq 0 then begin 
              plothist, fluxcolor, bin=cbin, $
                        /ylog, yrange=yrange, ystyle=3, ytickf='loglabels',$
                        xrange=xrange, xstyle=3, $
                        xtitle=xtitle, ytitle=ytitle, $
                        fill=fill, fcolor=!p.color
          endif else begin 
              plothist, fluxcolor, bin=cbin, /overplot, color=clr, $
                        fill=fill, fcolor=clr
          endelse 
          
          pj = pj+1
      endif 

      if bin eq 1 then begin 
          if tn eq 1 then fill = fillold
      endif 

  endfor 

  mess = 'ngals200 = '+nstr
  legend, mess, line=0, colors=pcolors, /right, box=0, charsize=lcharsize

  if keyword_set(dops) then endplot, /landfix

end 









pro maxbcg::plot_ngals200_vs_ngals, st, dops=dops

  if keyword_set(dops) then begin 
      dir = self->lensdir('plot', /base)
      file = $
        'maxbcg_'+self->catalog()+'_ngals200_vs_ngals_scatterplot.ps'
      file = concat_dir(dir, file)
      begplot, file
  endif 

  n_n200 = ntostr(st.ngals)+'-'+ntostr(st.ngals200)
  rmd = rem_dup(n_n200)

  xtitle = 'Ngals'
  ytitle = 'Ngals200'
  pplot, st[rmd].ngals, st[rmd].ngals200, $
    /xlog, /ylog, xrange=[0.4, 300], yrange=[0.4, 300], xstyle=3, ystyle=3,$
    xtitle=xtitle, ytitle=ytitle, aspect=1, psym=8, symsize=0.2
;  oplot, [0.01,1000],[0.01,1000]

  if keyword_set(dops) then endplot
  
end 








; The cut on ngals200 must work for the other variables too
function maxbcg::_convert_tagname, tagname

  p1 = stregex(tagname, '\[')
  if p1 ne -1 then begin 
      p2 = stregex(tagname, '\]')
      el = strmid(tagname,p1+1, p2-p1-1)
      tname = strmid(tagname, 0, p1)
      tname = tname + '_'+el
      return,tname
  endif else begin 
      return,tagname
  endelse 

end 
function maxbcg::tag_vs_tag, st, tagname1, tagname2, ind=ind, nperbin=nperbin, xinteger=xinteger

  on_error, 2
  if n_params() LT 3 then begin 
      message,'-Syntax: comp = mb->tag_vs_tag(struct, tagname1, tagname2, nperbin=, ind, /xinteger)'
  endif 

  if n_elements(ind) ne 0 then begin   
      cmd1 = 'v1 = st[ind].'+tagname1
      cmd2 = 'v2 = st[ind].'+tagname2
  endif else begin 
      cmd1 = 'v1 = st.'+tagname1
      cmd2 = 'v2 = st.'+tagname2
  endelse 

  if not execute(cmd1) then begin 
      message,'Could not get tag: '+tagname1
  endif 
  if not execute(cmd2) then begin 
      message,'Could not get tag: '+tagname2
  endif 


  if keyword_set(xinteger) then begin 
      bs = binner(v1, v2, binsize=1, sigma_clip=1)
  endif else begin 
      if n_elements(nperbin) eq 0 then nperbin=1000
      bs = binner(v1, v2, nperbin=nperbin, /sigma_clip)
  endelse 
  min_tag1 = min(v1, max=max_tag1)
  min_tag2 = min(v2, max=max_tag2)

  outst = create_struct('tag1_range', [min_tag1, max_tag1], $
                        'tag2_range', [min_tag2, max_tag2], $
                        'tag1_mean', bs.xbinned, $
                        'tag2_mean', bs.ybinned, $
                        'tag2_err',  bs.ybinned_err, $
                        'tag2_sdev', bs.ybinned_sdev)

  return,outst

end 
pro maxbcg::plot_tag_vs_tag, st, tagname1, tagname2, dops=dops, $
          xinteger=xinteger, $
          nperbin=nperbin, $
          ind=ind, $
          error=error, scatter=scatter, $
          xlog=xlog, ylog=ylog, $
          xrange=xrange, yrange=yrange, $
          xtitle=xtitle, ytitle=ytitle, comp=comp, $
          overplot=overplot, $
          _extra=_extra

;  on_error, 2
  nst = n_elements(st)
  nt1 = n_elements(tagname1)
  nt2 = n_elements(tagname2)
  ntot = (nst GT 0) + nt1 + nt2
  if ntot LT 3 then begin 
      print,'-Syntax: mb->plot_tag_vs_tag, bcg, xtag, ytag, '+$
        '/xinteger, /scatter, /dops, '+$
        'xlog=, ylog=, xtitle=, ytitle=, comp='
      print
      message,'Halting'
  endif 


  comp = self->tag_vs_tag(st, tagname1, tagname2, $
                          xinteger=xinteger, nperbin=nperbin, ind=ind)
  comptags = tag_names(comp)



  ;; xrange sent?
  if n_elements(xrange) eq 0 then begin 
      xrange = [ 0.5*min(comp.tag1_mean), max(comp.tag1_mean)*1.5 ]
  endif 

  ;; deal with errors/scatter
  if keyword_set(error) then begin 
      yerror = comp.tag2_err
      cyrange = [ 0.5*min( comp.tag2_mean  - yerror ), $
                  max( comp.tag2_mean+yerror )*1.5 ]
      addstr = '_error'
  endif else if keyword_set(scatter) then begin 
      yerror = comp.tag2_sdev
      cyrange = [ 0.5*min( comp.tag2_mean  - yerror ), $
                  max( comp.tag2_mean+yerror )*1.5 ]
      addstr = '_scatter'
  endif else begin 
      cyrange = [ 0.5*min(comp.tag2_mean), max(comp.tag2_mean)*1.5 ]
      addstr = ''
  endelse 


  if n_elements(yrange) eq 0 then yrange = cyrange

  if not keyword_set(xtitle) then xtitle = tagname1
  if not keyword_set(ytitle) then ytitle = tagname2
  if n_elements(xlog) eq 0 then xlog=1
  if n_elements(ylog) eq 0 then ylog=1  



  if keyword_set(dops) then begin 

      tname1 = self->_convert_tagname(tagname1)
      tname2 = self->_convert_tagname(tagname2)
      dir = self->lensdir('plot', /base)
      file = 'maxbcg_'+self->catalog()+'_'+$
             strlowcase(tname2)+'_vs_'+strlowcase(tname1)
      file = file + addstr + '.ps'
      file = concat_dir(dir, file)
      print,'File = ',file
      begplot, file, /color
  endif 


  pplot, comp.tag1_mean, comp.tag2_mean, $
    yerror=yerror, $
    aspect = 1, $
    xlog=xlog, ylog=ylog, $
    xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
    xtitle=xtitle, ytitle = ytitle, overplot=overplot, $
    _extra=_extra

;  if keyword_set(xinteger) AND n_elements(yerror) eq 0 then begin 
  if n_elements(yerror) eq 0 then begin 
      oplot, comp.tag1_mean, comp.tag2_mean, $
             psym=8, symsize=1, _extra=_extra
  endif 

  if keyword_set(dops) then endplot

end 

pro maxbcg::photoz_bias, st

  begplot,$
    name='~/plots/maxbcg/maxbcg_'+self->catalog()+'_photoz_bias.ps',$
    /color

  ngalsvals = [1,2,3,4,5,6,7]

  xrange=[0.0,0.4]
  yrange=[0.0,0.4]

  num=n_elements(ngalsvals)
  for i=0L, num-1 do begin 

      erase & multiplot, [1,2], /square

      ngalsval = ngalsvals[i]

      w1 = where(st.ngals200 eq ngalsval AND st.bcg_spec_z GT 0)
      help,w1
      
      comp = self->tag_vs_tag(st[w1], 'bcg_spec_z', 'photoz_cts')
      pplot, $
        st[w1].bcg_spec_z, st[w1].photoz_cts, $
        psym=3,xrange=xrange,yrange=xrange, xstyle=3, ystyle=3, $
        ytitle = 'maxBCG z'
      


      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean,$
        color=!darkgreen
      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean,$
        psym=8, $
        color=!darkgreen
      
      
      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean - comp.tag2_sdev,$
        color=!darkgreen
      
      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean + comp.tag2_sdev,$
        color=!darkgreen
      
      pplot, [0,1],[0,1], /overplot, color=c2i('red')
      
      
      legend,'ngals!D200!N = '+ntostr(ngalsval),/left,box=0, charsize=1
      
      key = prompt_kbrd('hit a key')
      
      
      multiplot
      
      comp = self->tag_vs_tag(st[w1], 'bcg_spec_z', 'photoz_z')
      pplot, $
        st[w1].bcg_spec_z, st[w1].photoz_z, $
        psym=3,xrange=xrange,yrange=xrange, xstyle=3, ystyle=3, $
        xtitle='specz', ytitle = 'NN photoz'
      
      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean,$
        color=!darkgreen
      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean,$
        psym=8,$
        color=!darkgreen
      
      
      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean - comp.tag2_sdev,$
        color=!darkgreen
      
      pplot, /overplot, $
        comp.tag1_mean,$
        comp.tag2_mean + comp.tag2_sdev,$
        color=!darkgreen
      
      pplot, [0,1],[0,1], /overplot, color=c2i('red')
      
      legend,'ngals!D200!N = '+ntostr(ngalsval),/left,box=0, charsize=1

      multiplot,/reset
      key = prompt_kbrd('hit a key')

  endfor 

  endplot

end 


function maxbcg::mb_plotfile, desc

  dir = self->lensdir('plot', /base)  
  file = 'maxbcg_'+self->catalog() + '_'+desc+'.ps'
  file = concat_dir(dir, file)
  return, file
end 

;; For objects with specz, Compare:
;    maxbcg z to photoz of all neihbors
pro maxbcg::compare_zbcg_zneigh, bcg, neigh, bcg_z, zdiff, dops=dops

  if keyword_set(dops) then begin 
      desc = 'compare_zbcg_zneigh'
      psfile = self->mb_plotfile(desc)
      begplot, psfile, /color
  endif 

  if n_elements(bcg) eq 0 then bcg=self->get()
  if n_elements(neigh) eq 0 then neigh=self->get(/neigh)

;  if n_elements(min_ngals) eq 0 then min_ngals = 1

  for min_ngals = 2, 30 do begin 

      w = where(bcg.bcg_spec_z GT 0.0 AND bcg.ngals200 GE min_ngals, nw)

      maxid = max(bcg[w].bcg_id)
      h = histogram( neigh.bcg_id, min=0, max=maxid, rev=rev)
      
      zbcg_ptrlist = ptrarr(nw)
      zdiff_ptrlist = ptrarr(nw)
      
      for i=0L, nw-1 do begin 
          
          bcg_id = bcg[w[i]].bcg_id
          
          if rev[bcg_id] ne rev[bcg_id+1] then begin 
              wneigh = rev[ rev[bcg_id]:rev[bcg_id+1]-1 ]
              
              wgood = where(neigh[wneigh].spec_z GT 0.0, ngood)
              
              if ngood ne 0 then begin 
                  wgood = wneigh[wgood]
                  
                  zbcg_ptrlist[i] = $
                    ptr_new( replicate(bcg[w[i]].bcg_spec_z, ngood), /no_copy )
                  zdiff_ptrlist[i] = $
                    ptr_new( bcg[w[i]].bcg_spec_z - neigh[wgood].spec_z, /no_copy)
                  
              endif 
          endif 
          
      endfor 
      
      zbcg = combine_ptrlist(zbcg_ptrlist)
      zdiff = combine_ptrlist(zdiff_ptrlist)
      
      pplot, zbcg, zdiff, psym=3, $
        yrange=[-0.2,0.2], ystyle=3, xrange=[0.0, 0.3], xstyle=1, $
        xtitle='z(BCG)', $
		ytitle='z(BCG) - z(neigh)', $
        /iso, $
		title=textoidl('ngals200 \geq ')+ntostr(min_ngals)
      
      ntot = n_elements(zdiff)
      good_vdiff = 3000         ;km/s
      good_zdiff = good_vdiff/3.e5
      oplot, [0,1],[-good_zdiff, -good_zdiff], color=c2i('red')
      oplot, [0,1],[good_zdiff, good_zdiff], color=c2i('red')
      
      wfrac = where( zdiff LT -good_zdiff OR zdiff GT good_zdiff, nfrac)
      
      frac = float(nfrac)/ntot
      
      add_arrval, frac, fracbad

      legend, 'fracbad = '+ntostr(frac,4,/round), /right, box=0
      legend, textoidl('\pm')+' '+ntostr(good_vdiff,4)+' km/s', $
        line=0, color=c2i('red'), /left, box=0

      key = prompt_kbrd('hit a key')
      if key eq 'q' then return

      add_arrval, min_ngals, min_ngals_array

  endfor 

  pplot, min_ngals_array, fracbad, aspect=1, $
    xtitle = 'min_ngals', ytitle='fracbad'

  if keyword_set(dops) then begin 
      endplot
  endif 

end 



pro maxbcg::coolplots, st

  if n_elements(st) eq 0 then begin
      st = self->get()
  endif 

  xtitle = 'ngals!D200!N'

  self->plot_zdensity, 9, bcg=st, /dops
  self->plot_zdensity, 18, bcg=st, /dops

  self->plot_tag_vs_tag, $
    st, 'ngals', 'ngals200', $
    /xinteger, /scatter, xrange=[0.4, 300], yrange=[0.4, 300], $
    xtitle='Ngals', ytitle='Ngals200', $
    /dops


  self->plot_tag_vs_tag, st, 'ngals200', 'ilum200', /xinteger, /dops,$
    xtitle=xtitle, ytitle='i-band L!D200!N', /scatter
;  self->plot_tag_vs_tag, st, 'ngals200', 'rlum200', /xinteger, /dops,$
;    xtitle=xtitle, ytitle='r-band L!D200!N', /scatter

;  self->plot_tag_vs_tag, st, 'ngals200', 'med_kgmr200', /xinteger, /dops,$
;    ylog=0, yrange=[1.1,1.6], xtitle=xtitle, ytitle='median(g-r)', $
;    /scatter
;  self->plot_tag_vs_tag, st, 'ngals200', 'med_krmi200', /xinteger, /dops,$
;    ylog=0, yrange=[0.4,0.7], xtitle=xtitle, ytitle='median(r-i)', $
;    /scatter


  ;; Ben's measurements
;  for set=1,4 do self->plothist_ngals200_ilum200, st, set, /dops
;  for set=1,5 do self->plothist_ngals200_color, st, set, 'gr', /dops
;  for set=1,5 do self->plothist_ngals200_color, st, set, 'ri', /dops

  ;; My measurements
  for set=1,4 do self->plothist_ngals200_ilum200, st, set, /dops
;  for set=1,5 do self->plothist_ngals200_totcolor, st, set, 'gr', /dops
;  for set=1,5 do self->plothist_ngals200_totcolor, st, set, 'ri', /dops

  ;; Comparisons with specz
  self->photoz_bias, st

  for set=1,4 do begin 
      self->plothist_ngals200_ilum200_compare_specz, st, set, /dops
  endfor 

  self->plot_ngals200_ilum200_compare_specz, st, /dops

;  self->compare_bcg_ilum_bcg_lumsolar, st, /dops

  self->plothist_ngals200, st, /dops
;  self->plothist_ngals200, st, /paper, /dops

  
  w=where(st.ngals200 GE 10 AND st.ngals200 LE 15)

  ;; Ben's measurements
;  self->plot_tag_vs_tag, $
;    st[w], 'photoz_cts', 'ilum200', $
;    nperbin=500,/err, xlog=0, ylog=0,xrange=[0.04, 0.35], yrange=[15,22], $
;    xtitle='z!Dc!N', ytitle='i-band L!D200!N [10!U10!N Lsun]', $
;    /dops

;  self->plot_tag_vs_tag, $
;    st[w], 'photoz_cts', 'rlum200', $
;    nperbin=500,/err, xlog=0, ylog=0,xrange=[0.04, 0.35], yrange=[10,17], $
;    xtitle='z!Dc!N', ytitle='r-band L!D200!N [10!U10!N Lsun]', $
;    /dops

;  self->plot_tag_vs_tag, $
;    st[w], 'photoz_cts', 'med_kgmr200', $
;    nperbin=500, /err, $
;    xlog=0, ylog=0,xrange=[0.08,0.35], yrange=[1.37,1.45], $
;    xtitle='z!Dc!N', ytitle='g-r', $
;    /dops

;  self->plot_tag_vs_tag, $
;    st[w], 'photoz_cts', 'med_krmi200', $
;    nperbin=500, /err, $
;    xlog=0, ylog=0,xrange=[0.08,0.35], yrange=[0.51,0.55], $
;    xtitle='z!Dc!N', ytitle='r-i', $
;    /dops

  ;; My measurements
  self->plot_tag_vs_tag, $
    st[w], 'photoz_cts', 'ilum200', $
    nperbin=500,/err, xlog=0, ylog=0,xrange=[0.04, 0.35], yrange=[16,28], $
    xtitle='z!Dc!N', ytitle='i-band L!D200!N [10!U10!N Lsun]', $
    /dops

  legend,'10 <= N200 <= 15',/right, box=0

end 




















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Binning.  For each bin type, you need to be able to return a where
;    string for each bin through the where_string method.  This is
;    usually done with a series of pro/funcs.  For example:
;
;      ::ilum200_bins, nbin, lowlim, highlim
;      ::ilum200_where_string(nbin, labels=labels)
;
;    Then add a CASE to 
;      ::subtype_nbin(subtype)
;      where_string(subtype)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




pro maxbcg::ngals_bins, nsub, lowlim, highlim

  CASE self->catalog() OF 
      'dr3': begin 
          maxngals = 343

          CASE nsub OF
              6: begin 
                  lowlim =  [1, 4, 6, 9,  16, 37] 
                  highlim = [3, 5, 8, 15, 36, maxngals] 
              end 
              9: begin 
                  lowlim =  [1, 4, 5, 6, 8, 10, 16, 31, 51] 
                  highlim = [3, 4, 5, 7, 9, 15, 30, 50, maxngals] 
              end 
              9: begin 
                  lowlim =  [1, 4, 5, 6, 8, 10, 16, 31, 51] 
                  highlim = [3, 4, 5, 7, 9, 15, 30, 50, maxngals] 
              end 
              12: begin 
                  lowlim =  [1, 4, 5, 6, 8, 10, 16, 21, 31, 41, 51, 76] 
                  highlim = [3, 4, 5, 7, 9, 15, 20, 30, 40, 50, 75, maxngals] 
              end 
              ELSE: message,'unknown number of bins: '+ntostr(nsub)
          endcase 
      end 
      'dr4plus': begin 
          CASE nsub OF 
              12: begin 
                  lowlim =  [0, 3, 4, 5, 7,  9, 15, 20, 30, 40, 50, 75] 
                  highlim = [3, 4, 5, 7, 9, 15, 20, 30, 40, 50, 75, 261] 
              end 
              ELSE: message,'unknown number of bins: '+ntostr(nsub)
          endcase 
      end 
      'dr406': begin 
          CASE nsub OF 
              12: begin 
                  lowlim =  [0, 3, 4, 5, 7,  9, 15, 20, 30, 40, 50, 75] 
                  highlim = [3, 4, 5, 7, 9, 15, 20, 30, 40, 50, 75, 261] 
              end 
              ELSE: message,'unknown number of bins: '+ntostr(nsub)
          endcase 
      end 
  endcase 

end 

function maxbcg::ngals_where_string

  self->ngals_bins, nbin, lowlim, highlim  
  for i=0L, nbin-1 do begin 

      tstring = $
        'struct.ngals GE '+ntostr(lowlim[i])+' AND '+$
        'struct.ngals LT '+ntostr(highlim[i])

      add_arrval, tstring, where_string

  endfor 
  return,where_string

end 











pro maxbcg::lambda_bins, nbin, lowlim, highlim

	case nbin of
		12: begin 
			; ngals bins
			;lowlim =  [3, 4, 5, 6, 7, 8,  9, 12, 18,  26, 41,  71]
			;highlim = [3, 4, 5, 6, 7, 8, 11, 17, 25,  40, 70, 200]
			
			lowlim  = [3, 4, 5, 6, 7, 8, 9,  11, 17, 25, 40, 70]
			highlim = [4, 5, 6, 7, 8, 9, 11, 17, 25, 40, 70, 200]
		end 
		else: message,"Unsupported nbin = "+ntostr(nbin)
	endcase 

end

function maxbcg::lambda_string, lowlim, highlim
	lamstr = textoidl('\lambda')
	ltstr = ' < '

	lowstr = ntostr(lowlim)
	highstr = ntostr(highlim)

	label = lowstr + ltstr + lamstr + ltstr + highstr
	return, label
end

function maxbcg::lambda_where_string, nbin, labels=labels, nodisplay=nodisplay

    self->lambda_bins, nbin, lowlim, highlim

    delvarx, labels
    if not keyword_set(nodisplay) then begin 
        nstr = textoidl('\lambda')
    endif else begin
        nstr = 'lambda'
    endelse
    for i=0L, nbin-1 do begin 

        lowstr = ntostr(lowlim[i])
        highstr = ntostr(highlim[i])

        tlabel = lowstr+' < '+nstr+' < '+highstr
        add_arrval, tlabel, labels

        tstring = $
            '(struct.photoz_cts GE 0.1) AND '+$
			'(struct.photoz_cts LE 0.3) AND '+$
            'struct.lambda GE '+ntostr(lowlim[i])+' AND '+$
            'struct.lambda LT '+ntostr(highlim[i])
        add_arrval, tstring, where_string
     endfor 
    return,where_string

end 

function maxbcg::lambda_z_where_string, nbin, labels=labels

  on_error, 2
  if n_elements(nbin) eq 0 then begin 
      print,'-Syntax: ws = mb->lambda_z_where_string(nbin)'
      print
      message,'Halting'
  endif 

  delvarx, labels
  lambda_nbin = nbin/2
  self->lambda_bins, lambda_nbin, lowlim, highlim

;  nstr = 'N200'
;  lstr = 'L200'

  nstr = 'N!D200!N'
  lstr = 'L!D200!N'

  ltstr = ' '+textoidl('\leq')+' '
  gtstr = ' '+textoidl('\geq')+' '

  delvarx, labels
  for i=0L, lambda_nbin-1 do begin 

      tstring = $
        '(struct.lambda  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.lambda  LE '+ntostr(highlim[i])+') AND '+$
        '(struct.photoz_cts GE 0.1 AND struct.photoz_cts LT 0.25)'
      add_arrval, tstring, where_string


      tlabel = self->lambda_string(lowlim[i], highlim[i])

      tlabel = tlabel + ' z < 0.25'
      add_arrval, tlabel, labels

      tstring = $
        '(struct.lambda  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.lambda  LE '+ntostr(highlim[i])+') AND '+$
        '(struct.photoz_cts GE 0.25) and (struct.photoz_cts lt 0.3)'
      add_arrval, tstring, where_string

      tlabel = 'z > 0.25'
      add_arrval, tlabel, labels


  endfor 
  return,where_string


end 






pro maxbcg::ngals200_bins, nbin, lowlim, highlim

  CASE nbin OF
      12: begin 
          lowlim =  [3, 4, 5, 6, 7, 8,  9, 12, 18,  26, 41,  71]
          highlim = [3, 4, 5, 6, 7, 8, 11, 17, 25,  40, 70, 220]
      end 
      ELSE: message,"Unsupported nbin = "+ntostr(nbin)
  endcase 

end 

function maxbcg::ngals200_string, lowlim, highlim

  nstr = 'N!D200!N'
  ltstr = ' '+textoidl('\leq')+' '

  ;; We know these are integers
  if lowlim eq highlim then begin 
      tlabel = nstr + ' = '+ntostr(lowlim)
  endif else begin 
      lowstr = ntostr(lowlim)
      highstr = ntostr(highlim)
      
      tlabel = lowstr+ltstr+nstr+ltstr+highstr
  endelse 
  return, tlabel
end 
function maxbcg::ngals200_where_string, nbin, labels=labels, nodisplay=nodisplay

    self->ngals200_bins, nbin, lowlim, highlim
    if not keyword_set(nodisplay) then begin 
        nstr = 'N!D200!N'
    endif else begin
        nstr = 'N200'
    endelse
 
    if self->catalog() eq 'hv1halo' then begin
        ztag='z'
        ntag='ngals'
    endif else begin
        ztag='photoz_cts'
        ntag='ngals200'
    endelse

  delvarx, labels
  ltstr = ' '+textoidl('\leq')+' '
  for i=0L, nbin-1 do begin 

      ;; We know these are integers
      if lowlim[i] eq highlim[i] then begin 
          tlabel = nstr + ' = '+ntostr(lowlim[i])
      endif else begin 
          lowstr = ntostr(lowlim[i])
          highstr = ntostr(highlim[i])

          tlabel = lowstr+ltstr+nstr+ltstr+highstr
      endelse 
      add_arrval, tlabel, labels

      tstring = $
        'struct.'+ztag+' GE 0.1 AND '+$
        'struct.'+ntag+' GE '+ntostr(lowlim[i])+' AND '+$
        'struct.'+ntag+' LE '+ntostr(highlim[i])
      add_arrval, tstring, where_string
  endfor 
  return,where_string

end 


pro maxbcg::sngals_bins, nsub, lowlim, highlim

	case self->catalog() of
		'gmbcg10': begin
			max_sngals = 278
			case nsub of
				6: begin
					lowlim =  [10, 13, 20,  31, 46,  81]
					highlim = [12, 19, 30,  45, 80, max_sngals]
				end
					else: message,string('Unkown nsub: ',nsub)
			endcase
		end
		else: message,'unknown catalog: '+self->catalog()
	endcase

end


function maxbcg::sngals_where_string, nbin, labels=labels, nodisplay=nodisplay

    self->sngals_bins, nbin, lowlim, highlim
    if not keyword_set(nodisplay) then begin 
        nstr = textoidl('N_{s}')
		ltstr = ' '+textoidl('\leq')+' '
    endif else begin
        nstr = 'Ns'
		ltstr = ' <= '
    endelse
 
	ztag = 'photoz'
	ntag = 'sngals'

	delvarx, labels
	for i=0L, nbin-1 do begin 

		;; We know these are integers
		if lowlim[i] eq highlim[i] then begin 
			tlabel = nstr + ' = '+ntostr(lowlim[i])
		endif else begin 
			lowstr = ntostr(lowlim[i])
			highstr = ntostr(highlim[i])

			tlabel = lowstr+ltstr+nstr+ltstr+highstr
		endelse 
		add_arrval, tlabel, labels

		tstring = $
			'struct.'+ntag+' GE '+ntostr(lowlim[i])+' AND '+$
			'struct.'+ntag+' LE '+ntostr(highlim[i])
		add_arrval, tstring, where_string
	endfor 
	return,where_string

end 









;;
;; calculate the quantiles  ilum200 in ngals200 bins; used for splits
;;
function maxbcg::ngals200_ilum200_quantile_file, ngals200_nbin
  on_error, 2
  if n_elements(ngals200_nbin) eq 0 then begin 
      print,'-Syntax: file = mb->ngals200_ilum200_quantile_file(nbin)'
      print
      message,'Halting'
  endif 

  subtype = 'ngals200_ilum200_'+strn(ngals200_nbin,len=2,padchar='0')+'_2'
  dir = self->objshear::lensdir('combined', subtype=subtype, /createdir)
  dir = repstr(dir, 'combined', '')
  file = subtype+'_quantiles.st'
  file = concat_dir(dir, file)
  return, file
end 
pro maxbcg::calc_ngals200_ilum200_quantiles, ngals200_nbin, struct=struct

  on_error, 2
  if n_elements(ngals200_nbin) eq 0 then begin 
      print,'-Syntax: mb->calc_ngals200_ilum200_quantiles,nbin'
      on_error, 2
      message,'Halting'
  endif 

  if n_elements(struct) eq 0 then begin 
      struct = self->get()
  endif 

  self->ngals200_bins, ngals200_nbin, lowlim, highlim

  quant = 2.0/3.0

  for i=0L, ngals200_nbin-1 do begin
      w=where(struct.ngals200 GE lowlim[i] AND $
              struct.ngals200 LE highlim[i])

      tmp_ilum200_quant = weighted_quantile(struct[w].ilum200, quant=quant)
      
      print,lowlim[i],highlim[i],tmp_ilum200_quant

      add_arrval, tmp_ilum200_quant, ilum200_quant

  endfor 

  outstruct = $
    { $
      quantile:         quant,        $
      ngals200_lowlim:  lowlim,       $
      ngals200_highlim: highlim,      $
      ilum200_quantiles: ilum200_quant $
    }

  file= self->ngals200_ilum200_quantile_file(ngals200_nbin)
  print
  print,'Writing to file: ',file
  write_idlstruct, outstruct, file, /ascii

end 

;; used calc_ngals200_ilum200_median to get this
function maxbcg::ngals200_ilum200_quantiles, ngals200_nbin

  on_error, 2
  if n_elements(ngals200_nbin) eq 0 then begin 
      print,'-Syntax: quant = mb->ngals200_ilum200_quantiles(nbin)'
      print
      message,'Halting'
  endif 

  CASE ngals200_nbin OF
      12: begin 
          file = self->ngals200_ilum200_quantile_file(ngals200_nbin)
          struct = read_idlstruct(file)
          return, struct.ilum200_quantiles
      end 
      ELSE: message,'unsupported number of bins: '+ntostr(nsub)
  endcase 

end 
function maxbcg::ngals200_ilum200_where_string, nbin, labels=labels

  on_error, 2
  if n_elements(nbin) eq 0 then begin 
      print,'-Syntax: ws = mb->ngals200_ilum200_where_string(nbin)'
      print
      message,'Halting'
  endif 

  delvarx, labels
  ngals200_nbin = nbin/2
  self->ngals200_bins, ngals200_nbin, lowlim, highlim
  qilum = self->ngals200_ilum200_quantiles(ngals200_nbin)

;  nstr = 'N200'
;  lstr = 'L200'

  nstr = 'N!D200!N'
  lstr = 'L!D200!N'

  ltstr = ' '+textoidl('\leq')+' '
  gtstr = ' '+textoidl('\geq')+' '

  delvarx, labels
  for i=0L, ngals200_nbin-1 do begin 

      tstring = $
        '(struct.photoz_cts GE 0.1) AND '+$
        '(struct.ngals200 GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ngals200 LE '+ntostr(highlim[i])+') AND '+$
        '(struct.ilum200  LT '+ntostr(qilum[i])+')'
      add_arrval, tstring, where_string


      if qilum[i] GE 10 then begin 
          format='(F10.1)'
      endif else begin 
          format='(F10.2)'
      endelse 

      qstr = ntostr(qilum[i], format=format)

      ;; We know these are integers
      if lowlim[i] eq highlim[i] then begin 
          nlabel = nstr + ' = '+ntostr(lowlim[i])
      endif else begin 
          lowstr = ntostr(lowlim[i])
          highstr = ntostr(highlim[i])

          nlabel = $
            lowstr+ltstr+nstr+ltstr+highstr+'!c'
      endelse       

      tlabel = lstr+' < '+qstr
      add_arrval, tlabel, labels

      tstring = $
        '(struct.photoz_cts GE 0.1) AND '+$
        '(struct.ngals200 GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ngals200 LE '+ntostr(highlim[i])+') AND '+$
        '(struct.ilum200  GE '+ntostr(qilum[i])+')'
      add_arrval, tstring, where_string


      tlabel = lstr+' > '+qstr+'!c'+nlabel
      add_arrval, tlabel, labels


  endfor 
  return,where_string


end 



function maxbcg::ngals200_z_where_string, nbin, labels=labels

  on_error, 2
  if n_elements(nbin) eq 0 then begin 
      print,'-Syntax: ws = mb->ngals200_z_where_string(nbin)'
      print
      message,'Halting'
  endif 

  delvarx, labels
  ngals200_nbin = nbin/2
  self->ngals200_bins, ngals200_nbin, lowlim, highlim

;  nstr = 'N200'
;  lstr = 'L200'

  nstr = 'N!D200!N'
  lstr = 'L!D200!N'

  ltstr = ' '+textoidl('\leq')+' '
  gtstr = ' '+textoidl('\geq')+' '

  delvarx, labels
  for i=0L, ngals200_nbin-1 do begin 

      tstring = $
        '(struct.ngals200  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ngals200  LE '+ntostr(highlim[i])+') AND '+$
        '(struct.photoz_cts GE 0.1 AND struct.photoz_cts LT 0.25)'
      add_arrval, tstring, where_string


      tlabel = self->ngals200_string(lowlim[i], highlim[i])

      tlabel = tlabel + ' z < 0.25'
      add_arrval, tlabel, labels

      tstring = $
        '(struct.ngals200  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ngals200  LE '+ntostr(highlim[i])+') AND '+$
        '(struct.photoz_cts GE 0.25)'
      add_arrval, tstring, where_string

      tlabel = 'z > 0.25'
      add_arrval, tlabel, labels


  endfor 
  return,where_string


end 













pro maxbcg::ilum200_bins, nbin, lowlim, highlim

  CASE nbin OF
      12: begin 
          lowlim =  $
            [0.60,1.16,1.35,1.78,2.53,3.42,6.07,8.91,15.5,23.5,33.8,63.6]
          highlim = $
            [1.16,1.35,1.78,2.53,3.42,6.07,8.91,15.5,23.5,33.8,63.6,1587.0]
      end 
      16: begin 
          ;; Because of the exponential cutoff, we will bin
          ;; logarithmically up to L of 141*10^10, then throw
          ;; all the rest into the last bin
;          llmin = alog10( 0.49 ) ; 10^10

          ;; Note this must be coupled with ngals200 >= 3

          ;; A first bin from min to 2.49

          
          lmin = 5.00
          lmax = 140.0

          llmin = alog10( lmin )
          llmax = alog10( lmax )

          nedge = 16
          bin_spacing = (llmax - llmin)/(nedge-1.0)
          bin_edges = llmin + bin_spacing*findgen(nedge)

          lowlim  = 10.0^bin_edges[0:nedge-2]
          highlim = 10.0^bin_edges[1:nedge-1]

          lowlim = [lowlim, highlim[nedge-2]]
          highlim = [highlim, 450.0]

          break

          tmp_min = 2.49
          logtmin = alog10( 2.49 )

          tnbin = 15
          bin_spacing = (llmax-logtmin)/(tnbin-1.0)
          bin_edges = logtmin + bin_spacing*findgen(tnbin)

          lowlim  = [lmin,    10.0^bin_edges[0:tnbin-2] ]
          highlim = [tmp_min, 10.0^bin_edges[1:tnbin-1] ]
          
          lowlim = [lowlim, highlim[tnbin-1]]
          highlim = [highlim, 450.0]

      end 
      ELSE: message,"Unsupported nbin = "+ntostr(nbin)
  endcase 

end 
function maxbcg::keepfigs, val
  if val GE 1000 then begin 
      nc=6
  endif else if val GE 100 then begin 
      nc=5
  endif else begin 
      nc=4
  endelse 
  return, nc
end 

function maxbcg::ilum200_where_string, nbin, labels=labels, nodisplay=nodisplay

    self->ilum200_bins, nbin, lowlim, highlim

    delvarx, labels
    if not keyword_set(nodisplay) then begin 
        nstr = 'L!D200!N'
    endif else begin
        nstr = 'L200'
    endelse
    for i=0L, nbin-1 do begin 

        nc = self->keepfigs(lowlim[i])
        lowstr = ntostr(lowlim[i],nc,/round)
        nc = self->keepfigs(highlim[i])
        highstr = ntostr(highlim[i],nc,/round)

        tlabel = lowstr+' < '+nstr+' < '+highstr
        add_arrval, tlabel, labels

        tstring = $
            '(struct.photoz_cts GE 0.1) AND '+$
            'struct.ngals200 ge 3 AND '+$
            'struct.ilum200 GE '+ntostr(lowlim[i])+' AND '+$
            'struct.ilum200 LT '+ntostr(highlim[i])
        add_arrval, tstring, where_string
     endfor 
    return,where_string

end 











;;
;; calculate the quantiles  ngals200 in ilum200 bins; used for splits
;;
function maxbcg::ilum200_ngals200_quantile_file, ilum200_nbin
  on_error, 2
  if n_elements(ilum200_nbin) eq 0 then begin 
      print,'-Syntax: file = mb->ilum200_ngals200_quantile_file(nbin)'
      print
      message,'Halting'
  endif 

  subtype = 'ilum200_ngals200_'+strn(ilum200_nbin,len=2,padchar='0')+'_2'
  dir = self->objshear::lensdir('combined', subtype=subtype, /createdir)
  dir = repstr(dir, 'combined', '')
  file = subtype+'_quantiles.st'
  file = concat_dir(dir, file)
  return, file
end 
pro maxbcg::calc_ilum200_ngals200_quantiles, ilum200_nbin, struct=struct

  on_error, 2
  if n_elements(ilum200_nbin) eq 0 then begin 
      print,'-Syntax: file = mb->calc_ilum200_ngals200_quantiles(nbin)'
      print
      message,'Halting'
  endif 

  if n_elements(struct) eq 0 then begin 
      struct = self->get()
  endif 

  ws = self->ilum200_where_string(ilum200_nbin)

  quant = 2.0/3.0

  keep = self->objshear::struct_select(struct, ws)

  for i=0L, ilum200_nbin-1 do begin

      w = *keep[i]

      tmp_ngals200_quant = $
        long( weighted_quantile(struct[w].ngals200, quant=quant) )
      
      print,ws[i],tmp_ngals200_quant

      add_arrval, tmp_ngals200_quant, ngals200_quant

  endfor 

  outstruct = $
    { $
      quantile:         quant,        $
      ilum200_where_string: ws, $
      ngals200_quantiles: ngals200_quant $
    }

  file= self->ilum200_ngals200_quantile_file(ilum200_nbin)
  print
  print,'Writing to file: ',file
  write_idlstruct, outstruct, file

end 

;; used calc_ilum200_ngals200_median to get this
function maxbcg::ilum200_ngals200_quantiles, ilum200_nbin

  on_error, 2
  if n_elements(ilum200_nbin) eq 0 then begin 
      print,'-Syntax: quant = mb->ilum200_ngals200_quantiles(nbin)'
      print
      message,'Halting'
  endif 

  CASE ilum200_nbin OF
      16: begin 
          file = self->ilum200_ngals200_quantile_file(ilum200_nbin)
          struct = read_idlstruct(file)
          return, struct.ngals200_quantiles
      end 
      ELSE: message,'unsupported number of bins: '+ntostr(nsub)
  endcase 

end 
function maxbcg::ilum200_ngals200_where_string, nbin, labels=labels

  on_error, 2
  if n_elements(nbin) eq 0 then begin 
      print,'-Syntax: ws = mb->ilum200_ngals200_where_string(nbin)'
      print
      message,'Halting'
  endif 

  ilum200_nbin = nbin/2
  self->ilum200_bins, ilum200_nbin, lowlim, highlim
  qngals = self->ilum200_ngals200_quantiles(ilum200_nbin)


;  nstr = 'N200'
;  lstr = 'L200'

  nstr = 'N!D200!N'
  lstr = 'L!D200!N'

  ltstr = ' '+textoidl('\leq')+' '
  gtstr = ' '+textoidl('\geq')+' '

  delvarx, labels
  for i=0L, ilum200_nbin-1 do begin 

      tstring = $
        '(struct.photoz_cts GE 0.1) AND '+$
        '(struct.ngals200 GE 3) AND '+$
        '(struct.ilum200  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ilum200  LT '+ntostr(highlim[i])+') AND '+$
        '(struct.ngals200 LT '+ntostr(qngals[i])+')'
      add_arrval, tstring, where_string


      nc = self->keepfigs(lowlim[i])
      lowstr = ntostr(lowlim[i],nc,/round)
      nc = self->keepfigs(highlim[i])
      highstr = ntostr(highlim[i],nc,/round)

      tlabel = lowstr+' < '+nstr+' < '+highstr
      qstr = ntostr(qngals[i])

      tlabel = tlabel + ' ' +nstr + ' < '+qstr
      add_arrval, tlabel, labels



      tstring = $
        '(struct.photoz_cts GE 0.1) AND '+$
        '(struct.ngals200 GE 3) AND '+$
        '(struct.ilum200  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ilum200  LT '+ntostr(highlim[i])+') AND '+$
        '(struct.ngals200 GE '+ntostr(qngals[i])+')'
      add_arrval, tstring, where_string

      tlabel = nstr + gtstr+qstr
      add_arrval, tlabel, labels


  endfor 
  return,where_string


end 













function maxbcg::ilum200_z_where_string, nbin, labels=labels

  on_error, 2
  if n_elements(nbin) eq 0 then begin 
      print,'-Syntax: ws = mb->ilum200_ngals200_where_string(nbin)'
      print
      message,'Halting'
  endif 

  ilum200_nbin = nbin/2
  self->ilum200_bins, ilum200_nbin, lowlim, highlim

;  nstr = 'N200'
;  lstr = 'L200'

  nstr = 'N!D200!N'
  lstr = 'L!D200!N'

  ltstr = ' '+textoidl('\leq')+' '
  gtstr = ' '+textoidl('\geq')+' '

  delvarx, labels
  for i=0L, ilum200_nbin-1 do begin 

      tstring = $
        '(struct.ngals200 GE 3) AND '+$
        '(struct.ilum200  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ilum200  LT '+ntostr(highlim[i])+') AND '+$
        '(struct.photoz_cts GE 0.1 AND struct.photoz_cts LT 0.25)'
      add_arrval, tstring, where_string


      nc = self->keepfigs(lowlim[i])
      lowstr = ntostr(lowlim[i],nc,/round)
      nc = self->keepfigs(highlim[i])
      highstr = ntostr(highlim[i],nc,/round)

      tlabel = lowstr+' < '+lstr+' < '+highstr

      tlabel = tlabel + ' z < 0.25'
      add_arrval, tlabel, labels



      tstring = $
        '(struct.ngals200 GE 3) AND '+$
        '(struct.ilum200  GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ilum200  LT '+ntostr(highlim[i])+') AND '+$
        '(struct.photoz_cts GE 0.25)'
      add_arrval, tstring, where_string

      tlabel = 'z > 0.25'
      add_arrval, tlabel, labels


  endfor 
  return,where_string


end 


pro maxbcg::bcg_kgmr_calculate_medians, struct, subtype

  keep = self->where_select(struct, subtype)

  allind = combine_ptrlist(keep, /nofree)

  nbin = n_elements(keep)

  binsize = 0.02
  simpctable, colorlist=clist
  plothist, struct[allind].bcg_kgmr, bin=binsize, $
    /ylog, yrange=[0.9,5.e4], ystyle=3, ytickf='loglabels', $
    xtitle='g-r', ytitle='N'

  print,'[',format='($,a)'
  for i=0L, nbin-1 do begin 
      ind = *keep[i]
      plothist, struct[ind].bcg_kgmr, bin=binsize, color=clist[i], /overplot

      print,ntostr(median(struct[ind].bcg_kgmr)), format='($,a)'

      if i lt nbin-1 then print,', ',format='($,a)'

  endfor 
  print,']'
  ptr_free, keep

end 
function maxbcg::ilum200_bcg_kgmr_medians, ilum200_nbin

  case ilum200_nbin of
      16: begin 
          meds = [1.39531, 1.39975, 1.40431, 1.40879, 1.41337, 1.41670, $
                  1.42038, 1.42541, 1.42903, 1.42924, 1.43234, 1.43177, $
                  1.42830, 1.43166, 1.42850, 1.44497]
      end 
      else: message,'Unsupported subtype: '+ntostr(subtype)
  endcase 

  return, meds
end 

function maxbcg::ilum200_bcg_kgmr_where_string, nbin, labels=labels

  on_error, 2
  if n_elements(nbin) eq 0 then begin 
      print,'-Syntax: ws = mb->ilum200_bcg_kgmr_where_string(nbin)'
      print
      message,'Halting'
  endif 

  ilum200_nbin = nbin/2
  self->ilum200_bins, ilum200_nbin, lowlim, highlim
  meds = self->ilum200_bcg_kgmr_medians(ilum200_nbin)

;  nstr = 'N200'
;  lstr = 'L200'

  nstr = 'N!D200!N'
  lstr = 'L!D200!N'

  ltstr = ' '+textoidl('\leq')+' '
  gtstr = ' '+textoidl('\geq')+' '

  delvarx, labels
  for i=0L, ilum200_nbin-1 do begin 

      medstr = ntostr(meds[i])
      medstr_print = ntostr(meds[i],4)
      tstring = $
        '(struct.ngals200 GE 3) AND '+$
        '(struct.ilum200 GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ilum200 LT '+ntostr(highlim[i])+') AND '+$
        '(struct.bcg_kgmr LT '+medstr+')'
      add_arrval, tstring, where_string


      nc = self->keepfigs(lowlim[i])
      lowstr = ntostr(lowlim[i],nc,/round)
      nc = self->keepfigs(highlim[i])
      highstr = ntostr(highlim[i],nc,/round)

      tlabel = lowstr+' < '+lstr+' < '+highstr

      tlabel = tlabel + ' BCG g-r < '+medstr_print
      add_arrval, tlabel, labels



      tstring = $
        '(struct.ngals200 GE 3) AND '+$
        '(struct.ilum200 GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.ilum200 LT '+ntostr(highlim[i])+') AND '+$
        '(struct.bcg_kgmr GE '+medstr+')'
      add_arrval, tstring, where_string

      tlabel = 'BCG g-r > '+medstr_print
      add_arrval, tlabel, labels


  endfor 
  return,where_string


end 






pro maxbcg::zbins, nbin, zlowlim, zhighlim

    delvarx, zlowlim, zhighlim
    case nbin of
        4: begin
            ; equal number per bin
            zlowlim  = [0.10, 0.15, 0.20, 0.25]
            zhighlim = [0.15, 0.20, 0.25, 0.30]
        end
        else: message,'Unsupported number of z bins: '+ntostr(nbin)
    endcase 
   
end

function maxbcg::zbin_where_string, nbin, labels=labels

    delvarx, labels
    self->zbins, nbin, lowlim, highlim

    f='(F0.2)'
    for i=0L, nbin-1 do begin
        ts = $
            '(struct.ngals200 ge 3) and '+$
            '(struct.photoz_cts ge '+ntostr(lowlim[i])+') and '+$
            '(struct.photoz_cts lt '+ntostr(highlim[i])+')'

        label = string(lowlim[i],f=f)+' < z < '+string(highlim[i],f=f)
        add_arrval, ts, where_string
        add_arrval, label, labels
    endfor

    return, where_string

end












pro maxbcg::lx_bins, nbin, lowlim, highlim

  CASE nbin OF
      2: begin 
          lowlim = [0.0, 4.58]
          highlim = [4.58, 33]
      end 
      ELSE: message,'Unsupported number of bins: '+ntostr(nbin)
  endcase 

end 
function maxbcg::lx_where_string, nbin

  self->lx_bins, nbin, lowlim, highlim
  for i=0L, nbin-1 do begin 

      tstring = $
        'struct.lx GE '+ntostr(lowlim[i])+' AND '+$
        'struct.lx LT '+ntostr(highlim[i])
      add_arrval, tstring, where_string
  endfor 
  return,where_string

end 




function maxbcg::centerclass_kdeclass_binvals, subtype
    case strlowcase(subtype) of
        'centerclass5_alt3': return, [0, 1, 2]
        'centerclass5': return, [-9999, 0, 1, 2, 3]
        'centerclass6': return, [-9999, 0, 1, 2, 3, 4]
        else:message,'Unmatched subtype: '+subtype
    endcase
end

function maxbcg::centerclass_where_string, nbin, labels=labels, binvals=binvals

    binvals = self->centerclass_kdeclass_binvals(nbin)
    for i=0L, nbin-1 do begin
        ws = 'struct.ngals200 ge 5 and struct.kde_class eq '+ntostr(binvals[i])
        lab = textoidl('N_{200} \geq 5 class = '+ntostr(binvals[i]))

        add_arrval, ws, wstrings
        add_arrval, lab, labels
    endfor

    return, wstrings
end

; Obsolete
function maxbcg::centerclass_select_match, struct, nkeep, labels=labels

    w0=where(struct.ngals200 ge 5 and struct.kde_class eq 0)
    w1=where(struct.ngals200 ge 5 and struct.kde_class eq 1)

    ; now match histograms of ilum200
    binsize = 0.1
    logL1 = alog10(struct[w1].ilum200)
    logL0 = alog10(struct[w0].ilum200)
    ;!p.multi=[0,0,2]

    lmin = 0.45
    lmax = 2.0
    matchind0 = match_hist(logL1, logL0, binsize, min=lmin, max=lmax, indices1=matchind1)

    ;plothist, logL0[matchind0], binsize=binsize, min=lmin, max=lmax
    ;plothist, logL1[matchind1], binsize=binsize, min=lmin, max=lmax, $
    ;    /over, color=!green
    ;!p.multi=0

    keep = ptrarr(2)
    keep[0] = ptr_new(w0[matchind0], /no_copy)
    keep[1] = ptr_new(w1[matchind1], /no_copy)
   
    nkeep = lonarr(2)
    nkeep[0] = n_elements(matchind0)
    nkeep[1] = n_elements(matchind1)
    return, keep 
end

; This matches to the first value in matchvals
function maxbcg::centerclass_select_nmatch, $
        subtype, struct, nkeep, labels=labels

    matchvals = self->centerclass_kdeclass_binvals(subtype)
    matchval0 = matchvals[0]
    w0=where(struct.ngals200 ge 5 and struct.kde_class eq matchval0)

    binsize = 0.1
    lmin = 0.45
    lmax = 2.0
    logL0 = alog10(struct[w0].ilum200)

    bs = binner(logL0, binsize=binsize, min=lmin, max=lmax, rev=rev, $
        indices=ind0)

    nmatch=n_elements(matchvals)
    keep = ptrarr(nmatch)
    nkeep = lonarr(nmatch)

    keep[0] = ptr_new(w0[ind0], /no_copy)
    nkeep[0] = n_elements(ind0)

    !p.multi = [0,2,2]
    for i=1L, nmatch-1 do begin

        matchval=matchvals[i]
        wi=where(struct.ngals200 ge 5 and struct.kde_class eq matchval)



        logLi = alog10(struct[wi].ilum200)
        matchind0 = match_hist(logLi, logL0, binsize, min=lmin, max=lmax, indices1=matchindi)

        if 1 then begin
            plothist, logL0[matchind0], binsize=binsize, min=lmin, max=lmax, $
                /norm
            plothist, logLi[matchindi], binsize=binsize, min=lmin, max=lmax, $
                /norm, $
                /over, color=!green
            key=prompt_kbrd('hit a key')
        endif

        keep[i] = ptr_new(wi[matchindi], /no_copy)
        nkeep[i] = n_elements(matchindi)
    endfor
    !p.multi=0
    return, keep 
end

function maxbcg::lxid_dir
    dir = self->catdir()
    dir = path_join(dir, 'lxid')
    return,dir
end
function maxbcg::lxid_read
    dir=self->lxid_dir()
    f = path_join(dir, 'lx_ids_for_erin.fit.gz')
    print,'Reading lxid file: ',f
    lxid = mrdfits(f,1)
    return, lxid
end
function maxbcg::lxid_match, struct, subtype, nkeep
    common lxid_match_block, lxid
    ; will have to generalize when nbin != 1

    if  n_elements(lxid) eq 0 then begin
        lxid=self->lxid_read()
    endif

    case subtype of
        'lxid-all': mem_match_id = lxid.all
        'lxid-sig3': mem_match_id = lxid.sig3
        'lxid-sig5': mem_match_id = lxid.sig5
        'lxid-noras': mem_match_id = lxid.noras
        else: message,'bad subtype: '+subtype
    endcase

    match, struct.mem_match_id, mem_match_id, mcat, mlx, /sort
    if mcat[0] eq -1 then nkeep=0 else nkeep=n_elements(mcat)

    splog,'Found ',nkeep,' matches for ',subtype,format='(a,i0,a,a)'

    return, ptr_new(mcat)
end



function maxbcg::subtype_nbin, subtype
  if n_elements(subtype) eq 0 then begin 
      message,'-Syntax: nbin=mb->nbin(subtype)'
  endif 
  CASE subtype OF
	  'sngals6': return, 6
      'ngals12': return,12

      'ngals200_8': return, 8
      'ngals200_12': return, 12
      'ngals200_ilum200_12_2': return,24
      'ngals200_z_12_2': return, 24
	  
	  'lambda12': return, 12
      'lambda_z_12_2': return, 24

      'ilum200_12': return, 12
      'ilum200_16': return, 16
      'ilum200_ngals200_16_2': return,32
      'ilum200_z_16_2': return,32
      'ilum200_bcg_kgmr_16_2': return, 32

      'lx2': return,2

      'zbin4': return, 4

      'centerclass5': return, 5
      'centerclass5_alt3': return, 3
      'centerclass6': return, 6
      'centerclass_alt2': return, 2

      'lxid-all': return, 1
      'lxid-sig3': return, 1
      'lxid-sig5': return, 1
      'lxid-noras': return, 1
      ELSE: message,'Unknown subtype: '+strn(subtype)
  endcase 
end 

function maxbcg::select_type, subtype
    case strlowcase(subtype) of
        'centerclass5_alt3': return, 'alt'
        'lxid-all': return, 'alt'
        'lxid-sig3': return, 'alt'
        'lxid-sig5': return, 'alt'
        'lxid-noras': return, 'alt'
        else: return,'where_string'
    endcase
end

function maxbcg::average_tags, subtype=subtype
  deftags = ['z','photoz_cts','photoz','sngals','wngals','ngals200','ilum200','bcg_ilum','bcg_kgmr']
  return, deftags
end 

function maxbcg::where_string, subtype, nbin=nbin, labels=labels, nodisplay=nodisplay, average_tags=average_tags

	if n_elements(subtype) eq 0 then begin 
		on_error, 2
		print,'-syntax: ws=mb->where_string(subtype, labels=, average_tags=, nbin=)'
		message,'halting'
	endif 

	nbin = self->subtype_nbin(subtype)

	average_tags = self->average_tags(subtype=subtype)

	cat = self->catalog()
	case subtype of 
		'sngals6': begin
			return, self->sngals_where_string(nbin, labels=labels, $
				nodisplay=nodisplay)
		end
		'ngals200_8': begin 
			return, self->ngals200_where_string(nbin, labels=labels)
		end 
		'ngals200_12': begin 
			return, self->ngals200_where_string(nbin, labels=labels,$
				nodisplay=nodisplay)
		end 
		'ngals200_ilum200_12_2': begin 
			return, self->ngals200_ilum200_where_string(nbin, labels=labels)
		end 
		'ngals200_z_12_2': begin 
			return, self->ngals200_z_where_string(nbin, labels=labels)
		end 

		'lambda12': begin
			return, self->lambda_where_string(nbin, labels=labels)
		end
		'lambda_z_12_2': begin 
			return, self->lambda_z_where_string(nbin, labels=labels)
		end 


		'ilum200_12': begin 
			return, self->ilum200_where_string(nbin, labels=labels)
		end 
		'ilum200_16': begin 
			return, self->ilum200_where_string(nbin, labels=labels, $
				nodisplay=nodisplay)
		end 
		'ilum200_ngals200_16_2': begin 
			return, self->ilum200_ngals200_where_string(nbin, labels=labels)
		end 
		'ilum200_z_16_2': begin 
			return, self->ilum200_z_where_string(nbin, labels=labels)
		end 
		'ilum200_bcg_kgmr_16_2': begin 
			return, self->ilum200_bcg_kgmr_where_string(nbin, labels=labels)
		end 


		'lx2': begin 
			return, self->lx_where_string(nbin)
		end 

		'zbin4': begin
			return, self->zbin_where_string(nbin, labels=labels)
		end

		'centerclass5': begin
			return, self->centerclass_where_string(nbin, labels=labels)
		end 
		'centerclass6': begin
			return, self->centerclass_where_string(nbin, labels=labels)
		end 
		else: message,'unknown type: '+ntostr(subtype)
	endcase 

end 
function maxbcg::labels, subtype, bin=bin
  t=self->where_string(subtype, labels=labels)
  if n_elements(labels) eq 0 then message,'labels undefined for subtype: '+$
    strlowcase(subtype)

  if n_elements(bin) ne 0 then begin 
      minb=min(bin, max=maxb)
      nl = n_elements(labels)
      if minb LT 0 OR minb GE nl then message,'bin keyword out of range'
      labels = labels[bin]
  endif 
  return,labels
end 

function maxbcg::subtypes
  subtypes = $
    ['ilum200_16','ilum200_ngals200_16_2','ilum200_z_16_2',$
     'ngals200_12','ngals200_ilum200_12_2','ngals200_z_12_2']
  return,subtypes
end 


function maxbcg::altselect, struct, subtype, nkeep=nkeep, nbin=nbin, $
    average_tags=average_tags

    case subtype of 
        'centerclass_alt2': begin
            keep = self->centerclass_select_match(struct, nkeep)
        end
        'centerclass5_alt3': begin
            keep = self->centerclass_select_nmatch(subtype, struct, nkeep)
        end
        'lxid-all': begin
            keep = self->lxid_match(struct, subtype, nkeep)
        end
        'lxid-sig3': begin
            keep = self->lxid_match(struct, subtype, nkeep)
        end
        'lxid-sig5': begin
            keep = self->lxid_match(struct, subtype, nkeep)
        end
        'lxid-noras': begin
            keep = self->lxid_match(struct, subtype, nkeep)
        end
        else: message,'Only support centerclass_alt2 now'
    endcase
    average_tags = self->average_tags(subtype=subtype)
    nbin=n_elements(keep)

    return, keep
end



; Print some info and statistics for a binning
pro maxbcg::describe_binning, type, bcg, meantags=meantags

    if n_elements(bcg) eq 0 then bcg=self->get()
    w=self->where_string(type, labels=labels, /nodisplay)
    k = self->where_select(bcg,type)
    nbin=n_elements(labels)
    for i=0L,nbin-1 do begin 
      
        s = labels[i]+' '+ntostr(n_elements(*k[i])) 
        if n_elements(meantags) ne 0 then begin 
            ms=struct_stats(bcg, meantags, index=*k[i])
            ; Build up a string to print
            nt = n_elements(ms)
            s=s + ' '
            for it=0L, nt-1 do begin 
                s = s+ms[it].tagname+': '+$
                    ntostr(ms[it].mean)+' '+$
                    ntostr(ms[it].sdev)+' '+$
                    ntostr(ms[it].err)
            endfor 
        endif
     
        print,s
    endfor
    ptr_free, k
end




pro maxbcg::run_some_stuff  
  subtypes = self->subtypes()
  nsub = n_elements(subtypes)
;  for i=0L, nsub-1 do begin 
;      subtype = subtypes[i]
;      
;      self->sub,subtype
;  endfor 

;  for i=0L, nsub-1 do begin 
;      subtype = subtypes[i]
;      
;      self->corr,subtype=subtype
;  endfor 

  for i=0L, nsub-1 do begin 
      subtype = subtypes[i]
      
      self->jackknife,subtype=subtype
  endfor 


end 












pro maxbcg::ngals_vs_z, bcg, zbin, s, low, high, num

  if n_elements(bcg) eq 0 then bcg = self->get()

  plot_dir = '~/plots/MaxBCG/'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Bin by ngals and look at redshift histogram
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  begplot, name=plot_dir+'zhist_ngalsbins.ps',/color, xsize=7, ysize=7

  simpctable, colorlist=colors

  nbin = 6
  self->ngals_bins, nbin, lowlim, highlim

  binsize_z = 0.01  
  w = where(bcg.ngals GE lowlim[0] AND bcg.ngals LE highlim[0], n)
  plothist, bcg[w].z, bin=binsize_z, /norm, $
    xtitle = 'z', ytitle = 'P(z)'
  for i=1L, nbin-1 do begin 
      w = where(bcg.ngals GE lowlim[i] AND bcg.ngals LE highlim[i], n)
      plothist, bcg[w].z, bin=binsize_z, /norm, /overplot, color=colors[i]
  endfor 

  mess = ntostr(lowlim)+' <= Ngals <= '+ntostr(highlim)
  lclr = colors[0:nbin-1]
  legend, mess, line=replicate(0, nbin), color=lclr, box=0

  endplot

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; bin by number in z and plot mean ngals
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  begplot, name=plot_dir + 'ngals_vs_z.ps'

  yrange = [0, 17]
  xtitle = 'z'
  ytitle = '<N!Dgals!N>'
  psym = 8
  symsize = 0.2
  xticklen=0.04
  errthick=1


  erase & multiplot, [1,4]
  
  
  binner_bynum, bcg.z, bcg.ngals, 500, zb, ngb, ngberr
  ploterror, zb, ngb, ngberr, psym=psym, yrange=yrange, ystyle=1, $
    ytitle=ytitle, symsize=symsize, xticklen=xticklen, $
    /nohat, errthick=errthick
  legend,'nperbin = 500',box=0, charsize=1

  multiplot
  binner_bynum, bcg.z, bcg.ngals, 1000, zb, ngb, ngberr
  ploterror, zb, ngb, ngberr, psym=psym, yrange=yrange, ystyle=1, $
    ytitle=ytitle, symsize=symsize, xticklen=xticklen, $
    /nohat, errthick=errthick
  legend,'nperbin = 1000',box=0, charsize=1

  multiplot
  binner_bynum, bcg.z, bcg.ngals, 5000, zb, ngb, ngberr
  ploterror, zb, ngb, ngberr, psym=psym, yrange=yrange, ystyle=1, $
    ytitle=ytitle, symsize=symsize, xticklen=xticklen, $
    /nohat, errthick=errthick
  legend,'nperbin = 5000',box=0, charsize=1

  multiplot
  binner_bynum, bcg.z, bcg.ngals, 10000, zb, ngb, ngberr
  ploterror, zb, ngb, ngberr, psym=psym, yrange=yrange, ystyle=1, $
    xtitle=xtitle, ytitle=ytitle, symsize=symsize, xticklen=xticklen, $
    /nohat, errthick=errthick
  legend,'nperbin = 10000',box=0, charsize=1
  multiplot,/reset

  endplot


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now look at histograms in the fixed bins of z
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  begplot, name=plot_dir + 'zhist_byzbin.ps', /color
  
  if n_elements(s) eq 0 then begin 
      self->zbin, bcg, zbin, s, low, high, num
  endif 

  !p.multi=[0,0,2]
  plot, zbin, num, psym=10, $
    thick=1, $
    xtitle = 'z', ytitle='Ncluster'
  oplot, zbin, num, psym=10, color=!blue, thick=1

  key = ''
  key = prompt_kbrd()
  plot, zbin, num, psym=10, $
    thick=1, $
    xtitle = 'z', ytitle='Ncluster', xrange=[0, 0.2]
  oplot, zbin, num, psym=10, color=!blue, thick=1
  key = prompt_kbrd()

  endplot
  !p.multi=0

  plot_dir = '~/plots/MaxBCG/ngals_hist_eachbin/'

  setupplot,'z'
  simpctable, rct, gct, bct

  device, set_res=[880,780]
  
  !p.multi = [0,4,4]

  !p.charsize = 1.3
  ocolor = !green
  ocolorold = c2i('red')
  nz = n_elements(low)

  iplot = 0
  iimage = 0

  for i=0L, nz-1 do begin

      ind = s[ low[i]:high[i] ]

      ngals = bcg[ s[ind] ].ngals
      title = 'z = '+ntostr(zbin[i])+' Number = '+ntostr(num[i])

      self->plothist_ngals, ngals, num[i], ocolor, title

      if i GT 0 then begin 
          self->plothist_ngals, ngalsold, num[i-1], ocolorold, /old

;          legend,$
;            ['This bin', 'Last bin'], line=[0,0], color=[ocolor,ocolorold],$
;            /right, box=0, charsize=1

      endif 
      ngalsold = ngals

      iplot = iplot + 1
      if (iplot MOD 16) eq 0 then begin 
          file = plot_dir + 'plot'+ntostr(iimage)+'.png'
          print,'Writing png file: ',file
          write_png, file, tvrd(), rct, gct, bct
          print,'ok'
          iimage = iimage+1
      endif 

      key = prompt_kbrd()
      if key eq 'q' then return

  endfor 

  !p.multi=0
  !p.charsize = 1.0

  setupplot, 'X'
  
end 


function maxbcg::schechter_func, x, nstar, xstar, alpha

  func = nstar * (x/xstar)^alpha * exp(-(x/xstar)) / xstar
  return,func

end 







pro maxbcg::calculate_zphot, bcg, neigh, bcg_photoz, bcg_photoz_err, $
          verbose=verbose


  if n_elements(bcg) eq 0 then begin 
      bcg = self->bcg_get()
  endif 
  if n_elements(neigh) eq 0 then begin 
      neigh = self->neighbor_get()
  endif 
  
  nbcg = n_elements(bcg)
  bcg_photoz = fltarr(nbcg)
  bcg_photoz_err = fltarr(nbcg)
  bcg_match = lonarr(nbcg)


  ustripes = bcg[ rem_dup(bcg.stripe) ].stripe
  nstripe = n_elements(ustripes)
  

  for ist = 0L, nstripe-1 do begin 

      stripe = ustripes[ist]

      print,'//////////////////////////////////////////////////////////////'
      print,'Doing stripe: ',stripe


      wneighst = where( neigh.stripe eq stripe, nneighst )
      wbcgst   = where( bcg.stripe   eq stripe, nbcgst)

      ;; histogram the neighbor bcgid
      print
      print,'Histogramming neighbor bcg_id'
      h = histogram(neigh[wneighst].bcg_id, min=0, rev=rev)

      print
      print,'Mathing to bcg and calculating photozs'
      print

      if keyword_set(verbose) then begin 
          print,$
            'bcg_id','nmatch','Ngals','bcgmatch','z','photoz','photozerr',$
            format='(A10,A10,A10,A10,A10,A10,A12)'
      endif 

      for iibcg=0L, nbcgst-1 do begin 

          ibcg = wbcgst[iibcg]

          bcg_id = bcg[ibcg].id
     
          ;; get neighbors
          if rev[bcg_id] ne rev[bcg_id+1] then begin 


              wmatch = rev[ rev[bcg_id]:rev[bcg_id+1]-1 ]
              wmatch = wneighst[wmatch]
              nmatch = n_elements(wmatch)


              ;; The bcg
              bcgmatch = $
                where(neigh[wmatch].photoid eq bcg[ibcg].photoid, nbcgmatch)
              if nbcgmatch ne 0 then bcgmatch = wmatch[bcgmatch]


              ;; Good photozs
              wkeep = where(neigh[wmatch].photoz GT 0.0, nkeep)
              

              
              if nkeep ne 0 then begin 
                  
                  wkeep = wmatch[wkeep]
                  
                  
                  wmom, $
                    neigh[wkeep].photoz, neigh[wkeep].photoz_err,$
                    photoz_mean, photoz_sig, photoz_err
                  
                  
                  
                  

                  xi = arrscl(findgen(100), $
                              photoz_mean-3.5*photoz_err, $
                              photoz_mean+3.5*photoz_err)
                  gp = gaussprob(xi, photoz_mean, photoz_err)
                  yrange = [0, max(gp)]
                  
                  pz = [photoz_mean, neigh[wkeep].photoz]
                  pzerr = [photoz_err, neigh[wkeep].photoz_err]
                  plotgauss, $
                    neigh[wkeep].photoz, neigh[wkeep].photoz_err, $
                    yrange=yrange, comb_xi=comb_xi, comb_gauss=comb_gauss
                  
                  oplot, comb_xi, comb_gauss, color=!dodgerBlue, thick=2
                  
                  plotgauss, photoz_mean, photoz_err,$
                             color=!green,thick=2, /overplot
                  
                  if nbcgmatch ne 0 then begin 
                      if neigh[bcgmatch].photoz GT 0 then begin 
                          plotgauss, $
                            neigh[bcgmatch].photoz, $
                            neigh[bcgmatch].photoz_err,$
                            color=c2i('red'),thick=2, line=2, /overplot
                      endif 
                  endif 
                  
                  oplot,$
                    [bcg[ibcg].z, bcg[ibcg].z],[0,1000],color=!yellow,thick=5
                  
;                  xi = arrscl(findgen(100), $
;                              photoz_mean-3.5*photoz_err, $
;                              photoz_mean+3.5*photoz_err)
;                  gp = gaussprob(xi, photoz_mean, photoz_err)
;                  gp = gp/max(gp)
;                  oplot, xi, gp, color=!green,thick=2
                  
              endif else begin 
                  erase
                  photoz_mean = -9999.0
                  photoz_err = 9999.0
              endelse 
                  
              if keyword_set(verbose) then begin 



                  print,$
                    strn(bcg_id),$
                    strn(nmatch),$
                    strn(bcg[ibcg].ngals),$
                    strn(bcgmatch),$
                    strn(bcg[ibcg].z),$
                    strn(photoz_mean),$
                    strn(photoz_err),$
                    format='(A10,A10,A10,A10,A10,A10,A12)'
                  
                  if nmatch ne bcg[ibcg].ngals then begin 
                      print,' ----  Ngals different'
                  endif 
                  
                  if nbcgmatch eq 0 then begin 
                      w2 = where(neigh.photoid eq bcg[ibcg].photoid)
                      print,'    w2 = ',w2
                  endif 
              
              
                  key = prompt_kbrd()
                  if key eq 'q' then return
              endif 

              bcg_match[ibcg] = bcgmatch
              bcg_photoz[ibcg] = photoz_mean
              bcg_photoz_err[ibcg] = photoz_err

              
          endif 
      endfor 

  end 



end 






function maxbcg::calc_loglike, gmr, rmi, emgauss, loglike_all=loglike_all

  COMMON calc_loglike_block, onesz

  if n_elements(onesz) eq 0 then onesz = replicate(1d, n_elements(emgauss))

  nclr = n_elements(gmr)
  if nclr eq 1 then begin 

      x = gmr - emgauss.cen[0]
      y = rmi - emgauss.cen[1]

      loglike = $
        x^2*emgauss.invcov[0,0]   +$
        2*x*y*emgauss.invcov[0,1] + $
        y^2*emgauss.invcov[1,1] 

  endif else begin 

      ;; matrix operations to improve speed emmensely
      onesc = replicate(1d, nclr)
      x = onesz#gmr - emgauss.cen[0]#onesc
      y = onesz#rmi - emgauss.cen[1]#onesc
      invcov00 = emgauss.invcov[0,0]#onesc
      invcov01 = emgauss.invcov[0,1]#onesc
      invcov11 = emgauss.invcov[1,1]#onesc
      
      loglike_all = x^2*invcov00  + 2*x*y*invcov01 + y^2*invcov11

      loglike = total(loglike_all, 2)

  endelse 

  return,loglike

end 


function maxbcg::zphot_em_file
  return,'~/benMaxBCG/new_photoz.st'
end 
pro maxbcg::zphot_em_write, struct
  file = self->zphot_em_file()
  print,'Writing to file: ',file
  write_idlstruct, struct, file
end 
function maxbcg::zphot_em_read
  return,read_idlstruct(self->zphot_em_file())
end 
pro maxbcg::calculate_zphot_em, bcg, neigh, spec, phz_struct, $
  verbose=verbose

  fm = obj_new('fastem')

  pg=obj_new('postgres')
  if n_elements(bcg) eq 0 then begin 
;      bcg = self->bcg_get()
      query = 'select photoid, bcg_id, ngals, z, gmr, rmi from bcg'
      print,query
      bcg = pg->query(query)
  endif 
  if n_elements(neigh) eq 0 then begin 
;      neigh = self->neighbor_get()
      query = 'select photoid, bcg_id, gmr, rmi from bcg_neighbor'
      print,query
      neigh = pg->query(query)
  endif 
  if n_elements(spec) eq 0 then begin 
      query = 'select match_photoid as photoid, z from specgal'
      print,query
      spec = pg->query(query)

      rmd = rem_dup(spec.photoid)
      spec = spec[rmd]
  endif 

  ;; matches between objects and spectroscopy
  neigh_specmatch = matchint(neigh.photoid, spec.photoid) 
  bcg_specmatch   = matchint(bcg.photoid, spec.photoid) 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set EM algorithm gaussians.  Crude spacing is 0.00434344
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ncrudez = 100
  nfinez = 100

  crude_zmin = 0.02
  crude_zmax = 0.45
  crudez = arrscl(findgen(ncrudez), crude_zmin, crude_zmax)

  emgauss = fm->redspec_interpolate(crudez)

  lognorm = -0.5*alog(2*!pi*sqrt(emgauss.detcov))


  random_width = 0.004*2        ; 0.004 both ways

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up outputs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nbcg = n_elements(bcg)
  defval = -9999.0
  phz_struct = {ben_photoz: defval, $
                photoz: defval, bcg_photoz: defval, $
                mean_specz:defval, bcg_specz:defval}
  phz_struct = replicate(phz_struct, nbcg)
  phz_struct.ben_photoz = bcg.z


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; histogram the neighbor bcgid
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Histogramming neighbor bcg_id'
  h = histogram(neigh.bcg_id, min=0, rev=rev)

  print
  print,'Matching to bcg and calculating photozs'
  print
  
;  if keyword_set(verbose) then begin 
;      print,$
;        'bcg_id','nmatch','Ngals','bcgmatch','z','photoz','photozerr',$
;        format='(A10,A10,A10,A10,A10,A10,A12)'
;  endif 

  for ibcg=0L, nbcg-1 do begin 

      bcg_id = bcg[ibcg].bcg_id
     
      ;; get neighbors
      if rev[bcg_id] ne rev[bcg_id+1] then begin 

          wmatch = rev[ rev[bcg_id]:rev[bcg_id+1]-1 ]
          nmatch = n_elements(wmatch)
          
          ;; Now find the maximum likelihood redshift.  Start on 
          ;; A crude grid, then zoom in.


          tloglike = $
            self->calc_loglike(neigh[wmatch].gmr, neigh[wmatch].rmi, emgauss,$
                               loglike_all=loglike_all)
          tbcg_loglike = $
            self->calc_loglike(bcg[ibcg].gmr, bcg[ibcg].rmi, emgauss)


          ;; normalization
          tloglike = -0.5*tloglike + nmatch*lognorm
          tbcg_loglike = -0.5*tbcg_loglike + lognorm

          ;; linear for plotting only
          like = exp( tloglike )
          bcg_like = exp( tbcg_loglike )

          like = like/qgauss(like, crudez, 100)
          bcg_like = bcg_like/qgauss(bcg_like, crudez, 100)



          maxlike = max(tloglike, wmax)
          tphotoz = crudez[wmax]
          maxlike = max(tbcg_loglike, wmax)
          tbcg_photoz = crudez[wmax]
          

          ;; refine overall photoz
          tphotoz = tphotoz + random_width*(randomu(seed) - 0.5)
          fine_zmin = tphotoz - 0.02 
          fine_zmax = tphotoz + 0.02
          finez = arrscl(findgen(nfinez), fine_zmin, fine_zmax)
          emgauss2 = fm->redspec_interpolate(finez)

          loglike = $
            self->calc_loglike(neigh[wmatch].gmr, neigh[wmatch].rmi, emgauss2)

          loglike = -0.5*loglike - nmatch*0.5*alog(2*!pi*sqrt(emgauss2.detcov))

          maxlike = max(loglike, wmax)
          photoz = finez[wmax]
          
          ;; refine bcg photoz
          tbcg_photoz = tbcg_photoz + random_width*(randomu(seed) - 0.5)
          fine_zmin = tbcg_photoz - 0.02
          fine_zmax = tbcg_photoz + 0.02
          finez = arrscl(findgen(nfinez), fine_zmin, fine_zmax)
          emgauss2 = fm->redspec_interpolate(finez)

          bcg_loglike = $
            self->calc_loglike(bcg[ibcg].gmr, bcg[ibcg].rmi, emgauss2)
          bcg_loglike = -0.5*bcg_loglike -0.5*alog(2*!pi*sqrt(emgauss2.detcov))

          maxlike = max(bcg_loglike, wmax)
          bcg_photoz = finez[wmax]


          phz_struct[ibcg].photoz = photoz
          phz_struct[ibcg].bcg_photoz = bcg_photoz


          emgauss2 = 0


          ;; spectroscopic redshifts
          if bcg_specmatch[ibcg] ne -1 then begin 
              phz_struct[ibcg].bcg_specz = $
                spec[bcg_specmatch[ibcg]].z
          endif 
          w=where(neigh_specmatch[wmatch] ne -1, nw)
          if nw ne 0 then begin 
              phz_struct[ibcg].mean_specz = $
                mean(spec[neigh_specmatch[wmatch[w]]].z)
          endif 
          if keyword_set(verbose) then begin 

              print,'BCG z: '+ntostr(bcg[ibcg].z)+' nmatch: '+ntostr(nmatch)+$
                    ' bcgspec: '+ntostr(phz_struct[ibcg].bcg_specz)+$
                    ' meanspec: '+ntostr(phz_struct[ibcg].mean_specz)+ $
                    ' nspec: '+ntostr(nw)+$
                    ' ngals: '+ntostr(nmatch)

              like = like/qgauss(like, crudez, 100)
              bcg_like = bcg_like/qgauss(bcg_like, crudez, 100)

              !p.multi = [0,0,2]
              
              plot,crudez, tloglike, psym=3
              
              yrange = [0, max([max(like),max(bcg_like)])]
              plot,crudez, like, yrange=yrange
              oplot, crudez, bcg_like, color=!green
              oplot, [bcg[ibcg].z,bcg[ibcg].z],[-1.e6,1.e6], color=c2i('red')

              ;; photoz/spec from neighbors
              oplot, [photoz,photoz],[-1.e6,1.e6]
              oplot, $
                [phz_struct[ibcg].mean_specz, phz_struct[ibcg].mean_specz], $
                [-1.e6,1.e6], line=2

;              for i=0L, nmatch-1 do begin 
                  
;                  tloglike = -0.5*loglike_all[*,i] + lognorm
;                  tlike = reform(tloglike)
;                  oplot, crudez, exp(tlike)

;                  if i eq 0 then totlike = tloglike else totlike=totlike + tloglike
;              endfor 
;              totlike2 = exp(totlike)
;              totlike2 = totlike2/qgauss(totlike2, crudez, 100)
;              oplot, crudez, totlike2, color=!cyan

              ;; photoz/specz from bcg
              oplot, [bcg_photoz,bcg_photoz],[-1.e6,1.e6], color=!green
              oplot, $
                [phz_struct[ibcg].bcg_specz, phz_struct[ibcg].bcg_specz], $
                [-1.e6,1.e6], line=1, color=!green

              legend, $
                ['ben z','mean photoz','mean specz','bcg photoz', 'bcg specz'],$
                line = [0,0,2,0,1], color=c2i(['red','white','white','green','green']),$
                /right, box=0, charsize=1

              key =prompt_kbrd()
              if key eq 'q' then begin 
                  obj_destroy,fm
                  return
              endif 
          endif else begin 
              if ((ibcg+1) MOD 1000) eq 0 then $
                print,ntostr(ibcg+1)+'/'+ntostr(nbcg)
          endelse 

    
      endif 
  endfor 
  print

  self->zphot_em_write, phz_struct
  obj_destroy,fm

end 







pro maxbcg::compare_spec_photoz_plot, z1, z2, range, _extra=_extra

  w = where(z1 GT range[0] AND $
            z1 LT range[1] AND $
            z2 GT range[0] AND $
            z2 LT range[1],nw)

  ploth, z1[w], z2[w], $
         xrange=range, yrange=range, _extra=_extra
  oplot, [0,1],[0,1],color=c2i('red')

  nperbin = nw/10
  binner_bynum, z1[w], z2[w], nperbin, $
                xb, yb, yberr, ybinned_sdev=ybsdev
  oploterror, xb, yb, ybsdev, color=!dodgerBlue, errc=!dodgerBlue


end 
pro maxbcg::compare_spec_photoz, bcg, phz

  pg=obj_new("postgres")
  if n_elements(bcg) eq 0 then begin 
      query = 'select photoid, bcg_id, ngals, z, gmr, rmi from bcg'
      print,query
      bcg = pg->query(query)
  endif 
  if n_elements(phz) eq 0 then begin 
      phz = self->zphot_em_read()
  endif 

  !p.multi=[0,2,2]

  range = [0, 0.4]

  ;; Ben's photozs compared to bcg_specz
  self->compare_spec_photoz_plot, phz.bcg_specz, phz.ben_photoz, range, $
    xtitle = 'BCG z', ytitle = 'Ben photoz'

  ;; Ben's photozs compared to mean_specz
  self->compare_spec_photoz_plot, phz.mean_specz, phz.ben_photoz, range, $
    xtitle = 'Mean z', ytitle = 'Ben photoz'

  ;; Ben's photozs compared to bcg_photoz
  self->compare_spec_photoz_plot, phz.ben_photoz, phz.bcg_photoz, range,$
    xtitle = 'Ben photoz', ytitle = 'BCG photoz'

  ;; Ben's photozs compared to mean_photoz
  self->compare_spec_photoz_plot, phz.ben_photoz, phz.mean_photoz, range,$
    xtitle = 'Ben photoz', ytitle = 'Mean photoz'

  write_png, '~/plots/MaxBCG/new_photoz/ben_compare.png',tvrd(/true)
  
  key = prompt_kbrd("hit a key")


  ;; BCG photozs compared to bcg_specz
  self->compare_spec_photoz_plot, phz.bcg_specz, phz.bcg_photoz, range, $
    xtitle = 'BCG z', ytitle = 'BCG photoz'

  ;; BCG photozs compared to mean_specz
  self->compare_spec_photoz_plot, phz.mean_specz, phz.bcg_photoz, range,$
    xtitle = 'Mean z', ytitle = 'BCG photoz'


  ;; Mean photozs compared to bcg_specz
  self->compare_spec_photoz_plot, phz.bcg_specz, phz.mean_photoz, range, $
    xtitle = 'BCG z', ytitle = 'Mean photoz'

  ;; Mean photozs compared to mean_specz
  self->compare_spec_photoz_plot, phz.mean_specz, phz.mean_photoz, range, $
    xtitle = 'Mean z', ytitle = 'Mean photoz'

  write_png, '~/plots/MaxBCG/new_photoz/new_compare.png',tvrd(/true)

  !p.multi = 0

  

end 





;; Add random numbers to z and plot histogram
pro maxbcg::plotz, z, _extra=_extra

  sig = 0.001
  print,'Adding random noise of: ',sig
  nz = n_elements(z)
  z2 = z + sig*randomn(seed, nz)

  plothist, z2, _extra=_extra

end 

















pro maxbcg::find_dup, neigh

  print,'Sorting'
  s = sort(neigh.photoid)
  photoid_old = -1LL

  print,'Finding duplicates'
  nn = n_elements(neigh)
  for i=0L, nn-1 do begin 
      photoid = neigh[s[i]].photoid
      if photoid eq photoid_old then begin 
          print,'photoid dup = ',photoid
      endif 
      photoid_old = photoid
  endfor 

end 



;; I never finished this code.  This was an attempt to do the matching, since
;; the table join was going to take two days.

function maxbcg::match_neigh2adatc_extra_query, run
  
  photoid_range = sdss_photoid_range(run)
  query = 'SELECT photoid, counts_model, reddening FROM adatc_extra WHERE photoid BETWEEN '+ntostr(photoid_range[0])+' AND '+ntostr(photoid_range[1])
  return,query
end 
pro maxbcg::match_neigh2adatc_extra, neigh

  if n_elements(neigh) eq 0 then neigh = self->neighbor_get()

  print
  print,'Extracting id info from photoids'
  photoid_extract, neigh.photoid, runs, reruns, camcols, fields, ids

;  print
;  print,'Histogramming runs'
;  hrun = histogram(runs, rev=revrun)
;  nh = n_elements(h)

  ;; getting sdss_hist
  print,'Getting sdss_hist'
  idstruct = sdss_histid(runs,reruns,camcols,fields)

  stop

;  pruns = idstruct.runs
;  for ri=0L, idstruct.nruns-1 do begin 
;      runStr = ntostr( (*pruns)[ri].run )
;      preruns = (*pruns)[ri].reruns
      




  ;; Deal with each run
  pg=obj_new('postgres')
  for i=0L, nh-1 do begin 

      if revrun[i] ne revrun[i+1] then begin 

          wrun = revrun[ revrun[i]:revrun[i+1]-1 ]
          run = runs[wrun[0]]
          print,'Working on run: ',run

          query = self->match_neigh2adatc_extra_query(run)
          print,query

          struct = pg->query(query)

          ;; The problem here is there are many duplicates in the
          ;; neighbor list

          print,'Matching'



          stop

      endif 
 
  endfor 
  

end 















pro maxbcg::simple_ngals200_bins, nbin, lowlim, highlim
  CASE nbin OF
      9: begin 
          lowlim =[10,11,12,14,20,29,40,55,80]
          highlim=[10,11,13,19,28,39,54,79,200]
      end 
      16: begin 
          lowlim =[3,4,5,6,7,8,9,10,11,12,14,20,29,40,55,80]
          highlim=[3,4,5,6,7,8,9,10,11,13,19,28,39,54,79,200]
      end 
      ELSE: message,'only have nbin=9,16 for now'
  endcase 
end 

; area of this catalog. Run stand alone
; calculate_area code in sdsspixIDL/app
; note, dr406 now uses the "BOUND" mask
function maxbcg::area, radians=radians
  CASE self->catalog() OF
      'dr3': area = 6467.284785d
      'dr4plus': area = 7541.709773d
      'dr406': area = 7398.233579d
      ELSE: message,'Do not have an area for catalog '+self->catalog()
  end 
  if keyword_set(radians) then begin 
      area = area * (!dpi/180d)^2
  endif 
  return, area
end 

; Take the input xhist, yhist for a redshift
; histogram and return the number density in z bins assuming the
; concordance model
function maxbcg::zdensity, xhist, yhist, comoving=comoving
 

  area = self->area(/rad) ;steradians

  factor = area/(4d*!dpi)

  c = obj_new('cosmology')

  nhist = n_elements(xhist)

  ;; Assume equal sized bins and xhist is
  ;; the center of the bin
  dz = xhist[1]-xhist[0]

  for i=0L, nhist-1 do begin 

      tvolume = c->volume(xhist[i]-dz/2.0, xhist[i]+dz/2.0,$
                              omega_m=0.27, comoving=comoving)

      add_arrval, tvolume, volume

  endfor 


  density = yhist/volume*factor

  obj_destroy, c

  return, density

end 

;; in bins of ngals200.
pro maxbcg::plot_zdensity, nbin, bcg=bcg, color=color

  self->simple_ngals200_bins, nbin, lowlim, highlim

  dir = '~/plots/maxbcg/'
  file = 'maxbcg_'+self->catalog()+'_number_density_z_'+$
      'ngals200_nbin'+strn(nbin,len=2,padchar='0')

  file=file+'.eps'
  if keyword_set(color) then file=repstr(file, '.eps', '_color.eps')

  file = concat_dir(dir, file)
  begplot, name=file, color=color, /encapsulated

;      !p.charsize=1

  ;; read in the objects
  if n_elements(bcg) eq 0 then bcg = self->get()

  ngals_min = min([lowlim,highlim], max=ngals_max)
  use = where(bcg.ngals200 GE ngals_min AND bcg.ngals200 LE ngals_max, nbcg)

  if keyword_set(color) then begin
      simpctable, colorlist=clist
  endif else begin
      loadct, 0
      clist = reverse(arrscl(findgen(nbin), 30, 200))
      clist = [0, clist]
  endelse

  zmin = 0.1
  zmax = 0.3
  xrange = [0.08, 0.5]

  bin = 0.015

  plothist, bcg[use].photoz_cts, xhist, yhist, min=zmin, max=zmax, bin=bin, /noplot

  comoving=0
  density = self->zdensity(xhist, yhist, comoving=comoving)

  yhisterr = sqrt(yhist > 0)
  densityerr = density/yhist*yhisterr

  ;delta = textoidl('\Delta')
  ytitle = textoidl('\DeltaN/\DeltaV [h^3 Mpc^{-3}]')
  ;ytitle = delta+'N/'+delta+'V [ h!U3!N Mpc!U-3!N ]'


  
  yrange = [1.e-8, max(density)]

  aploterror, !gratio, xhist, density, densityerr, psym=10, $
    xtitle='z', ytitle=ytitle, hat=0, $
    xrange=xrange,xstyle=3, yrange=yrange, ystyle=3, /ylog, ytickf='loglabels'

  mess = 'all'
  colors = !p.color

  for i=0L, nbin-1 do begin 
      
      w = where(bcg[use].ngals200 GE lowlim[i] AND $
                bcg[use].ngals200 LE highlim[i], nw)

      if nw GT 10 then begin 
          
          w = use[w]
          plothist, $
            bcg[w].photoz_cts, xhist, yhist, min=zmin, max=zmax, bin=bin, /noplot 
          
          density = self->zdensity(xhist, yhist, comoving=comoving)
          
          yhisterr = sqrt(yhist > 0)
          densityerr = density/yhist*yhisterr
          oploterror, xhist, density, densityerr, psym=10, $
            color=clist[i+1], errc=clist[i+1], hat=0
;              oplot, xhist, density, psym=10, color=clist[i+1]
          
          add_arrval, self->ngals200_string(lowlim[i], highlim[i]), mess
          add_arrval, clist[i+1], colors
          
;              key = prompt_kbrd('hit a key')
          
      endif 
      
  endfor 

  legend, mess, /right, box=0, line=0, colors=colors, charsize=1.0

  endplot, /trim_bbox


end 









;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Input catalogs for Juan's auto-correlation code.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function maxbcg::juan_dir, dtype

  CASE strlowcase(dtype) OF
      'input': begin 
          dir = esheldon_config('maxbcg_dir')
          dir = concat_dir(dir, 'juan_input')
      end 
      'output': begin 
          dir = '~mm1330/smallscale/maxbcg'
      end 
  endcase 

  return, dir
end 

function maxbcg::juan_file, type, randnum=randnum

  dir = self->juan_dir(type)

  CASE strlowcase(type) OF 
      'output': begin 
          file = 'maxbcg-20-'
          if n_elements(randnum) ne 0 then begin 
              file = file + 'random.fits'
          endif else begin 
              file = file + 'data.fits'
          endelse 
      end 
      'input': begin 
          file = 'maxbcg_lum_input'

          nrand = n_elements(randnum)
          if nrand ne 0 then begin 
              for i=0L, nrand-1 do begin
                  taddstr = '_random_'+strn(randnum[i],length=2,padchar='0')
                  add_arrval, taddstr, addstr
              endfor 
              file = file + addstr
          endif 
          file = file + '.fit'
      end 
      ELSE: message,'Unsupported type: '+ntostr(type)
  endcase 
  file = concat_dir(dir, file)
  return,file
end 

function maxbcg::juan_read, type, randnum=randnum

  on_error, 2
  if n_elements(type) eq 0 then begin 
      print,'-Syntax: st = mb->juan_read(type, randnum=)'
      print,'type=input|output'
      message,'Halting'
  endif 
  file = self->juan_file(type, randnum=randnum)
  print
  print,'Reading file: ',file
  struct = mrdfits(file, 1)
  return, struct
end 









function maxbcg::juan_struct, num, random=random
  if keyword_set(random) then begin 
      struct = { $
                 bcg_id: 0L, $
                 ra: 0d, $
                 dec: 0d, $
                 z: 0.0 $
                 }
  endif else begin 
      struct = { $
                 bcg_id: 0L, $
                 ra:0d, $
                 dec:0d, $
                 z:0.0,$
                 bcg_spec_z: 0.0, $
                 ngals: 0, $
                 ngals200: 0,$
                 bcg_ilum: 0.0, $
                 ilum200: 0.0 $
               }
  endelse 
  if n_elements(num) ne 0 then begin 
      struct = replicate(struct, num[0])
  endif 
  return,struct
end 












function maxbcg::juan_cut, clambda, ceta, z, nkeep

  ;; run through 1Mpc cut that Ben uses. This is slightly
  ;; different than Ben's, so we should do it on both

  rMpc = 1.0
  maxangle = rMpc/angdist_lambda(z, omega=0.27)*180.0/!pi


  apply_pixel_mask, $
    clambda, ceta, masked, unmasked, $
    maxangle = maxangle, /basic

  if unmasked[0] eq -1 then nkeep = 0 else nkeep = n_elements(unmasked)
  
  return,unmasked


end 

function maxbcg::juan_genrandz, rlam, reta, bcg

  nrand = n_elements(rlam)
  nbcg = n_elements(bcg)

  rbcgid = long( nbcg*randomu(seed, nrand) )
  zrand = bcg[rbcgid].photoz_cts
  bcg_id = bcg[rbcgid].bcg_id

  print
  print,'Making edge cut'
  rkeep = self->juan_cut(rlam, reta, zrand, nrkeep)

  print,'Kept: ',nrkeep
  struct = self->juan_struct(nrkeep, /random)

  csurvey2eq, rlam, reta, rra, rdec
  rra = rra[rkeep]
  rdec = rdec[rkeep]
  zrand = zrand[rkeep]
  bcg_id = bcg_id[rkeep]
  
  struct.bcg_id = bcg_id
  struct.ra = rra
  struct.dec = rdec
  struct.z = zrand

  return,struct

end 

pro maxbcg::juan_input, bcg=bcg, keep=keep, randnum=randnum

  if n_elements(bcg) eq 0 then begin 
      bcg=self->get(columns=['bcg_id', $
                             'ra','dec','photoz_cts','bcg_spec_z',$
                             'ngals','ngals200','bcg_ilum','ilum200'])
  endif 

  nbcg = n_elements(bcg)

  nrand = n_elements(randnum)
  if nrand eq 0 then begin 

      eq2csurvey, bcg.ra, bcg.dec, clambda, ceta
      print
      print,'Making edge cut'
      keep = self->juan_cut(clambda, ceta, bcg.photoz_cts)


      nkeep = n_elements(keep)
      print,'Kept: ',nkeep

      outst = self->juan_struct(nkeep)

      outst.bcg_id = bcg[keep].bcg_id

      outst.ra = bcg[keep].ra
      outst.dec = bcg[keep].dec

      outst.z = bcg[keep].photoz_cts
      outst.bcg_spec_z = bcg[keep].bcg_spec_z

      outst.ngals    = bcg[keep].ngals
      outst.ngals200 = bcg[keep].ngals200
      outst.bcg_ilum = bcg[keep].bcg_ilum
      outst.ilum200  = bcg[keep].ilum200

      outfile = self->juan_file('input')

      hdr=['END']
      sxaddpar, hdr, 'catalog', self->catalog()
      print
      print,'Writing to file: ',outfile
      mwrfits, outst, outfile, hdr, /create

  endif else begin 

      for i=0L, nrand-1 do begin

          rand = self->random_read(randnum[i])
          rstruct = self->juan_genrandz(rand.clambda, rand.ceta, bcg)
          routfile = self->juan_file('input',randnum=randnum[i])

          hdr=['END']
          sxaddpar, hdr, 'catalog', self->catalog()
 
          print
          print,'Writing random file: ',routfile
          mwrfits, rstruct, routfile, hdr, /create

      endfor 

  endelse 


end 








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Stuff for Matt Becker
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function maxbcg::lssmatch_dir
  dir = self->catdir()
  dir = concat_dir(dir, 'lssmatch')
  return,dir
end 
function maxbcg::lssmatch_file, random=random
  dir = self->lssmatch_dir()
  file = 'bcg-match-lss-dr4plus'
  
  if keyword_set(random) then begin 
      file = file + '-random'
  endif 
  file = file + '.fits'
  file = concat_dir(dir, file)
  return, file
end 

function maxbcg::lss_spec_file, random=random
  dir = self->lssmatch_dir()
  file = 'lss-dr4plus'
  
  if keyword_set(random) then begin 
      file = file + '-random'
  endif 
  file = file + '.fits'
  file = concat_dir(dir, file)
  return, file
end 


function maxbcg::lssmatch_read, random=random, rows=rows, columns=columns
  file = self->lssmatch_file(random=random)
  print,'Reading file: ',file
  struct = mrdfits(file, 1, rows=rows, columns=columns)
  return, struct
end 
function maxbcg::lss_spec_read, random=random, rows=rows, columns=columns
  file = self->lss_spec_file(random=random)
  print,'Reading file: ',file
  struct = mrdfits(file, 1, rows=rows, columns=columns)
  return, struct
end 



function maxbcg::lssmatch_fgotcut
  return, 0.5
end 
function maxbcg::lss_zrange
  return, [0.002, 0.6]
end 
function maxbcg::match2lss, bcg

  v = obj_new('vagc')
  lss = v->read('lss', 'lss_index', columns=['ra','dec','z','fgotmain'])

  zrange = self->lss_zrange()
  w = where(lss.fgotmain GT self->lssmatch_fgotcut() AND $
            lss.z GT zrange[0] AND  lss.z LT zrange[1])
  lss = lss[w]

  angle = 5d/3600d*!dpi/180d
  htm_match, $
    bcg.ra, bcg.dec, lss.ra, lss.dec, angle, $
    mbcg, mlss, d12, $
    maxmatch=1

  match_struct = {mbcg: 0L, mlss: 0L, z: 0.0}

  nmatch = n_elements(mbcg)
  match_struct = replicate(match_struct, nmatch)

  match_struct.mbcg = mbcg
  match_struct.mlss = mlss
  match_struct.z = lss[mlss].z

  obj_destroy, v
  return, match_struct

end 

pro maxbcg::write_lssmatch

  bcg = self->get()
  nbcg = n_elements(bcg)
  print
  print,'Matching to lss'
  match_struct = self->match2lss(bcg)

  bcgmatch = bcg[ match_struct.mbcg ]
  bcgmatch.bcg_spec_z = match_struct.z

  print,'Matched '+ntostr(n_elements(bcgmatch))+'/'+ntostr(nbcg)

  ;; mask cut
  maskfile = self->maxbcg_maskfile()
  print
  print,'Running through mask: ',maskfile
  apply_pixel_mask, bcgmatch.clambda, bcgmatch.ceta, m, um, $
    maskfile=maskfile

  bcgmatch = bcgmatch[um]

  print,'Kept '+ntostr(n_elements(bcgmatch))+'/'+ntostr(nbcg)

  ;; Now Ben's 1Mpc edge cut.  We must use the *wrong* redshifts
  ;; for the edgecut!  The photoz.
  print
  print,'Running through edge mask'

  print
  print,'Doing edge cuts with maxbcg photoz'

  rMpc = 1.0
  
  maxAngle = rMpc/angdist_lambda(bcgmatch.z, omega=0.27)*180.0/!pi
  apply_pixel_mask, $
    bcgmatch.clambda, bcgmatch.ceta, m, um, maxangle = maxangle, /basic

  bcgmatch = bcgmatch[um]

  print,'Finally kept '+ntostr(n_elements(bcgmatch))+'/'+ntostr(nbcg)

  outfile = self->lssmatch_file()
  print
  print,'Writing file: ',outfile
  mwrfits, bcgmatch, outfile, /create

end 

pro maxbcg::write_lss_spec

  v = obj_new('vagc')

  lss = v->read('lss', 'lss_index', columns=['ra','dec','z','fgotmain'])

  zrange = self->lss_zrange()
  w = where(lss.fgotmain GT self->lssmatch_fgotcut() AND $
            lss.z GT 0.0 AND lss.z GT zrange[0] AND lss.z LT zrange[1], nlss)
  lss = lss[w]

  file = self->lss_spec_file()
  print
  print,'Writing lss spec file: ',file
  mwrfits, lss, file, /create

end 

;; This is a set of random points that have the same spatial and
;; redshift distribution as the spec
pro maxbcg::write_lss_spec_random

  v = obj_new('vagc')

  
  ;; Read as many points as we need
  randf = v->file('lss', 'lss_random', randnum=lindgen(10))
  rand = mrdfits_multi(randf)

  geom = v->read('lss', 'lss_geometry')

  nrand = n_elements(rand)

  fgotmain = geom[rand.ilss].fgotmain

  ;; completeness cut
  print
  print,'Cutting on fgot'
  w = where(fgotmain GT self->lssmatch_fgotcut(), nkeep)
  rand = rand[w]
  fgotmain = fgotmain[w]

  print,'Kept '+ntostr(nkeep)+'/'+ntostr(nrand)+' after fgot cuts'

  ;; assign random redshifts from actual catalog
  print
  print,'Assigning redshifts'
  lss = self->lss_spec_read()
  
  nlss = n_elements(lss)
  rand_index = long( nlss*randomu(seed, nkeep) )

  z = lss[rand_index].z

  ;; add some noise with width that depends on redshift
  zwidth = 0.01*z
  z = z + zwidth*( randomu(seed, nkeep)-0.5 ) > 0.0

  outst = replicate({ra: 0d, dec:0d, z: 0.0, fgotmain: 0.0}, nkeep)
  
  outst.ra = rand.ra
  outst.dec = rand.dec
  outst.z = z
  outst.fgotmain = fgotmain

  file = self->lss_spec_file(/random)
  print
  print,'Writing to file: ',file
  mwrfits, outst, file, /create
end 

;; This is random bcgs; matt won't actually use these
pro maxbcg::write_lssmatch_random

  message,'We used all randoms in spec, so consider this'

  v = obj_new('vagc')
  
  ;; Read as many points as we need
  randf = v->file('lss', 'lss_random', randnum=[0,1,2,3])
  rand = mrdfits_multi(randf)

  geom = v->read('lss', 'lss_geometry')

  nrand = n_elements(rand)

  fgotmain = geom[rand.ilss].fgotmain

  ;; completeness cut
  print
  print,'Cutting on fgot'
  w = where(fgotmain GT self->lssmatch_fgotcut(), nkeep)
  rand = rand[w]

  print,'Kept '+ntostr(nkeep)+'/'+ntostr(nrand)+' after fgot cuts'

  ;; run through the bound mask
  maskfile = self->maxbcg_maskfile()
  eq2csurvey, rand.ra, rand.dec, rlam, reta
  print
  print,'Running through mask: ',maskfile
  apply_pixel_mask, rlam, reta, m, um, $
    maskfile=maskfile

  rand = rand[um]
  rlam = rlam[um]
  reta = reta[um]
  nkeep = n_elements(rand)
  print,'Kept '+ntostr(nkeep)+'/'+ntostr(nrand)+' after pixel mask'

  ;; assign random redshifts from actual catalog
  print
  print,'Assigning redshifts'
  bcg = self->lssmatch_read()


  nbcg = n_elements(bcg)
  rand_index = long( nbcg*randomu(seed, nkeep) )

  outst = replicate({ra: 0d, dec:0d, z: 0.0}, nkeep)
  
  outst.ra = rand.ra
  outst.dec = rand.dec
  outst.z = bcg[rand_index].bcg_spec_z
  zphot = bcg[rand_index].z

  ;; Now Ben's 1Mpc edge cut.  We must use the *wrong* redshifts
  ;; for the edgecut!  The photoz.
  print
  print,'Doing edge cuts with maxbcg photoz'

  rMpc = 1.0
  
  maxAngle = rMpc/angdist_lambda(zphot, omega=0.27)*180.0/!pi
  apply_pixel_mask, $
    rlam, reta, m, um, maxangle = maxangle, /basic

  outst = outst[um]
  nkeep = n_elements(outst)
  print,'Kept '+ntostr(nkeep)+'/'+ntostr(nrand)+' after edge cuts'

  outfile = self->lssmatch_file(/random)
  print
  print,'Writing to file: ',outfile
  mwrfits, outst, outfile, /create

end 

pro maxbcg::plot_lssmatch, dops=dops, dopng=dopng, spec=spec

  plot_dir = self->lensdir('plot',/base)
  if keyword_set(spec) then begin 
      real = self->lss_spec_read()
      rand = self->lss_spec_read(/random, rows=lindgen(1000000))
      pngfile = concat_dir(plot_dir,'lss_spec_positions.png')
      psfile  = concat_dir(plot_dir,'lss_spec_positions.ps')
  endif else begin 
      real = self->lssmatch_read()
      rand = self->lssmatch_read(/random, rows=lindgen(1000000))
      pngfile = concat_dir(plot_dir,'lssmatch_positions.png')
      psfile  = concat_dir(plot_dir,'lssmatch_positions.ps')
  endelse 

  eq2csurvey, real.ra, real.dec, clam, ceta
  eq2csurvey, rand.ra, rand.dec, rlam, reta

  delvarx, real, rand




  if keyword_set(dops) then begin 
      begplot, name=psfile, /color, /landscape
  endif 

  ytitle = textoidl('\eta_c')
  xtitle = textoidl('\lambda_c')
  xrange = [-70, 70]
  botyrange = [-40,50]
  topyrange = [120,170]
  oplotcolor = c2i('red')

  regions = self->ceta2region(ceta)
  rregions = self->ceta2region(reta)

  w1 = where(regions eq 1, nw1)
  w2 = where(regions eq 2, nw2)
  rw1 = where(rregions eq 1, nrw1)
  rw2 = where(rregions eq 2, nrw2)


  if n_elements(psym_lens) eq 0 then psym_lens=3
  if nrw1 ne 0 AND nrw2 ne 0 then begin 

      plottwo, $
        rlam[rw2], reta[rw2], $
        rlam[rw1], reta[rw1], $
        toppsym=3, botpsym=3, $
        botyrange=botyrange, topyrange=topyrange, $
        xrange=xrange, xstyle=1, $
        xtitle = xtitle, topytitle=ytitle, botytitle=ytitle, $
        frac1 = 0.25, $
        $
        xoplot = clam[w2], yoplot = ceta[w2], $
        oplotsym=3, oplotcolor=oplotcolor, $
        /ynozero, ystyle=1+2
      
      oplot,clam[w1],ceta[w1],psym=psym_lens,$
        symsize=symsize_lens,color=oplotcolor
  endif else begin 

      plot,rlam,reta,$
        /iso, $
        psym=3,$
        /ynozero, $
        title=title, xtitle=xtitle, ytitle=ytitle, xrange=xrange, xstyle=1
      
      oplot,clam,ceta,color=oplotcolor, psym=psym_lens, symsize=symsize_lens
  endelse 

  legend,$
    ['real','random'],$
    psym=[8,8],color=[oplotcolor,!p.color],/right,box=0,charsize=1.5

  if keyword_set(dops) then begin 
      endplot, /landfix
  endif else if keyword_set(dopng) then begin 
      print
      print,'Writing to file: ',pngfile
      write_png, pngfile, tvrd(/true)
  endif 
  

end 

pro maxbcg::plot_lssmatch_zdist, dopng=dopng

  bcg = self->lssmatch_read()
  rbcg = self->lssmatch_read(/random)

  charsize=2
  plothist, bcg.bcg_spec_z, bin=0.01,xstyle=3,/norm, $
    xtitle = 'z', ytitle='Number', charsize=charsize, $
    position=aspect(1.0/!gratio)
  plothist, rbcg.z, bin=0.01, color=!green, /norm, /overplot

  legend, ['Clusters', 'Random'], line=[0,0], color=[!p.color, !green], $
    box=0, /right

  plot_dir = self->lensdir('plot',/base)
  pngfile = concat_dir(plot_dir,'lssmatch_zhist.png')

  print
  print,'Writing to file: ',pngfile
  write_png, pngfile, tvrd(/true)


end 


pro maxbcg::icl_match, bcg=bcg

    if n_elements(bcg) eq 0 then bcg=self->get()
    dir = esheldon_config('maxbcg_dir')
    outfile=concat_dir(dir, 'icl_sample.fits')

    w=where(bcg.photoz_cts gt 0.1 and bcg.photoz_cts lt 0.3 and bcg.ngals200 ge 10,nw)

    sp=obj_new('sdss_postgres')
    t = sp->read_photoids($
        'adatc',$
        bcg[w].photoid,$
        columns=['m_e1_corr_h','m_e2_corr_h'])

    obj_destroy, sp

    newstruct=create_struct($
        bcg[0], $
        'e1', fltarr(5), 'e2', fltarr(5), $
        'aratio', fltarr(5), 'posangle', fltarr(5))
    newstruct = replicate(newstruct, nw)
    struct_assign, bcg[w], newstruct
    newstruct.e1 = t.m_e1_corr_h
    newstruct.e2 = t.m_e2_corr_h

    for i=0,4 do begin
        findabtheta, newstruct.e1[i], newstruct.e2[i], aratio, posangle
        newstruct.aratio[i]=aratio
        newstruct.posangle[i]=posangle
    endfor

    print,'Writing to file: ',outfile
    mwrfits, newstruct, outfile, /create
end


; For sarah bridle
function maxbcg::bridle_file
    if self->catalog() ne 'public' then message,'Public catalog only'

    outfile = self->catname()
    outfile = repstr(outfile, '.fit', '_evans.fit')
    return,outfile
end
pro maxbcg::make_bridle_maskflags, struct

    t=systime(1)

    if self->catalog() ne 'public' then message,'Public catalog only'

    outfile = self->bridle_file()
    print,'Writing file: ',outfile

    if n_elements(struct) eq 0 then begin
        struct = self->get()
    endif

    eq2csurvey, struct.ra, struct.dec, clambda, ceta

    maskfile = sdssidl_config('pixel_mask_princeton_basic')

    rmpcvals = 1+lindgen(25)*0.5
    nr=n_elements(rmpcvals)

    nrstr=ntostr(nr)
    add_tags, struct, ['rmax','maskflags'], ['fltarr('+nrstr+')','intarr('+nrstr+')'], newstruct
  
    maskflagvals = intarr(nr)
    for i=0L, nr-1 do begin

        rmpc = rmpcvals[i]
        print,'r = '+ntostr(rmpc)+' Mpc'

        ;; convert to degrees for each lens
        maxangle = rmpc/angdist(0.0, struct.z, omega_m=0.3)*180.0/!pi

        apply_pixel_mask, $
            clambda, ceta, masked, unmasked, maskflags, $
            maxangle=maxangle, maskfile=maskfile, $
            status=status

        if status ne 0 then message,'Failed'
        
        newstruct.rmax[i]=rmpc
        newstruct.maskflags[i] = maskflags
    endfor
  
    print,'Writing file: ',outfile
    mwrfits, newstruct, outfile, /create
    ptime,systime(1)-t

end

pro maxbcg::plot_bridle_stats

    if self->catalog() ne 'public' then message,'Public catalog only'

    file = self->bridle_file()
    psfile=repstr(file,'.fit','.eps')
    begplot,psfile,/encap
    t=mrdfits(file,1)

    ; plot fraction as a function of radius
    ntot=float(n_elements(t))
    nr=n_elements(t[0].rmax)
    for i=0L, nr-1 do begin
        w=where(t.maskflags[i] eq 0, ngood)
        tfrac = ngood/ntot
        print,tfrac
        add_arrval, tfrac, frac
    endfor

    pplot, t[0].rmax, frac, $
        xtitle=textoidl('r_{max} [h^{-1} Mpc]'), $
        ytitle='Unmasked Fraction', aspect=1

    endplot,/trim
end



;+
; NAME:
;   sdss_kcorrect
; PURPOSE:
;   calculate K-corrections for standard SDSS input
; CALLING SEQUENCE:
;   kcorrect= sdss_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                           calibobj=, tsobj=, flux=, band_shift=,$
;                           chi2=, rmaggies=, omaggies=, vname=, $
;                           oivar=, mass=, mtol=, absmag=, amivar=, $
;                           omega0=, omegal0= ])
; INPUTS:
;   redshift - [N] redshifts
;   calibobj - [N] photoop-style structure, containing:
;                  .PETROFLUX[5]
;                  .PETROFLUX_IVAR[5]
;                  .MODELFLUX[5]
;                  .MODELFLUX_IVAR[5]
;                  .PSFFLUX[5]
;                  .PSFFLUX_IVAR[5]
;                  .EXTINCTION[5]
;   tsobj - [N] opdb-style structure, containing:
;                  .PETROCOUNTS[5]
;                  .PETROCOUNTSERR[5]
;                  .COUNTS_MODEL[5]
;                  .COUNTS_MODELERR[5]
;                  .PSFCOUNTS[5]
;                  .PSFCOUNTSERR[5]
;                  .REDDENING[5]
;   nmgy, ivar - [5, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [5, N] asinh magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   flux - use this version of the fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'PETRO'] if tsobj or calibobj keywords are
;          used 
;   band_shift    - blueshift of bandpasses to apply (to get ^{z}b
;                   type bands) [default 0.]
;   vname - name of fit to use (defaults to 'default')
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
; OPTIONAL KEYWORDS:
;   lrg - do "luminous red galaxy" fit; this means changing vname to
;         'lrg1', using "model" fluxes, and ignoring the u-band
;         (setting ivar[0,*]=0); this uses a single template
;         appropriate for the SDSS Luminous Red Galaxy sample
; OUTPUTS:
;   kcorrect - [5, ngals] K-corrections in ugriz satisfying
;                m = M + DM(z) + K(z)
;              based on the best fit sum of templates
;   mtol - [5, ngals] current stellar mass-to-light ratios from model
;          in each band
;   mass - [ngals] total current stellar mass from model 
;   mets - [ngals] average metallicity in current stars 
;   intsfh - [ngals] total integrated star formation history
;   absmag - [5, ngals] absolute magnitude (for missing data, substitutes
;            model fit). (evolution correction *not* applied)
;   amivar - [5, ngals] inverse variance of absolute magnitude (for
;            missing data = 0)
; OPTIONAL OUTPUTS:
;   coeffs - [Nt, ngals] coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [5, ngals] reconstructed maggies from the fit (ugriz)
;   omaggies, oivar - [5, ngals] maggies and inverse variances used for fit
;                           (after extinction, AB correction, etc)  (ugriz)
;   b300 - [ngals] star-formation within last 300Myrs relative to average
;          star-formation rate
;   b1000 - [ngals] star-formation within last 1Gyrs relative to average
;           star-formation rate
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro which is almost always
;   just what you want. It keeps a version of rmatrix and zvals in
;   memory to save time, recalculating them each time you change
;   band_shift.
;
;   You must specify nmgy,ivar OR mag,err OR calibobj OR tsobj. If
;   nmgy or mag, make sure they are AB calibrated and Galactic
;   extinction corrected.
;
;   Uses sdss_to_maggies to convert tsobj or calibobj structure to
;   AB, Galactic extinction corrected maggies. Passes optional
;   argument "flux" to sdss_to_maggies to indicate which type of flux
;   to use. Note that if you get magnitudes like petroMag or modelMag
;   from the Catalog Archive Servers, these numbers are exactly like
;   the petroCounts and counts_model numbers in the tsObj structures.
;
;   For v4_0b templates and later, coefficients are in units of:
;     1 solar mass / (D/10pc)^2
;   That is, sum the coefficients and multiply by (D/10pc)^2 to get
;   masses. (In fact, for Omega0=0.3 and OmegaL0=0.7, this is what the
;   "mass" keyword returns).
; EXAMPLE:
;   For using with photoop system:
; 
;    ra=136.
;    dec=20.
;    obj= sdss_findobj(ra, dec, rerun=137, childobj=calibobj)
;    findspec, ra, dec, slist=slist
;    readspec, slist.plate, slist.fiberid, mjd=slist.mjd, zans=zans
;    kc= sdss_kcorrect(zans.z, calibobj=calibobj) 
;  
;   For reading tsobj structures:
;
;    tsobj=mrdfits('tsObj-01336-3-0456.fit',1,row=100)
;    findspec, tsobj.ra, tsobj.dec, slist=slist
;    readspec, slist.plate, slist.fiberid, mjd=slist.mjd, zans=zans
;    kc= sdss_kcorrect(zans.z, tsobj=tsobj) 
; 
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function maxbcg::sdss_kcorrect, $
						redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                        calibobj=calibobj, tsobj=tsobj, flux=flux, $
                        band_shift=in_band_shift, chi2=chi2, coeffs=coeffs, $
                        rmaggies=rmaggies, omaggies=omaggies, $
                        oivar=oivar, vname=in_vname, mass=mass, mtol=mtol, $
                        absmag=absmag, amivar=amivar, omega0=omega0, $
                        omegal0=omegal0, lrg=lrg, mets=mets, b300=b300, $
                        b1000=b1000, intsfh=intsfh, $
						clear=clear

common com_maxbcg_sdss_kcorrect, rmatrix, zvals, band_shift, vname, ermatrix

if keyword_set(clear) then begin
	delvarx, rmatrix, zvals, band_shift, vname, ermatrix
endif

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_elements(coeffs) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library, 'sdss_kcorrect'
    return, -1
endif 

if(keyword_set(lrg)) then $
  flux='model'

if(n_elements(in_vname) gt 0) then begin
    use_vname=in_vname
endif else begin
    if(keyword_set(lrg)) then $
      use_vname='lrg1' $
    else $
      use_vname='default'
endelse
if(n_elements(vname) gt 0) then begin
    if(vname ne use_vname) then begin
        rmatrix=0
        ermatrix=0
        zvals=0
    endif
endif
vname=use_vname

;; need to reset rmatrix if band_shift changes
if(n_elements(in_band_shift) gt 0) then $
  use_band_shift=in_band_shift $
else $
  use_band_shift=0. 
if(n_elements(band_shift) gt 0) then begin
    if(band_shift ne use_band_shift) then begin
        rmatrix=0
        ermatrix=0
        zvals=0
    endif
endif
band_shift=use_band_shift

if(keyword_set(mag) AND keyword_set(err)) then begin
    mgy=(10.D)^(-(0.4D)*(mag))
    mags_ivar=1./err^2
    mgy_ivar= mags_ivar/(0.4*alog(10.)*mgy)^2.
endif
if(keyword_set(nmgy) AND keyword_set(ivar)) then begin
    mgy=1.e-9*nmgy
    mgy_ivar=1.e+18*ivar
endif
if(n_tags(tsobj) gt 0 OR n_tags(calibobj) gt 0) then $
  sdss_to_maggies, mgy, mgy_ivar, calibobj=calibobj, tsobj=tsobj, flux=flux

if(keyword_set(lrg)) then $
  if(keyword_set(mgy_ivar)) then $
  mgy_ivar[0,*]=0.

;; call kcorrect
kcorrect, mgy, mgy_ivar, redshift, kcorrect, band_shift=band_shift, $
  rmatrix=rmatrix, zvals=zvals, coeffs=coeffs, rmaggies=rmaggies, $
  vname=vname, mass=mass, mtol=mtol, absmag=absmag, amivar=amivar, $
  omega0=omega0, omegal0=omegal0, chi2=chi2, mets=mets, b300=b300, $
  intsfh=intsfh, b1000=b1000

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, kcorrect

end








function maxbcg::cleanup
  ;; Nothing to clean up
  return,1
end 

pro maxbcg__define

  struct = { maxbcg, $
             catalog:'', $
             INHERITS objshear, $
             INHERITS postgres}

end 
