

FUNCTION vagc::init
  return,1
END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; File names, directories, and reading
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION vagc::dir, type, release=release, version=version, letter=letter, post=post

  on_error, 2
  IF n_elements(type) EQ 0 THEN BEGIN 
      print,'-Syntax:  dir = v->dir(type, version=, letter=, post=)'
      print,"  type=('vagc' | 'lss')"
      print," Default letter='safe', post='0'"
      print
      message,'Halting'
  ENDIF 

  CASE strlowcase(type) OF
      'vagc': BEGIN 
          dir = sdssidl_config('vagc_dir')
          IF n_elements(release) EQ 0 THEN BEGIN 
              release = sdssidl_config('vagc_release')
          ENDIF 
          IF n_elements(version) EQ 0 THEN BEGIN 
              version = sdssidl_config('vagc_vers')
          ENDIF 

          dir = concat_dir(dir, release)
          dir = concat_dir(dir, version)
      END 
      'lss': BEGIN 
          dir = sdssidl_config('lss_dir')
          IF n_elements(version) EQ 0 THEN BEGIN 
              version = sdssidl_config('lss_vers')
          ENDIF 
          dir = concat_dir(dir, version)

          IF n_elements(letter) NE 0 THEN BEGIN 
              dir = concat_dir(dir, strlowcase(letter))
          ENDIF 
          IF n_elements(post) NE 0 THEN BEGIN 
              dir = concat_dir(dir, strlowcase(post))
          ENDIF 
      END 
      ELSE: message,'Unknown type: '+ntostr(type)
  ENDCASE 

  return, dir
END 

FUNCTION vagc::vagc_file, subtype, run=run, camcol=camcol, release=release, version=version, collisiontype=collisiontype, kcorrmag=kcorrmag

	;; This will return version if not sent
	dir = self->dir('vagc', release=release, version=version)

	CASE strlowcase(subtype) OF
		; These are under the /sdss subdirectory
		'sdss_spectro_catalog': begin
			file = path_join('sdss', 'sdss_spectro_catalog.fits')
		end
		'sdss_spectro_objects': begin
			file = path_join('sdss', 'sdss_spectro_objects.fits')
		end

		'sdss_tiling_catalog': begin
			file = path_join('sdss', 'sdss_tiling_catalog.fits')
		end
		'sdss_tiling_geometry': begin
			file = path_join('sdss', 'sdss_tiling_geometry.fits')
		end
		'sdss_tiling_objects': begin
			file = path_join('sdss', 'sdss_tiling_objects.fits')
		end
		
		'sdss_target_geometry': begin
			file = path_join('sdss', 'sdss_target_geometry.fits')
		end


		'sdss_imaging_catalog': begin
			file = path_join('sdss', 'sdss_imaging_catalog.fits')
		end
		'sdss_imaging_geometry': begin
			file = path_join('sdss', 'sdss_imaging_geometry.fits')
		end
		'sdss_imaging_objects': begin
			file = path_join('sdss', 'sdss_imaging_objects.fits')
		end
		

		; extra imaging info under /sdss/parameters
		'calibobj': BEGIN 
			file = $
				'calibObj-'+run2string(run)+'-'+ntostr(camcol)+'.fits'
			file=path_join(['sdss','parameters'], file)
		END 


		; main directory.  This is only a subset, add more as needed
		; These all get a object_ prefix added
		'catalog': BEGIN 
			file = 'object_catalog.fits'
		END 
		'sdss_spectro': BEGIN 
			file = 'object_sdss_spectro.fits'
		END 
		'sdss_imaging': BEGIN 
			file = 'object_sdss_imaging.fits'
		END 
		'sdss_tiling': BEGIN 
			file = 'object_sdss_tiling.fits'
		END 
		'twodf': BEGIN 
			file = 'object_twodf.fits'
		END 
		'twomass': BEGIN 
			file = 'object_twomass.fits'
		END 
		'rc3': BEGIN 
			file = 'object_rc3.fits'
		END 
		'pscz': BEGIN 
			file = 'object_pscz.fits'
		END 

		'radec': BEGIN 
			file = 'object_radec.fits'
		END 


		;; kcorrect subdir
		'kcorrect': BEGIN 
			IF n_elements(collisiontype) EQ 0 THEN BEGIN 
				message,'You must send collisiontype= for "kcorrect" files'
			ENDIF 
			IF n_elements(kcorrmag) EQ 0 THEN BEGIN 
				message,'You must send kcorrmag= for "kcorrect" files'
			ENDIF 
			file = 'kcorrect.'+collisiontype+'.'+kcorrmag+'.z0.10.fits'
			file = path_join('kcorrect', file)
		END 
		ELSE: message,'Unknown subtype: '+ntostr(subtype)
	ENDCASE 

	file = path_join(dir, file)
	return, file
      
END 

;; Currently only returns a few types of catalogs
FUNCTION vagc::lss_file, type, version=version, letter=letter, post=post, randnum=randnum

  CASE strlowcase(type) OF 
      ;; Main directory
      'lss_index': BEGIN 
          dir = self->dir('lss', version=version)
          file = 'lss_index.'+strlowcase(version)+'.fits'
          file = concat_dir(dir, file)
      END 
      'lss_geometry': BEGIN 
          dir = self->dir('lss', version=version)
          file = 'lss_geometry.'+strlowcase(version)+'.fits'
          file = concat_dir(dir, file)
      END 
      'lss_bsmask': BEGIN 
          dir = self->dir('lss', version=version)
          file = 'lss_bsmask.'+strlowcase(version)+'.fits'
          file = concat_dir(dir, file)
      END 
      'lss_combmask': BEGIN 
          dir = self->dir('lss', version=version)
          file = 'lss_combmask.'+strlowcase(version)+'.fits'
          file = concat_dir(dir, file)
      END 
      'lss_imta': BEGIN 
          dir = self->dir('lss', version=version)
          file = 'lss_imta.'+strlowcase(version)+'.fits'
          file = concat_dir(dir, file)
      END 
      'lss_random': BEGIN 
          IF n_elements(randnum) EQ 0 THEN  BEGIN 
              message,'You must send randnum= for "lss_random" files'
          ENDIF 
          dir = self->dir('lss', version=version)
          dir = concat_dir(dir, 'random')
          file = 'lss_random-'+ntostr(randnum)+'.'+version+'.fits'
          file = concat_dir(dir, file)
      END 


      ;; Letter directory
      'window': BEGIN 
          IF n_elements(letter) NE 0 THEN BEGIN 
              dir = self->dir('lss', version=version, letter=letter)
          ENDIF ELSE BEGIN 
              message,'You must send letter= for "window" files'
          ENDELSE 
          file = 'window.'+strlowcase(version)+strlowcase(letter)+'.fits'
          file = concat_dir(dir, file)
      END 
      'mask':BEGIN 
          IF n_elements(letter) NE 0 THEN BEGIN 
              dir = self->dir('lss', version=version, letter=letter)
          ENDIF ELSE BEGIN 
              message,'You must send letter= for "mask" files'
          ENDELSE 
          file = 'mask.'+strlowcase(version)+strlowcase(letter)+'.fits'
          file = concat_dir(dir, file)
      END 
      'combmask':BEGIN 
          IF n_elements(letter) NE 0 THEN BEGIN 
              dir = self->dir('lss', version=version, letter=letter)
          ENDIF ELSE BEGIN 
              message,'You must send letter= for "combmask" files'
          ENDELSE 
          file = 'combmask.'+strlowcase(version)+strlowcase(letter)+'.fits'
          file = concat_dir(dir, file)
      END 
      'lss_to_window':BEGIN 
          IF n_elements(letter) NE 0 THEN BEGIN 
              dir = self->dir('lss', version=version, letter=letter)
          ENDIF ELSE BEGIN 
              message,'You must send letter= for "lss_to_window" files'
          ENDELSE 
          file = 'lss_to_window.'+strlowcase(version)+strlowcase(letter)+'.fits'
          file = concat_dir(dir, file)
      END 



      'letter_catalog': BEGIN 
          IF n_elements(letter) NE 0 THEN BEGIN 
              dir = self->dir('lss', version=version, letter=letter)
              file = $
                'letter_catalog.'+$
                strlowcase(version)+strlowcase(letter)+'.fits'
              file = concat_dir(dir, file)
          ENDIF ELSE BEGIN 
              message,'You must send letter= for "letter_catalog" files'
          ENDELSE 
      END 

      ;; randoms for subsamples
      'letter_random': BEGIN 
          IF n_elements(randnum) EQ 0 THEN  BEGIN 
              message,'You must send randnum= for "random" files'
          ENDIF 
          IF n_elements(letter) NE 0 THEN BEGIN 
              dir = self->dir('lss', version=version, letter=letter)
              dir = concat_dir(dir, 'random')
              file = $
                'random-'+ntostr(randnum)+'.'+$
                strlowcase(version)+strlowcase(letter)+'.fits.gz'
              file = concat_dir(dir, file)
          ENDIF ELSE BEGIN 
              message,'You must send letter= for "random" files'
          ENDELSE           
      END 


      ;; post directory
      'post_catalog': BEGIN 
          IF n_elements(letter) NE 0 AND n_elements(post) NE 0 THEN BEGIN 
              dir = self->dir('lss',version=version,letter=letter,post=post)
              file = $
                'post_catalog.'+$
                strlowcase(version)+$
                strlowcase(letter)+$
                strlowcase(post)+$
                '.fits'
              file = concat_dir(dir, file)
          ENDIF ELSE BEGIN 
              message,'You must send letter= and post= for "post_catalog" files'
          ENDELSE 
      END 

      ELSE: message,'Unknown type: '+ntostr(type)
  ENDCASE 

  return, file


END 

FUNCTION vagc::file, type, subtype, run=run, camcol=camcol, version=version, letter=letter, post=post, randnum=randnum, collisiontype=collisiontype, kcorrmag=kcorrmag

  on_error, 2
  IF n_elements(type) EQ 0 OR n_elements(subtype) EQ 0 THEN BEGIN 
      print,'-syntax: f = v->file(type, subtype, version=, letter=, post=, randnum=, run=, camcol=, collisiontype=, kcorrmag=)'
      print
      message,'Halting'
  ENDIF 
  CASE strlowcase(type) OF 
      'vagc': BEGIN 
          file = self->vagc_file(subtype, version=version, $
                                 run=run, camcol=camcol, $
                                 collisiontype=collisiontype, kcorrmag=kcorrmag)
      END 
      'lss': BEGIN 
          file = self->lss_file(subtype, version=version, $
                                letter=letter, post=post, randnum=randnum)
      END 
  ENDCASE 
  return, file
END 
FUNCTION vagc::read, type, subtype, run=run, camcol=camcol, version=version, letter=letter, post=post, randnum=randnum, collisiontype=collisiontype, kcorrmag=kcorrmag, rows=rows, columns=columns

    file = self->file(type, subtype, $
                    run=run, camcol=camcol, $
                    version=version, $
                    letter=letter, post=post, randnum=randnum, $
                    collisiontype=collisiontype, kcorrmag=kcorrmag)

    nf = n_elements(file)
    if nf gt 1 then begin
        comm = 'struct=mrdfits_multi(file, columns=columns, status=status)'
    endif else begin
        comm = 'struct=mrdfits(file, 1, rows=rows, columns=columns, status=status)'
    endelse

    print,'Reading file: ',file
    ;struct = mrdfits(file, 1, rows=rows, columns=columns, status=status)
    if not execute(comm) then message,'Failed to run mrdfits'
    IF status NE 0 THEN BEGIN 
        IF strmatch(file, '*.gz') THEN BEGIN 
            message,'Failed to read .gz file, trying without .gz',/inf
            file = (strsplit(file, '\.gz$',/extract,/regex))[0]
            if not execute(comm) then message,'Failed to run mrdfits'
            ;struct = mrdfits(ff, 1, rows=rows, columns=columns, status=status)          
        ENDIF ELSE BEGIN 
            message,'Failed to read file, trying with .gz',/inf
            file = file + '.gz'
            ;struct = mrdfits(ff, 1, rows=rows, columns=columns, status=status)          
            if not execute(comm) then message,'Failed to run mrdfits'
        ENDELSE 
        IF status NE 0 THEN message,'Failed to read file'
    ENDIF 
  return,struct
END 















function vagc::jackie_defaults
    return, {vagc_release:'vagc-dr6',vagc_version:'vagc0', $
             version:'dr6', letter:'full', post:'0'}
end
function vagc::jackie_dir
    dir = '/global/early2/esheldon/jackiecat'
    return, dir
end 
function vagc::jackie_file, version=version, letter=letter, post=post, randnum=randnum

    outdir = self->jackie_dir()

    def = self->jackie_defaults()
    if n_elements(version) eq 0 then version=def.version
    if n_elements(letter) eq 0 then letter=def.letter
    if n_elements(post) eq 0 then post=def.post

    outfile = 'jackiecat-'+version+letter+post
    if n_elements(randnum) ne 0 then begin 
        outfile = outfile + '-random'+strn(randnum, len=2, padchar='0')
    endif 
    outfile = outfile + '.st'
    outfile = concat_dir(outdir, outfile)
    return, outfile
end 

function vagc::jackie_read, version=version, letter=letter, post=post, randnum=randnum
    file = self->jackie_file(version=version, letter=letter, post=post, $
                            randnum=randnum)
    struct = read_idlstruct(file)
    return, struct
end 

pro vagc::make_jackiecat_geometry, version=version

    def = self->jackie_defaults()
    if n_elements(version) eq 0 then version=def.version
    if n_elements(letter) eq 0 then letter=def.letter

    dir = self->jackie_dir()
;  file = self->lss_file('lss_geometry', version=version)

    letter = 'full'
    file = self->lss_file('window', version=version, letter=letter)

    dirsep, file, tdir, tfile
    file = repstr(tfile, '.fits', '.st')
    file = concat_dir(dir, file)
    print,file

    g = self->read('lss','window', version=version, letter=letter)
    ng = n_elements(g)
    outst = {area_str: 0d, sector: 0L, fgot: 0.0}
    outst = replicate(outst, ng)
    outst.area_str = g.str
    outst.sector = g.sector
    outst.fgot = g.fgot

    print,'Writing file: ',file
    write_idlstruct, outst, file, /ascii

end 
pro vagc::make_jackiecat, version=version, random=random

    def = self->jackie_defaults()
    if n_elements(version) eq 0 then version=def.version
 
    vagc_release = def.vagc_release
    vagc_version = def.vagc_version
    letter = def.letter
    post = def.post
    collisiontype='nearest'
    fgotlim = 0.0

  

    if not keyword_set(random) then begin 

        outfile = self->jackie_file(version=version, letter=letter, post=post)

        slss = self->read('lss', 'post_catalog', $
            version=version, letter=letter, post=post)
        rows = slss.object_position
        nlss = n_elements(slss)
        print,"Only Using "+ntostr(nlss)+" rows from the spec file"

        print
        slss_index = self->read('lss', 'lss_index', version=version)
        window = self->read('lss','window', version=version, letter=letter)

        spcolumns = ['mjd', 'plate', 'fiberid', 'objid', 'primtarget', $
            'sectarget', 'class', 'subclass', 'z', 'z_err', 'rchi2', $
            'dof', 'zwarning', 'vdisp', 'vdisp_err', $
            'vdispchi2', 'vdispdof']
        print

;      sp = self->read('vagc', 'spectro', rows=rows, columns=spcolumns)
;      
;      kpetro = self->read('vagc','kcorrect', $
;                          collisiontype=collisiontype, kcorrmag='petro', $
;                          rows=rows)
      
      
        print
        print,'Making completeness cut'
        fgot = window[ slss.iwindow ].fgot
        ws = where(fgot GT fgotlim, nws)
        print,'Kept '+ntostr(nws)+'/'+ntostr(nlss)

        outStruct = create_struct($
            'ra',             0d, $
            'dec',            0d, $
            'abs_petro_mag',  fltarr(5), $
            'z',             0.0, $
            'fgot',          0.0, $
            'plate',         -1L, $
            'sector',        -1L, $
            'ilss',          -1L, $
            'mmax',          0.0)

        outStruct = replicate(outStruct, nws)

        outStruct.ra = slss[ws].ra
        outStruct.dec = slss[ws].dec

        outStruct.abs_petro_mag[0] = slss[ws].absm[0]
        outStruct.abs_petro_mag[1] = slss[ws].absm[1]
        outStruct.abs_petro_mag[2] = slss[ws].absm[2]
        outStruct.abs_petro_mag[3] = slss[ws].absm[3]
        outStruct.abs_petro_mag[4] = slss[ws].absm[4]

        outStruct.z = slss[ws].z
        outStruct.fgot = fgot[ws]
        outStruct.mmax = slss[ws].mmax

        wsind = slss[ws].object_position
        outStruct.sector = slss_index[wsind].sector
        outStruct.ilss = slss_index[wsind].ilss

        hdr = {$
            vagc_release:vagc_release, $
            vagc_version:vagc_version, $
            version: version, $
            letter: letter, $
            post:post, $
            collisiontype: collisiontype, $
            fgotlim: fgotlim}
        print
        print,'Writing outFile = ',outFile
        write_idlstruct, outStruct, outFile, /ascii, hdr=hdr

        delvarx, outStruct, sp, kPetro, slss, window, slss_index

    endif else begin 

        hdr = { $
            vagc_release:vagc_release, $
            vagc_version:vagc_version, $
            version: version, $
            letter: letter, $
            fgotlim: fgotlim}

        rcolumns = ['ra','dec','fgot','sector','mmax']

        irand = 0L
        infile = self->file('lss', 'letter_random', $
            version=version, letter=letter, randnum=irand)
        while fexist(infile) do begin 

            print,'-----------------------------------------------------'
            outfile = $
                self->jackie_file(version=version, letter=letter, post=post, $
                                    randnum=irand)

            print
            print,'Reading random file: ',infile
            r=mrdfits(infile, 1, columns=rcolumns)
            
            print
            print,'Making completeness cut'
            w = where(r.fgot GT fgotlim, nkeep)
            print,'Keeping '+ntostr(nkeep)

            r = r[w]

            print
            print,'Writing file: ',outFile
            write_idlstruct, r, outFile, /ascii, hdr=hdr

            irand = irand + 1
            infile = self->file('lss', 'letter_random', $
                version=version,letter=letter,randnum=irand)
        endwhile 

    endelse 
    return
END 










PRO vagc::make_photoz_train, sp, im

  IF n_elements(sp) EQ 0 THEN BEGIN 
      sp = self->read('vagc','spectro', columns=['z','primtarget'])
  ENDIF 
  IF n_elements(im) EQ 0 THEN BEGIN 
      im = self->read('vagc','imaging', columns=['run','rerun','camcol','field','id','modelflux','modelflux_ivar'])
  ENDIF 

  w = where(sp.z GT 0.01 AND sp.z LT 0.6)
  ww = !sdss->flag_select(sp.primtarget, $
                          'primtarget', $
                          {galaxy:'y',galaxy_red:'y'}, $
                          nkeep, $
                          input_index=w, $
                          /orflags)
  
  arrval = replicate(-9999.0, 5)
  outstruct = { $
                id: 0LL, $
                modelmag: arrval, modelmag_err: arrval, $
                z:0.0 $
              }
  outstruct = replicate(outstruct,nkeep)

  id = !sdss->photoid(im)
  outstruct.id = id[ww]

  FOR band=0L, 4 DO BEGIN 

      wm = where( im[ww].modelflux_ivar[band] GT 0 AND $
                  im[ww].modelflux[band] GT 0, nwm)
      
      flux = im[ww[wm]].modelflux[band]
      fluxerr = 1.0/sqrt( im[ww[wm]].modelflux_ivar[band] )
      outstruct[wm].modelmag[band] = 22.5 - 2.5*alog10( flux )
      outstruct[wm].modelmag_err[band] = 2.5*fluxerr/(flux*alog(10.0))

  ENDFOR 

  outstruct.z = sp[ww].z

  file = '/global/evolve1/esheldon/photoz_trainingset/photoz_trainingset.st'
  print
  print,'Writing to file: ',file
  write_idlstruct, outstruct, file, /ascii
END 









PRO vagc::mag2color, mags, color, fratio, abs=abs

  nobj = n_elements(mags)/5
  defval = -9999.0
  color = replicate(defval, 4, nobj)
  fratio = color
  
  IF keyword_set(abs) THEN BEGIN 
      minmag = -25.0
      maxmag = -5.0
  ENDIF ELSE BEGIN 
      minmag = 10.0
      maxmag = 27.5
  ENDELSE 

  FOR i=0,3 DO BEGIN 

      w = where(mags[i,*] GT minmag AND mags[i,*] LT maxmag AND $
                mags[i+1,*] GT minmag AND mags[i+1,*] LT maxmag, ngood)

      IF ngood NE 0 THEN BEGIN 
          color[i,w] = reform( mags[i,w] - mags[i+1,w] )
          fratio[i,w] = 10.0^(-color[i,w]/2.5)
      ENDIF 
      

  ENDFOR 


END 

PRO vagc::flux2fratio, flux, flux_ivar, fratio, fratio_ivar

  nobj = n_elements(flux)/5

  defval = -9999.0

  fratio = replicate(-9999.0, 4, nobj)
  fratio_ivar = replicate(0.0, 4, nobj)


  minflux = 0.01                ; 27.5 mag
  maxflux = 1.e5                ; 10th mag


  FOR i=0,3 DO BEGIN 

      w = where(flux_ivar[i,  *] GT 0 AND $
                flux_ivar[i+1,*] GT 0 AND $
                flux[i,  *] GT minflux AND $
                flux[i+1,*] GT minflux AND $
                flux[i,  *] LT maxflux AND $
                flux[i+1,*] LT maxflux, ngood)
      
      IF ngood NE 0 THEN BEGIN 
          fratio[i, w] = reform( flux[i,w]/flux[i+1,w] )
          sigma2 = fratio[i, w]^2*( 1.0/flux_ivar[i,w]/flux[i,w]^2 + $
                                    1.0/flux_ivar[i+1,w]/flux[i+1,w]^2 )

          fratio_ivar[i, w] = reform( 1.0/sigma2 )
      ENDIF 

  ENDFOR 

END 


; The spec tags to keep for main table
FUNCTION vagc::vagc_stuff_spec_tags
  tags = $
    ['plate','fiberid','mjd',$
	 'progname',$
     'objtype','class','subclass',$
     'z','z_err','zwarning',$
     'vdisp','vdisp_err', $
     'primtarget','sectarget','plug_ra','plug_dec']
  return,tags
END 

FUNCTION vagc::vagc_stuff_structdef, n

  iarrval = replicate(-9999, 5)
  farrval = replicate(-9999.0, 5)
  larrval = replicate(-9999L, 5)
  st = $
    {vagc_index: 0L, $
     photoid: 0LL, $
     run: 0, $
     rerun:0, $
     camcol:0, $
     field:0, $
     id:0, $
     stripe: 0, $
     ra: 0d, $
     dec: 0d, $
     $
     plate: 0, $
     fiberid: 0, $
     mjd: 0L, $
     primtarget: 0L, $
     sectarget: 0L, $
     $
     objtype: '', $
     class: '', $
     subclass: '', $
	 progname:'',$
     $
     z: 0.0, $
     z_err: 0.0, $
     zwarning: 0.0, $
     $
     vdisp: 0.0, $
     vdisp_err: 0.0, $
     $
     mjd_imaging: 0L, $
     nchild: 0, $
     parent:0, $
     objc_type: 0, $
     resolve_status: 0, $
     flags: larrval, $
     flags2: larrval, $
     objc_flags: 0L, $
     objc_flags2: 0L, $
     rowc: farrval, $
     colc: farrval, $
     $
     modelflux: farrval, $
     modelflux_ivar: farrval, $
     modelfratio: replicate(-9999.0, 4), $
     modelfratio_ivar: replicate(-9999.0, 4), $
     petroflux: farrval, $
     petroflux_ivar: farrval, $
     cmodelflux: farrval, $
     cmodelflux_ivar: farrval, $
     psfflux: farrval, $
     psfflux_ivar: farrval, $
	 fiberflux: farrval, $
	 fiberflux_ivar: farrval, $
	 $
	 r_dev: farrval, $
	 r_exp: farrval, $
     $
     kcorrect_nearest: farrval, $
     modelabsmag: farrval, $
     kmodelfratio: replicate(-9999.0, 4), $
     extinction: farrval, $
     $
     petror50: farrval, $
     petror90: farrval, $
     $
     m_e1_corr: farrval, $
     m_e2_corr: farrval, $
     m_e1e1err: farrval, $
     m_e1e2err: farrval, $
     m_e2e2err: farrval, $
     m_rr_cc: farrval, $
     m_r: farrval, $
     m_seeing: farrval, $
     $
     vagc_select: 0L, $
     score: 0.0, $
     calib_status: iarrval, $
     ifield: 0L, $
     balkan_id: 0L}

  IF n_elements(n) NE 0 THEN BEGIN 
      st = replicate(st, n)
  ENDIF 
  return, st

END 


FUNCTION vagc::getsmear, pstruct
  
  nst = n_elements(pstruct)
  corr = replicate(-9999.0, 5, nst)

  FOR clr=0,4 DO BEGIN 
      w = where(pstruct.m_rr_cc_psf GT 0 AND $
                pstruct.m_rr_cc GT 0 AND $
                pstruct.m_cr4 NE 0.25, nw)
      
      IF nw NE 0 THEN BEGIN 
          corr[clr,w] = pstruct[w].m_rr_cc_psf[clr]/pstruct[w].m_rr_cc[clr]*(4./pstruct[w].m_cr4_psf[clr] -1.)/(4./pstruct[w].m_cr4[clr]-1.)
      ENDIF 
  ENDFOR 
  return, corr

END 

FUNCTION vagc::getseeing, struct

  nst = n_elements(struct)
  seeing = replicate(9999.0, 5, nst)

  FOR clr=0,4 DO BEGIN 
      spsf = struct.m_rr_cc_psf[clr] 
      w = where(spsf GT 0.0, nw)
      IF nw NE 0 THEN seeing[clr,w] = 2.35*0.4*sqrt( spsf[w]/2.0 )
  ENDFOR 
  return, seeing

END 

PRO vagc::add_shape_correct, instruct, outstruct
  
  corr = self->getsmear(instruct)
  outstruct.m_r = corr

  FOR clr=0,4 DO BEGIN 
      w = where(corr[clr,*] GT 0 AND $
                instruct.m_e1[clr] NE -1000 AND $
                instruct.m_e1[clr] NE -9999, nw)
      IF nw NE 0 THEN BEGIN 
          outstruct[w].m_e1_corr[clr] = instruct[w].m_e1[clr] - corr[clr,w]*instruct[w].m_e1_psf[clr]
          outstruct[w].m_e2_corr[clr] = instruct[w].m_e2[clr] - corr[clr,w]*instruct[w].m_e2_psf[clr]
      ENDIF 
  ENDFOR 

END 

PRO vagc::add_cmodelflux, instruct, outstruct

  FOR clr=0,4 DO BEGIN 
      
      w = where(instruct.expflux_ivar[clr] GT 0.0 AND $
                instruct.devflux_ivar[clr] GT 0.0, nw)


      fracdev = instruct[w].fracpsf[clr] > 0.0 < 1.0

      outstruct[w].cmodelflux[clr] = $
        fracdev*instruct[w].devflux[clr] + (1.0-fracdev)*instruct[w].expflux[clr]

      cvar = $
        fracdev^2/instruct[w].devflux_ivar[clr] + $
        (1.0-fracdev)^2/instruct[w].expflux_ivar[clr]

      outstruct[w].cmodelflux_ivar[clr] = 1.0/cvar

  ENDFOR 


END 

pro vagc::vagc_stuff, im=im, sp=sp, kcorr=kcorr
  
    table = 'vagc'
    metatable = 'vagc_meta'

    if self->postgres::table_exists(table) then begin
        message,'Please drop old table '+table+' first'
    endif

    ; metadata table
    meta = {release: sdssidl_config('vagc_release'), $
            version: sdssidl_config('vagc_vers'), $
            collisiontype: 'nearest'}


    if self->postgres::table_exists(metatable) then begin
        message,'Dropping table '+metatable+' first',/inf
        self->postgres::query, 'drop table '+metatable, conn='user=postgres' 
    endif
    self->postgres::struct2table, meta, metatable
    self->postgres::query, 'grant select on '+metatable+' to sdss', $
        conn='user=postgres'



    sptags = self->vagc_stuff_spec_tags()

    if n_elements(sp) eq 0 then begin 
        sp = self->read('vagc','spectro',columns=sptags)
		; only actual spectra
		wkeep=where(sp.plate ne -1)
		sp=sp[wkeep]

        im = self->read('vagc', 'imaging', rows=wkeep, $
                        columns=['run','rerun','camcol','field','id', $
                                 'calibobj_position'])
        kcorr = self->read('vagc', 'kcorrect', rows=wkeep, $
                            collisiontype=meta.collisiontype, $
                            kcorrmag = 'model', $
                            columns = ['kcorrect','absmag'])
    endif 

    vagc_index = lindgen(n_elements(sp))

    ;; Make tree structure of run,rerun,camcol
    print,'Running histid'
    idstruct = self->sdss_util::histid(im.run, fix(im.rerun), im.camcol)
      

    sdssidl_setup
    pruns = idStruct.runs
    for ri=0l, idstruct.nruns-1 do begin 
        print,'-----------------------------------------------------------'
        print
        run = (*pruns)[ri].run
        stripe = !sdss->stripe(run)

        if run ge 0 then begin 

            preruns = (*pruns)[ri].reruns
            for rri=0l, (*pruns)[ri].nreruns-1 do begin 
                rerun = (*preruns)[rri].rerun 
                pcamcols = (*preruns)[rri].camcols
                for ci=0l, (*preruns)[rri].ncamcols-1 do begin 
                  
                    camcol = (*pcamcols)[ci].camcol 
                    w= *(*pcamcols)[ci].indices

                    calib=self->read('vagc', 'calibObj', run=run, camcol=camcol)

                    ;; position in the calibobj file
                    pos = im[w].calibobj_position

                    nw = n_elements(w)
                    nc = n_elements(calib)

                    if nw gt nc then begin 
                        message,'Mismatch in number'
                    endif 
                  
                    help,w,calib

                    outst = self->vagc_stuff_structdef(nw)

                    ;; index
                    outst.vagc_index = vagc_index[w]

                    ;; copy imaging

                    print,'Copying imaging'
                    copy_struct, calib[pos], outst
                    outst.mjd_imaging = calib[pos].mjd
                    outst.rerun = fix(calib[pos].rerun)

                    ;; copy spectro
                    print,'Copying spectro'
                    copy_struct, sp[w], outst

                    ;; copy kcorreccted stuff
                    print,'Copying kcorr stuff'
                    outst.kcorrect_nearest = kcorr[w].kcorrect[0:4]
                    outst.modelabsmag = kcorr[w].absmag[0:4]

                    print,'Getting flux ratios (model)'
                    self->flux2fratio, $
                        calib[pos].modelflux, calib[pos].modelflux_ivar, $
                        fratio, fratio_ivar

                    print,'Getting kfratio (model)'
                    self->mag2color, kcorr[w].absmag[0:4], kcolor, kfratio, /abs

                    outst.modelfratio = fratio
                    outst.modelfratio_ivar = fratio_ivar
                    outst.kmodelfratio = kfratio

                    print,'Getting corrected shapes'
                    self->add_shape_correct, calib[pos], outst
                    outst.m_seeing = self->getseeing(calib[pos])
                    print,'Getting cmodelflux'
                    self->add_cmodelflux, calib[pos], outst


                    outst.stripe = stripe
                    outst.photoid = !sdss->photoid(outst)


                    self->postgres::struct2table, outst, 'vagc', $
                        primary_key='vagc_index', $
                        /varchar, $
                        conn='user=postgres', $
                        tmpdir = '/mount/early2/esheldon/tmp', $
                        status=status

                    if status ne 0 then begin 
                        message,'Error stuffing'
                    endif 
                  
                    outst = 0
                endfor ;; camcols
            endfor ;; reruns
        endif ;; this if can help with crashes
    endfor   ;; runs

    self->sdss_util::histid_destroy, idstruct

    self->postgres::query, 'grant select on vagc to sdss', $
        conn='user=postgres'

end 




;; for now, this is just a single "letter" and "post"
;; I may be able to combine all letters and posts into
;; a single table with flags set, but I'm not sure yet
;; about the window.

FUNCTION vagc::lss_structdef, num
  st = $
    {vagc_index:0L, $
     ra: 0d, $
     dec: 0d, $
     z: 0.0, $
     absmag: fltarr(5), $
     fgot: 0.0, $
     sector: 0L, $
     ilss: 0L $
    }
  IF n_elements(num) NE 0 THEN BEGIN 
      st = replicate(st, num)
  ENDIF 
  return, st
END 
PRO vagc::lss_stuff, letter=letter, post=post



    vagc_release = sdssidl_config('vagc_release')
    vagc_version = sdssidl_config('vagc_vers')
    if n_elements(letter) eq 0 then letter = 'full'
    if n_elements(post) eq 0 then post='0'

    ;; e.g. v->lss_stuff, 'full','0'
    table = 'lss'+letter+post
    metatable = table + '_meta'

    if self->postgres::table_exists(table) then begin
        message,'Please drop old table '+table+' first'
    endif



    if self->postgres::table_exists(metatable) then begin
        message,'Dropping old table '+metatable+' first', /inf
        self->postgres::query, 'drop table '+metatable, conn='user=postgres' 
    endif

    meta = {release: sdssidl_config('vagc_release'), $
            version: sdssidl_config('vagc_vers'), $
            letter: letter, $
            post: post}
    self->postgres::struct2table, meta, metatable
    self->postgres::query, 'grant select on '+metatable+' to sdss', $
        conn='user=postgres'


    print
    print,'Will create table: ',table
    print,'---------------------------------------------------------'


    post = self->read('lss','post_catalog',letter=letter,post=post)
    window = self->read('lss', 'window', letter=letter)
    index = self->read('lss','lss_index')

    fgot = window[ post.iwindow ].fgot
    sector = index[ post.object_position ].sector

    num = n_elements(post)
    outst = self->lss_structdef(num)

    outst.vagc_index = post.object_position
    outst.ra = post.ra
    outst.dec = post.dec
    outst.z = post.z
    outst.absmag = post.absm[0:4]
    outst.fgot = fgot
    outst.sector = sector
    outst.ilss = post.ilss

    self->postgres::struct2table, $
        outst, table, primary_key='vagc_index', $
        tmpdir = '/global/early2/esheldon/tmp', status=status

    if status ne 0 then begin 
        message,'table was not stuffed'
    endif 

    self->postgres::query, 'grant select on '+table+' to sdss', $
        conn='user=postgres'
end 



PRO vagc::lowz_stuff

  table = 'lowz'

  dir = concat_dir(esheldon_config('lensinput_dir'), 'lowz')
  file = concat_dir(dir, 'lowz_comb.dr4.fits')
  print
  print,'Reading file: ',file
  lz = mrdfits(file,1)

  lz = rename_tags(lz, 'object_position', 'vagc_index')

  query = 'select vagc_index,ra,dec from vagc'
  print,query
  vagc = self->postgres::query(query)

  ;; match by ra/dec
  match_angle = 3               ;arcsec
  print
  print,'Matching by ra/dec within '+ntostr(match_angle)+' arcsec'

  match_angle = match_angle/3600.0*!dpi/180.0
  htm_match, lz.ra, lz.dec, vagc.ra, vagc.dec, match_angle, mlz, mv, d12, $
    maxmatch=1

  nlz = n_elements(lz)
  nm = n_elements(mlz)
  print,ntostr( (nlz-nm)/float(nlz)) + ' did not match'

  IF nlz NE nm THEN message,'Not all matched'

  lz[mlz].vagc_index = vagc[mv].vagc_index


  print
  print,'Sorting by vagc_index'
  s = sort(lz.vagc_index)
  lz = lz[s]

  self->postgres::struct2table, lz, table, primary_key='vagc_index', $
    conn='user=postgres'

END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; A match table between samples and vagc, to simplify things
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO vagc::make_matchtable

  table = 'vagc_match'

  LZBIT = 0
  FULL0BIT = 1
  ALL0BIT = 2

  ;; The lowz sample
  query = 'select vagc_index from lowz'
  print,query
  lz = self->postgres::query(query)
  
  ;; full0
  query = 'select vagc_index from lssfull0'
  print,query
  full0 = self->postgres::query(query)

  ;; all0
  query = 'select vagc_index from lssall0'
  print,query
  all0 = self->postgres::query(query)

  ;; now full vagc
  query = 'select vagc_index from vagc'
  print,query
  vagc = self->postgres::query(query)

  nv = n_elements(vagc)
  matchflags = intarr(nv)

  
  
  match, lz.vagc_index, vagc.vagc_index, mlz, mvagc
  nmatch = n_elements(mlz)
  nlz = n_elements(lz)
  print,'Matched '+ntostr(nmatch)+'/'+ntostr(nlz)+' from lowz'
  matchflags[mvagc] = matchflags[mvagc] + 2^LZBIT


  match, full0.vagc_index, vagc.vagc_index, mfull0, mvagc
  nmatch = n_elements(mfull0)
  nfull0 = n_elements(full0)
  print,'Matched '+ntostr(nmatch)+'/'+ntostr(nfull0)+' from full0'
  matchflags[mvagc] = matchflags[mvagc] + 2^FULL0BIT


  match, all0.vagc_index, vagc.vagc_index, mall0, mvagc
  nmatch = n_elements(mall0)
  nall0 = n_elements(all0)
  print,'Matched '+ntostr(nmatch)+'/'+ntostr(nall0)+' from all0'
  matchflags[mvagc] = matchflags[mvagc] + 2^ALL0BIT

  outst = {vagc_index:0L, matchflags:0}
  outst = replicate(outst, nv)

  outst.vagc_index = vagc.vagc_index
  outst.matchflags = matchflags

  self->postgres::struct2table, outst, table, primary_key='vagc_index',$
    conn='user=postgres'

END 




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





FUNCTION vagc::select_main, struct
  main = $
    sdss_flag_select(struct.primtarget, 'primtarget', {galaxy:'Y'})
  return, main
END 
FUNCTION vagc::select_lrg_old, struct
  w=where(struct.z GT 0.15)
  gred = $
    sdss_flag_select(struct[w].primtarget,'primtarget',{galaxy_red:'Y'})
  gred = w[gred]
  return, gred
END 

function vagc::select_lrg, struct, nw, zcut=zcut
	; this struct might be the result of reading
	;  ->file('vagc','spectro')
	; or ->read('vagc','spectro')
	logic = (struct.primtarget AND 2L^5 + 2L^26) NE 0

	; make sure we actually got a spectra!
	if tag_exist(struct,'plate') then begin
		logic = logic and (struct.plate ne -1)
	endif
	if tag_exist(struct,'progname') then begin
		logic = logic and strmatch(struct.progname,'main*')
	endif
	if tag_exist(struct,'class') then begin
		logic = logic and strmatch(struct.class,'GALAXY*')
	endif
	if tag_exist(struct,'zwarning') then begin
		logic = logic and (struct.zwarning eq 0)
	endif

	w = where(logic, nw)

	if keyword_set(zcut) then begin
		w2=where(struct.z gt 0.15,nw)
		if nw ne 0 then begin
			w=w[w2]
		endif else begin
			w=-1
		endelse
	endif
	return, w
end



FUNCTION vagc::select_qso, struct

  w=where(struct.z GT 0.03)

  st = {qso_hiz: 'Y', qso_cap:'Y', qso_skirt:'Y', $
        qso_first_cap:'Y', qso_first_skirt:'Y', qso_mag_outlier:'Y'}
  qso = !sdss->flag_select(struct.primtarget, 'primtarget', st, $
                            input=w,/orflags)
  return,qso
END 


PRO vagc::plot_zdist, struct, ylog=ylog

  IF n_elements(struct) NE 0 THEN BEGIN 
      struct = $
        self->query('select z,primtarget from vagc')
  ENDIF 

  main = self->select_main(struct)
  lrg = self->select_lrg(struct)
  qso = self->select_qso(struct)

  mclr = !darkGreen
  rclr = !red
  qclr = !steelBlue


  zmin = 0.01
  zmax = 0.6
  IF keyword_set(ylog) THEN BEGIN 
      yrange = [0.9, 5.e4]
      ystyle=3
  ENDIF 

  plothist, struct.z, bin=0.01, min=zmin, max=zmax, $
    xrange=[zmin,zmax], xstyle=3, $
    yrange=yrange, ystyle=ystyle, ylog=ylog, $
    xtitle='Z', ytitle='Number'
  plothist, struct[main].z, bin=0.01, /overplot, color=mclr
  plothist, struct[lrg].z, bin=0.01, /overplot, color=rclr
  plothist, struct[qso].z,  bin=0.01, /overplot, color=qclr


  legend, $
    ['All', 'Main', 'LRG', 'QSO'], $
    line=[0,0,0,0], $
    color=[!p.color, mclr, rclr, qclr], $
    /right, box=0

END 


function vagc::read_good_sdss_galaxies, columns=columns
	if n_elements(columns) ne 0 then begin
		cols=[columns, 'specprimary','zwarning', 'objtype']
		cols=cols[rem_dup(cols)]
	endif
	st=self->read('vagc', 'sdss_spectro_catalog', columns=cols)
	logic_specprimary = st.specprimary gt 0
	print,'specprimary: ',n_elements(where(logic_specprimary))
	logic_zwarning = st.zwarning eq 0
	print,'zwarning eq 0: ',n_elements(where(logic_zwarning))
	logic_zrange = st.z ge 1.e-3 and st.z le 0.6
	print,'zrange [1.e-3,0.6]: ',n_elements(where(logic_zrange))
	logic_galaxy = strmatch(st.objtype, 'GAL*') ne 0
	print,'galaxy: ',n_elements(where(logic_galaxy))

	wgood=where(logic_specprimary $
		and logic_zwarning $
		and logic_zrange $
		and logic_galaxy, ngood)
	print,'combined: ',ngood
	st=st[wgood]
	if n_elements(columns) ne 0 then begin
		st = extract_tags(st, columns)
	endif
	return, st
end

function vagc::match_spectro_radec, ra, dec, spectro=spectro
	; match against the sdss ones only
	if n_elements(spectro) eq 0 then begin
		spectro = self->read('vagc','sdss_spectro_catalog')
	endif

	; this gives almost unique objects
	wsp=where(spectro.specprimary eq 1 $
				and strmatch(spectro.objtype, 'QA*') eq 0)

	; match by ra/dec
	match_radius_deg = 5.0/3600.0

	tm1=systime(1)
	spherematch, ra, dec, spectro[wsp].plug_ra, spectro[wsp].plug_dec,$
		match_radius_deg, minput, mspectro, maxmatch=1
	tm2=systime(1)
	ptime, tm2-tm1

	if mspectro[0] eq -1 then return, -1

	mspectro = wsp[mspectro]
	outstruct = {plate:0L, fiberid:0L, mjd:0L}
	outstruct =replicate(outstruct, n_elements(mspectro))

	outstruct.plate = spectro[mspectro].plate
	outstruct.fiberid = spectro[mspectro].fiberid
	outstruct.mjd = spectro[mspectro].mjd
	
	return, outstruct

end











PRO vagc__define 

  struct = { vagc, $
             vagc_dummy: 0,$
             INHERITS sdss_util, $
             INHERITS postgres}

END 
