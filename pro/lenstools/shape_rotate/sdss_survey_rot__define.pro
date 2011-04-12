function sdss_survey_rot::init
  return,1
end 


PRO sdss_survey_rot::run_calcrot, system

	minscore = 0.1+0.001
	print, 'getting window_runlist, ..,  minscore=',minscore,$
		format='(a,f)'
	window_runlist, runs, rerun=reruns, minscore=minscore

	nruns = n_elements(runs)
	for i=0l, nruns-1 do begin 
		self->calcrot, system, runs[i], rerun=reruns[i]
	endfor 

END 


pro sdss_survey_rot::calcrot, system, run, rerun=rerun, $
                    indir=indir, outdir=outdir, status=status, $
					pgsql=pgsql

    if system eq 'survey' then begin
        self->calcrot_survey, run, rerun=rerun, $
              indir=indir, outdir=outdir, status=status, $
              pgsql=pgsql
    endif else if system eq 'eq' then begin
        self->calcrot_eq, run, rerun=rerun, $
              indir=indir, outdir=outdir, status=status, $
              pgsql=pgsql
    endif else message,'Bad system: '+system

end
pro sdss_survey_rot::calcrot_survey, run, rerun=rerun, $
                    indir=indir, outdir=outdir, status=status, $
					pgsql=pgsql

    ;; status is 1 unless we reach end
    status=1

    if n_params() lt 1 then begin 
        print,'-Syntax: obj->calcrot, run, rerun=, status=, indir=, outdir=,/pgsql,status='
        return
    endif 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Output file  names
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if n_elements(rerun) eq 0 then rerun=sdss_rerun(run)


    if n_elements(outdir) eq 0 then begin 
        outdir = self->dir('survey')
    endif 

    outfile = self->file('survey', run,rerun=rerun, dir=outdir)
    pgsql_file = self->file('survey', run,rerun=rerun,/postgres, dir=outdir)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; astrometry structure
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    radeg = 180./!dpi
    val = .4d/3600.
    astr={cd: double(identity(2)),$
        cdelt: [val,val], $
        crpix: dblarr(2), $
        crval: [0., 0.]}

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; define set of pixels to rotate
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    lrow = 1489-1
    lcol = 2048-1
    nnr=20.
    vert=dblarr(2,nnr)
    horz=dblarr(2,nnr)
    ccol = 1023.5
    crow = 744.0
    FOR i=0L, nnr-1 DO BEGIN
        vert[0,i] = ccol          ;columns
        vert[1,i] = lrow-i        ;rows
        IF i EQ 10 THEN vert[1,i] = 1000
        horz[0,i] = lcol - i      ;columns
        horz[1,i] = crow          ;rows
    ENDFOR 

    ; this is worse?
    ;vert[0,*] = ccol; columns
    ;vert[1,*] = arrscl(findgen(nnr), crow, lrow) ; rows

    ;horz[0,*] = arrscl(findgen(nnr), ccol, lcol) ; columns
    ;horz[1,*] = crow ; rows

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Output structure
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    clr = [0,1,2,3,4]
    colors = ['u','g','r','i','z']
    nclr = n_elements(clr)

    rotst = self->structdef(clr)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; loop over bandpasses
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;; get number of fields

    ;run_status = sdss_runstatus()
    ;IF n_elements(indir) EQ 0 THEN BEGIN 
    ;    w=where(run_status.run EQ run AND run_status.rerun EQ rerun, nw)
    ;    IF nw EQ 0 THEN BEGIN 
    ;        print,'Error: run/rerun not found: ',run,rerun
    ;        return
    ;    ENDIF 
    ;ENDIF 

    ttrans = sdss_read('astrans', run, 1, band=clr[0], rerun=rerun, $
        node=tnode, inc=tinc)

    nf = n_elements(ttrans.field)

    ptrlist = ptrarr(6)
    FOR camcol=1,6 DO BEGIN 

        rotstruct = replicate(rotst, nf)
        rotstruct.run = run
        rotstruct.rerun = rerun
        rotstruct.camcol = camcol
        rotstruct.field = ttrans.field

        rotstruct.fieldid = sdss_photoid(rotstruct.run,   $
            rotstruct.rerun, $
            rotstruct.camcol,$
            rotstruct.field)

        FOR ic=0L, nclr-1 DO BEGIN 

            ;; read in asTrans file
            trans=sdss_read('astrans',run,camcol,band=clr[ic],rerun=rerun,$
                node=node, inc=inc, status=astatus, /silent)
            IF astatus ne 0 THEN BEGIN
                message,'Could not read astrans'
            ENDIF 

            fields = trans.field

            mvert = fltarr(nf)
            mhorz = fltarr(nf)
            angle = fltarr(nf)
            FOR fi=0L, nf-1 DO BEGIN 
                field = fields[fi]

                ;print,'field:',field

                ;; convert center (row,col) to (mu,nu) great circle coords
                rowcol2munu, trans, field, crow, ccol, cmu, cnu
                ;; convert (mu,nu) to (clambda,ceta) corrected survey coords.
                gc2csurvey, cmu, cnu, node, inc, clam, ceta

                ;; same for the vertical/horizontal points
                rowcol2munu, trans, field, vert[1,*], vert[0,*], vmu, vnu
                rowcol2munu, trans, field, horz[1,*], horz[0,*], hmu, hnu
                gc2csurvey, vmu, vnu, node, inc, vlam, veta
                gc2csurvey, hmu, hnu, node, inc, hlam, heta

                ;; tangent project
                astr.crval = [ceta, clam]
                rd2xy,ceta,clam,astr,cxx,cyy
                rd2xy, veta, vlam, astr, vxx, vyy
                rd2xy, heta, hlam, astr, hxx, hyy

                ;; angle of points in (lambda,eta)
                lameta_a_ang = atan(vyy, vxx)
                lameta_b_ang = atan(hyy, hxx)

                ;mvert[fi] = mean(!dpi/2. - lameta_a_ang)
                ;mhorz[fi] = mean(0. - lameta_b_ang)
                ;angle[fi] = mean( [!dpi/2. - lameta_a_ang, 0. - lameta_b_ang] )

                ; can get outliers sometimes, near the poles?
                mvert[fi] = median(!dpi/2. - lameta_a_ang)
                mhorz[fi] = median(0. - lameta_b_ang)
                angle[fi] = median( [!dpi/2. - lameta_a_ang, 0. - lameta_b_ang] )

                ;if angle[fi]*180/!dpi gt 3 then begin
                ;if abs(sin(2*angle[fi])) gt 0.2 then begin
                ;	if fi gt 0 then begin
                ;		angle[fi] = angle[fi-1]
                ;	endif else begin
                ;		message,'cannot deal with bad angle at beginning'
                ;	endelse
                ;endif
            ENDFOR ; over fields


            rotstruct.hangle[ic] = mhorz
            rotstruct.vangle[ic] = mvert
            rotstruct.angle[ic] = angle
            rotstruct.cos2angle[ic] = cos(2*angle)
            rotstruct.sin2angle[ic] = sin(2*angle)


        ENDFOR ;; over bandpasses

        ptrlist[camcol-1] = ptr_new(rotstruct)

    ENDFOR ;; over camcols
    status=0

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; output file
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
    rotstruct = combine_ptrlist(ptrlist)
    print,'Writing to output file: ',outfile, f='(a,a)'
    mwrfits, rotstruct, outfile, /create
    if keyword_set(pgsql) then begin
        print,'Writing  to pgsql file: ',pgsql_file, f='(a,a)'
        ascii_write, rotstruct, pgsql_file, /bracket_arrays, append=append
    endif


	; make a plot

    data = mrdfits(outfile, 1)
	sin2angle = data.sin2angle[2]

	;data.angle = data.angle*180d/!dpi

	;w=where(data.angle[2] gt 20)
	;data[w].angle[2] = data[w].angle[2] - 180d
	

	psfile=self->psfile('survey',run,rerun=rerun)
    dir = file_dirname(psfile)
    if not file_test(dir) then file_mkdir, dir
	begplot,psfile,/color,xsize=15,ysize=8.5,/encap
	;!p.thick=1

	xrange=[min(data.field), max(data.field)]
	yrange=[min(sin2angle),max(sin2angle)]

	pcolors=make_rainbow(6)
	for camcol=1,6 do begin
		w=where(data.camcol eq camcol)
		if camcol eq 1 then begin
			xtitle='field'
			ytitle='sin(2*theta)'
			pplot,[0],/nodata,xrange=xrange,yrange=yrange, $
				xtitle=xtitle, ytitle=ytitle
		endif
		pdata = sin2angle
		pplot, /over, $
			data[w].field, pdata[w], color=pcolors[camcol-1], $
			psym=-8,symsize=0.25
	endfor
	plegend,'camcol: '+string([1,2,3,4,5,6],f='(i0)'), $
		line=0, $
		color=pcolors
	plegend,run2string(run)+'-'+strn(rerun),/right

	endplot,/trim,/png

  return
END 


pro sdss_survey_rot::calcrot_eq, run, rerun=rerun, $
                    indir=indir, outdir=outdir, status=status, $
					pgsql=pgsql

    ;; status is 1 unless we reach end
    status=1

    if n_params() lt 1 then begin 
        on_error, 2
        print,'-Syntax: obj->calcrot_eq, run, rerun=, status=, indir=, outdir=,/pgsql,status='
        message,'Halting'
    endif 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Output file  names
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if n_elements(rerun) eq 0 then rerun=301


    if n_elements(outdir) eq 0 then begin 
        outdir = self->dir('eq')
    endif 

    outfile = self->file('eq', run,rerun=rerun, dir=outdir)
    pgsql_file = self->file('eq', run,rerun=rerun,/postgres, dir=outdir)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; astrometry structure
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    radeg = 180./!dpi
    val = .4d/3600.
    astr={cd: double(identity(2)),$
        cdelt: [val,val], $
        crpix: dblarr(2), $
        crval: [0., 0.]}

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; define set of pixels to rotate
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    lrow = 1489-1
    lcol = 2048-1
    nnr=20.
    vert=dblarr(2,nnr)
    horz=dblarr(2,nnr)
    ccol = 1023.5
    crow = 744.0
    FOR i=0L, nnr-1 DO BEGIN
        vert[0,i] = ccol          ;columns
        vert[1,i] = lrow-i        ;rows
        IF i EQ 10 THEN vert[1,i] = 1000
        horz[0,i] = lcol - i      ;columns
        horz[1,i] = crow          ;rows
    ENDFOR 

    ; this is worse?
    ;vert[0,*] = ccol; columns
    ;vert[1,*] = arrscl(findgen(nnr), crow, lrow) ; rows

    ;horz[0,*] = arrscl(findgen(nnr), ccol, lcol) ; columns
    ;horz[1,*] = crow ; rows

    ;; Output structure

    clr = [0,1,2,3,4]
    colors = ['u','g','r','i','z']
    nclr = n_elements(clr)

    rotst = self->structdef(clr)


    ttrans = sdss_read('astrans', run, 1, band=clr[0], rerun=rerun, $
        node=tnode, inc=tinc)

    nf = n_elements(ttrans.field)

    ptrlist = ptrarr(6)
    FOR camcol=1,6 DO BEGIN 

        rotstruct = replicate(rotst, nf)
        rotstruct.run = run
        rotstruct.rerun = rerun
        rotstruct.camcol = camcol
        rotstruct.field = ttrans.field

        rotstruct.fieldid = sdss_photoid(rotstruct.run,   $
            rotstruct.rerun, $
            rotstruct.camcol,$
            rotstruct.field)

        FOR ic=0L, nclr-1 DO BEGIN 

            ;; read in asTrans file
            trans=sdss_read('astrans',run,camcol,band=clr[ic],rerun=rerun,$
                node=node, inc=inc, status=astatus, /silent)
            IF astatus ne 0 THEN BEGIN
                message,'Could not read astrans'
            ENDIF 

            fields = trans.field

            mvert = fltarr(nf)
            mhorz = fltarr(nf)
            angle = fltarr(nf)
            FOR fi=0L, nf-1 DO BEGIN 
                field = fields[fi]

                ;print,'field:',field

                ;; convert center (row,col) to (ra,dec)
                rowcol2munu, trans, field, crow, ccol, cmu, cnu
                rowcol2eq, trans, node, inc, field, crow, ccol, cenra, cendec


                ;; same for the vertical/horizontal points
                rowcol2eq, trans, node, inc, field, vert[1,*], vert[0,*], vra, vdec
                rowcol2eq, trans, node, inc, field, horz[1,*], horz[0,*], hra, hdec

                ;; tangent project
                ; eta is longitude
                ;astr.crval = [ceta, clam]
                ;rd2xy, ceta,clam,astr,cxx,cyy
                ;rd2xy, veta, vlam, astr, vxx, vyy
                ;rd2xy, heta, hlam, astr, hxx, hyy
                astr.crval = [cenra, cendec]
                rd2xy, cenra, cendec, astr, cenxx, cenyy
                rd2xy, vra,   vdec,   astr, vxx,   vyy
                rd2xy, hra,   hdec,   astr, hxx,   hyy

                ;; angle of points in (ra,dec)
                eq_a_ang = atan(vyy, vxx)
                eq_b_ang = atan(hyy, hxx)

                ; can get outliers sometimes, near the poles?
                mvert[fi] = median(!dpi/2. - eq_a_ang)
                mhorz[fi] = median(0. - eq_b_ang)
                angle[fi] = median( [!dpi/2. - eq_a_ang, 0. - eq_b_ang] )


            ENDFOR ; over fields


            rotstruct.hangle[ic] = mhorz
            rotstruct.vangle[ic] = mvert
            rotstruct.angle[ic] = angle
            rotstruct.cos2angle[ic] = cos(2*angle)
            rotstruct.sin2angle[ic] = sin(2*angle)


        ENDFOR ;; over bandpasses

        ptrlist[camcol-1] = ptr_new(rotstruct)

    ENDFOR ;; over camcols
    status=0

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; output file
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
    rotstruct = combine_ptrlist(ptrlist)
    print,'Writing to output file: ',outfile, f='(a,a)'
    mwrfits, rotstruct, outfile, /create
    if keyword_set(pgsql) then begin
        print,'Writing  to pgsql file: ',pgsql_file, f='(a,a)'
        ascii_write, rotstruct, pgsql_file, /bracket_arrays, append=append
    endif


	; make a plot

    data = mrdfits(outfile, 1)
	sin2angle = data.sin2angle[2]

	;data.angle = data.angle*180d/!dpi

	;w=where(data.angle[2] gt 20)
	;data[w].angle[2] = data[w].angle[2] - 180d
	

	psfile=self->psfile('eq',run,rerun=rerun)
    dir = file_dirname(psfile)
    if not file_test(dir) then file_mkdir, dir
	begplot,psfile,/color,xsize=15,ysize=8.5,/encap
	;!p.thick=1

	xrange=[min(data.field), max(data.field)]
	yrange=[min(sin2angle),max(sin2angle)]

	pcolors=make_rainbow(6)
	for camcol=1,6 do begin
		w=where(data.camcol eq camcol)
		if camcol eq 1 then begin
			xtitle='field'
			ytitle='sin(2*theta)'
			pplot,[0],/nodata,xrange=xrange,yrange=yrange, $
				xtitle=xtitle, ytitle=ytitle
		endif
		pdata = sin2angle
		pplot, /over, $
			data[w].field, pdata[w], color=pcolors[camcol-1], $
			psym=-8,symsize=0.25
	endfor
	plegend,'camcol: '+string([1,2,3,4,5,6],f='(i0)'), $
		line=0, $
		color=pcolors
	plegend,run2string(run)+'-'+strn(rerun),/right

	endplot,/trim,/png

  return
END 




function sdss_survey_rot::structdef, clr

  nclr = n_elements(clr)
  arrval = fltarr(nclr)
  rotst=create_struct('fieldid', 0ULL,$
                      'run',0L,$
                      'rerun',0,$
                      'camcol',0b,$
                      'field',0,$
                      'hangle', arrval, $
                      'vangle', arrval, $
                      'angle', arrval, $
					  'cos2angle', arrval, $
					  'sin2angle', arrval)

  return,rotst

end 



PRO sdss_survey_rot::run_calcrot_old, num

  run_status = sdss_runstatus()
  w = where(run_status.tsobj_photo_v GE 5.4, nw)

  eachdef = nw/4
  remainder = nw - (eachdef*4)
  CASE num OF
      1: w = w[0:eachdef-1]
      2: w = w[eachdef:2*eachdef-1]
      3: w = w[2*eachdef:3*eachdef-1]
      4: w = w[3*eachdef:nw-1]
      ELSE: message,'unknown number: '+ntostr(num)
  ENDCASE 

  nruns = n_elements(w)
  FOR i=0L, nruns-1 DO BEGIN 

      run = run_status[w[i]].run
      rerun = run_status[w[i]].rerun

      self->calcrot, run, rerun=rerun

  ENDFOR 

END 




PRO sdss_survey_rot::princeton_run_calcrot

    message,'need to fix paths'
  indir = '/net/cheops1/data7/imaging.local/princeton_astrom/'
  outdir = '/net/cheops1/data7/imaging.local/princeton_surveyrot/'
  files = file_search(indir+'*', count=nf)
  FOR i=0L, nf-1 DO BEGIN 
      dirsep, files[i], tdir, tf
      run = long((strsplit( (strsplit(tf,'-',/extract))[1], '.', /ext))[0])
      rerun = 137

      self->calcrot, run, rerun=rerun, indir=indir, outdir=outdir
  ENDFOR 

END 

pro sdss_survey_rot::stuff, princeton=princeton

    table = 'field_rotation'
    dir = self->dir(system, princeton=princeton)
    pattern = 'surveyrot_*_*.st'

    pattern = concat_dir(dir, pattern)
    files = file_search(pattern, count=nf)

    pg=obj_new('postgres')
    for i=0L, nf-1 do begin
        f=files[i]

        t=read_idlstruct(f)
        pg->struct2table, t, table, primary_key='fieldid', status=status
        if status ne 0 then message,'Failed'
    endfor

    query='create index '+table+'_rrcf_index on '+table+' (run,rerun,camcol,field)' 
    pg->query, query, status=status
    if status ne pg->statusval('no_result') then begin
        print,'Failed to create index'
    endif
    pg->query, 'analyze '+table

    obj_destroy, pg

end


PRO sdss_survey_rot::tabledef, sqlfile

  IF n_elements(sqlfile) EQ 0 THEN BEGIN 
      lun = -1
  ENDIF ELSE BEGIN 
      openw,lun,sqlfile,/get_lun
  ENDELSE 

  struct = self->structdef([1,2,3,4,5])

  coldefs = self->struct2coldefs(struct)

  ;; Add the primary key
  coldefs = [coldefs, 'PRIMARY KEY (fieldid)']

  ncoldefs = n_elements(coldefs)
  printf,lun,'CREATE TABLE field_rotation'
  printf,lun,'('
  FOR i=0L, ncoldefs-2 DO BEGIN 
      printf,lun,coldefs[i]+', '
  ENDFOR 
  printf,lun,coldefs[i]
  printf,lun,');'

  IF n_elements(sqlfile) NE 0 THEN free_lun,lun

END 


function sdss_survey_rot::dir, system, princeton=princeton
    dir = '~/lensing/sdss-shape-rot/'+system
	dir = expand_tilde(dir)
    return,dir
end 
function sdss_survey_rot::psdir, system
	dir=self->dir(system)
	psdir=path_join(dir,'plots')
	return,psdir
end

function sdss_survey_rot::file, system, run, dir=dir, rerun=rerun, postgres=postgres, princeton=princeton

	if n_elements(rerun) eq 0 then rerun = sdss_rerun(run)
	if n_elements(dir) eq 0 then begin 
		outdir = self->dir(system, princeton=princeton)
	endif else begin 
		outdir = dir
	endelse 
	outfile = system+'rot-'+run2string(run)+'-'+ ntostr(rerun)

	outfile=path_join(outdir, outfile)
	if keyword_set(postgres) then begin
		ext = '.pgsql'
	endif else begin
		;ext = '.st'
		ext = '.fits'
	endelse

	outfile = outfile + ext
	return,outfile
end 
function sdss_survey_rot::psfile, system, run, rerun=rerun
	psdir=self->psdir(system)
	file=file_basename(self->file(system, run,rerun=rerun))
	file=repstr(file,'.fits','.eps')
	file=path_join(psdir, file)
	return, file
end

function sdss_survey_rot::read, system, run, rerun=rerun, dir=dir, status=status,$
  princeton=princeton, silent=silent
  
  status = 1
  if n_elements(system) eq 0 or n_elements(run) eq 0 then begin 
      print,'-Syntax: rot = ss->read(system, run, rerun=, /princeton, dir=, /silent, status=)'
      message,'stop'
  endif 
  file = self->file(system, run, rerun=rerun, princeton=princeton, dir=dir)
  if not fexist(file) then begin 
      print,'File not found: ',file
      return,-1
  endif 
  if not keyword_set(silent) then begin 
      print
      print,'Reading rotation struct: ',file
  endif 
  ;struct = read_idlstruct(file, status=status, silent=silent)
  struct = mrdfits(file, 1, status=status, silent=silent)
  return,struct
end 






PRO sdss_survey_rot__define

  struct = {$
             sdss_survey_rot, $
             dummy: 0, $
             INHERITS postgres $
           }

END 
