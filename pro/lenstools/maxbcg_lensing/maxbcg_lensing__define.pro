;+
; NAME:
;  MAXBCG_LENSING__DEFINE
;
;
; PURPOSE:
;  Class file containing routines to make the lensing inputs from the maxbcg
;  catalog or random points.  Also, will contain routines for running the
;  lensing code.
;
;  sub-sampling:
;       See maxbcg_lensing_runsubs
;
;    SIGMACRIT_STYLES:
;       1: treat photozs as truth
;       2: use the mean inverse critical density calculated from the 
;          deconvolved overall distribution. 
;       3: Interpolate to get the mean given zl, zs
;       4: Use pre-calculated siginv vs zl on grid for each source gal
;
; INHERITED CLASSES:
;   MAXBCG
;   OBJSHEAR
;
; MODIFICATION HISTORY:
;  Created: Early May, 2005 from separate programs.  
;  Author: Erin Sheldon, UofChicago
;
;-

function maxbcg_lensing::init, sample

	if n_elements(sample) eq 0 then begin 
		message,'You must send the sample number',/inf
		message,"ml = obj_new('maxbcg_lensing', sample)"
	endif 

	;; The other defaults are mostly set in objshear::init.  See call below
	rmin = 20.0
	logbin = 1

	sample = fix(sample)

	type = 'maxbcg'
	par = $
		{maskfile: '',      $
		sample: sample,    $
		catalog: '',       $
		source_sample:  1, $
		random_sample: -1, $
		scinv_sample:  -1, $
		dopairs: 0}


	case sample of 
		1: begin 
			;; 20kpc -10 Mpc, using photozs as truth
			;; dr3_final_bcgs_2_16_05_fix.st
			rmax = 11500.
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 2 ; princeton analytic
			par.catalog = 'dr3'
		end 
		2: begin 
			;; 20kpc -30 Mpc, using photozs as truth
			;; dr3_final_bcgs_2_16_05_fix.st
			rmax = 36567.0
			nbin = 21
			sigmacrit_style = 1
			shape_correction_style = 2 ; princeton analytic
			par.catalog = 'dr3'
		end 
		3: begin 
			;; 20kpc -10 Mpc, using overall photoz distribution
			;; also uses newer config files
			;; dr3_final_bcgs_2_16_05_fix.st
			rmax = 11500.
			nbin = 18
			sigmacrit_style = 2
			shape_correction_style = 2 ; princeton analytic
			par.catalog = 'dr3'
		end
		4: begin 
			;; using dr3_final_bcgs_2_16_05_fix.st and
			;; PRINCETON source shapes. Another update of config file
			;; to include this
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton re-gaussianization
			par.catalog = 'dr3'
		end 

		5: begin
			;; same spacing, but smaller radius, 14 bins.  Better for comparisons
			;; with princeton
			rmax = 2801.8093
			nbin = 14
			sigmacrit_style = 1
			shape_correction_style = 2 ; princeton analytic
			;; Maskfile for princeton is actually taken care of in
			;; objshear through knowing shape_correction_style.  But for
			;; others you will have to put it here.
			par.maskfile = sdssidl_config('PIXEL_MASK_PRINCETON_BASIC')
			par.catalog = 'dr3'
		end 
		6: begin
			;; same as 5 but with new princeton corrections
			rmax = 2801.8093
			nbin = 14
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton re-gaussianization
			par.catalog = 'dr3'
		end 
		7: begin
			;; test catalog with still discrete photozs
			rmax = 2801.8093
			nbin = 14
			sigmacrit_style = 1
			shape_correction_style = 2 ; princeton analytic
			par.catalog = 'dr4plus'
		END 
		8: BEGIN
			;; test catalog with still discrete photozs
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 2 ; princeton analytic
			par.catalog = 'dr4plus'
		END 
		9: BEGIN
			;; test catalog with still discrete photozs
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr4plus'
		END 
		10: BEGIN
			;; Still discrete
			rmax = 36567.0
			nbin = 21
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr4plus'
		END 



		;; Only those below here are on new naming/directory scheme
		11: BEGIN 
			;; Continuous photozs
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 9 ; The 10 Mpc sample
		END 
		12: BEGIN 
			;; Continuous photozs
			rmax = 36567.0
			nbin = 21
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 10 ; The 10 Mpc sample
		END 
		13: BEGIN 
			;; Continuous photozs, cut source catalog
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 9 ; The 10 Mpc sample
			par.source_sample = 2 ; r > 0.585
		END 
		14: BEGIN 
			;; Continuous photozs, cut source catalog
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 9 ; The 10 Mpc sample
			par.source_sample = 6 ; r > 1/3
		END 
		15: BEGIN 
			;; Continuous photozs, cut source catalog
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 9 ; The 10 Mpc sample
			par.source_sample = 8 ; r > 0.4
		END 
		16: BEGIN 
			;; Continuous photozs, cut source catalog
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 9 ; The 10 Mpc sample
			par.source_sample = 7 ; r > 2/3
		END 

		17: BEGIN 
			;; noref matched to maxbcg (noref positions, redshifts)
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'n2m'
			par.random_sample = 9 ; The 10 Mpc sample
			par.source_sample = 1 ; old cuts just for a check
		END 
		18: BEGIN 
			;; noref matched to maxbcg (noref positions, redshifts)
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'm2n'
			par.random_sample = 9 ; The 10 Mpc sample
			par.source_sample = 1 ; old cuts just for a check
		END 

		;; Ben's new catalog
		19: BEGIN 
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 4 ; The 10 Mpc sample
			par.source_sample = 6 ; r > 1/3
		END 
		20: BEGIN 
			rmax = 36567.0
			nbin = 21
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 5 ; The 30 Mpc sample
			par.source_sample = 6 ; r > 1/3
		END 


        ;; 21-22 were used in the papers
		;; new photozs, same numbers as before note
		21: BEGIN 
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 11 ; The 10 Mpc sample
			par.source_sample = 6 ; r > 1/3
		END 
		22: BEGIN 
			rmax = 36567.0
			nbin = 21
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 12 ; The 10 Mpc sample
			par.source_sample = 6 ; r > 1/3
		END 


		; For doing pairs
		23: BEGIN 
			rmax = 2801.8093
			nbin = 14
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 13 ; The 10 Mpc sample
			par.source_sample = 6 ; r > 1/3
			par.dopairs = 1
		END 

		; testing new mean inverse critical density stuff
		24: BEGIN 
			rmax = 2801.8093
			nbin = 14
			sigmacrit_style = 3
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 14 ; The 10 Mpc sample
			par.source_sample = 6 ; r > 1/3
			par.dopairs = 1
			par.scinv_sample = 2
		END 

		;; new photozs "dr6cc2".
		25: BEGIN 
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'dr406'
			par.random_sample = 11 ; The 10 Mpc sample
			par.source_sample = 9 ; r > 1/3
		END 

		; The new maxbcg catalog.  Will re-do with new source catalog later
		26: begin
			rmax = 11500.0
			nbin = 18
			sigmacrit_style = 1
			shape_correction_style = 3 ; princeton regauss
			par.catalog = 'gmbcg10'
			par.random_sample = 15 ; The 10 Mpc sample z=[0.1,0.35]
			par.source_sample = 6 ; r > 1/3
		end

		ELSE: message,'Unknown sample: '+ntostr(sample)
	ENDCASE 

	self.sample = sample

	print,'type   =                 ',type
	print,'sample =                 ',sample,format='(A,I7)'
	print,'source_sample =          ',par.source_sample,format='(A,I7)'
	print,'random_sample =          ',par.random_sample,format='(A,I7)'
	print,'rmin =                   ',rmin
	print,'rmax =                   ',rmax
	print,'nbin =                   ',nbin
	print,'sigmacrit_style =        ',sigmacrit_style,format='(A,I7)'
	print,'shape_correction_style = ',shape_correction_style,format='(A,I7)'
	print,'scinv_sample =           ',par.scinv_sample
	maxbcg_retval = self->maxbcg::init(par.catalog)
	objshear_retval = $
		self->objshear::init(type, rmin, rmax, nbin, $
		sigmacrit_style, shape_correction_style, $
		logbin=logbin, par_struct=par)

	return,( (maxbcg_retval + objshear_retval) EQ 2 )


END 













































;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; this is not a lensing edge cut, but ben's 1Mpc edge cut
;; This cut is currently more harsh than ben's because he ignored
;; the large holes in the survey
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function maxbcg_lensing::ben_edgecuts, clambda, ceta, z, nkeep

    par = self->par_struct()
  
    ;; this will return a defined maskfile in certain situations. if undefined,
    ;; then default /basic will be used below
    self->maskfile, maskfile, /edgecuts

    rmpc = 1.0
  
    ;; convert to degrees for each lens
    maxangle = rmpc/angdist_lambda(z, omega=par.omega_m)*180.0/!pi

    apply_pixel_mask, $
        clambda, ceta, masked, unmasked, $
        maxangle = maxangle, /basic, maskfile=maskfile

    if unmasked[0] eq -1 then nkeep = 0 else nkeep = n_elements(unmasked)
  
    return,unmasked

end 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generate the basic random points.  No lensing cuts here.  Mask
;; depends upon the source catalog
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;FUNCTION maxbcg_lensing::randname, randnum
;  IF n_elements(randnum) EQ 0 THEN message,'n=mb->randname(randnum)'
;  dir = self->lensinput_dir(/random)
;  return,dir + 'random'+self->objshear::randstring(randnum)+'.st'
;END 











;; Run all the code in order to set up this sample
PRO maxbcg_lensing::setupall, rand=rand

  IF keyword_set(rand) THEN BEGIN 
      ;; Run this and lenses AND randomn will be done
      self->setuplens, randnum=self->randnum()

      self->plot_lensrand, /dops
  ENDIF ELSE BEGIN 

      self->setuplens
  
      self->write_scripts
      self->write_parfiles
  ENDELSE 

END 

;; Because of the ben_edgecut, is lens sample dependent
PRO maxbcg_lensing::setuplens
  
	tm = systime(1)

	; Parameters
	par_struct = self->objshear::par_struct()

	; output file name(s)
	outfile = self->lensfile('input', /createdir)

	;; make sure dir exists
	tfile = self->lensfile('output', /createdir)
	print
	print,'Setting up file: ',outfile

	; get the data.  Note: zindex will be the same for all samples
	; created from this lens catalog
	bcg = self->get()
	nbcg = n_elements(bcg)
	zindex = lindgen(nbcg)

	nlens_init = n_elements(bcg)      

	print,'Creating lens struct'
	lstruct = self->objshear::lensinput_structdef()
	lstruct = replicate(lstruct, nlens_init)

	print,'Copying....'
	copy_struct, bcg, lstruct


	if not tag_exist(bcg, 'clambda') then begin 
		print,'Adding csurvey coords'
		eq2csurvey, bcg.ra, bcg.dec, clam, ceta
		lstruct.clambda = clam
		lstruct.ceta = ceta
	endif 
 


	;; use continuous z if available
	if tag_exist(bcg, 'photoz_cts') then begin 
		print,'Using tag "photoz_cts" for redshift'
		lstruct.z=bcg.photoz_cts
	endif else if tag_exist(bcg, 'photoz') then begin
		print,'Using tag "photoz" for redshift'
		lstruct.z = bcg.photoz
	endif else if tag_exist(bcg, 'z') then begin 
		; this already got copied in above in the copy_struct
		print,'Using tag "z" for redshift'
	endif else begin 
		message,'Neither PHOTOZ_CTS or PHOTOZ or Z found'
	endelse 



	;; Keep track of the objects and their redshifts
	lstruct.zindex = zindex

	;; Calculate some cosmology-dependent stuff and copy into struct
	self->objshear::calc_cosmo, lstruct

	;; Make some generic lens cuts
	wlens = self->objshear::lenscuts(lstruct, nkeep)

	;; remove the unwanted lenses
	lstruct = lstruct[wlens]

	print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'


	print
	print,'Applying maxbcg 1Mpc cut'
	keep = self->ben_edgecuts(lstruct.clambda,lstruct.ceta,lstruct.z,nkeep)
	lstruct = lstruct[keep]
	print,'Finally using '+ntostr(nkeep)+'/'+ntostr(nlens_init)

	print
	print,'Writing to file: ',outfile
	write_idlstruct, lstruct, outfile

	ptime,systime(1)-tm

end 


;; Because of the ben_edgecut, is lens sample dependent
PRO maxbcg_lensing::setup_paper_catalog
  
  tm = systime(1)

  ;; Parameters
  par_struct = self->objshear::par_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file name(s)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  ;; get the data.  Note: zindex will be the same for all samples
  ;; created from this lens catalog
  dir=self->catdir()
  dir=concat_dir(dir,'paper_catalog')
  infile=concat_dir(dir, 'maxbcg_cluster_catalog_no_bs.fit')
  outfile=concat_dir(dir, 'maxbcg_cluster_catalog_no_bs_cut.fit')
  bcg = mrdfits(infile,1)
  nbcg = n_elements(bcg)
  zindex = lindgen(nbcg)

  print
  print,'Setting up file: ',outfile

  nlens_init = n_elements(bcg)      

  print,'Creating lens struct'
  lstruct = self->objshear::lensinput_structdef()
  lstruct = replicate(lstruct, nlens_init)
  
  print,'Copying....'
  copy_struct, bcg, lstruct
  
  ;; use continuous z if available
  IF tag_exist(bcg[0], 'photoz_cts') THEN BEGIN 
      print,'Copying continuous redshift'
      lstruct.z=bcg.photoz_cts
  ENDIF ELSE IF tag_exist(bcg[0], 'z') THEN BEGIN 
      print,'PHOTOZ_CTS not found, using z'
  ENDIF ELSE BEGIN 
      message,'Neither PHOTOZ_CTS or Z found'
  ENDELSE 
  
  ;; Keep track of the objects and their redshifts
  lstruct.zindex = zindex
  
  ;; Calculate some cosmology-dependent stuff and copy into struct
  self->objshear::calc_cosmo, lstruct
  
  ;; Make some generic lens cuts
  wlens = self->objshear::lenscuts(lstruct, nkeep)
  
  ;; remove the unwanted lenses
  help,wlens,bcg
  bcg=bcg[wlens]
  lstruct = lstruct[wlens]

  print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'
  
  print
  print,'Applying maxbcg 1Mpc cut'
  keep = self->ben_edgecuts(lstruct.clambda,lstruct.ceta,lstruct.z,nkeep)
  lstruct = lstruct[keep]
  bcg=bcg[keep]
  print,'Finally using '+ntostr(nkeep)+'/'+ntostr(nlens_init)
  
  print
  print,'Writing to file: ',outfile
  mwrfits, bcg, outfile, /create

  ptime,systime(1)-tm

END 




















;; Because of the ben_edgecut, is lens sample dependent
PRO maxbcg_lensing::setuplens_old, randnum=randnum
  
  tm = systime(1)

  ;; Parameters
  par_struct = self->objshear::par_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file name(s)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outfile = self->lensfile('input', /createdir)

  print
  print,'Setting up file: ',outfile

  ;; get the data.  Note: zindex will be the same for all samples
  ;; created from this lens catalog
  bcg = self->get()
  nbcg = n_elements(bcg)
  zindex = lindgen(nbcg)

  nrand = n_elements(randnum)
  IF nrand EQ 0 THEN BEGIN 

      nlens_init = n_elements(bcg)      

      print,'Creating lens struct'
      lstruct = self->objshear::lensinput_structdef()
      lstruct = replicate(lstruct, nlens_init)
      
      print,'Copying....'
      copy_struct, bcg, lstruct

      ;; use continuous z if available
      IF tag_exist(bcg[0], 'photoz_cts') THEN BEGIN 
          print,'Copying continuous redshift'
          lstruct.z=bcg.photoz_cts
      ENDIF ELSE BEGIN 
          print,'No continuous redshift found'
      ENDELSE 

      ;; Keep track of the objects and their redshifts
      lstruct.zindex = zindex

      ;; Calculate some cosmology-dependent stuff and copy into struct
      self->objshear::calc_cosmo, lstruct
      
      ;; Make some generic lens cuts
      wlens = self->objshear::lenscuts(lstruct, nkeep)

      ;; remove the unwanted lenses
      lstruct = lstruct[wlens]

      print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'

      print
      print,'Applying maxbcg 1Mpc cut'
      keep = self->ben_edgecuts(lstruct.clambda,lstruct.ceta,lstruct.z,nkeep)
      lstruct = lstruct[keep]
      print,'Finally using '+ntostr(nkeep)+'/'+ntostr(nlens_init)

      print
      print,'Writing to file: ',outfile
      write_idlstruct, lstruct, outfile

  ENDIF ELSE BEGIN 

      bcgz = bcg.z
      nbcgz = n_elements(bcg)

      ;; clear bcg memory
      bcg = 0

      FOR i=0L, nrand-1 DO BEGIN 

          tmi = systime(1)

          print
          rand = self->genrand()

          outfile = self->lensinput_file(randnum=randnum[i])
          print,'*********************************************************'
          print,'Output file: ',outFile

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

          bcg_index = long( nbcgz*randomu(seed, nlens_init) )

          lstruct.z = bcgz[bcg_index]
          lstruct.zindex = zindex[bcg_index]

          ;; Calculate some stuff and copy int
          self->objshear::calc_cosmo, lstruct

          ;; Make some generic lens cuts
          wlens = self->objshear::lenscuts(lstruct, nkeep)

          ;; remove the unwanted lenses
          lstruct = lstruct[wlens]

          print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'
          
          print
          print,'Applying maxbcg 1Mpc cut'
          keep = self->ben_edgecuts(lstruct.clambda, lstruct.ceta, lstruct.z, nkeep)
          lstruct = lstruct[keep]
          print,'Finally using '+ntostr(nkeep)+'/'+ntostr(nlens_init)

          print
          print,'Writing to file: ',outfile
          write_idlstruct, lstruct, outfile

          print,'Time for this random'
          ptime,systime(1)-tmi

      ENDFOR 
      
  ENDELSE 

  ptime,systime(1)-tm

END 










  












PRO maxbcg_lensing::setuplens_region, region, randnum=randnum

  ;; Parameters
  par_struct = self->objshear::par_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file name(s)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  outfile = self->lensinput_file(region, randnum=randnum)

  bcg = self->get()
  nbcg = n_elements(bcg)
  zindex = lindgen(nbcg)
  regions = self->objshear::ceta2region(bcg.ceta)

  nrand = n_elements(randnum)
  IF nrand EQ 0 THEN BEGIN 

      nlens_init = n_elements(bcg)      

      print,'Creating lens struct'
      lstruct = self->objshear::lensinput_structdef()
      lstruct = replicate(lstruct, nlens_init)
      
      print,'Copying....'
      copy_struct, bcg, lstruct

      ;; Keep track of the objects and their redshifts
      lstruct.zindex = zindex


      ;; Now select subset in the region of interest
      w = where(regions EQ region)
      lstruct = lstruct[w]

      ;; Calculate some cosmology-dependent stuff and copy into struct
      self->objshear::calc_cosmo, lstruct
      
      ;; Make some generic lens cuts
      wlens = self->lenscuts(lstruct, nkeep)

      ;; remove the unwanted lenses
      lstruct = lstruct[wlens]

      print,'Finally using '+ntostr(nkeep)+'/'+ntostr(nlens_init)

      print
      print,'Writing to file: ',outfile
      write_idlstruct, lstruct, outfile

  ENDIF ELSE BEGIN 

      ;; Since randoms generated in same region, cut out
      ;; all but input region
      w = where(regions EQ region)
      bcg = bcg[w]
      zindex = zindex[w]

      bcgz = bcg.z
      nbcgz = n_elements(bcg)

      ;; clear bcg memory
      bcg = 0

      FOR i=0L, nrand-1 DO BEGIN 

          print
          rand = self->get(randnum=randnum[i], region=region)

          outfile = self->lensinput_file(region,randnum=randnum[i])
          print,'Output file: ',outFile

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

          bcg_index = long( nbcgz*randomu(seed, nlens_init) )

          lstruct.z = bcgz[bcg_index]
          lstruct.zindex = zindex[bcg_index]

          ;; Calculate some stuff and copy int
          self->objshear::calc_cosmo, lstruct

          ;; Make some generic lens cuts
          wlens = self->lenscuts(lstruct, nkeep)

          ;; remove the unwanted lenses
          lstruct = lstruct[wlens]

          print,'Finally using '+ntostr(nkeep)+'/'+ntostr(nlens_init)

          print
          print,'Writing to file: ',outfile
          write_idlstruct, lstruct, outfile

      ENDFOR 
      
  ENDELSE 

END 














; Plothist ngals

PRO maxbcg_lensing::plothist_ngals, struct

  begplot, name='~/plots/MaxBCG/tngals_function_'+self->maxbcg::catalog()+'.eps',/encap, /color
  
  data_clr = c2i('blue')
  phi_clr = c2i('red')

  xtitle = 'N!S!Dgals!N!R!Ut!N'
  ytitle = 'Number'
  xrange = [1, 500]
  yrange = [1, 6.e5]
  aplot, 1, [5], xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
    xtitle=xtitle, ytitle=ytitle, $
    /nodata, /xlog, /ylog, ytickf='loglabels'
  plothist, struct.tngals, Ngals, number, min=3, color=data_clr, /overplot

;  ysmooth = smooth(yhist, 5)

;  w = where(xhist GT 50 AND xhist LT 175)
;  yhist[w] = smooth(yhist[w], 5)
;  w = where(xhist LT 175)
;  Ngals = xhist[w]
;  number = yhist[w]
;  oplot, Ngals, number, color=c2i('green')
  
  Nstar = 35
  norm = 4.5e3
  alpha = -2
  phi = norm*(Ngals/Nstar)^alpha*exp(-(Ngals/Nstar)^1.3)
  oplot, Ngals, phi, color=phi_clr

  mess = ['Data',$
          'Fit: '+!csym.alpha+' = '+ntostr(alpha)+' N* = '+ntostr(Nstar)]
  legend, mess, color=[data_clr, phi_clr], line=[0,0], $
    /right, box=0, charsize=1

;, yrange=[0.1, 6.e5], ystyle=1, $
;            /ylog, /xlog, min=3;, xrange=[0.01, max(struct.tngals)], xstyle=1+2

  endplot, /trim_bbox

END 












;; add this to the subsample code
PRO maxbcg_lensing::ngals_calcmean, nbin, bcg=bcg, lensout=lensout, $
  tagname=tagname

  IF n_elements(nbin) EQ 0 THEN BEGIN 
      print,'-Syntax: mb->ngals_calcmean, nbin, bcg=, lensout=, tagname='
      return
  ENDIF 

  IF n_elements(lensout) EQ 0 THEN lensout = self->lensoutput_get()

  nlens = n_elements(lensout)

  zindex = lensout.zindex

  IF n_elements(bcg) EQ 0 THEN BEGIN 
      bcg = self->get()
      bcg = bcg[zindex]
  ENDIF 

  self->maxbcg::ngals_bins, nbin, lowlim, highlim
  zcuts = self->zcuts()


  ;; Get arrays of indices for each ngals cut
  keeparray = self->ngals_select(bcg, nbin, nkeeparray, tagname=tagname)

  IF n_elements(tagname) NE 0 THEN BEGIN 
      tag=where(tag_names(bcg) EQ strupcase(tagname), nw)
      IF nw EQ 0 THEN message,'tag not found: '+tagname
  ENDIF 

  print,'ngals_range','ncluster','ngals_mean','ngals_sdev','ngals_err',$
        format = '(A15,A15,A15,A15,A15)'
  print, '-------------------------------------------------------------------------------'
  FOR i=0L, nbin-1 DO BEGIN 
      
      IF nkeeparray[i] NE 0 THEN BEGIN 
          
          ngals_str = '['+ntostr(lowlim[i])+','+ntostr(highlim[i])+']'

          IF lowlim[i] EQ highlim[i] THEN BEGIN 
              ngals_mean = float(lowlim[i])
              ngals_sdev = 0.0
              ngals_err = 0.0

              nclust = n_elements(*keeparray[i])
          ENDIF ELSE BEGIN 
          
              w = *keeparray[i]
              nclust = n_elements(w)
              
              wmom, bcg[w].(tag), blah, ngals_mean, ngals_sdev, ngals_err,$
                    inputweight=lensout[w].weight
          ENDELSE 
          print,ngals_str,nclust,ngals_mean,ngals_sdev,ngals_err,$
                format='(A15,I15,F15,F15,F15)'

;          plothist, bcg[w].ngals, xtitle='Ngals'
;          key = prompt_kbrd(1)
          
;          IF i EQ 2 THEN stop
      ENDIF 

  ENDFOR 

  ptr_free,keeparray

END 



PRO maxbcg_lensing::zquantiles, bcg, lensum, nbin, zlowlim, zhighlim

  num = n_elements(lensum)

  zindlow = lonarr(nbin)
  zindhigh = lonarr(nbin)

  FOR i=0L, nbin-1 DO BEGIN 
      zindlow[i] = long(i*num/float(nbin))
      zindhigh[i] = long((i+1)*num/float(nbin) -1)
  ENDFOR 
  zindhigh[nbin-1] = num-1


  ;; print some stuff
  s = sort(bcg[lensum.zindex].z)
  print,'zrange   zmean    mean_tNgals'
  print,'------------------------------------------'
  FOR i=0L, nbin-1 DO BEGIN 
      
      w = s[zindlow[i]:zindhigh[i]]

      minz = min(bcg[lensum[w].zindex].z, max=maxz)


      wtot = total( lensum[w].weight )
      mean_ngals = total( lensum[w].weight*bcg[lensum[w].zindex].tNgals )/wtot
      mean_z = total( lensum[w].weight*bcg[lensum[w].zindex].z )/wtot

      print,'['+ntostr(minz)+', '+ntostr(maxz)+']' + $
            '  '+ntostr(mean_z)+'   '+ntostr(mean_ngals)

  ENDFOR 

END 

















PRO maxbcg_lensing::process_outputs, ngals_nbin

  ngals_nbin = 12
  tlum_nbin = 12

  self->tlum_sub, tlum_nbin
  self->combine_allrand, tlum_nbin=tlum_nbin
  self->runcorr, tlum_nbin=tlum_nbin
  self->jackknife, tlum_nbin=tlum_nbin

  return
  ;; The order is important. E.g. some routines try to
  ;; deal with ngals bins 

  ;; This assumes the randoms are done.

  ;; Sub-sample. Note, this already does  "combine"
  self->ngals_sub, ngals_nbin, tagname='tngals'
  self->tlum_sub, tlum_nbin

  self->combine_lensum
  self->combine_lensum, randnum=self->randnum()

  self->combine_allrand
  self->combine_allrand, ngals_nbin=ngals_nbin
  self->combine_allrand, tlum_nbin=tlum_nbin
;  self->combine_allrand, ngals_nbin=ngals_nbin, zbin_nbin=3


  ;; Correct
  self->corr
  self->runcorr, ngals_nbin=ngals_nbin
  self->runcorr, tlum_nbin=tlum_nbin

  ;; Jackknife
  self->jackknife, ngals_nbin=ngals_nbin, tagname='tngals'
  self->jackknife, tlum_nbin=tlum_nbin


END 


PRO  maxbcg_lensing::makeplots, ngals_nbin

  IF n_params() LT 1 THEN BEGIN  
      print,'-Syntax: mb->makeplots, ngals_nbin'
      return
  ENDIF 

  ;; plots
  self->ngals_plot, ngals_nbin, /dops
  self->ngals_plot, ngals_nbin, /dops, /add_uncorr
  self->ngals_plot, ngals_nbin, /dops, /uncorr
  self->ngals_plot, ngals_nbin, /dops, /jack
  self->ngals_plot, ngals_nbin, /dops, /jack, /add_uncorr

  self->ngals_plot_corr, ngals_nbin, /dops

  self->plotrand_dispersion, ngals_nbin, 8, /dops
  self->ngals_rand_plot,ngals_nbin,/dops
  self->ngals_rand_plot,ngals_nbin,/dops,/sigma

  self->ngals_plot_covariance, ngals_nbin, /dops

  ;; png hard to do from here since we need to adjust window size
  self->ngals_overplot, 'deltasig', ngals_nbin, /dops
;  self->ngals_overplot, 'mass', ngals_nbin, /dops
;  self->ngals_overplot, 'rho', ngals_nbin, /dops

  FOR ngals_bin=0,ngals_nbin-1 DO BEGIN 
      self->ngals_plot_linear,ngals_nbin,ngals_bin,/jack,/dops
      self->ngals_plot_linear,ngals_nbin,ngals_bin,/jack,/dops,/addrand
      self->ngals_plot_linear,ngals_nbin,ngals_bin,/jack,/dops, /restrict_range
      self->ngals_plot_linear,ngals_nbin,ngals_bin,/jack,/dops,/addrand, /restrict_range
  ENDFOR 


END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plotting routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO maxbcg_lensing::plot_profile, $
		type, subtype=subtype, bin=bin, sample=sample, $
		nsplit=nsplit, $
		xlog=xlog, xmin=xmin, xmax=xmax, $
		ylog=ylog, ymin=ymin, ymax=ymax, $
		incolor=incolor, dops=dops
                  
  if n_elements(xmin) eq 0 then xmin = 0.015
  IF n_elements(subtype) NE 0 THEN BEGIN 
      CASE subtype OF
          'lambda12': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.7
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
          END 
          'lambda_z_12_2': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.5
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1

              IF n_elements(ymin) EQ 0 THEN BEGIN 
                  ymin = 0.03
              ENDIF 

              nsplit = 2
          END 


          'sngals6': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.7
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
          END 



          'ngals200_8': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=1
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
          END 
          'ngals200_12': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.7
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
          END 
          'ngals200_ilum200_12_2': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.5
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1

              IF n_elements(ymin) EQ 0 THEN BEGIN 
                  ymin = 0.03
              ENDIF 

              nsplit = 2
          END 
          'ngals200_z_12_2': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.5
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1

              IF n_elements(ymin) EQ 0 THEN BEGIN 
                  ymin = 0.03
              ENDIF 

              nsplit = 2
          END 

          'ilum200_12': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=1
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
          END 
          'ilum200_16': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.7
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
              IF type EQ 'clustcorr' AND n_elements(ymin) EQ 0 THEN BEGIN 
                  ymin = 1.e-3
              ENDIF ELSE IF n_elements(ymin) EQ 0 THEN BEGIN 
                  ymin = 0.03
              ENDIF 
          END 
          'ilum200_ngals200_16_2': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.5
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1

              IF n_elements(ymin) EQ 0 THEN BEGIN 
                  ymin = 0.03
              ENDIF 

              nsplit = 2
          END 

          'ilum200_z_16_2': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.5
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1

              IF n_elements(ymin) EQ 0 THEN BEGIN 
                  ymin = 0.03
              ENDIF 

              nsplit = 2
          END 

          'centerclass6': begin
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  label_charsize=0.7
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
          end


          ELSE: message,'subtype: '+subtype+' not yet supported'
      ENDCASE  
  ENDIF 


  self->objshear::plot_profile, $
	  type, subtype=subtype, bin=bin, sample=sample, $
	  xlog=xlog, xmin=xmin, xmax=xmax, $
	  ylog=ylog, ymin=ymin, ymax=ymax, $
	  dops=dops, $
	  charsize=charsize, label_charsize=label_charsize, $
	  landscape=landscape, /encapsulated, $
	  incolor=incolor, $
	  $
	  nsplit=nsplit

END 






PRO maxbcg_lensing::ngals_plot_linear, nbin, bin, $
                  corr=corr, jack=jack, dops=dops, addrand=addrand, $
                  inset=inset, $
                  restrict_range=restrict_range, _extra=_extra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: mb->ngals_plot_bin1, ngals_nbin, ngals_bin, /corr, /jack, /dops, /addrand, /restrict_range, _extra=_extra'
      return
  ENDIF 

  plotdir = self->plotdir()+'linear/'

  restricted_yrange = [-2,2]
  IF keyword_set(restrict_range) THEN yrange=restricted_yrange

  if keyword_set(dops) then begin 
      nbinstr = strn(nbin, len=2,padchar='0')
      binstr = strn(bin, len=2, padchar='0')

      psfile = $
        concat_dir(plotdir,'ngals'+nbinstr+'_bin'+binstr+'_deltasig_linear_')

      if keyword_set(corr) then begin 
          psfile=psfile+'corr_'
      endif else if keyword_set(jack) then begin 
      	  psfile=psfile+'jack_'
      ENDIF
      
      IF keyword_set(restrict_range) THEN BEGIN 
          psfile=psfile + 'yrange_'
      ENDIF 

      IF keyword_set(addrand) THEN psfile = psfile + 'addrand_'
      IF keyword_set(inset) THEN psfile = psfile + 'inset_'
      
      psfile = psfile + self->objshear::sample_string()+'.eps'
      begplot, name=psfile, /encapsulated, /color

      rand_color = c2i('red')
  ENDIF ELSE BEGIN 
      rand_color = c2i('red')
  ENDELSE 

  sh = self->objshear::combined_get(ngals_nbin=nbin, ngals_bin=bin, corr=corr,jack=jack)

  esheldon_setup
  xtitle=estitle('mpcxtitle')
  ytitle = estitle('deltaytitle')
  aploterror, !gratio, sh.meanr/1000.0, sh.sigma, sh.sigmaerr,$ 
    ytitle=ytitle, xtitle=xtitle, $
    xtickf='loglabels', $
    psym=8, hat=0, $
    /xlog, yrange=yrange, _extra=_extra
 
;  plot_density_contrast, sh, $
;    /xlog, xrange=[10, 15000.], xstyle=1, $
;    yrange = [-5, 40], ystyle=1

  oplot, [.005, 300], [0,0]

  IF keyword_set(addrand) THEN BEGIN 
;      rsh = self->objshear::combined_get(/rand_combined)
      rsh = self->objshear::combined_get(/rand_combined, ngals_nbin=nbin, ngals_bin=bin)
      oploterror, rsh.meanr/1000, rsh.sigma, rsh.sigmaerr, $
                  psym=4, symsize=0.7, hat=0, $
                  color=rand_color, errc=rand_color
  ENDIF 
  self->maxbcg::ngals_bins, nbin, lowlim, highlim
  legend, self->maxbcg::ngals_string(lowlim[bin], highlim[bin]), $
          /right, box=0

  IF keyword_set(inset) THEN BEGIN 

      IF !d.name EQ 'PS' THEN BEGIN 

          xsize = 0.4
          ysize = xsize/!gratio/1.7
          
          xmin = 0.5
          ymin = 0.27
          xmax = xmin + xsize
          ymax = ymin + ysize


      ENDIF ELSE BEGIN 
          
          
          xsize = 0.4
          ysize = xsize/!gratio
          
          xmin = 0.5
          ymin = 0.4
          xmax = xmin + xsize
          ymax = ymin + ysize

      ENDELSE 

      pos = [xmin, ymin, xmax, ymax]

	  ytitle2=estitle('deltaytitle')
	  xtitle2=estitle('mpcxtitle')
      ploterror, sh.meanr/1000.0, sh.sigma, sh.sigmaerr,$ 
        ytitle=ytitle2, xtitle=xtitle2, $
        xtickf='loglabels', $
        psym=8, hat=0, symsize=0.7, xticklen=0.04, $
        /xlog, yrange=restricted_yrange, _extra=_extra, $
        position = pos, /noerase, charsize=0.9

      oplot, [.005, 300], [0,0]
  ENDIF 


  if keyword_set(dops) then endplot, /trim_bbox

END ;; ngals_plot_linear


;; See how the random pair density changes with ngals bin
;; This is primarily due to the *redshift cut* on the steep
;; number counts funcion, not the difference in mean lens redshift, 
;  which is relatively small.
;; send /sigma to plot sigma instead

PRO maxbcg_lensing::ngals_rand_plot, $
  ngals_nbin, sigma=sigma, dops=dops, dopng=dopng


  ;; build up output file names
  IF keyword_set(sigma) THEN BEGIN 
      addstr = 'deltasig' 
  ENDIF ELSE BEGIN 
      addstr='npair'
  ENDELSE 

  plotdir = self->plotdir()
  sampstring = self->objshear::sample_string()
  nbinstr = strn(ngals_nbin, len=2, padchar='0')
  psfile = concat_dir(plotdir,$
                      'ngals'+nbinstr+'_rand_'+addstr+'_'+sampstring+'.eps')
  pngfile = concat_dir(plotdir,$
                       'ngals'+nbinstr+'_rand_'+addstr+'_'+sampstring+'.png')
  IF keyword_set(dops) THEN BEGIN 
      charsize = 0.7
      begplot, name=psfile, /encapsulated, xsize = 8.5, ysize=11, /color
      othick = 5
  ENDIF ELSE BEGIN 
      charsize=1
      othick = 2
  ENDELSE 

  esheldon_setup


  r = self->objshear::combined_get(/rand_combined)

  nrand = float(r.nlenses)


  simpctable, colorlist=colorlist
  ncolor = n_elements(colorlist)


  self->maxbcg::ngals_bins, ngals_nbin, lowlim, highlim


  FOR ngals_bin=0L,ngals_nbin-1 DO BEGIN 
 
      ri = self->objshear::combined_get(/rand_combined, $
                              ngals_nbin=ngals_nbin, ngals_bin=ngals_bin)

      ri_nrand = float(ri.nlenses)

      pclr = colorlist[ngals_bin+1]

      x = ri.meanr/1000
      IF keyword_set(sigma) THEN BEGIN  

		  ytitle = textoidl('\Delta\Sigma_i/\Delta\Sigma_{all}')
          y = ri.sigma
          yerr = ri.sigmaerr
          yrange = [-1,1]
;          yerr = ratio*sqrt( (ri.sigmaerr/ri.sigma)^2 + $
;                                  (r.sigmaerr/r.sigma)^2 )
          iso = 1
          xticklen=0.04
      ENDIF ELSE BEGIN 
          y = (ri.npair/ri_nrand) / (r.npair/nrand)
          yerr = y*sqrt( 1.0/ri.npair + 1.0/r.npair )
          ytitle='Npair_i/Npair_all'
          yrange = [0.8, 1.1]

          iso = 0
      ENDELSE 

      IF ngals_bin EQ 0 THEN BEGIN 
          plot, x, y, $
                xrange=[0.0005, 100], xstyle=3, $
                yrange = yrange, $
                /xlog, $
                xtitle=estitle('mpcxtitle'), $
                ytitle=ytitle, $
                xtickf='loglabels', iso=iso, $
                xticklen=xticklen
      ENDIF  
      oploterror, x, y, yerr, $
                  color=pclr, errc=pclr

      add_arrval, self->maxbcg::ngals_string(lowlim[ngals_bin], $
                                             highlim[ngals_bin]),mess
      add_arrval, pclr, pcolors

  ENDFOR 
  
  lines = intarr(ngals_nbin)
  legend, mess, box=0, colors=pcolors, line=lines, charsize=charsize

  IF keyword_set(sigma) THEN BEGIN 
      oplot, r.meanr/1000, r.sigma, thick=othick
      oplot, [2.e-2, 1.e8], [0,0]
  ENDIF ELSE BEGIN 
      oplot, [2.e-2, 1.e8], [1,1]
  ENDELSE 

  IF keyword_set(dops) THEN BEGIN 
      endplot,/trim_bbox
  ENDIF ELSE IF keyword_set(dopng) THEN BEGIN 
      print
      print,'Writing png: ',pngfile
      write_png, pngfile, tvrd(/true)
  ENDIF 


END 







PRO maxbcg_lensing::ngals_plot_corr, nbin, dops=dops

;  IF n_params() LT 1 THEN BEGIN 
;      print,'-Syntax: ml->ngals_plot_corr, nbin'
;      return
;  ENDIF 

  subtype='ngals12'
  nbin = self->subtype_nbin(subtype)
    

  plotdir = self->plot_dir()
  nbinstr = strn(nbin, len=2, padchar='0')

  IF keyword_set(dops) THEN BEGIN 
      psfile = $
        concat_dir(plotdir,$
               'ngals'+nbinstr+'_corr_'+self->objshear::sample_string()+'.ps')
      begplot, name=psfile, /landscape

      symsize = 0.7
      !p.charsize = 1.0
  ENDIF 
  
  esheldon_setup
  self->maxbcg::ngals_bins, nbin, lowlim, highlim
  mplot_value = self->mplot_value(nbin)
  psym = 8
  
  xtickf='loglabels'
  ytickf='loglabels'
  xticklen=0.04
  yticklen=0.04
  yrange = [0.001, 10]
      

  sharr = self->objshear::lensread('corr',subtype=subtype)
  
  
  xrange = [min(sharr[0].meanr), max(sharr[0].meanr)]/1000.
  xrange[0] = 0.5*xrange[0]
  xrange[1] = 1.5*xrange[1]

  erase & multiplot, mplot_value, /square
  FOR i=0L, nbin-1 DO BEGIN 
      
      shc = sharr[i]

      self->ngals_setupplot, nbin, i, xtickf, ytickf, xtitle, ytitle,/corr
      ploterror, shc.meanr/1000.0, shc.corr-1, shc.corr_err, psym=psym,$
        /xlog, /ylog, $
        xtitle=xtitle, ytitle=ytitle, $
        yrange=yrange, ystyle=1+2, $
        xrange=xrange, xstyle=1+2, $
        xtickf=xtickf,ytickf=ytickf,xticklen=xticklen,yticklen=yticklen,$
        /nohat, symsize=0.7
      
      str = self->maxbcg::ngals_legend(lowlim[i], highlim[i])
      legend,str,/right,box=0,charsize=1
      IF i NE nbin-1 THEN multiplot
      
  ENDFOR 
  multiplot, /reset


  IF keyword_set(dops) THEN endplot,/landfix

END 



PRO maxbcg_lensing::lx_plot, subtype, cumulative=cumulative, dops=dops

  IF n_params() LT 1 THEN BEGIN 
      message,'-syntax: mb->lx_plot, subtype, /cumulative, /dops'
  ENDIF 

  IF keyword_set(dops) THEN BEGIN 
      plotdir = self->plot_dir(subtype=subtype, /createdir)

      file = 'maxbcg_'+self->sample_string()+'_'+subtype

      IF keyword_set(cumulative) THEN BEGIN 
          file = file + '_cumulative'
      ENDIF 
      file = file + '.eps'

      file = concat_dir(plotdir, file)

      begplot, file, /color, /encapsulated
  ENDIF 

  esheldon_setup

  sh =self->lensread('corr')
  lx = self->lensread('corr', subtype=subtype)

  clr0 = c2i('darkgreen')
  clr1 = c2i('red')

 
  ytitle = estitle('deltaytitle')
  xtitle = estitle('mpcxtitle')

  IF keyword_set(cumulative) THEN BEGIN 
      shtag = where(tag_names(sh) EQ 'TSIGMA')
      shetag = where(tag_names(sh) EQ 'TSIGMAERR')
      lxtag = where(tag_names(lx) EQ 'TSIGMA')
      lxetag = where(tag_names(lx) EQ 'TSIGMAERR')
      
      addstr = 'Cumulative '
  ENDIF ELSE BEGIN 
      shtag = where(tag_names(sh) EQ 'SIGMA')
      shetag = where(tag_names(sh) EQ 'SIGMAERR')
      lxtag = where(tag_names(lx) EQ 'SIGMA')
      lxetag = where(tag_names(lx) EQ 'SIGMAERR')
      addstr = ''
  ENDELSE 
  ytitle = addstr + ytitle
  erase & multiplot, [1,2]

  ploterror, sh.meanr/1000, sh.(shtag), sh.(shetag), $
    /xlog, /ylog,yrange=[1,1.e3],ystyle=3,xrange=[0.02,20],xstyle=3,$
    ytitle=ytitle, ytickf='loglabels'
  oploterror, lx[1].meanr/1000, lx[1].(lxtag), lx[1].(lxetag), $
    color=clr1
  oploterror, lx[0].meanr/1000, lx[0].(lxtag), lx[0].(lxetag), $
    color=clr0

  units = '4.58'+!csym.times+'10!U43!N ergs/cm!U2!N/s'
  legend, $
    ['all', 'L!DX!N < '+units, 'L!DX!N > '+units], $
    line=0, color=[!p.color, clr0, clr1], /right, box=0, charsize=1

  multiplot

  ratio0 = lx[0].tsigma/sh.tsigma
  ratio0err = ratio0*sqrt( (sh.(shetag)/sh.(shtag))^2 + $
                           (lx[0].(lxetag)/lx[0].(lxtag))^2 )

  ratio1 = lx[1].tsigma/sh.tsigma
  ratio1err = ratio1*sqrt( (sh.(shetag)/sh.(shtag))^2 + $
                           (lx[1].(lxetag)/lx[1].(lxtag))^2 )


  plot, [0], /nodata, $
    /xlog, yrange=[0,2],ystyle=3,xrange=[0.02,20],xstyle=3, $
    xtitle=xtitle, ytitle='L!DX!N!Ubin!N / all', $
    xtickf='loglabels'

  oplot, [0.001, 10000], [1,1]
  oploterror, sh.meanr/1000, ratio0, ratio0err, $
    color=clr0
  oploterror, sh.meanr/1000, ratio1, ratio1err, $
    color=clr1
    


  multiplot, /reset

  IF keyword_set(dops) THEN endplot, /trim_bbox


END 











;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plot all the separate randoms
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO maxbcg_lensing::plotrand_dispersion, ngals_nbin, ngals_bin, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      plotdir=self->plotdir()
      rmstr = self->objshear::sample_string()
      file = concat_dir(plotdir,'rand_dispersion_'+rmstr+'.ps')
      begplot,name=file,/color,/landscape
  ENDIF 

  esheldon_setup

  randnum = self->randnum()
  sh = self->objshear::combined_get(randnum=randnum)
  lsh = self->objshear::combined_get(ngals_nbin=ngals_nbin, $
                           ngals_bin=ngals_bin, $
                           /corr)


  self->maxbcg::ngals_bins, ngals_nbin, lowlim, highlim
  leg = self->maxbcg::ngals_legend(lowlim[ngals_bin], highlim[ngals_bin])

  nrand = n_elements(sh)

  shave = combine_shstructs(sh)

  sh.meanr = sh.meanr/1000.0
  lsh.meanr = lsh.meanr/1000.0
  shave.meanr = shave.meanr/1000.0

  xtickf='loglabels'
  xticklen = 0.04
  randline = 0
  lensline = 3


  simpctable, colorlist=colors
  w = where(colors NE c2i('red'))
  colors = colors[w]

  !p.multi=[0,0,2]
  
  erase & multiplot, [1,2]

  cc=c2i('red')
  cthick=!p.thick*1.5

  yrange = [-40, 40]

  xtitle = estitle('mpcxtitle')
  ytitle = estitle('deltaytitle')

  plot, sh[0].meanr, sh[0].sigma, /xlog, $
    yrange=yrange, ystyle=1+2, xstyle=2, $
    ytitle=ytitle, xticklen=xticklen

  FOR i=1L, nrand-1 DO BEGIN 

      oplot, sh[i].meanr, sh[i].sigma, color=colors[i]

  ENDFOR 



  oploterror, lsh.meanr, lsh.sigma, lsh.sigmaerr, psym=4
  oplot, lsh.meanr, lsh.sigma, line=lensline

  oplot,  [0.001, 100000], [0,0]
  oploterror,shave.meanr,shave.sigma,shave.sigmaerr,$
    color=cc,errc=cc,psym=8, thick=cthick
  oplot,shave.meanr,shave.sigma, color=cc, thick=cthick

  legend,['rand', leg], line=[randline, lensline], /right, /bottom, box=0,$
    charsize=1,color=[cc, !p.color]

;  key = prompt_kbrd()

  multiplot

  yrange = [-3, 3]
  plot, sh[0].meanr, sh[0].sigma, /xlog, $
    yrange=yrange, ystyle=1+2, xstyle=2, $
    xtitle=xtitle, ytitle=ytitle, xticklen=xticklen, xtickf=xtickf

  FOR i=1L, nrand-1 DO BEGIN 

      oplot, sh[i].meanr, sh[i].sigma, color=colors[i]

  ENDFOR 
  oplot,  [0.001, 100000], [0,0]
  oploterror,shave.meanr,shave.sigma,shave.sigmaerr,$
    color=cc,errc=cc,psym=8, thick=cthick
  oplot,shave.meanr,shave.sigma, color=cc, thick=cthick

  oplot, lsh.meanr, lsh.sigma, line=lensline
  oploterror, lsh.meanr, lsh.sigma, lsh.sigmaerr, psym=4

  multiplot,/reset

  IF keyword_set(dops) THEN BEGIN 
      endplot,/landfix
  ENDIF 

END 







PRO maxbcg_lensing::compare_rand_regauss

  mb8 = obj_new('maxbcg_lensing',8)
  mb9 = obj_new('maxbcg_lensing',9)

  r8 = mb8->objshear::combined_get(/rand_combined)
  r9 = mb9->objshear::combined_get(/rand_combined)
  obj_destroy,mb8,mb9

  begplot,name='~/plots/MaxBCG/compare_regauss_rand.eps',/color,/encap,$
          xsize=11,ysize=8.5

  esheldon_setup

  yrange = [-0.3, 0.3]
  ytitle = estitle('randDeltaYtitle')
  xtitle = estitle('mpcxtitle')

  aploterror,!gratio,r8.meanr/1000,r8.sigma,r8.sigmaerr,$
             psym=8,/xlog,$
             yrange=[-.3,.3],xstyle=2,ystyle=3,ytit=ytitle, xtit=xtitle
  oploterror,r9.meanr/1000,r9.sigma,r9.sigmaerr,psym=4,color=c2i('blue'),errc=c2i('blue')
  oplot,[.001,30],[0,0]
  legend,['Old','Regauss'],$
         psym=[8,4],$
         color=[!p.color,c2i('blue')],/right,box=0,charsize=1

  endplot;, /trim_bbox



END 








function maxbcg_lensing::subtypes, nsubs
    subs = $
        ['ilum200_16',$
        'ilum200_ngals200_16_2',$
        'ilum200_z_16_2',$
        'ngals200_12',$
        'ngals200_ilum200_12_2']
    nsubs = n_elements(subs)
    return, subs
end
pro maxbcg_lensing::jackknife_run

    subs = self->subtypes(nsubs)
    for i=0L, nsubs-1 do begin
        self->objshear::jackknife, subtype=subs[i]
    endfor

end
pro maxbcg_lensing::combine_samples_run, sample1, sample2, dtype

    subs = self->subtypes(nsubs)
    for i=0L, nsubs-1 do begin
        self->objshear::combine_samples, sample1, sample2, dtype, subtype=subs[i]
    endfor

end

pro objshear::lensfile_convert_run, dtype, sample=sample, overwrite=overwrite

    subs = self->subtypes(nsubs)
    for i=0L, nsubs-1 do begin
        self->objshear::lensfile_convert, dtype, subtype=subs[i], sample=sample, overwrite=overwrite
    endfor

end


PRO maxbcg_lensing::compare_jack

  ngals_nbin = 12
  FOR ngals_bin = 0L, ngals_Nbin-1 DO BEGIN 
      t = self->objshear::combined_get(ngals_nbin=ngals_nbin, ngals_bin=ngals_bin, /corr)
      tj = self->objshear::jackknife_get(ngals_nbin=ngals_nbin, ngals_bin=ngals_bin)
      
      plot_density_contrast, t, /log
      oploterror, tj.meanr, tj.sigma, tj.sigmaerr, psym=4, color=c2i('green'), errc=c2i('green')
      
      colprint,t.sigma,tj.sigma,t.sigma/tj.sigma

      key = prompt_kbrd()
  ENDFOR 

END 







PRO maxbcg_lensing::ngals_plot_covariance, nbin, dops=dops

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: ml->ngals_plot_covariance, nbin'
      return
  ENDIF 
  

  nbinstr = strn(nbin, len=2, padchar='0')
  IF keyword_set(dops) THEN BEGIN 

      plotdir=self->plotdir()
      psfile = $
        concat_dir(plotdir,$
                   'ngals'+nbinstr+'_covariance_'+$
                   self->objshear::sample_string()+'.ps')

      begplot, name=psfile, /landscape
;      !p.charsize = 1
      !p.thick = 2

      print,!x.margin
      print,!y.margin

      !x.margin = [5, 1]
      !y.margin = [2,1]

  ENDIF ELSE BEGIN 
      !p.charsize = 1.5
  ENDELSE 

  esheldon_setup

  mplot_value = self->mplot_value(nbin)
  self->maxbcg::ngals_bins, nbin, lowlim, highlim


  symsize = 0.7

;  erase & multiplot,mplot_value, /square
  !p.multi=0
  !p.multi[1:2] = mplot_value


  FOR i=0L, nbin-1 DO BEGIN 

      shc = self->objshear::jackknife_get(ngals_nbin=nbin, ngals_bin=i)

      self->ngals_setupplot, nbin, i, xtickf, ytickf, xtitle, ytitle

      str = self->maxbcg::ngals_legend(lowlim[i], highlim[i])
      tvim2, shc.correlation, title=str, /scale

;      ploterror, shc.meanr/1000.0, shc.sigma, shc.sigmaerr, $
;        psym=psym, symsize=symsize, $
;        /xlog, /ylog, $
;        xtitle=xtitle, ytitle=ytitle, $
;        xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, $
;        xtickf=xtickf, ytickf=ytickf, $
;        xticklen=ticklen, yticklen=ticklen



;      legend,str,/right,box=0,charsize=1
;      IF i NE nbin-1 THEN multiplot


  ENDFOR 

;  multiplot,/reset
  !p.multi=0

  IF keyword_set(dops) THEN BEGIN 
      endplot,/landfix

      !x.margin = [10,4]
      !y.margin = [4,2]

  ENDIF 

END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Dave's mass profiles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION maxbcg_lensing::mass_dir,$
  ngals_nbin=ngals_nbin, tlum_nbin=tlum_nbin
  dir = $
    self->lensoutput_dir(ngals_nbin=ngals_nbin,tlum_nbin=tlum_nbin)+'masses/'
  return,dir
  
END 

FUNCTION maxbcg_lensing::mass_file,$
  ngals_nbin=ngals_nbin, ngals_bin=ngals_bin, $
  tlum_nbin=tlum_nbin, tlum_bin=tlum_bin, $
  zbin_nbin=zbin_nbin, zbin_bin=zbin_bin

  dir = self->mass_dir(ngals_nbin=ngals_nbin, tlum_nbin=tlum_nbin)

  tfile = self->combined_file(ngals_nbin=ngals_nbin, ngals_bin=ngals_bin,$
                              tlum_nbin=tlum_nbin, tlum_bin=tlum_bin, $
                              zbin_nbin=zbin_nbin, zbin_bin=zbin_bin)
  
  dirsep, tfile, tdir, file
  file = repstr(file, 'MaxBCG_combined', 'MaxBCG_jackknife')
  file = repstr(file, '.st', '_xm.fit')

  return,dir + file

END 

FUNCTION maxbcg_lensing::add_rhocum, struct

  nrad = n_elements(struct[0].mass)
  arrval = dblarr(nrad)
  newstruct = create_struct(struct[0], $
                            'massrad', arrval, $
                            'drhocum',arrval, $
                            'drhocum_err', arrval, $
                            'rhocrit', 0.0)
  nstruct = n_elements(struct)
  newstruct = replicate(newstruct, nstruct)
  copy_struct, struct, newstruct


  rmpc = struct.meanr[0:nrad-1]/1000.0
  volume = 4.0/3.0*!dpi*rmpc^3

  ;; mass in units of 10^2/Mpc^3
  drhocum = struct.mass/volume
  drhocum_err = struct.mass_err/volume
  ;; Units solar masses/Mpc^3
  rhocrit = !rhocrit * (1.0 + struct.mean_z)^3
  ;; Now 10^12 Msolar/Mpc^3
  rhocrit = rhocrit/1.e12

  newstruct.massrad = rmpc
  newstruct.drhocum = drhocum
  newstruct.drhocum_err = drhocum_err
  newstruct.rhocrit = rhocrit
  return, newstruct

END 

FUNCTION maxbcg_lensing::add_means, struct, jstruct
  newstruct = create_struct(struct[0], $
                            'mean_z', 0.0, $
                            'mean_tngals', 0.0, $
                            'mean_tlum', 0.0)
  newstruct = replicate(newstruct, n_elements(struct))
  copy_struct, struct, newstruct
  newstruct.mean_z = jstruct.mean_z
  newstruct.mean_tngals = jstruct.mean_tngals
  newstruct.mean_tlum = jstruct.mean_tlum
  return,newstruct
END 

FUNCTION maxbcg_lensing::mass_get, $
  ngals_nbin=ngals_nbin, ngals_bin=ngals_bin, $
  tlum_nbin=tlum_nbin, tlum_bin=tlum_bin, $
  zbin_nbin=zbin_nbin, zbin_bin=zbin_bin
  combine=combine

  files = self->mass_file(ngals_nbin=ngals_nbin, ngals_bin=ngals_bin, $
                          tlum_nbin=tlum_nbin, tlum_bin=tlum_bin)


  struct = mrdfits_multi(files)
  IF NOT tag_exist(struct[0], 'mean_z') THEN BEGIN  
      jstruct = self->jackknife_get(ngals_nbin=ngals_nbin, $
                                    ngals_bin=ngals_bin, $
                                    tlum_nbin=tlum_nbin, $
                                    tlum_bin=tlum_bin)
      
      struct = self->add_means(struct, jstruct)
  ENDIF

  struct = self->add_rhocum(struct)

  IF keyword_set(combine) AND n_elements(struct) GT 1 THEN BEGIN 
      print,'Combining'
      struct = combine_shstructs(struct)
  ENDIF 

  return,struct

END 


FUNCTION maxbcg_lensing::m200_file, $
  ngals_nbin=ngals_nbin, tlum_nbin=tlum_nbin

  dir = self->mass_dir(ngals_nbin=ngals_nbin, tlum_nbin=tlum_nbin)

  front = ''
  IF n_elements(ngals_nbin) NE 0 THEN BEGIN 
      front = 'ngals'+strn(ngals_nbin, len=2, padchar='0')+'_'
  ENDIF ELSE IF n_elements(tlum_nbin) NE 0 THEN BEGIN 
      front = 'tlum'+strn(tlum_nbin, len=2, padchar='0')+'_'
  ENDIF 
  file = dir + front + 'all_masses.st'
  return,file
END 
FUNCTION maxbcg_lensing::m200_get, $
  ngals_nbin=ngals_nbin, tlum_nbin=tlum_nbin

  file = self->m200_file(  ngals_nbin=ngals_nbin, tlum_nbin=tlum_nbin)
  print,'Reading file: ',file
  t=read_idlstruct(file)
  return,t
END 



PRO maxbcg_lensing::run_m200, ngals_struct, tlum_struct, dops=dops

  nbin = 12

  plotdir = self->plotdir()
  psfile = $
    concat_dir(plotdir,'tngals'+strn(nbin,len=2,padchar='0')+$
               '_fit_m200_'+self->objshear::sample_string()+'.ps')
  IF keyword_set(dops) THEN BEGIN 
      begplot,name=psfile,/color, /encapsulated
  ENDIF 

  simpctable, colorlist = colorlist
  ncolor = n_elements(colorlist)


  xtitle = estitle('mpcxtitle3d')
  ytitle = textoidl('\rho(<r)/\rho_{crit}')

  

  print
  print,'Doing ngals bins'
  self->maxbcg::ngals_bins, nbin, lowlim, highlim
  FOR i=0L, nbin-1 DO BEGIN 
      t = self->mass_get(ngals_nbin=nbin, ngals_bin=i)
      t = self->m200(t)

      rhomult = t.drhocum/t.rhocrit
      rhomult_err = t.drhocum_err/t.rhocrit

      pclr = colorlist[i MOD ncolor]
      add_arrval, self->maxbcg::ngals_string(lowlim[i], highlim[i]),mess
      add_arrval, pclr, pcolors
      add_arrval, 0, lines

      IF i EQ 0 THEN BEGIN 
          ngals_struct = replicate(t, nbin)
          yrange = [1, 1.e5]
          xrange = [min(t.massrad), max(t.massrad)]
          xrange[0] = xrange[0] * 0.8
          xrange[1] = xrange[1] * 1.1
          
          aploterror, 1, t.massrad, rhomult, rhomult_err, psym=8, $
            yrange=yrange, ystyle=3, $
            xrange=xrange, xstyle=3, /xlog, /ylog, $
            xtitle=xtitle, ytitle=ytitle, color=pclr, errc=pclr, $
            xtickf='loglabels', ytickf='loglabels', /nohat

          oplot, [0.001, 1.e5], [200,200]
          
          oploterror, t.massrad, rhomult, rhomult_err, $
            color=pclr, errc=pclr, /nohat
      ENDIF  ELSE BEGIN 
          oploterror, t.massrad, rhomult, rhomult_err, psym=8, $
            color=pclr, errc=pclr, /nohat
          oploterror, t.massrad, rhomult, rhomult_err, $
            color=pclr, errc=pclr,$
            /nohat
      ENDELSE 

      ngals_struct[i] = t
      key = prompt_kbrd('hit a key')
      IF key EQ 'q' THEN return
  ENDFOR 

  legend, mess, box=0, colors=pcolors, lines=lines, charsize=1, /right

  IF keyword_set(dops) THEN endplot, /trim_bbox



  file = self->m200_file(ngals_nbin=nbin)
  print
  print,'Writing to file: ',file
  write_idlstruct, ngals_struct, file

return
  psfile = $
    concat_dir(plotdir,'tlum'+strn(nbin,len=2,padchar='0')+$
               '_fit_m200_'+self->objshear::sample_string()+'.ps')

  IF keyword_set(dops) THEN BEGIN 
      begplot,name=psfile,/color
  ENDIF 


  print
  print,'Doing tlum bins'
  FOR i=0L, nbin-1 DO BEGIN 
      t = self->m200(tlum_nbin=nbin, tlum_bin=i, /doplot)
      
      IF i EQ 0 THEN BEGIN 
          tlum_struct = replicate(t, nbin)
      ENDIF  

      tlum_struct[i] = t
      key = prompt_kbrd('hit a key')
      IF key EQ 'q' THEN return

  ENDFOR 

  file = self->m200_file(tlum_nbin=nbin)
  print
  print,'Writing to file: ',file
  write_idlstruct, tlum_struct, file

  IF keyword_set(dops) THEN endplot
;  ploterror, tlum_struct.mean_tlum, tlum_struct.m200, tlum_struct.m200err, $
;             psym=8, /xlog, /ylog


END 


function maxbcg_lensing::massvals, struct
    if tag_exist(struct, 'massout', ind=mind) then begin
        if tag_exist(struct, 'massout_err', ind=eind) then begin
            if tag_exist(struct, 'ir', ind=rind) then begin
                mass=struct.massout
                mass_err = struct.massout_err
                radius = struct.ir
                return,{r:radius, mass:mass, mass_err:mass_err}
            endif
        endif
    endif
    if tag_exist(struct, 'mass_out', ind=mind) then begin
        if tag_exist(struct, 'mass_out_err', ind=eind) then begin
            if tag_exist(struct, 'meanr') then begin
                mass = struct.mass_out
                mass_err = struct.mass_out_err
                radius = struct.meanr[0:n_elements(mass)-1]/1000
                return,{r:radius, mass:mass, mass_err:mass_err}
            endif
        endif
    endif
    on_error, 2
    message,'mass and mass err tags not found'

end
FUNCTION maxbcg_lensing::m200, struct, doplot=doplot

    mst = self->massvals(struct)
    mfac = 1.e12
    mass = mst.mass*mfac
    mass_err = mst.mass_err*mfac
    rmpc = mst.r
    ;rhoc = struct.rhocrit;*mfac
    rhoc = rhocrit(z=0.25)

    volume = 4.0/3.0*!dpi*rmpc^3
    ;; Msun/Mpc^3
    drhocum = mass/volume
    drhocum_err = mass_err/volume
 
    rhomult = drhocum/rhoc
    rhomult_err = drhocum_err/rhoc
    nrad = n_elements(rmpc)

    if 0 then begin
        i = 0
        WHILE (rhomult[i] GT 200) AND (i LT nrad) DO BEGIN 
            i = i+1
        ENDWHILE 

        IF i EQ i+1 THEN BEGIN 
            print,'r200 is greater than max radius'
            r200 = -9999.0
            m200 = -9999.0
            m200err = -9999.0
        ENDIF ELSE BEGIN 
            ;; simply take midpoint as r200
            r200 = (rmpc[i] + rmpc[i-1])/2.0
            M200 = interpol(mass, rmpc, r200)
            M200err = interpol(mass_err, rmpc, r200)
 
       ENDELSE 
    endif else begin
        w=where(rhomult gt 0)
        r200 = interpol(rmpc[w], rhomult[w], 200.0)
        M200 = interpol(mass[w], rmpc[w], r200)
        M200err = interpol(mass_err[w], rmpc[w], r200)
    endelse
    
    ;print,struct.r200,r200,struct.m200,m200

    if keyword_set(doplot) then begin

        !p.multi=[0,0,2]

        yrange = prange(rhomult, rhomult_err, /slack)
        if yrange[0] lt 0 then yrange[0] = 1.e-2
        pplot, rmpc, rhomult, yerr=rhomult_err, yrange=yrange, ystyle=3, $
            psym=8, /xlog, /ylog
        oplot, [1.e-3, 1.e16], [200,200], color=c2i('red')
        legend,'r200 = '+ntostr(r200,f='(F5.2)'), /right

        yrange = prange(mass, mass_err, /slack)
        if yrange[0] lt 0 then yrange[0] = 1.e10
        pplot, rmpc, mass, yerr=mass_err, psym=8, yrange=yrange, ystyle=3, $
            /xlog, /ylog
        pplot, [r200, r200], [1.e-3,1.e16], line=2, color=c2i('red'),/over
        pplot, [0.001, 1.e6], [m200, m200], line=2, color=c2i('red'),/over
        legend,'m200 = '+ntostr(m200/1.e12,f='(F7.2)')

        !p.multi=0
    endif

    sh = create_struct('r200', r200, $
                        'm200', m200, $
                        'm200err', m200err)


    return,sh

END 

PRO maxbcg_lensing::m200_plot, type, dops=dops

    nbin = 12
    ngals_struct = self->m200_get(ngals_nbin=nbin)
    tlum_struct = self->m200_get(tlum_nbin=nbin)

    self->titles, 'mass', xtitle, ytitle

    plotdir = self->plotdir()
    psfile = $
        concat_dir(plotdir,$
        type+strn(nbin,len=2,padchar='0')+'_m200_'+$
        self->objshear::sample_string()+'.eps')
    IF keyword_set(dops) THEN BEGIN 
        begplot,name=psfile,/encap,/color
        oclr = c2i('blue')
    ENDIF ELSE BEGIN 
        oclr = c2i('green')
    ENDELSE 

    esheldon_setup
    IF type EQ 'tngals'THEN BEGIN 

        xrange = [1,200]
        yrange = [2.e11, 5.e14]
        xtitle = 'N!S!UT!N!R!Dgals!N'
        ytitle = 'M!D200!N [h!U-1!N M'+sunsymbol()+']'

        m200 = ngals_struct.m200*1.e14
        m200err = ngals_struct.m200err*1.e14
        ploterror, $
            ngals_struct.mean_tngals, m200, m200err, $
            psym=8, $
            xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
            /xlog, /ylog, /iso, $
            xtickf='loglabels', ytickf='loglabels', $
            xtitle=xtitle, ytitle=ytitle

        ;      oploterror, $
            ;        tlum_struct.mean_tngals, tlum_struct.m200, tlum_struct.m200err, $
            ;        psym=4, color=c2i('green'), errc=c2i('green')

        guess = [1.e10, 2.0]
        fitpower, $
            ngals_struct.mean_tngals, m200, m200err, $
            guess, yfit, aout, aout_err


        minn = min(ngals_struct.mean_tngals, max=maxn)
        xp = arrscl(findgen(100), minn, maxn)
        yp = aout[0]*xp^aout[1]
        oplot, xp, yp, color=oclr


        aout[0] = aout[0]/1.e10
        aout_err[0] = aout_err[0]/1.e10

        legend,['M = A Ngals!U'+!csym.alpha+'!N',$
            'A = '+ntostr(aout[0],4,/round) + !csym.plusminus+$
            ntostr(aout_err[0],4,/round)+' 10!U10!N M'+sunsymbol(), $
            !csym.alpha +' = '+ntostr(aout[1],4,/round)+$
            !csym.plusminus+ntostr(aout_err[1],4,/round) ], $
            /left, box=0, spacing=2.5

    ENDIF ELSE IF type EQ 'tlum'THEN BEGIN 

        xrange = [3.e10, 1.e13]
        yrange = [1.e11, 2.e14]
        xtitle = 'i-band Lum!UT!N [L'+sunsymbol()+']'
        ytitle = 'M!D200!N [h!U-1!N M'+sunsymbol()+']'

        m200 = tlum_struct.m200*1.e14
        m200err = tlum_struct.m200err*1.e14
        ploterror, $
            tlum_struct.mean_tlum, m200, m200err, $
            psym=8, /xlog, /ylog, /iso, $
            xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
            xtickf='loglabels', ytickf='loglabels', $
            xtitle=xtitle, ytitle=ytitle


        guess = [4., 1.5]
        fitpower, $
            tlum_struct.mean_tlum/1.e13, tlum_struct.m200, tlum_struct.m200err, $
            guess, yfit, aout, aout_err

        minn = min(tlum_struct.mean_tlum, max=maxn)
        xp = arrscl(findgen(100), minn, maxn)
        yp = 1.e14*aout[0]*(xp/1.e13)^aout[1]
        oplot, xp, yp, color=oclr


        ;      aout[0] = aout[1]/1.e14
        ;      aout_err[0] = aout_err[1]/1.e14
        legend,['M = A (L/10!U10!N L'+sunsymbol()+')!U'+!csym.alpha+'!N',$
            'A = '+ntostr(aout[0],4,/round) + !csym.plusminus+$
            ntostr(aout_err[0],4,/round)+' 10!U14!N M'+sunsymbol(), $
            !csym.alpha +' = '+ntostr(aout[1],4,/round)+$
            !csym.plusminus+ntostr(aout_err[1],4,/round) ], $
            /left, box=0, spacing=2.5




    ENDIF 

    IF keyword_set(dops) THEN endplot,/trim_bbox

END 




PRO maxbcg_lensing::compare_mass_sis, struct, dops=dops
  
    dir = '~/tmp/'
    files = findfile(dir+'maxbcg_ilum200_16/*')

    ws = self->where_string('ilum200_16', labels=labels)

    print,files

    t=mrdfits_multi(files)

    nt = n_elements(t)

    nrad = n_elements(t[0].meanr)

    xtitle = estitle('mpcxtitle')
    rhocrit = !rhocrit * (1.0 + 0.25)^3 ; Msolar/Mpc^3
    rhocrit = rhocrit/1.e14       ; 10^14 Msolar/Mpc^3


    FOR i=0L, nt-1 DO BEGIN 

        !p.multi = [0,2,2]
        tsigma_sis_fit, t[i], sismass, sismass_err, sissigma2, sissigma2_err
        sissigma = sqrt(sissigma2)
        sissigma_err = sqrt(sissigma2_err)/2.0

        ;; correct for projection
        sismass = sismass[0:nrad-2]*2d/!pi/1.e14
        sismass_err = sismass_err[0:nrad-2]*2d/!pi/1.e14

        mass = t[i].mass_out/100.0
        mass_err = t[i].mass_out_err/100.0

        r = t[i].meanr[0:nrad-2]/1000

        yrange = [0.001,1.e2]
        xrange = [0.01, 30]
        pplot, r, mass, yerr=mass_err, hat=0, psym=8, $
            /xlog, /ylog, yrange=yrange, ystyle=3, xrange=xrange, xstyle=3, $
            aspect=1,/center, $
            xtitle=xtitle, ytitle='M(<r) [10!U14!N M'+sunsymbol()+']', $
            ytickf='loglabels', xtickf='loglabels'

        pplot, r, sismass, yerr=sismass_err, /overplot, psym=4, color=c2i('darkgreen')

        legend, ['Inversion','SIS'],psym=[8,4],color=[!p.color, c2i('darkgreen')], box=0, charsize=1
        legend, labels[i], /right, box=0,charsize=1


        ;; overplot a line corresponding to 200*rhocrit*4/3*!pi*r^3

        volume = 4.0/3.0*!pi*r^3
        rho200 = 200*rhocrit

        M200crit = rho200*volume
        oplot, r, M200crit, color=c2i('red')

        ;; get r200,M200
        mean_density = mass/volume
        mean_density_err = mass_err/volume

        w = max( where(mean_density GT rho200) )

        slope = (mean_density[w+1] - mean_density[w])/(r[w+1]-r[w])
        offset = mean_density[w]-slope*r[w]

        tr200 = (rho200 - offset)/slope

        tM200 = interpol(mass, r, tr200)
        tM200err = interpol(mass_err, r, tr200)

        ;; Get estimated sis sigmav 
        tsigmav = interpol(sissigma[0:nrad-2], r, tr200)
        tsigmav_err = interpol(sissigma_err[0:nrad-2], r, tr200)


        add_arrval, tr200, r200
        add_arrval, tM200, M200
        add_arrval, tM200err, M200err

        add_arrval, tsigmav, sigmav
        add_arrval, tsigmav_err, sigmav_err

        oplot, [tr200], [tM200], psym=7,color=c2i('red'), symsize=2



        pplot, r, mean_density, yerr=mean_density_err, psym=8, /xlog, /ylog, $
            aspect=1, /center, yrange=[1.e-4,1.e4],ystyle=3,xrange=xrange, xstyle=3
        oplot, [1.e-3,1.e3], [rho200,rho200], color=c2i('red')
        oplot, [tr200,tr200],[1.e-5,1.e5],color=c2i('red')

        ratio = mass/sismass
        ratio_err = ratio*abs(mass_err/mass)
        pplot, r, ratio, yrange=[0,2], yerr=ratio_err, /xlog, $
            aspect=!gratio,/center, xtitle=xtitle,ytitle='M/M(SIS)'

        oplot, [1.e-5,1.e5],[1,1]


        plot_density_contrast, t[i], /log, /mpc

        ;; overplot the sigmav model at r200
        print,'sigmav (r200) = ',tsigmav
        w = where(t[i].meanr/1000 LE tr200)
        smod = 3.36e3*(tsigmav/170.)^2/t[i].meanr[w]

        oplot, r[w], smod, color=c2i('red')
        fitsis, $
            t[i].meanr[w], t[i].sigma[w], t[i].sigmaerr[w], 700.0, $
            yfit, tfit_sigmav, tfit_sigmav_err

        oplot, t[i].meanr[w]/1000, yfit, color=c2i('green')


        add_arrval, tfit_sigmav, fit_sigmav
        add_arrval, tfit_sigmav_err, fit_sigmav_err

        ;      key = prompt_kbrd('hit a key')
        ;      IF key EQ 'q' THEN return
    ENDFOR 


    !p.multi=0
    pplot, $
        t.mean_ilum200*1.e10, M200*1.e14, yerr=M200err*1.e14, $
        psym=8, /xlog, /ylog, $
        aspect=1,/center, $
        xtitle='L!D200!N [ L'+sunsymbol()+' ]', $
        ytitle='M!D200!N [ M'+sunsymbol()+' ]'

    format = '(8A16)'
    print,$
        'L200','r200','M200','M200err','sigmav','sigmav_err','fit_sigmav','fit_sigmav_err',$
        format=format
    colprint,$
        double(t.mean_ilum200),$
        r200, M200, M200err, $
        sigmav, sigmav_err, $
        double(fit_sigmav),double(fit_sigmav_err)

    key = get_kbrd(1)

    pplot, t.mean_ilum200*1.e10, sigmav, yerr=sigmav_err, /xlog, /ylog, $
        yrange = [50,3000], ystyle=3

    pplot, t.mean_ilum200*1.e10, fit_sigmav, yerr=fit_sigmav_err, color=c2i('green'),psym=8, $
        /overplot

    struct = $
        { $
        iLum:double(t.mean_ilum200*1.e10), $
        r200:r200,$
        m200:m200,$
        m200_err:m200err, $
        sigmav:sigmav, $
        sigmav_err:sigmav_err, $
        fit_sigmav:double( fit_sigmav ), $
        fit_sigmav_err:double( fit_sigmav_err ) $
        }


END 



PRO maxbcg_lensing::powerlaw_fits, sample=sample, color=color

  subtype='ngals200_12'
  pdir = self->plot_dir(sub=subtype, sample=sample, /createdir)
  pfile = subtype+'_powerlaw_fits.ps'

  pfile = concat_dir(pdir, pfile)
  begplot,pfile,color=color


  t=self->lensread('jackknife', sub=subtype, sample=sample)
  nbin=n_elements(t)
  nrad = n_elements(t[0].meanr)

  IF keyword_set(color) THEN BEGIN 
      simpctable, colorlist=colors
  ENDIF ELSE BEGIN 
      colors = replicate(!p.color, 100)
  ENDELSE 

  FOR i=0L, nbin-1 DO BEGIN 


      r = t[i].meanr/1000
      sig = t[i].sigma
      sigerr = t[i].sigmaerr
      cov = t[i].covariance


      fitpower, r, sig, sigerr, [10d,-1d], yfit, afit, asig

      powrange = afit[1] + [-4.0, 4.0]*asig[1]
      normrange = afit[0] + [-4.0, 4.0]*asig[0]

      if i eq nbin-1 then begin
          cov = t[i].sigmaerr
      endif

      pow_chisq_conf_gen, r, sig, cov, powrange, normrange, 100, 100, $
        chisq_surf, bp, bn, $
        perrlow=perrlow, perrhigh=perrhigh, $
        nerrlow=nerrlow, nerrhigh=nerrhigh, $
        minchisq=minch, degfree=degf, /noprompt

      print,'bestn = ',ntostr(bn)+' '+!plusminus+' '+ntostr(nerrhigh[0])
      print,'bestp = ',ntostr(bp)+' '+!plusminus+' '+ntostr(perrhigh[0])
      print,'minchisq = ',minch
      print,'degfree = ',degf

      add_arrval, bn, bestn
      add_arrval, nerrhigh[0], bestnerr
      add_arrval, bp, bestp
      add_arrval, perrhigh[0], bestperr

      add_arrval, minch, minchisq
      add_arrval, degf, degfree

;      key=prompt_kbrd('hit a key')
;      IF key EQ 'q' THEN return

  ENDFOR 

  endplot

  IF keyword_set(color) THEN BEGIN 
      pfile = subtype+'_powerlaw_fits_oplot_color.eps'
  ENDIF ELSE BEGIN 
      pfile = subtype+'_powerlaw_fits_oplot.eps'
  ENDELSE 
  pfile = concat_dir(pdir, pfile)
  begplot, pfile,/color,/encap

  lines = [0,1,2,3,4,5]
  thick = [1,1,1,1,1,1,5,5,5,5,5,5]

  xrange = [0.01, 1.e4]
  yrange = [0.1,1.e9]

  offsets = 10.0^(0.7*lindgen(nbin))
  offsets = 10.0^(0.7*lindgen(nbin+1))

  ws = self->where_string(subtype, labels=labels)

  self->ngals200_bins, 12, lowlim, highlim

  FOR i=0L, nbin-1 DO BEGIN 

      r = t[i].meanr/1000
      sig = t[i].sigma
      sigerr = t[i].sigmaerr

      IF i GT 9 THEN BEGIN 
          w=where(sig/sigerr LE 1, nw)
          IF nw NE 0 THEN sigerr[w] = 0.01*sig[w]
      ENDIF 
;      w = where(sig/sigerr GT 1, nw)
      
      yfit = bestn[i]*r^bestp[i]

      sig = sig/yfit*offsets[i]
      sigerr = sigerr/yfit*offsets[i]

      ;; for plotting purposes
      


      IF i EQ 0 THEN overplot=0 ELSE overplot=1
      pplot, r, sig, yerr=sigerr, /nohat, $
        color=colors[i], overplot=overplot, $
        /xlog, /ylog, $
        xrange=xrange, xstyle=3, $
        yrange=yrange, ystyle=3, $
        xtitle=estitle('mpcxtitle'), $
        xtickf='loglabels', ytickf='loglabels',$
		ytitle=textoidl('\Delta\Sigma/power law [Arbitrary Units]'); $
        ;line=lines[i MOD 6], thick=thick[i], errthick=thick[i], $

      pplot, r, sig, psym=8, symsize=0.5, /overplot, color=colors[i]

      oplot, r, replicate(offsets[i],nrad), color=colors[i]
      mstr=ntostr(minchisq[i], 4, /round)
      dfstr = ntostr(degfree[i])
      chpstr = ntostr(minchisq[i]/degfree[i],4,/round)

      ;mess = labels[i] + '   '+!csym.chi+'!U2!N/'+!csym.nu+' = '+chpstr
      mess = labels[i] + '   '+chpstr
      xyouts, 60.0, offsets[i], mess, $
        charsize=1


      ;; Print the table
      ;; ngals bin, norm, power, chi2/deg
      IF lowlim[i] EQ highlim[i] THEN BEGIN 
          label = 'N_{200} = '+ntostr(lowlim[i])
      ENDIF ELSE BEGIN 
          label = ntostr(lowlim[i])+' \le N_{200} \le '+ntostr(highlim[i])
      ENDELSE 
      
      IF i NE nbin-1 THEN cont = '\\\\' ELSE cont = ''
      ;;         bindef            norm/err                pow/err
      ;;         chi2/deg

      IF abs(bestp[i]) GE 1 THEN bpf='%10.3g' ELSE bpf='%10.2g'
      IF bestn[i] GE 10 THEN BEGIN 
          IF bestnerr[i] LT 1 THEN BEGIN 
              bn = bestn[i]
              bnf = '%10.3g'
          ENDIF ELSE BEGIN 
              bn = round(bestn[i])
              bnf = '%10d' 
          ENDELSE 
      ENDIF ELSE BEGIN 
          bn = bestn[i]
          bnf = '%10.2g'
      ENDELSE 
      IF bestnerr[i] GE 1 THEN BEGIN 
          bnerr = round(bestnerr[i])
          bnerrf = '%10d' 
      ENDIF ELSE BEGIN 
          bnerr = bestnerr[i]
          bnerrf = '%10.1g'
      ENDELSE 


      format='(%"$%s$ & '+bnf+' $\\pm$ '+bnerrf+' & '+bpf+' $\\pm$ %10.1g & %10.3g/%d = %10.3g '+cont+'")'
      print, format=format, $
        label, $
        bn, bnerr, abs(bestp[i]), bestperr[i], $
        minchisq[i], degfree[i], minchisq[i]/degfree[i]


;      print, $
;        label, $
;        ntostr(bestn[i],4,/round), ntostr(bestnerr[i], 4, /round), $
;        ntostr(abs(bestp[i]),4,/round), ntostr(bestperr[i],4,/round), $
;        ntostr(minchisq[i],4,/round), ntostr(degfree[i]), $
;        ntostr(minchisq[i]/degfree[i],4,/round)

  ENDFOR 
    xyouts, 60.0, offsets[nbin], '     range        '+textoidl('\chi^2/\nu'), $
        charsize=1
  
  endplot, /trim_bbox

END 





pro maxbcg_lensing::redshift_regression_combined, type, subtype=subtype, sample=sample, indices=indices, dops=dops

  pdir = self->plot_dir(sub=subtype, sample=sample, /createdir)
  pfile = subtype+'_zregress_combined.ps'
  pfile = concat_dir(pdir, pfile)



    t=self->lensread(type, sub=subtype, sample=sample)

	if n_elements(indices) ne 0 then begin
		t=t[indices]
		newstr='_sub'+ntostr(indices[0],f='(i02)')+$
			'_'+ntostr(indices[n_elements(indices)-1], f='(i02)')+'.ps'
		pfile = repstr(pfile, '.ps', newstr)
	endif

  IF keyword_set(dops) THEN begin
	  begplot,pfile, /color
	  !p.thick=3
  endif



    ntot = n_elements(t)
    nbin=n_elements(t)/2
    nrad = n_elements(t[0].meanr)
   
    meansig = dblarr(nrad)
    meansig1 = dblarr(nrad)
    meansig2 = dblarr(nrad)
    meansig_err = dblarr(nrad)
    meansig1_err = dblarr(nrad)
    meansig2_err = dblarr(nrad)

    ind1 = lindgen(ntot/2)*2
    ind2 = lindgen(ntot/2)*2 + 1
    j=0
    for irad=0L, nrad-1 do begin
        wmom, t.sigma[irad], t.sigmaerr[irad], wm, ws, we
        meansig[irad] = wm
        meansig_err[irad] = we

        wmom, t[ind1].sigma[irad], t[ind1].sigmaerr[irad], wm, ws, we
        meansig1[irad] = wm
        meansig1_err[irad] = we

        wmom, t[ind2].sigma[irad], t[ind2].sigmaerr[irad], wm, ws, we
        meansig2[irad] = wm
        meansig2_err[irad] = we
    endfor

    ytitle = textoidl('\Delta\Sigma(R) [h M_\odot pc^{-2}]')
    xtitle = textoidl('R [h^{-1} Mpc]')
    xrange=[0.015, max(t[0].meanr/1000)*1.5]
    yrange=[0.3, 1.e3]
    pplot, t[0].meanr/1000, meansig, yerr=meansig_err, psym=-8, $
        yrange=yrange, ystyle=3, xrange=xrange, xstyle=3, $
        /xlog, /ylog, aspect=1, $
        xtitle=xtitle,ytitle=ytitle, $
		xtickf='loglabels', ytickf='loglabels'
    pplot, t[0].meanr/1000, meansig1, yerr=meansig1_err, psym=-4, $
        /over, color=c2i('blue')
    pplot, t[0].meanr/1000, meansig2, yerr=meansig2_err, psym=-7, $
        /over, color=c2i('darkgreen')

	legend, ['all','lowz','highz'], $
		color=[!p.color, c2i('blue'), c2i('darkgreen')], $
		psym=[-8, -4, -7], $
		/right

    key=prompt_kbrd('hit a key')

    diff = (meansig2-meansig1)
    differr = sqrt(meansig2^2 + meansig1^2)
    chi = diff/differr
    chierr = replicate(1.0, nrad)

    pplot, t[0].meanr/1000, chi, yerr=chierr, psym=8, $
        xrange=xrange, xstyle=3, yrange=[-2,2], ystyle=3, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=textoidl('\chi')

	if keyword_set(dops) then endplot

end

PRO maxbcg_lensing::redshift_regression, type, subtype=subtype, sample=sample, indices=indices, dops=dops

  ;subtype='ilum200_z_16_2'
  pdir = self->plot_dir(sub=subtype, sample=sample, /createdir)
  pfile = subtype+'_zregress.eps'
  pfile = concat_dir(pdir, pfile)

  t=self->lensread(type, sub=subtype, sample=sample)
	if n_elements(indices) ne 0 then begin
		t=t[indices]
		newstr='_sub'+ntostr(indices[0],f='(i02)')+$
			'_'+ntostr(indices[n_elements(indices)-1], f='(i02)')+'.eps'
		pfile = repstr(pfile, '.eps', newstr)
	endif



  IF keyword_set(dops) THEN begplot,pfile, /encap


  nbin=n_elements(t)/2
  nrad = n_elements(t[0].meanr)

  
  if nbin eq 12 then begin
      mp=[4,3]
  endif else begin
      mp=[4,4]
  endelse

  xtitle = estitle('mpcxtitle')
  ytitle = textoidl('\chi')
  erase & multiplot, mp, /square, mXtitle=xtitle, mYtitle=ytitle

  xrange = [0.012,50]
  yrange = [-3.5,5]


  IF !d.name EQ 'PS' THEN BEGIN 
      !p.charsize=1
      !x.thick=2
      !y.thick=2
      !p.thick=2
      symsize=0.7
  ENDIF 

  running_chi2 = 0.0
  numpoints = 0LL

  i=0
  FOR il=0L,nbin-1 DO BEGIN 

      r = t[i].meanr/1000.0

      IF ( (iL MOD 4) EQ 0) THEN BEGIN 
      ENDIF 
      
      diff = t[i+1].sigma - t[i].sigma
      differr = sqrt( t[i+1].sigmaerr^2 + t[i].sigmaerr^2 )

      chi = (diff/differr)
      chierr = replicate(1.0, nrad)

      pplot, r, chi, yerr=chierr, psym=8, symsize=symsize, $
        /xlog, xrange=xrange,xstyle=3, xticklen=0.04, $
        yrange=yrange, ystyle=3
      oplot, [1.e-5, 1.e5], [0,0]

      
      wmom, chi, chierr, wmean, wsig, werr
      chi2 = total(chi^2)
      redchi2 = chi2/(nrad-1.0)
;      oplot, [1.e-5,1.e5],[wmean,wmean],color=c2i('red'), line=2

      mess = textoidl('\chi^2 = ')
      IF chi2 GE 10 THEN cstr = ntostr(chi2, 4, /round) $
      ELSE cstr = ntostr(chi2, 3, /round)
      mess = mess + cstr
      legend, mess, /right, box=0,charsize=1

      running_chi2 = running_chi2 + chi2
	  numpoints = numpoints + n_elements(t[i].meanr)
      add_arrval, chi, allpoints
      add_arrval, chierr, allpointserr
      i=i+2
      multiplot

  ENDFOR 

  wmom, allpoints, allpointserr, wmean, wsig, werr
  print,'chi mean/err: ',wmean,werr

  print,'total chi^2: ',running_chi2
  print,'reduced chi^2: '+ntostr(running_chi2)+'/'+ntostr(numpoints)+' = '+ntostr(running_chi2/numpoints)
  print,'probability happening by chance: ',prob_chisq(numpoints, running_chi2)

  multiplot,/default


  IF keyword_set(dops) THEN endplot, /trim_bbox


END 


PRO maxbcg_lensing::plot_ngals_lumsplit, sample=sample, color=color, dops=dops

  subtype='ngals200_ilum200_12_2'
  pdir = self->plot_dir(sub=subtype, sample=sample, /createdir)

  IF keyword_set(color) THEN BEGIN 
      pfile = subtype+'_ratios_color.eps'
  ENDIF ELSE BEGIN 
      pfile = subtype+'_ratios.eps'
  ENDELSE 
  pfile = concat_dir(pdir, pfile)

  IF keyword_set(dops) THEN BEGIN 
      begplot,pfile,color=color, /encap, xsize=8, ysize=8
      thick=2
      !p.thick=thick
      !p.charthick=thick
      !x.thick=thick
      !y.thick=thick
  ENDIF 

  IF keyword_set(color) THEN BEGIN 
      pcolor = c2i('blue')
  ENDIF ELSE BEGIN 
      pcolor = !p.color
  ENDELSE 

  t1 = self->lensread('jackknife', sub='ngals200_12', sample=sample)
  t=self->lensread('jackknife', sub=subtype, sample=sample)
  nbin=n_elements(t)
  nrad = n_elements(t[0].meanr)

  ws = self->where_string(subtype, labels=labels)
  
  xtitle = textoidl('R [h^{-1} Mpc]')
  erase & multiplot, [4,3], /square, $
      mxtitle=xtitle, xtickf='loglabels'

  xrange = [0.012,50]
  yrange = [-3.9,6.4]


  IF !d.name EQ 'PS' THEN BEGIN 
      !p.charsize=1
      !x.thick=2
      !y.thick=2
      !p.thick=2
      symsize=0.7

      lcharsize=0.6
  ENDIF 

  running_chi2 = 0.0
  i=0
  FOR il=0L,11 DO BEGIN 

      ;; Fit sis to the unsplit sample
;      fitsis, t1[iL].meanr/1000, t1[iL].sigma, t1[iL].sigmaerr, $
;        400.0, yfit, sigmav, sigmaverr

      ;; Average the two, fit to it
;      allsigma = t[i].sigma
;      allsigmaerr = allsigma
;      FOR irad=0L, nrad-1 DO BEGIN 
;          wmom, $
;            [ t[i].sigma[irad], t[i+1].sigma[irad] ], $
;            [ t[i].sigmaerr[irad], t[i+1].sigmaerr[irad] ], $
;            wmean, wsig, werr
;          allsigma[irad] = wmean
;          allsigmaerr[irad] = werr
;      ENDFOR 
;      fitpower, t[i].meanr/1000, allsigma, allsigmaerr, [10.0,-1.0], yfit

      fitpower, t1[il].meanr/1000, t1[il].sigma, t1[il].sigmaerr, [10.0, -1.0], yfit

      r = t[i].meanr/1000.0

      lsigma = t[i].sigma/yfit
      lsigmaerr = t[i].sigmaerr/yfit
      hsigma = t[i+1].sigma/yfit
      hsigmaerr = t[i+1].sigmaerr/yfit

      ratio = t[i+1].sigma/t[i].sigma
      ratioerr = ratio*sqrt( (t[i].sigmaerr/t[i].sigma)^2 + $
                             (t[i+1].sigmaerr/t[i+1].sigma)^2 )

      wmom, ratio, ratioerr, ratmean, ratsig, raterr

      pplot, r, lsigma, yerr=lsigmaerr, /nohat, psym=8, symsize=symsize, $
        /xlog, xrange=xrange,xstyle=3, xticklen=0.04, $
        yrange=yrange, ystyle=3

      pplot, r, hsigma, yerr=hsigmaerr, /nohat, psym=4, symsize=symsize, $
        /overplot,color=pcolor
      oplot, [1.e-5, 1.e5], [1,1]


      legend, labels[i:i+1], psym=[8,4], color=[!p.color, pcolor], $
        /right, charsize=lcharsize, box=0, margin=0.

      mess = $
        'ratio = '+ntostr(ratmean,4,/round)+!csym.plusminus+ntostr(raterr,4,/round)
      legend, mess, /right, /bottom, box=0,charsize=lcharsize, margin=0.

      i=i+2
      multiplot

  ENDFOR 


  multiplot,/default


  IF keyword_set(dops) THEN endplot, /trim_bbox


END 



pro maxbcg_lensing::plot_inversions

  esheldon_setup
  sub='ilum200_16'
  t = self->lensread('jackknife_invert',sub=sub)
  nt = n_elements(t)

  ws=self->where_string(sub,labels=labels)

  erase & multiplot,[4,4], /square
  
  yrange = [1.e11,5.e15]

;  if self->sample() eq 19 then begin 
;      xrange = [0.015,10]
;  endif else begin 
      xrange = [0.015,30]
;  endelse 

  if !d.name eq 'PS' then begin 
      !p.thick=2
      !x.thick=2
      !y.thick=2
  endif 

  for i=0L, nt-1 do begin 

      delvarx, xt, yt, xtickf
      if (i mod 4) eq 0 then yt='M(<r) [M'+sunsymbol()+' ]' 
      if i gt 11 then begin 
          xt=estitle('mpcxtitle')
          xtickf='loglabels'
      endif 

      pplot, t[i].ir, t[i].massout*1.e12, yerr=t[i].massout_err*1.e12, $
        /nohat, $
        psym=8, symsize=0.5, $
        /xlog, /ylog, $
        xticklen=0.04, yticklen=0.04, $
        charsize=1, $
        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
        xtitle=xt, ytitle=yt, xtickf=xtickf

      oplot, t[i].ir, t[i].mass_yfit, color=c2i('red')
      oplot, t[i].ir, t[i].bias_yfit, color=c2i('blue')
      oplot, t[i].ir, t[i].nfw_yfit, color=c2i('darkgreen')

      legend, labels[i], /left,box=0,charsize=0.5, margin=0

      if i ne nt-1 then multiplot
  endfor 
  multiplot,/reset

  return


  for i=0L, nt-1 do begin 

      pplot, t[i].ir, t[i].drho, yerr=t[i].drho_err, $
        /xlog, /ylog, yrange=[1.e-2,1.e4], ystyle=3, $
        xrange=[0.015,20],xstyle=3, aspect=1

      key = prompt_kbrd('hit a key')

      m=t[i].massout*1.e12
      me=t[i].massout_err*1.e12
      pplot, t[i].ir, m, yerr=me, $
        /xlog, /ylog, yrange=[1.e11,5.e15], ystyle=3, $
        xrange=[0.015,20],xstyle=3, aspect=1

      key = prompt_kbrd('hit a key')

  endfor 

end 


;
; Analyzing the centerclass stuff
;

function maxbcg_lensing::centerclass_rescale, d, wgood, w
    ;; Rescale good centers to this sample mean at large radius

    sigmagood = d[wgood].sigma
    wrad = where(d[0].meanr/1000 gt 2)
    wrad2 = where(d[w].sigma[wrad] gt 0 and d[wgood].sigma[wrad] gt 0)
    wrad2 = wrad[wrad2]
    scale = mean(d[w].sigma[wrad2])/mean(d[wgood].sigma[wrad2])
    print,'Rescaling by ',scale
    sigmagood = sigmagood*scale

    print,'Compare: ',mean(sigmagood[wrad2]),mean(d[w].sigma[wrad2])
    return, sigmagood

end

pro maxbcg_lensing::plot_compare2_centerclass, dtype, sample=sample, $
        rescale=rescale

    subtype='centerclass_alt2'
    psfile = self->plotfile(dtype, sub=subtype, sample=sample,/color,/encap,$
        /create)
    if keyword_set(rescale) then begin
        psfile = repstr(psfile, '.eps', '_compare_rescaled.eps')
    endif else begin
        psfile = repstr(psfile, '.eps', '_compare.eps')
    endelse

    begplot, psfile, /color, /encap
    d=self->lensread(dtype, subtype=subtype, sample=sample)

    if keyword_set(rescale) then begin
        sigmagood = self->centerclass_rescale(d, 0, 1)
    endif else begin
        sigmagood = d[0].sigma
    endelse

    xrange = [0.015, 75.0]
    yrange = [0.01, 4000]

    xtitle=textoidl('R [h^{-1} Mpc]')
    ytitle=textoidl('\Delta\Sigma [h M_{\odot} pc^{-2}]')
    pplot, d[0].meanr/1000, sigmagood, yerr=d[0].sigmaerr, $
        /xlog, /ylog, xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
        xtitle=xtitle, ytitle=ytitle, $
        aspect=1, xtickf='loglabels', ytickf='loglabels'

    pplot, d[1].meanr/1000, d[1].sigma, yerr=d[1].sigmaerr, psym=8, $
        /overplot,color=c2i('blue')

    legend, ['class 0', 'class 1'], psym=[0,8], color=[!p.color,c2i('blue')],$
        /right
    endplot,/trim

    print,['class','nlens','N200','err(N200)','L200','err(L200)','Lbcg','err(Lbcg)'], $
        format='(2A7,6A10)'
    colprint,$
        [0,1], $
        d.nlens,$
        d.mean_ngals200, d.err_ngals200, $
        d.mean_ilum200, d.err_ilum200, $
        d.mean_bcg_ilum, d.err_bcg_ilum, $
        format='(2I7, 6F10.2)'



end

pro maxbcg_lensing::plot_compare_centerclass, dtype, $
        sample=sample, subtype=subtype, $
        rescale=rescale

    psfile = self->plotfile(dtype, sub=subtype, sample=sample,/color,/encap,$
                /createdir)
    if keyword_set(rescale) then begin
        psfile = repstr(psfile, '.eps', '_compare_rescaled.eps')
    endif else begin
        psfile = repstr(psfile, '.eps', '_compare.eps')
    endelse
    begplot, psfile, /color, /encap, xsize=8, ysize=8
    !p.thick=2
    !x.thick=2
    !y.thick=2
    !p.charthick=2

    d = self->lensread(dtype, sub=subtype, sample=sample)
    nbin = n_elements(d)

    binvals = self->centerclass_kdeclass_binvals(subtype)
    w99=where(binvals eq -9999, n99)
    if n99 ne 0 then binvals[w99] = -1


    ;wrad = where(d.meanr/1000 gt 2 and d.meanr/1000 lt 10)
    wrad = where(d[0].meanr/1000 gt 2)

    ;; The "good" centers
    wgood=where(binvals eq 0, nw, comp=wrest, ncomp=nrest)

    mplot_value=self->mplot_value(nrest)
    print,'mplot_value = ',mplot_value
    erase & multiplot, mplot_value, /square, $
        mxtitle=textoidl('R [h^{-1} Mpc]'), $
        mytitle=textoidl('\Delta\Sigma [h M_{\odot} pc^{-2}]'), $
        xtickformat='loglabels', ytickformat='loglabels', $
        mxtitoffset=1, mytitoffset=1, mxtitsize=1.5, mytitsize=1.5
    xrange = [0.015, 75.0]
    yrange = [0.01, 4000]
    for i=0L, nrest-1 do begin

        w=wrest[i]

        pplot, d[w].meanr/1000, d[w].sigma, yerr=d[w].sigmaerr, $
            psym=8, symsize=0.7, /nohat, $
            /xlog, /ylog, xrange=xrange, xstyle=3, yrange=yrange, ystyle=3,$
            xticklen=0.04, yticklen=0.04



        sigmagood = d[wgood].sigma
        gleg = 'class=0'
        if keyword_set(rescale) then begin
            ;; Rescale good centers to this sample mean at large radius

            wrad2 = where(d[w].sigma[wrad] gt 0 and d[wgood].sigma[wrad] gt 0)
            wrad2 = wrad[wrad2]
            ;wmom, d[w].sigma[wrad2], d[w].sigmaerr[wrad2], wm, wsig, werr
            ;wmom, d[wgood].sigma[wrad2], d[wgood].sigmaerr[wrad2], wmg, wsig, werr
            ;scale = wm/wmg
            scale = mean(d[w].sigma[wrad2])/mean(d[wgood].sigma[wrad2])
            print,'Rescaling by ',scale
            sigmagood = sigmagood*scale
            gleg = 'scaled '+gleg

            print,'Compare: ',mean(sigmagood[wrad2]),mean(d[w].sigma[wrad2])

            ;pplot, d[wgood].meanr[wrad2]/1000, sigmagood[wrad2],/over, color=c2i('blue'), psym=7
            ;pplot, d[w].meanr[wrad2]/1000, d[w].sigma[wrad2],/over, color=c2i('green'), psym=7
        endif

        pplot, d[wgood].meanr/1000, sigmagood, color=c2i('blue'),/over


        legend, ['class='+ntostr(binvals[w]),gleg], $
            psym=[8,0], color=[!p.color, c2i('blue')],/right, $
            charsize=0.7

        multiplot
    endfor

    multiplot, /reset

    endplot, /trim






    print,['class','nlens','N200','err(N200)','L200','err(L200)','Lbcg','err(Lbcg)'], $
        format='(2A7,6A10)'
    colprint,$
        binvals, $
        d.nlens,$
        d.mean_ngals200, d.err_ngals200, $
        d.mean_ilum200, d.err_ilum200, $
        d.mean_bcg_ilum, d.err_bcg_ilum, $
        format='(2I7, 6F10.2)'



    psfile = self->plotfile(dtype, sub=subtype, sample=sample,/encap)
    psfile = repstr(psfile, '.eps', '_stats.eps')
    begplot, psfile, /color, /encap, xsize=11
    !p.thick=2
    !x.thick=2
    !y.thick=2
    !p.charthick=2



    erase & multiplot, /square, [3,1], xgap=0.05
    ;erase & multiplot, [3,1], xgap=0.05

    xrange = [-2,5]
    yrange = [0.7, 1.7]
    calc_ratio_cov, $
        d.mean_ngals200, d.err_ngals200, $
        replicate(d[wgood].mean_ngals200,nrest+1), $
        replicate(d[wgood].err_ngals200,nrest+1), $
        rat, raterr
    raterr[wgood] = 0
    pplot, binvals, rat, yerr=raterr, psym=-8, $
        xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, $
        ytitle=textoidl('N_{200}^i/N_{200}^{0}'), $
        xtitle='Center Class'

    calc_ratio_cov, $
        d.mean_ilum200, d.err_ilum200, $
        replicate(d[wgood].mean_ilum200,nrest+1), $
        replicate(d[wgood].err_ilum200,nrest+1), $
        rat, raterr
    raterr[wgood] = 0
    multiplot, /doxaxis, /doyaxis
    pplot, binvals, rat, yerr=raterr, psym=-8, $
        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3, $
        xtitle='Center Class', $
        ytitle=textoidl('L_{200}^i/L_{200}^{0}')


    calc_ratio_cov, $
        d.mean_bcg_ilum, d.err_bcg_ilum, $
        replicate(d[wgood].mean_bcg_ilum,nrest+1), $
        replicate(d[wgood].err_bcg_ilum,nrest+1), $
        rat, raterr
    raterr[wgood] = 0
    multiplot, /doyaxis
    pplot, binvals, rat, yerr=raterr, psym=-8, $
        xrange=xrange, xstyle=3, yrange=yrange, ystyle=3,$
        xtitle='Center Class', $
        ytitle=textoidl('L_{bcg}^i/L_{bcg}^{0}')


    multiplot, /default

    endplot,/trim

end

; make a csv table and the first few lines of a tex table
pro maxbcg_lensing::create_paper_csv


    N = self->lensread('jackknife', sub='ngals200_12', sample=[21,22])
    L = self->lensread('jackknife', sub='ilum200_16', sample=[21,22])

    dir=self->lensdir('jackknife',sample=[21,22])
    dir = repstr(dir, 'jackknife', 'csv')

    if not fexist(dir) then file_mkdir, dir
    file = 'maxbcg_sample21-22_ngals200_12_ilum200_12_jackknife.csv'
    texfile = 'maxbcg_sample21-22_ngals200_12_ilum200_12_jackknife.tex'
    file = path_join(dir, file)
    texfile = path_join(dir, texfile)

    print,'Writing to file: ',file
    openw, lun, file, /get


    nbin = n_elements(N)
    nrad = n_elements(N[0].meanr)

    corrs = 'corr'+ntostr(lindgen(nrad))
    corrs = strjoin(corrs, ',')
    printf, lun, $
        'type,binval,rmpc,delta_sigma,delta_sigma_err,'+corrs

    for ibin=0L, nbin-1 do begin
        for irad=0L, nrad-1 do begin
            binval = N[ibin].mean_ngals200
            if (binval lt 9) then binval=rnd(binval)
            printf, lun, $
                'N',binval, N[ibin].meanr[irad]/1000, $
                N[ibin].sigma[irad], N[ibin].sigmaerr[irad], $
                f='(a,",",e0,",",e0,",",e0,",",e0,$)'
            for irad2=0L,nrad-1 do begin
                printf, lun, N[ibin].correlation[irad,irad2],$
                    f='(",",e0,$)'
            endfor
            printf, lun
        endfor
    endfor

    nbin = n_elements(L)
    nrad = n_elements(L[0].meanr)
    for ibin=0L, nbin-1 do begin
        for irad=0L, nrad-1 do begin
            printf, lun, $
                'L',L[ibin].mean_ngals200, L[ibin].meanr[irad]/1000, $
                L[ibin].sigma[irad], L[ibin].sigmaerr[irad], $
                f='(a,",",e0,",",e0,",",e0,",",e0,$)'
            for irad2=0L,nrad-1 do begin
                printf, lun, L[ibin].correlation[irad,irad2],$
                    f='(",",e0,$)'
            endfor
            printf, lun
        endfor
    endfor

    free_lun, lun

    print,'Writing to file: ',texfile
    openw, lun, texfile, /get

    printf, lun,'\begin{deluxetable}{ccccc}'
    printf, lun,'\tablecaption{\deltasig\ data for \maxbcg\ Clusters \label{tab:deltasig}}'
    printf, lun,'\tablewidth{0pt}'
    printf, lun,'\tablehead{'
    printf, lun,'    \colhead{$X$}       &'
    printf, lun,'    \colhead{$\langle X \rangle$}         &'
    printf, lun,'    \colhead{$r$}  &'
    printf, lun,'    \colhead{$\Delta \Sigma$}  &'
    printf, lun,'    \colhead{$\sigma(\Delta \Sigma)$}\\ '
    printf, lun,'    (N or L) &'
    printf, lun,'     &'
    printf, lun,'    ($h^{-1}$ Mpc) &'
    printf, lun,'    ($h M_{\odot}$ pc$^{-2}$) &'
    printf, lun,'    ($h M_{\odot}$ pc$^{-2}$)'
    printf, lun,'}'
    printf, lun,'\'
    printf, lun,'\startdata'

    ;format = '(a," & ",e0," & ",e0," & ",e0," & ",e0,$)'
    format = '(a," & ",e10.3," & ",e10.3," & ",e10.3," & ",e10.3,$)'
    nprint = 10
    count=0L
    for ibin=0L, nbin-1 do begin
        for irad=0L, nrad-1 do begin
            binval = N[ibin].mean_ngals200
            if (binval lt 9) then binval=rnd(binval)
            printf, lun, $
                'N',binval, N[ibin].meanr[irad]/1000, $
                N[ibin].sigma[irad], N[ibin].sigmaerr[irad], $
                f=format
            count = count+1
            if count gt nprint then begin
                printf, lun
                break
            endif
            printf, lun, ' \\'
        endfor
        if count gt nprint then break
    endfor
    printf, lun, '\enddata'
    printf, lun, "\tablecomments{\deltasig\ data corresponding to figures "
    printf, lun, "\ref{fig:deltasig_ngals200_12} and \ref{fig:deltasig_ilum200_16}."
    printf, lun, "The first column indicates binning on either \nvir, labeled ``N'', "
    printf, lun, "or \lvir, labeled ``L''. The second column is the mean value of "
    printf, lun, "\nvir\ or \lvir\ for the bin, with units of either number or  "
    printf, lun, "$10^{10} h^{-2} L_{\odot}$.  This table "
    printf, lun, "is presented in its entirety in the electronic edition of "
    printf, lun, "the {\it Astrophysical Journal}, including the full dimensionless "
    printf, lun, "correlation matrix for \deltasig.  A portion is shown here for "
    printf, lun, "guidance regarding its form and content.}"
    printf, lun, '\end{deluxetable}'


    free_lun, lun

end


pro maxbcg_lensing::print_jeremy, subtype

	outdir='~/tmp/ngals200_12_deltasigma_for_jeremy'
	if not fexist(outdir) then file_mkdir, outdir

    t = self->lensread('jackknife', sub=subtype, sample=[21,22])
	nt=n_elements(t)

	f='($,g)'
	nrad=n_elements(t[0].meanr)
	for i=0L, nt-1 do begin
		file=path_join(outdir, subtype+'_deltasigma_bin'+ntostr(i, format='(i02)')+'.dat')

		openw, lun, file, /get_lun

		corr = t[i].correlation
		for j=0L, nrad-1 do begin
			printf, lun, float(t[i].meanr[j]/1000), f=f
			printf, lun, float(t[i].sigma[j]), f=f
			printf, lun, float(t[i].sigmaerr[j]), f=f
			for k=0L, nrad-1 do begin
				printf, lun, float(corr[j,k]), f=f 
			endfor
			printf, lun
		endfor

		free_lun, lun
	endfor
	

end


pro maxbcg_lensing::compare_lxid
    common compare_lxid_block, cat, lxid
    if n_elements(cat) eq 0 then begin
        cat = self->get()
    endif

    dir=self->lxid_dir()

    ; first compare histograms of n200
    psfile=path_join(dir, 'compare-n200-lxid.eps')
    begplot, psfile, /color, /encap
    color2 = 'blue'
    color3='red'
    color4='darkgreen'
    plothist, alog10(cat.ngals200), $
        bin=0.1, min=alog10(1), /ylog, /norm, yrange=[1.e-6,1.e2], $
        xtitle=textoidl('log_{10}(N_{200})'), $
        aspect=1
    sig3match = self->lxid_match(cat, 'lxid-sig3')
    sig5match = self->lxid_match(cat, 'lxid-sig5')
    norasmatch = self->lxid_match(cat, 'lxid-noras')

    plothist, alog10(cat[*sig3match].ngals200), bin=0.1, min=alog10(1), /norm, /overplot, color=color2
    plothist, alog10(cat[*sig5match].ngals200), bin=0.1, min=alog10(1), /norm, /overplot, color=color3
    plothist, alog10(cat[*norasmatch].ngals200), bin=0.1, min=alog10(1), /norm, /overplot, color=color4

    ptr_free, sig3match, sig5match, norasmatch

    plegend, ['all','sig3','sig5','noras'], color=[!p.color, c2i(color2), c2i(color3),c2i(color4)], line=0, $
        /right

    endplot, /trim, /png, dpi=90

    ; now compare some stats
    psym2 = 7
    psym3 = 4
    psym4 = 6
    psfile=path_join(dir, 'compare-dsig-lxid.eps')
    begplot, psfile, /color, /encap

    n12 = self->lensread('jackknife', sub='ngals200_12')
    sig3 = self->lensread('jackknife', sub='lxid-sig3')
    sig5 = self->lensread('jackknife', sub='lxid-sig5')
    noras = self->lensread('jackknife', sub='lxid-noras')

    minr = 0.2
    maxr = 2.0
    minr_str = string(minr, f='(f0.1)')
    maxr_str = string(maxr, f='(f0.1)')
    w=where(n12[7].meanr ge minr*1000.0 and n12[7].meanr le maxr*1000)

    wmom, sig3.sigma[w], sig3.sigmaerr[w], mean_dsig_sig3, junk, err_dsig_sig3
    wmom, sig5.sigma[w], sig5.sigmaerr[w], mean_dsig_sig5, junk, err_dsig_sig5
    wmom, noras.sigma[w], noras.sigmaerr[w], mean_dsig_noras, junk, err_dsig_noras
    for i=0L, n_elements(n12)-1 do begin
        wmom, n12[i].sigma[w], n12[i].sigmaerr[w], wmean, wsig, werr
        add_arrval, wmean, mean_dsig12
        add_arrval, werr, err_dsig12
    endfor

    pplot, n12.mean_ngals200, mean_dsig12, yerr=err_dsig12, psym=8, $
        /xlog, ylog=0, $
        xtitle=textoidl('N_{200}'), $
        ytitle=textoidl('<\Delta\Sigma>('+minr_str+' < R < '+maxr_str+')'), $
        aspect=1, xstyle=3, ystyle=3
    pplot, [sig3.mean_ngals200], [mean_dsig_sig3], yerr=[err_dsig_sig3], psym=psym2, /over, color=color2
    pplot, [sig5.mean_ngals200], [mean_dsig_sig5], yerr=[err_dsig_sig5], psym=psym3, /over, color=color3
    pplot, [noras.mean_ngals200], [mean_dsig_noras], yerr=[err_dsig_noras], psym=psym4, /over, color=color4

    leg = ['all   <z>='+string(mean(n12.mean_z),f='(f0.2)'), $
           'sig3  <z>='+string(sig3.mean_z,f='(f0.2)'), $
           'sig5  <z>='+string(sig5.mean_z,f='(f0.2)'), $
           'noras <z>='+string(noras.mean_z,f='(f0.2)')]
    plegend, $
        leg, $
        color=[!p.color, c2i(color2), c2i(color3), c2i(color4)], $
        psym=[8,psym2,psym3,psym4]

    endplot, /trim, /png, dpi=90
end



PRO maxbcg_lensing__define

  struct = { $
             maxbcg_lensing, $
             sample: 0, $
             INHERITS maxbcg $
           }

END 
