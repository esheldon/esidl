function maxbcg_correlate::init, sample

    par = correlate_pardef()
    maxbcg_catalog = 'dr406'

    case sample of
        1: begin 
            ;; rmin = 0.02, rmax=11.5, nbin=18, ... 
            ;; see correlate_pardef for other defaults
            par.psample = 'p01'   
            par.DPsample = 'maxbcg01'
        end 
        2: begin 
            par.psample = 'p01'
            par.DPsample = 'maxbcg02'
        end 
        3: begin ;; This is a test sample to 10Mpc with just radbins lum. 
            par.psample = 'p01' 
            par.DPsample = 'maxbcg02'
            par.numerator_output_type='cl_r'
            par.nlum=0
            par.nkgmr=0
        end
        ; This is the old way where we set nlum=0 but that is not right
        ; we still need to set nlum and nkgmr because of the way the code
        ; works
        4: begin 
            par.psample = 'p02'  ;; new psample since parameters changed. 
            par.DPsample = 'maxbcg02'
            par.rmin = 0.014
            par.rmax =  36
            par.nrad = 22
            par.nlum = 0
            par.nkgmr = 0
            par.numerator_output_type='cl_r'
        end 


        5: begin 
            ;; This is sample 1 from corrsim. Only read this from here, 
            ;; don't try to create inputs or anything.
            maxbcg_catalog = 'hv1'
            par.psample = 'p03'   
            par.DPsample = 'hv1bcg'
            par.RPsample = 'rhv1bcg'
            par.DSsample = 'hv1gals'
            par.RSsample = 'rhv1gals'

            par.omega_m = 0.3
        end 
        6: begin 
            ;; This is sample 2 from corrsim. Only read this from here, 
            ;; don't try to create inputs or anything.

            maxbcg_catalog = 'hv1halo'
            par.psample = 'p03'   
            par.DPsample = 'hv1halo'
            ; We can re-use all these since we re-sample the redshifts
            ; anyway and the secondaries are the same.
            par.RPsample = 'rhv1bcg'
            par.DSsample = 'hv1gals'
            par.RSsample = 'rhv1gals'

            par.omega_m = 0.3
        end 



        ; We will do luminosity limits in r-band
        ; We don't need to worry about the secondary catalog I checked the
        ; flux limit is ok.  In other words, we can use all the usual 
        ; inputs
        7: begin
            par.psample = 'p04'  ;; new psample since parameters changed. 
            par.DPsample = 'maxbcg02'
            ; this gives a mean radius of 2 Mpc for the last bin
            par.rmax = 2.802
            par.nrad = 14

            ; r-band.  lmin corresponds to Mr = -19
            par.lumband = 2
            par.loglmin = 9.59
            par.loglmax = 11.83
            par.numerator_output_type='cl_r'
        end
        ; same as above but full cube
        8: begin
            par.psample = 'p04'  ;; new psample since parameters changed. 
            par.DPsample = 'maxbcg02'
            ; this gives a mean radius of 2 Mpc for the last bin
            par.rmax = 2.802
            par.nrad = 14

            ; r-band.  lmin corresponds to Mr = -19
            par.lumband = 2
            par.loglmin = 9.59
            par.loglmax = 11.83
        end


        else: message,'Unsupported sample number: '+ntostr(sample)
    endcase 

    cinit = self->correlate::init(par)
    mbinit = self->maxbcg::init(maxbcg_catalog)
    return, (cinit+mbinit) eq 2

end 








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; create inputs
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





function maxbcg_correlate::zrange
  return, [0.1, 0.3]
end 

function maxbcg_correlate::cuts, clambda, ceta, z, nkeep

  ;; run through 1mpc cut that ben uses. this is slightly
  ;; different than ben's, so we should do it on both real
  ;; and random

  ntot = n_elements(clambda)

  par = self->par_struct()

  zrange = self->zrange()

  print
  print,'applying redshift cuts'
  w = where(z ge zrange[0] and z le zrange[1], nkeep)

  print,'Kept '+ntostr(nkeep)+'/'+ntostr(ntot)

  rMpc = 1.0
  maxangle = rMpc/angdist_lambda(z[w], omega=par.omega_m, h=par.h)*180.0/!pi

  print
  print,'Applying edge cuts'

  apply_pixel_mask, $
    clambda[w], ceta[w], masked, unmasked, $
    maxangle = maxangle, /basic


  if unmasked[0] eq -1 then begin 
      nkeep = 0 
      keep = -1
  endif else begin 
      nkeep = n_elements(unmasked)
      keep = w[unmasked]
  endelse 
  
  print,'Kept '+ntostr(nkeep)+'/'+ntostr(ntot)
  return,keep


end 


function maxbcg_correlate::genrand, type

  if n_elements(type) eq 0 then begin  
      print,'-Syntax: randstruct = mb->genrand(type)'
      print,'type = "primary" or "secondary"'
      on_error, 2
      message,'Halting'
  endif 

  par = self->par_struct()


  if type eq 'primary' then begin 

      ;; Get some random clambda, ceta points
      print
      print,'Generating positions'
      rand = self->maxbcg::genrand()
      
      csurvey2eq, rand.clambda, rand.ceta, rra, rdec
      
      nrand = n_elements(rand)


      print
      print,'Generating redshifts'
      cosmo = obj_new('cosmology')
      
      zrange = self->zrange()
      zrand = cosmo->genrandz(nrand, $
                              zrange[0], zrange[1], $
                              omega_m=par.omega_m)
      obj_destroy, cosmo
      
      print
      print,'Making edge cut'
      rkeep = self->cuts(rand.clambda, rand.ceta, zrand, nrkeep)

      print,'Kept: ',nrkeep
      struct = self->primary_struct(nrkeep)

      rra = rra[rkeep]
      rdec = rdec[rkeep]
      zrand = zrand[rkeep]
      bcg_id = lindgen(nrkeep)
      
      struct.bcg_id = bcg_id
      struct.ra = rra
      struct.dec = rdec
      struct.z = zrand
      
  endif else begin 

      ;; Get some random clambda, ceta points
      print
      print,'Generating positions'
      self->genrand, rlam, reta, nperfile=40000000
      
      print,'Converting to ra/dec'
      csurvey2eq, temporary(rlam), temporary(reta), rra, rdec
      
      nrand = n_elements(rra)

      struct = self->secondary_struct(nrand, /random)

      print,'Getting htm_index'
      struct.htm_index = htm_index(rra, rdec, par.depth)

      print,'copying to struct'
      struct.ra = temporary(rra)
      struct.dec = temporary(rdec)
  endelse 
  return,struct

end 


pro maxbcg_correlate::primary_input, randnum=randnum


  nrand = n_elements(randnum)
  if nrand eq 0 then begin 

      bcg=self->get()

      ntot = n_elements(bcg)
      ngals_min = 3
      print
      print,'Selecting ngals >= '+ntostr(ngals_min)
      keep = where(bcg.ngals ge ngals_min, nkeep)

      print
      print,'Kept '+ntostr(nkeep)+'/'+ntostr(ntot)

      bcg = bcg[keep]

      keep = $
        self->cuts(bcg.clambda, bcg.ceta, bcg.photoz_cts, nkeep)

      outst = self->primary_struct(nkeep)

      outst.bcg_id = bcg[keep].bcg_id

      outst.ra = bcg[keep].ra
      outst.dec = bcg[keep].dec

      outst.z = bcg[keep].photoz_cts

      outfile = self->corrfile('primary','input')

      print
      print,'Writing to file: ',outfile
      write_idlstruct, outst, outfile

  endif else begin 

      for i=0l, nrand-1 do begin

          print,'------------------------------------------------------------'
          rstruct = self->genrand('primary')
          routfile = $
            self->corrfile('primary','input', randnum=randnum[i])
          print
          print,'Writing random file: ',routfile
          write_idlstruct, rstruct, routfile

      endfor 

  endelse 


end 


; These secondary methods should actually probably go into a separate class.  
; Not really tied to maxbcg
pro maxbcg_correlate::secondary_input, randnum=randnum

  par = self->par_struct()

  nrand = n_elements(randnum)
  if nrand eq 0 then begin 
      
      pg = obj_new('postgres')

      res = pg->query('select stripes from maxbcg_input_meta')
      stripes = res.stripes
      nstripe = n_elements(stripes)

      outfile = self->corrfile('secondary','input')      
      print,outfile

      for ist=0l, nstripe-1 do begin 

          print,'------------------------------------------------------------'
          stripe = stripes[ist]

          query = $
            'SELECT '+$
              'ra,dec,clambda,ceta,gflux,rflux,iflux '+$
            'FROM '+$
              'maxbcg_input '+$
            'WHERE '+$
              'stripe = '+ntostr(stripe)

          print,query
          stripe_struct = self->postgres::query(query)
          ntot = n_elements(stripe_struct)

          print
          print,ntostr(ntot)+' objects in stripe '+ntostr(stripe)

          outst = self->secondary_struct(ntot)

          outst.ra = stripe_struct.ra
          outst.dec = stripe_struct.dec

          outst.gflux = stripe_struct.gflux
          outst.rflux = stripe_struct.rflux
          outst.iflux = stripe_struct.iflux

          outst.htm_index = htm_index(outst.ra, outst.dec, par.depth)

          if ist eq 0 then begin 
              print
              print,'Writing to file: ',outfile
              write_idlstruct, outst, outfile
          endif else begin 
              print
              print,'Appending to file: ',outfile
              write_idlstruct, outst, outfile, /append
          endelse 

      endfor 

  endif else begin 

      for i=0l, nrand-1 do begin

          print,'------------------------------------------------------------'
          rstruct = self->genrand('secondary')
          routfile = $
            self->corrfile('secondary','input',randnum=randnum[i])
          print
          print,'Writing random file: ',routfile
          write_idlstruct, rstruct, routfile

      endfor 

  endelse 


end 



; also not directly tied to maxbcg
function maxbcg_correlate::kcorr_description
  nz = 80L
  ngmr = 21L
  nrmi = 21L
  nband = 5L

  zmin = 0.0
  zstep = 0.005
  z = zmin + zstep*findgen(nz)

  gmrmin = 0.0
  gmrstep = 0.1
  gmr = gmrmin + gmrstep*findgen(ngmr)

  rmimin = 0.0
  rmistep = 0.05
  rmi = rmimin + rmistep*findgen(nrmi)

  bands = lindgen(5)

  struct = $
    {nz: nz, $
     zstep: zstep, $
     z: z, $
     $
     ngmr: ngmr, $
     gmrstep: gmrstep, $
     gmr: gmr, $
     $
     nrmi: nrmi, $
     rmistep: rmistep, $
     rmi: rmi, $
     $
     nband: nband, $
     bands: bands $
    }
  return, struct
end 
pro maxbcg_correlate::make_kcorr

  desc = self->kcorr_description()

  ;; The original file
  dir = self->correlate_dir('input')
  dir = concat_dir(dir, 'original_kcorrtable')
  file = 'mm_kcortable-maxbcg.dat'
  file = concat_dir(dir, file)

  kcorr=fltarr(desc.nz,desc.ngmr,desc.nrmi,desc.nband)
  openr,1,file
  readf,1,kcorr
  close,1

  w = where(kcorr eq 0.0, nw)
  if nw ne 0 then kcorr[w] = -9999.0

  ;; Write it out so that it can be read and indexed the same
  ;; way in C++

  outfile = self->corrfile('kcorr','input')
  print
  print,'Writing to file: ',outfile

  openw, lun, outfile, /get_lun

  writeu, lun, desc.nz
  writeu, lun, desc.z

  writeu, lun, desc.ngmr
  writeu, lun, desc.gmr


  writeu, lun, desc.nrmi
  writeu, lun, desc.rmi

  writeu, lun, desc.nband
  writeu, lun, desc.bands

  for iz=0l, desc.nz-1 do begin 
      for igmr=0l, desc.ngmr-1 do begin 
          for irmi=0l, desc.nrmi-1 do begin 
              for ib=0l, desc.nband-1 do begin 

                  writeu, lun, kcorr[iz, igmr, irmi, ib]

              endfor
          endfor
      endfor
  endfor 

  free_lun, lun

end 






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Might make this more general
pro maxbcg_correlate::plot_relative_error

  sub = 'ilum200_16'
  t=self->corr_read('dd','corrected',sub=sub)
  j=self->corr_read('dd','jackknife',sub=sub)

  nsub=n_elements(j)

  ws=self->where_string('ilum200_16',labels=labels)

  begplot,'~/tmp/relative_error_ilum200_16.ps',/color

  n = 5

  xt=!mpcxtitle2 
  yt='fractional error'

  simpctable,color=clist
  pplot, j[n].r, j[n].radnumdens_err/j[n].radnumdens,$
    /xlog,/ylog,xrange=[0.01,1000],xstyle=3,xtit=xt,ytit=yt, $
    aspect=1,$
    xtickf='loglabels', ytickf='loglabels'
  for i=0,nsub-1 do begin 
      pplot, j[i].r, j[i].radnumdens_err/j[i].radnumdens,/over,color=clist[i]
  endfor 

  legend,labels,/right,line=0,color=clist[0:nsub-1],$
    box=0,charsize=0.7


  yrange=[5.e-3,0.2]
  tit = labels[n]
  pplot, j[n].r, j[n].radnumdens_err/j[n].radnumdens,$
    /xlog,/ylog,$
    xrange=[0.01,30],xstyle=3,yrange=yrange,ystyle=3,$
    xtit=xt,ytit=yt,tit=tit, aspect=1.0,$
    xtickf='loglabels', ytickf='loglabels'
  pplot, t[n].r, t[n].radnumdens_err/t[n].radnumdens,/overplot,color=!blue
  legend,['jackknife','poisson'],line=[0,0],color=[!p.color,!blue],/right,box=0


  endplot

end 


pro maxbcg_correlate::run_plot_profile_over

    esheldon_setup
    self->plot_profile_over, $
        'jackknife', 'ilum200_16', 'radilumdens',$
        ytitle = 'L!Di!N/Area [ Lsun Mpc!U-2!N]', $
        yrange=[1.e-2,400],xrange=[1.e-2,5000],/dops,aspect=1.5, /reverse_labels
    self->plot_profile_over, $
        'jackknife', 'ilum200_16', 'radnumdens',$
        ytitle='#/Area [ h!U2!N Mpc!U-2!N ]',$
        yrange=[1.e-2,400],xrange=[1.e-2,5000],/dops,aspect=1.5, /reverse_labels

end


pro maxbcg_correlate::print_jeremy, subtype

	outdir='~/tmp/ngals200_12_Mr19_for_jeremy'
	if not fexist(outdir) then file_mkdir, outdir

	t=self->corr_read('dd','jackknife', subtype=subtype)
	nt=n_elements(t)

	f='($,g)'
	nrad=n_elements(t[0].r)
	for i=0L, nt-1 do begin
		file=path_join(outdir, subtype+'_bin'+ntostr(i, format='(i02)')+'.dat')

		openw, lun, file, /get_lun

		corr = cov2corr(t[i].radnumdens_cov)
		for j=0L, nrad-1 do begin
			printf, lun, float(t[i].r[j]), f=f
			printf, lun, float(t[i].radnumdens[j]), f=f
			printf, lun, float(t[i].radnumdens_err[j]), f=f
			for k=0L, nrad-1 do begin
				printf, lun, float(corr[j,k]), f=f 
			endfor
			printf, lun
		endfor

		free_lun, lun
	endfor
	

end


pro maxbcg_correlate::print_jeremy_thresholds_readme, lun, subtype, t

	band=2
	band_shift=0.25
	solarmags = k_solar_magnitudes(band_shift=band_shift, /silent)
	lumthresh = t.loglbins_min
	magthresh = -2.5*lumthresh + solarmags[band]
	

	nbin=self->subtype_nbin(subtype)
    self->ngals200_bins, nbin, lowlim, highlim
	printf, lun, 'h=1'
	printf, lun, 'luminosity thresholds in r-band'
	printf, lun, 'k-corrected to z=0.25'
	printf, lun
	printf, lun, 'ngals200 binning explanation:'
	printf, lun, '    bin   constraint'
	for i=0L, n_elements(t)-1 do begin
		printf,lun,i,f='("    ",i02,$)'
		if lowlim[i] eq highlim[i] then begin
			printf,lun,'    n200 = '+ntostr(lowlim[i])
		endif else begin
			printf,lun,'    '+$
				ntostr(lowlim[i])+' <= n200 <= '+ntostr(highlim[i])
		endelse
	endfor


	printf,lun
	printf,lun,'luminosity threshold explanation'
	printf, lun, '    bin      Llim         Mlim'
	for i=0L, t[0].nlum-1 do begin
		printf, lun, i, f='("    ",i02,$)'
		printf, lun, float(lumthresh[i]), f='(g,$)'
		printf, lun, magthresh[i]
	endfor

end

pro maxbcg_correlate::print_jeremy_thresholds_header, lun, t, ilum
	band_shift = 0.25
	band=2
	lumthresh = t.loglbins_min[ilum]
	solarmags = k_solar_magnitudes(band_shift=band_shift, /silent)
	magthresh = -2.5*lumthresh + solarmags[band]
	
	printf, lun, '# Threshold bin: log(Lr) > '+ntostr(lumthresh)
	printf, lun, '#                Mr < '+ntostr(magthresh)
	printf, lun, '# h=1.0, k-corrections to z=0.25'
	printf, lun, '# Columns for row i are: r_i(Mpc) n_i(Mpc^{-2}) C(i,j)'
	printf, lun, '#    where i,j run from 0,nrad-1 -- WARNING C(i,j) not yet implemented!!'
end

pro maxbcg_correlate::print_jeremy_thresholds, subtype
	
	t=self->corr_read('dd','corrected', subtype=subtype)
	nrad=t[0].nrad
	nlum=t[0].nlum

	;t = extract_tags(t, ['r','radnumdens_cumlum','radnumdens_cumlum_err'])

	; same as above but luminosity thresholds.  No covariance
	; matrices yet

	maindir='~/tmp/ngals200_12_thresholds_'+ntostr(nlum)+'_for_jeremy'
	if not fexist(maindir) then file_mkdir, maindir

	file = path_join(maindir, 'README.txt')
	print,'Writing README: ',file
	openw, lun, file, /get_lun
	self->print_jeremy_thresholds_readme, lun, subtype, t
	free_lun, lun

	nt=n_elements(t)

	f='($,g)'

	for ilum=0L,nlum-1 do begin
		lstr = ntostr(ilum,f='(i02)')
		dirname = 'threshold'+lstr
		outdir = path_join(maindir, dirname)
		if not fexist(outdir) then file_mkdir, outdir

		for i=0L, nt-1 do begin

			file = dirname+'_'+subtype+'_bin'+ntostr(i,f='(i02)')+'.dat'
			file=path_join(outdir, file)

			openw, lun, file, /get_lun
			print,'Writing to file: ',file

			self->print_jeremy_thresholds_header,lun, t[i], ilum

			; scale the covariance matrix of the whole sample to the 
			; diagonal errors

			;corr = cov2corr(t[i].radnumdens_cov)
			for irad=0L,nrad-1 do begin
				printf, lun, float(t[i].r[irad]), f=f
				printf, lun, float(t[i].radnumdens_cumlum[irad,ilum]), f=f
				printf, lun, float(t[i].radnumdens_cumlum_err[irad,ilum]), f=f

				;for k=0L, nrad-1 do begin
					;	printf, lun, float(corr[j,k]), f=f 
				;endfor
				printf, lun
			endfor ; radius
			free_lun, lun

		endfor ; clusterbin

	endfor ; luminosity threshold

end



PRO maxbcg_correlate__define

  struct = {$
             maxbcg_correlate, $
             maxbcg_correlate_dummy_variable: 0, $
             INHERITS correlate, $
             INHERITS maxbcg $
           }

END 
