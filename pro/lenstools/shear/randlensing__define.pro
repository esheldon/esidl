function randlensing::init, sample

  on_error, 2
  if n_elements(sample) eq 0 then begin
      message,'You must initialize the sample',/inf
      return,0
  endif 

  sample = fix(sample)
  self.sample = sample
  rmin = 20.0
  logbin = 1

  par = { $
          sample: sample , $
          catalog: '', $
          source_sample: 1, $
          max_allowed_angle: 6.0, $
          scinv_sample: -1, $
          edgecuts: 1 $
        }

  case sample of
      1: begin 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style=1
          shape_correction_style = 3
      end 
      2: begin 
          rmax = 36567.0
          nbin = 21
          sigmacrit_style=1
          shape_correction_style = 3
      end 

      ;; lowz randoms
      3: begin 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'lowz'
          par.source_sample = 6 ; r > 1/3
          par.max_allowed_angle = 20.0
      end 

      ;; new maxbcg randoms
      4: begin 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'maxbcg_dr406'
          par.source_sample = 6 ; r > 1/3
      end 
      5: begin 
          rmax = 36567.0
          nbin = 21
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'maxbcg_dr406'
          par.source_sample = 6 ; r > 1/3
      end 

      ;; these are using the older catalogs
      9: begin 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style=1
          shape_correction_style = 3
      end 
      10: begin 
          rmax = 36567.0
          nbin = 21
          sigmacrit_style=1
          shape_correction_style = 3
      end 

      ;; new maxbcg randoms
      11: begin 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'maxbcg_dr406'
          par.source_sample = 6 ; r > 1/3
      end 
       12: begin 
          rmax = 36567.0
          nbin = 21
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'maxbcg_dr406'
          par.source_sample = 6 ; r > 1/3
      end 
      13: begin 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 1           ; photoz as truth
          shape_correction_style = 3    ; princeton regauss
          par.catalog = 'maxbcg_dr406'
          par.source_sample = 6 ; r > 1/3
      end 
      14: begin 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 3           ; using mean scinv interpolated
          shape_correction_style = 3    ; princeton regauss
          par.scinv_sample = 2
          par.catalog = 'maxbcg_dr406'
          par.source_sample = 6 ; r > 1/3
      end 
 

      ;; new maxbcg randoms
      15: begin 
		  ; same as sample 11 but to higher redshift
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'gmbcg10_dr4' ; new catalog but old mask
          par.source_sample = 6 ; r > 1/3
      end 


      ; new photozs
      else: message,'Unknown random_sample: '+ntostr(random_sample)
  endcase 


  objshear_retval = $
    self->objshear::init('random',rmin, rmax, nbin, $
                         sigmacrit_style, shape_correction_style, $
                         logbin=logbin, par_struct = par)


  return,objshear_retval

end 





; number of randoms
function randlensing::numrand
  sample = self->sample()
  case sample of
      1: return,48
      ;; only use the first 2/20
      3: return,2
      4: begin 
          mb = obj_new('maxbcg', 'dr406', /silent)
          return,mb->numrand()
      end 
      5: begin 
          mb = obj_new('maxbcg', 'dr406', /silent)
          return,mb->numrand()
      end 
      9: return,24
      10: return,24
      11: return,24
      12: return,24
      13: return,24
      14: return,24
	  15: return,24
      else: message,"Don't know about random sample: "+ntostr(sample)
  endcase 
end 

function randlensing::randfile, dtype, randnum, createdir=createdir

  on_error, 2
  ntype = n_elements(dtype)
  if ntype eq 0 then begin 
      message,'-Syntax: files = rl->randfile(dirtype [, randnum])',/inf
      message,'       dirtype is "input" or "output"',/inf
      message,'       by default returns all files, send randnum for specific'
  endif 

  dirtype=strlowcase(dtype)
  if dirtype ne 'input' and dirtype ne 'output' then begin 
      message,'type must be "input" or "output"'
  endif 

  nrand = n_elements(randnum)
  if nrand eq 0 then begin 
      nrand = self->numrand()
      randnum = lindgen(nrand)
  endif 
  lfile = self->objshear::lensfile(dirtype, createdir=createdir)

  ;; strip the extension
  lfile = repstr(lfile, '.st', '')

  for i=0l, nrand-1 do begin 
      file = lfile + '_' + strn(randnum[i], length=2, padchar='0')+'.st'
      add_arrval, file, files
  endfor 

  return,files
end 

function randlensing::randread, dtype, randnum, hdr=hdr, columns=columns

  on_error, 2
  ntype = n_elements(dtype)
  if ntype eq 0 then begin 
      message,'-Syntax: struct = rl->randread(dirtype [, randnum])',/inf
      message,'       dirtype is "input" or "output"',/inf
      message,'       by default returns all, send randnum for specific'
  endif 

 
  files = self->randfile(dtype, randnum)

  return,self->objshear::_getfiles(files, hdr=hdr, columns=columns)

end 




function randlensing::rand_parfile, randnum, nodir=nodir
  tparfile = self->objshear::lensfile('par', /createdir, nodir=nodir)
  par_file = repstr(tparfile, $
                    'par.conf', $
                    strn(randnum,len=2,padchar='0')+'_par.conf')
  return, par_file

end 
pro randlensing::write_parfiles, mafalda=mafalda
  
  
    ;; The user can give their own source file name
    ;; to override the defaults

    p = self->par_struct()




    ms = obj_new('make_scat')
    case p.shape_correction_style of
        3: begin 
            source_file = ms->princeton_source_file(p.source_sample)
            htmrev_file = ms->princeton_htmrev_file(p.source_sample)
        end 
        4: begin 
            source_file = ms->intrinsic_source_file(p.source_sample)
            htmrev_file = ms->intrinsic_htmrev_file(p.source_sample)
        end 
        else: begin 
            source_file = ms->source_file(p.source_sample)
            htmrev_file = ms->htmrev_file(p.source_sample)
        end 
    endcase 

    obj_destroy,ms


    if p.scinv_sample ne -1 then begin
        sc = obj_new('sdss_sigma_crit', p.scinv_sample)
        scinv_file = sc->file('output', 'meanscinv', project='scinv')
        obj_destroy, sc
    endif else begin
        scinv_file='None'
    endelse

    if keyword_set(mafalda) then begin
        source_file = repstr(source_file, '/global/early2', '/home')
        htmrev_file = repstr(htmrev_file, '/global/early2', '/home')
        scinv_file  = repstr(scinv_file,  '/global/early2', '/home')
    endif
 
    nrand = self->numrand()
    for i=0l, nrand-1 do begin 

        lensin_file = self->randfile('input', i, /createdir)
        lensout_file = self->randfile('output', i, /createdir)

        par_file = self->rand_parfile(i)
        if keyword_set(mafalda) then begin
            lensin_file = repstr(lensin_file, '/global/early2', '/home')
            lensout_file =    repstr(lensout_file, '/global/early2', '/scratch')
            par_file = repstr(par_file, '.conf', '_mafalda.conf')
        endif

        if p.dopairs then begin
            pair_file = self->randfile('pairs', i, /createdir)
        endif else begin
            pair_file = 'None'
        endelse

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
        printf, lun, 'shape_correction_style '+string(p.shape_correction_style,format='(i0)')

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


    endfor 

end 


function randlensing::rand_pbsfile, randnum, nodir=nodir
  tparfile = self->objshear::lensfile('pbs', /createdir, nodir=nodir)
  par_file = repstr(tparfile, $
                    'pbs.pbs', $
                    strn(randnum,len=2,padchar='0')+'_pbs.pbs')
  return, par_file

end 
function randlensing::rand_pbs_suball_file, nodir=nodir
  file = self->objshear::lensfile('pbs', /createdir, nodir=nodir)
  file = repstr(file, 'pbs.pbs', 'pbs_suball.sh')
  return, file
end 


pro randlensing::write_pbs_suball
    file = self->rand_pbs_suball_file()
    print,'Writing to file: ',file
    openw,lun,file,/get_lun

    pbsf=self->objshear::lensfile('pbs', /nodir)
    pbsf=repstr(pbsf, '_pbs.pbs', '')

    nrand = self->numrand()
    nstr = ntostr(nrand-1)
    printf, lun
    printf, lun, 'for i in `seq -w 0 '+nstr+'`; do' 
    printf, lun, '    file='+pbsf+'_${i}_pbs.pbs'
    printf, lun, '    echo "qsub $file"'
    printf, lun, '    qsub $file'
    printf, lun, '    sleep 240'
    printf, lun, 'done'

        
    free_lun, lun
end
pro randlensing::write_pbs

    nrand = self->numrand()

    for i=0L, nrand-1 do begin

        parfile = self->rand_parfile(i)


        parfile = self->rand_parfile(i)
        pbs = self->rand_pbsfile(i)

		; to make sure directory is there
        output_file = self->randfile('output', i, /createdir)

        print,'Writing to pbs: ',pbs

        ;; The job name
        jobname = 'lrand'+ntostr(i, format='(I02)')

        openw, lun, pbs, /get
		printf,lun,"#PBS -j oe"
		; not useful since it only writes the file after the job finished
		printf,lun,"#PBS -o "+pbs+'.pbslog'
		printf,lun,"#PBS -m a"
		printf,lun,"#PBS -V"
		printf,lun,"#PBS -r n"
		printf,lun,"#PBS -l nodes=1:ppn=1"
		printf,lun,"#PBS -l walltime=48:00:00"
		printf,lun,"#PBS -N "+string(jobname)
		printf,lun,"#PBS -W umask=0022"

		printf,lun,"echo Running on `hostname`"
        printf, lun, 'objshear '+parfile+' &> '+pbs+'.out'
        free_lun, lun

    endfor


end 



; This is a script wich takes random numbers for args
pro objshear::write_script

    par_dir = self->lensdir('par', /createdir)
    ss=self->sample_string()

    script_file = 'run_'+ss+'.sh'
    script_file = concat_dir(par_dir, script_file)
    print
    print,'Writing to file: ',script_file

    openw, lun, script_file, /get_lun

    printf,lun,'#!/bin/bash'
    printf,lun
    printf,lun,'if [ $# -lt 2 ]; then'
    printf,lun,'    echo "run_'+ss+'.sh 01 02 03 ..."'
    printf,lun,'    exit 0'
    printf,lun,'fi'
    printf,lun
    printf,lun,'rnums=$*'
    printf,lun,'mac=`uname -n`'
    printf,lun,'echo "Machine: ${mac}"'
    printf,lun,'echo "Working with randoms: ${rnums}"'
    printf,lun
    printf,lun,'for rand in ${rnums}; do'
    printf,lun
    printf,lun,'    /usr/bin/time -p objshear random_'+ss+'_${rand}_par.conf'
    printf,lun,'    dt=`date`'
    printf,lun,'    echo "finished objshear random_'+ss+'_${rand}_par.conf ${dt}" | mail erin.sheldon@gmail.com -s "finished objshear random_'+ss+'_${rand}_par.conf"'
    printf,lun
    printf,lun,'done'

    free_lun, lun

end
pro objshear::write_scripts

    par_dir = self->lensdir('par', /createdir)

    nrand = self->numrand()
    for i=0l, nrand-1 do begin 
        par_file = self->rand_parfile(i, /nodir)

        script_file = 'run_'+self->sample_string()+'_'+strn(i,len=2,padchar='0')+'.sh'
        script_file = concat_dir(par_dir, script_file)

        address = 'erin.sheldon@gmail.com'

        print,'Writing to file: ',script_file

        openw, lun, script_file, /get_lun
        printf, lun
        printf, lun, '/usr/bin/time -p time objshear '+par_file

        mess = 'finished objshear '+par_file
        printf,lun,'dt=`date`'
        printf, lun, $
            'echo "'+mess+' $dt" | mail '+address+' -s "'+mess+'"'

        free_lun, lun
    endfor 

    return
end 





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setuplens
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro randlensing::setuprand, randnum=randnum
  sample = self->sample()
  case sample of
      3: self->setuplowz
      4: self->setup_maxbcg, randnum=randnum
      5: self->setup_maxbcg, randnum=randnum
      11: self->setup_maxbcg, randnum=randnum
      12: self->setup_maxbcg, randnum=randnum
      13: self->setup_maxbcg, randnum=randnum
      14: self->setup_maxbcg, randnum=randnum
      15: self->setup_maxbcg, randnum=randnum
      else: message,"Unsupported random sample: "+ntostr(sample)
  endcase 
end 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; lowz catalog
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function randlensing::lowzdir
  return, '/global/bias5/vagc-dr4/vagc0/lowz'
end 
function randlensing::lowzfile, randnum
  dir = self->lowzdir()
  file = 'lowz_random-'+ntostr(randnum)+'.dr4.fits'
  file = concat_dir(dir, file)
  return, file
END 
function randlensing::getlowz, randnum, status=status
  file=self->lowzfile(randnum)
  print
  print,'Reading file: ',file
  return, mrdfits(file, 1,status=status)
end 
pro randlensing::setuplowz

  if self->catalog() ne 'lowz' then begin 
      message,'Must be using lowz randoms'
  endif 

  nrand = self->numrand()

  outfiles = self->randfile('input', /createdir)
  print
  print,'Output files: '
  print,outfiles
  print,'-------------------------------------------------------------'
  ;; make sure dir exists
  tfile = self->randfile('output', /createdir)

  for i=0l, nrand-1 do begin 

      tm = systime(1)

      ;; Parameters
      par_struct = self->objshear::par_struct()

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; output file name(s)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      outfile = outfiles[i]

      print
      print,'Setting up file: ',outfile

      cat = self->getlowz(i)
      ncat = n_elements(cat)
      zindex = lindgen(ncat)

      nlens_init = n_elements(cat)      

      print,'Creating lens struct'
      lstruct = self->objshear::lensinput_structdef()
      lstruct = replicate(lstruct, nlens_init)
      
      print,'Copying....'
      copy_struct, cat, lstruct

      if not tag_exist(cat, 'clambda') then begin 
          eq2csurvey, cat.ra, cat.dec, clam, ceta
          lstruct.clambda = clam
          lstruct.ceta = ceta
      endif 
      
      ;; Keep track of the objects and their redshifts
      lstruct.zindex = zindex
      
      ;; Calculate some cosmology-dependent stuff and copy into struct
      self->objshear::calc_cosmo, lstruct
      
      ;; Make some generic lens cuts
      wlens = self->objshear::lenscuts(lstruct, nkeep, /doplot)
      
      ;; remove the unwanted lenses
      lstruct = lstruct[wlens]
      
      print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'
      
      print
      print,'Writing to file: ',outfile
      write_idlstruct, lstruct, outfile


      ptime,systime(1)-tm


  endfor 
  
end 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; MaxBCG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; MAKE SURE THIS IS EXACTLY THE SAME AS IN MAXBCG_LENSING!!!!
function randlensing::maxbcg_edgecuts, clambda, ceta, z, nkeep

  par = self->par_struct()
  
  ;; This will return a defined maskfile in certain situations. If undefined,
  ;; then default /basic will be used below
  self->maskfile, maskfile, /edgecuts

  rMpc = 1.0
  
  ;; convert to degrees for each lens
  maxAngle = rMpc/angdist(0.0, z, omega_m=par.omega_m)*180.0/!pi

  if n_elements(maskfile) ne 0 then begin 
      print
      message,'Using maskfile = '+maskfile,/inf
  endif 
  apply_pixel_mask, $
    clambda, ceta, masked, unmasked, $
    maxangle = maxangle, /basic, maskfile=maskfile

  if unmasked[0] eq -1 then nkeep = 0 else nkeep = n_elements(unmasked)
  
  return,unmasked

end 

pro randlensing::setup_maxbcg, randnum=randnum

	case self->catalog() of
		'maxbcg_dr406': mbcat = 'dr406'
		'gmbcg10_dr4': mbcat = 'gmbcg10'
		else: message,'bad cat: '+self->catalog()
	endcase

  mb = obj_new('maxbcg', mbcat)
  cosmo = obj_new('cosmology')

  nrand = n_elements(randnum)
  if nrand eq 0 then begin 
      nrand = mb->numrand()
      randnum = lindgen(nrand)
  endif else begin 
      nrand = n_elements(randnum)
  endelse 

  zrange = mb->zcuts()

  for ii=0l, nrand-1 do begin 

      irand = randnum[ii]

      print,'-------------------------------------------------------------'
      tm = systime(1)

      ;; Parameters
      par_struct = self->objshear::par_struct()

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; output file name(s)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      outfile = self->randfile('input', irand, /createdir)

      print
      print,'Setting up file: ',outfile

      cat = mb->random_read(irand)

      ncat = n_elements(cat)

      nlens_init = n_elements(cat)      

      print,'Creating lens struct'
      lstruct = self->objshear::lensinput_structdef()
      lstruct = replicate(lstruct, nlens_init)
      
      print,'Copying....'
      copy_struct, cat, lstruct

      csurvey2eq, lstruct.clambda, lstruct.ceta, ra, dec
      lstruct.ra = ra
      lstruct.dec = dec
      
      ;; These should be used with histogram matching
      lstruct.zindex = -1
      
      ;; Generate redshifts
      print
      print,'Generating random redshifts in range: ',zrange
      lstruct.z = cosmo->genrandz(ncat,$
                                  zrange[0],zrange[1],$
                                  omega_m=par_struct.omega_m) 


      ;; Calculate some cosmology-dependent stuff and copy into struct
      self->objshear::calc_cosmo, lstruct
      
      ;; Make some generic lens cuts
      wlens = self->objshear::lenscuts(lstruct, nkeep, /doplot)

      ;; remove the unwanted lenses
      lstruct = lstruct[wlens]
      
      print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'
      

      print
      print,'Applying maxbcg 1Mpc cut'
      keep = self->maxbcg_edgecuts(lstruct.clambda,lstruct.ceta,lstruct.z,$
                                   nkeep)
      lstruct = lstruct[keep]
      print,'Finally using '+ntostr(nkeep)+'/'+ntostr(nlens_init)


      print
      print,'Writing to file: ',outfile
      write_idlstruct, lstruct, outfile


      ptime,systime(1)-tm


  ENDFOR 

  obj_destroy,mb,cosmo

  ;; make sure dir exists
  tfile = self->randfile('output', 0, /createdir)
  
END 


PRO randlensing::maxbcg_fix

  nrand = self->numrand()
  FOR i=0L, nrand-1 DO BEGIN 

      print,'-----------------------------------------------------------'
      inf = self->randfile('input', i)
      fix_inf = inf+'_fix'

      print,'Fixing file: ',inf
      print,'Will write to fixed file: ',fix_inf
      rand = self->randread('input', i)

      csurvey2eq, rand.clambda, rand.ceta, ra, dec
      rand.ra = ra
      rand.dec = dec

      write_idlstruct, rand, fix_inf
      
  ENDFOR 

END 


PRO randlensing__define

  struct = { $
             randlensing,        $
             sample: 0,          $
             INHERITS objshear   $
           }

END 
