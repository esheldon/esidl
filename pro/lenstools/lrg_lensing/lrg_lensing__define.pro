FUNCTION lrg_lensing::init, sample

  IF n_elements(sample) EQ 0 THEN BEGIN 
      message,'You must send the sample number',/inf
      message,"ll = obj_new('lrg_lensing', sample)",/inf
      message,''
  ENDIF 

  rmin = 20.0
  logbin = 1

  sample = fix(sample)

  type = 'lrg'
  par = $
    { $
      maskfile: '', $
      sample: sample, $
      source_sample: 1, $
      random_sample: -1, $
      catalog: '',$
      zbuffer: 0.0 $
    }

  CASE sample OF
      1: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton re-gaussianization
          par.catalog = 'dr4plus'    ; lrg catalog
      END 
      2: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton re-gaussianization
          par.catalog = 'dr4plus'    ; lrg catalog
          par.zbuffer = 0.1     ; a buffer behind lens
      END
      3: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton re-gaussianization
          par.catalog = 'dr4plus'    ; lrg catalog
          par.zbuffer = 0.2     ; a buffer behind lens
      END
      4: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 2 ; princeton analytic
          par.catalog = 'dr4plus'    ; lrg catalog
      END 

      ;; Added source samples
      5: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 2 ; more well resolved r > 0.585
      END 
      6: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 3 ; less well resolved r <= 0.585
      END 
      7: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 4 ; seeing < 1.5 
      END 

      ;; Introduced a minimum angle of 20 arcsec just to see what happens
      ;; into the C code.  Will make this a parameter if it makes a difference
      8: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 1 ; normal
      END 
      ;; 10 arcsec min
      9: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 1 ; normal
      END 

      ;; Not blended
      ;; no min radius
      10: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 5 ; not blended
      END 

      ;; Rachel's R > 1/3 cut
      11: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 6 ; more well resolved r > 1/3
      END 

      ;; the top 1/3 in r comes at r=2/3, just a coinc.
      12: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1   ; photozs as truth
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr4plus'    ; lrg catalog
          par.source_sample = 7 ; more well resolved r > 2/3
      END 

      ELSE: message,'Unknown sample: '+ntostr(sample)
  ENDCASE 

  self.sample = sample

  print,'type = '                 ',type
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





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Where things are
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; plots

FUNCTION lrg_lensing::plotdir

  plotdir = '~/plots/lrg/'

  IF NOT keyword_set(base) THEN BEGIN 
      plotdir = plotdir + self->objshear::sample_string()+'/'
  ENDIF 

  IF NOT file_test(plotdir, /directory) THEN BEGIN 
      file_mkdir, plotdir
  ENDIF 
  return,plotdir
END 

;; lens inputs

FUNCTION lrg_lensing::lensinput_dir, random=random

  dir = esheldon_config('lensinput_dir')+'lrg_lensinput/'
  dir = dir + self->objshear::sample_string()+'/'

  IF keyword_set(random) THEN dir=dir+'random/'

  ;; Will make parent dirs too
  IF NOT file_test(dir, /directory) THEN file_mkdir, dir
  return,dir

END 

FUNCTION lrg_lensing::lensinput_file, $
                    randnum=randnum, allrand=allrand

  nrand = n_elements(randnum)
  IF nrand GT 0 OR keyword_set(allrand) THEN random=1 ELSE random=0
  dir = self->lensinput_dir(random=random)
  
  IF keyword_set(allrand) THEN BEGIN 
      files = file_search(dir + 'lrg_random_*.st')
      return,files
  ENDIF 

  ext = self->objshear::sample_string() + '.st'

  name = 'lrg'

  IF random THEN BEGIN 
      randstr = self->randstring(randnum)
      name = name + '_random_' + randstr
  ENDIF 

  name = name + '_'+ext


  name = dir + name
  return,name

END 


;; Lens outputs
FUNCTION lrg_lensing::lensoutput_dir, random=random, $
                    lum_nbin=lum_nbin, $
                    z_nbin=z_nbin

  dir = esheldon_config('lensout_dir')+'lrg/'
  dir = dir + self->objshear::sample_string()+'/'

  IF keyword_set(random) THEN dir = dir + 'random/'


  IF keyword_set(lum_nbin) THEN BEGIN 
      dir = dir + 'lum'+strn(lum_nbin, len=2, padchar='0') + '/'
  ENDIF ELSE IF keyword_set(z_nbin) THEN BEGIN 
      dir = dir + 'zbin'+strn(z_nbin, len=2, padchar='0') + '/'
  ENDIF 

  ;; Will make parent dirs too
  IF NOT file_test(dir, /directory) THEN file_mkdir, dir
  return,dir

END 

FUNCTION lrg_lensing::lensoutput_file, nfiles, $
                    randnum=randnum, $
                    z_nbin=z_nbin, z_bin=z_bin, $
                    lum_nbin=lum_nbin, lum_bin=lum_bin

  ;; This doesn't work for ngals subs because we only keep
  ;; the combined, but the *output* is used by ::combined_file

  random=(n_elements(randnum) GT 0)

  files = self->lensinput_file(randnum=randnum)
  dirsep, temporary(files), dirs, files

  n_lum_nbin = n_elements(lum_nbin)
  n_lum_bin = n_elements(lum_bin)

  n_z_nbin = n_elements(z_nbin)
  n_z_bin = n_elements(z_bin)

  lum_both = (n_lum_nbin GT 0) + (n_lum_bin GT 0)
  IF (lum_both NE 0) AND (lum_both NE 2) THEN BEGIN 
      message,'You must send both lum_nbin and lum_bin'
  ENDIF  

  z_both = (n_z_nbin GT 0) + (n_z_bin GT 0)
  IF (z_both NE 0) AND (z_both NE 2) THEN BEGIN 
      message,'You must send both z_nbin and z_bin'
  ENDIF  

  IF n_lum_nbin NE 0 THEN BEGIN 
      lumbinstr = self->lum_namestring(lum_nbin, lum_bin) + '_'
  ENDIF 

  IF n_z_nbin NE 0 THEN BEGIN 
      IF n_lum_nbin GT 0 AND n_z_bin GT 1 THEN BEGIN 
          message,'Only one z allowed when doing lum bins'
      ENDIF 
      zstr = self->z_namestring(z_nbin, z_bin) + '_'
  ENDIF ELSE BEGIN 
      zstr = ''
  ENDELSE 


  IF n_lum_nbin NE 0 THEN BEGIN 
      
      numbin = n_elements(lum_bin)
      FOR i=0L, numbin-1 DO BEGIN 

          binstr = strn(lum_bin[i], len=2, padchar='0')
          front = zstr + lumbinstr 

          ttfiles = front + files

          add_arrval, ttfiles, tfiles

      ENDFOR 

      files = tfiles

      dirs = self->lensoutput_dir(random=random, lum_nbin=lum_nbin)

  ENDIF ELSE BEGIN 
      dirs = self->lensoutput_dir(random=random, z_nbin=z_nbin)
      files = zstr + files
  ENDELSE 

  nfiles = n_elements(files)
  files = dirs + files
  return,files

END 


FUNCTION lrg_lensing::jackknife_dir, lum_nbin=lum_nbin

  dir = self->lensoutput_dir(lum_nbin=lum_nbin)+'jackknife/'
  IF NOT file_test(dir, /directory) THEN file_mkdir, dir
  return,dir
  
END 

FUNCTION lrg_lensing::jackknife_file,$
  lum_nbin=lum_nbin, lum_bin=lum_bin

  dir = self->jackknife_dir(lum_nbin=lum_nbin)

  tfile = self->combined_file(lum_nbin=lum_nbin, lum_bin=lum_bin)
  
  dirsep, tfile, tdir, file
  file = repstr(file, '_combined', '_jackknife')

  return,dir + file

END 







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Random points
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION lrg_lensing::randname, randnum
  IF n_elements(randnum) EQ 0 THEN message,'n=mb->randname(randnum)'
  dir = self->lensinput_dir(/random)
  return,dir + 'random'+self->objshear::randstring(randnum)+'.st'
END 
FUNCTION lrg_lensing::numrand
  return,1
END 
FUNCTION lrg_lensing::randnum
  numrand = self->numrand()
  return,lindgen(numrand)
END 









;; lens sample dependent
PRO lrg_lensing::setuplens, random=random
  
  tm = systime(1)

  ;; Parameters
  par_struct = self->objshear::par_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file name(s)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(random) THEN randnum=0
  outfile = self->lensinput_file(randnum=randnum)

  print
  print,'Setting up file: ',outfile

  ;; get the data.  Note: zindex will be the same for all samples
  ;; created from this lens catalog
  cat = self->get()
  ncat = n_elements(cat)
  zindex = lindgen(ncat)

  nrand = n_elements(randnum)
  IF nrand EQ 0 THEN BEGIN 

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
      
      ;; Make some generic lens cuts
      wlens = self->objshear::lenscuts(lstruct, nkeep)

      ;; remove the unwanted lenses
      lstruct = lstruct[wlens]

      print,'Kept '+ntostr(nkeep)+'/'+ntostr(nlens_init)+' from lenscuts'

      print
      print,'Writing to file: ',outfile
      write_idlstruct, lstruct, outfile

  ENDIF ELSE BEGIN 

      catz = cat.z
      ncatz = n_elements(cat)

      ;; clear memory
      cat = 0

      FOR i=0L, nrand-1 DO BEGIN 

          tmi = systime(1)

          print
          rand = self->get(random=random)

          ;; remove randoms from regions not covered by this sample.
          completeness_cut = 0.20
          print,'Making completeness cut at fgot = '+ntostr(completeness_cut)
          w = where(rand.fgot GT completeness_cut)
          rand = rand[w]

          outfile = self->lensinput_file(randnum=randnum[i])
          print,'*********************************************************'
          print,'Output file: ',outFile

          print,'Creating lensum struct'
          lstruct = self->objshear::lensinput_structdef()

          nlens_init = n_elements(rand)
          lstruct = replicate(lstruct, nlens_init)

          print,'Copying....'
          lstruct.ra = rand.ra
          lstruct.dec = rand.dec

          eq2csurvey, rand.ra, rand.dec, clam, ceta
          lstruct.clambda = clam
          lstruct.ceta = ceta

          delvarx, rand

          print
          print,'Assigning redshifts'

          cat_index = long( ncatz*randomu(seed, nlens_init) )

          lstruct.z = catz[cat_index]
          lstruct.zindex = zindex[cat_index]

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

          IF nrand GT 1 THEN BEGIN 
              print,'Time for this random'
              ptime,systime(1)-tmi
          ENDIF 

      ENDFOR 
      
  ENDELSE 

  ptime,systime(1)-tm

END 


FUNCTION lrg_lensing::lum_namestring, nbin, bin
  str = self->z_namestring(nbin, bin)
  str = repstr(str, 'zbin', 'lumbin')
  return,str
END 


PRO lrg_lensing::z_bins, nbin, zlowlim, zhighlim

  CASE nbin OF 
      2: BEGIN 
          zlowlim = [0.0, 0.3]
          zhighlim = [0.3, 1000.0]
      END 
      ELSE: message,'unsupported nbins in z: '+ntostr(nbin)
  ENDCASE 

END 

FUNCTION lrg_lensing::z_namestring, nbin, bin

  num_bin = n_elements(bin)
  IF num_bin GT 1 THEN str = strarr(num_bin) ELSE str = ''

  FOR i=0L, num_bin-1 DO BEGIN 
      str[i] = 'zbin'+strn(nbin[0], len=2, padchar='0')+$
        '_'+strn(bin[i], len=2, padchar='0')
  ENDFOR 
  return,str

END 

FUNCTION lrg_lensing::z_where_string, nbin

  self->z_bins, nbin, zlowlim, zhighlim

  IF nbin EQ 1 THEN where_string = '' $
  ELSE where_string = strarr(nbin)
  
  FOR i=0L, nbin-1 DO BEGIN 

      where_string[i] = $
        'struct.z GE '+ntostr(zlowlim[i])+' AND '+$
        'struct.z LT '+ntostr(zhighlim[i])

  ENDFOR 

  return,where_string

END 








FUNCTION lrg_lensing::where_string, type, nbin

  CASE type OF
      'z': return, self->z_where_string(nbin)
      ELSE: message,'Unknown select type: '+ntostr(type)
  ENDCASE 

END 

FUNCTION lrg_lensing::z_select, struct, nbin, nkeep

  where_string = self->z_where_string(nbin)
  keep = self->objshear::struct_select(struct, where_string, nkeep)
  return,keep

END 


PRO lrg_lensing::sub, type, nbin, randnum=randnum

  average_tags = ['z',$
                  'petroflux[0]', 'petroflux[1]', 'petroflux[2]',$
                  'petroflux[3]', 'petroflux[4]']

  self->objshear::sub_sample, type, nbin, randnum=randnum, $
    average_tags=average_tags

END 



PRO lrg_lensing__define

  struct = { $
             lrg_lensing, $
             sample: 0, $
             INHERITS objshear, $
             INHERITS lrg_morad $
           }

END 
