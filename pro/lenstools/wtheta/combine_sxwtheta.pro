PRO combine_sxwtheta_combine, outstruct, runs, strips, stripe

  IF stripe GT 45 THEN issouth=1 ELSE issouth=0

  lambda=outstruct.lambda
  
  IF issouth THEN rotate_lambda, lambda

  nrun=n_elements(runs)

  FOR irun=0L, nrun-1 DO BEGIN 

      wr=where(outstruct.run EQ runs[irun])
      minl = min(lambda[wr])
      maxl = max(lambda[wr])

      add_arrval, minl, tminl
      add_arrval, maxl, tmaxl
      IF strips[irun] EQ 'N' THEN BEGIN 
          add_arrval, minl, n_minlam
          add_arrval, maxl, n_maxlam
      ENDIF ELSE BEGIN 
          add_arrval, minl, s_minlam
          add_arrval, maxl, s_maxlam
      ENDELSE 

  ENDFOR 

  print,'Run   Strip   min-lambda  max-lambda'
  forprint,runs, ' '+strips+' ',tminl,tmaxl

  n_minl = min(n_minlam)
  n_maxl = max(n_maxlam)
  
  s_minl = min(s_minlam)
  s_maxl = max(s_maxlam)
  
  minl = max([n_minl, s_minl])
  maxl = min([n_maxl, s_maxl])

  print
  print,'Min lambda in overlap: ',minl
  print,'Max lambda in overlap: ',maxl
  print

  wgood = where(lambda LT maxl AND lambda GT minl, ngood)

  s=sort(lambda[wgood])
  wgood = wgood[s]

  outstruct = temporary(outstruct[wgood])

  print,'Removing duplicates'
  rmclose_lameta, outstruct, good, bad
  outstruct = temporary(outstruct[good])
  ngood=n_elements(outstruct)

  print
  print,'# in overlap: ',ngood
  print

END 

PRO combine_sxwtheta_selectgal, struct, clr, minmag, maxmag, good, ngood
  
  good = -1
  pet = struct.petrocounts[clr]
  s1 = where(pet GT minmag[clr] AND pet LT maxmag[clr], ngood)
  IF ngood NE 0 THEN good = s1

  return

END 

PRO combine_sxwtheta, stripe, clr, indir=indir, outdir=outdir, front=front, $
                      combine_secondary=combine_secondary

  ;; Take output files from split_sxwtheta_byrerun and combine into stripe
  ;; can also include secondary objects with /combine_secondary. Must have
  ;; gotten them from sx!

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: combine_sxwtheta, stripe, clr, indir=indir, outdir=outdir, front=front'
      return
  ENDIF 

  maxmag = [21.0, 21.0, 21.0, 21.0, 21.0]
  maxmagstr = ['21','21','21','21','21']
  minmag = [15.5, 15.5, 15.5, 15.5, 15.5]

  ;;;;;;;; Files
  IF n_elements(indir) EQ 0 THEN indir = '/sdss6/data0/wtheta/sxoutput/'
  IF n_elements(front) EQ 0 THEN front = 'getPrimaryWthetaGals'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss6/data0/wtheta/'
  bigoutfile = outdir+'wtheta-stripe'+stripe2string(stripe)+$
        '-gal-allmag3-'+!colors[clr]+maxmagstr[clr]+'.fit'

  ;;;;;;;; Get stripe/run info
  sdss_goodstripe, stripe, runs, reruns, stripes, strips, $
                  psfieldreruns=psfieldreruns

  rmd=rem_dup(runs)
  runs = runs[rmd]
  strips = strips[rmd]

  nrun = n_elements(runs)

  ;; this may not be the reruns in the directory, find them
  ;; and pick latest
  reruns = lonarr(nrun)
  FOR irun=0L, nrun-1 DO BEGIN 
      
      lookf=indir+front+'-'+run2string(runs[irun])+'*'
      print,'Searching for run '+run2string(runs[irun])+' files'
      spawn, 'ls '+lookf, answer

      IF answer[0] EQ '' THEN add_arrval, irun, remove_ids $
      ELSE BEGIN 
          nans = n_elements(answer)

          rrs = -1
          
          FOR ians=0L, nans-1 DO BEGIN 
              thdr=headfits(answer[ians])
              
              tr = sxpar(thdr, 'RERUN')
              IF tr GT rrs THEN rrs=tr
              
          ENDFOR 
          
          reruns[irun] = rrs
      ENDELSE 
  ENDFOR 

  nremove = n_elements(remove_ids)
  IF nremove EQ nrun THEN message,'No files for stripe: '+ntostr(stripe)
  IF nremove NE 0 THEN remove, remove_ids, runs, reruns, strips

  print
  print,'Some Files found. Using these reruns: '
  colprint,runs,reruns,'   '+strips

  colmin=1
  colmax=6
  ncol = (colmax-colmin+1)
;runs=runs[0:1]
;  print,runs
  nrun = n_elements(runs)

  keepst = create_struct(name='keepst_struct',$
                         'run',0L,$
                         'strip','',$
                         'lambda',0d,$
                         'eta',0d,$
                         'petrocounts', fltarr(5),$
                         'gr',0.0)
;GOTO,jump

  ptrlist = ptrarr(nrun)
  numlist = lonarr(nrun)

  FOR irun=0L, nrun-1 DO BEGIN 

      run = runs[irun]
      rerun = reruns[irun]
      strip = strips[irun]

      infile = indir + front + '-'+run2string(run)+'-'+ntostr(rerun)+'.fit'
      print,'Reading Primary File: ',infile
      tmp = mrdfits(infile,1)

      ;; look for secondary object file
      IF keyword_set(combine_secondary) THEN BEGIN 
          sec_infile = repstr(infile, 'Primary','Secondary')
          IF fexist(sec_infile) THEN BEGIN 
              print,'Combining secondary file: ',sec_infile
              tmp2=mrdfits(sec_infile, 1)
              concat_structs, temporary(tmp), temporary(tmp2), tmp
          ENDIF 
      ENDIF 

      combine_sxwtheta_selectgal, tmp, clr, minmag, maxmag, good, ngood

      keep = replicate(keepst, ngood)

      keep.lambda = tmp[good].lambda
      keep.eta = tmp[good].eta
      keep.petrocounts = tmp[good].petrocounts
      keep.gr = tmp[good].gr
      keep.run = run
      keep.strip = strip

      ptrlist[irun] = ptr_new(keep, /no_copy)
      numlist[irun] = ngood
      
  ENDFOR 
jump:

  print
  print,'Combining all the runs'
  comb_ptrstruct, keepst, ptrlist, numlist, outstruct

  combine_sxwtheta_combine, outstruct, runs, strips, stripe

  ;; cut off ends
  print,'Outputting to file ',bigoutfile
  mwrfits2, outstruct, bigoutfile, /create, /destroy

END 
