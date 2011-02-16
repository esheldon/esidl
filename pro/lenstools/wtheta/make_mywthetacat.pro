PRO make_mywthetacat_combine, outstruct, runs, strips, stripe

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

  print
  print,'# in overlap: ',ngood
  print

  s=sort(lambda[wgood])
  wgood = wgood[s]

  outstruct = temporary(outstruct[wgood])

END 

PRO make_mywthetacat_selectgal, struct, clr, minmag, maxmag, good, ngood

  galtype = 3

  s = where(struct.objc_type EQ galtype, ngood)
  IF ngood EQ 0 THEN return

  make_flag_struct,fs
  make_status_struct,ss
                                ;THESE CUTS FOR ALL OBJECTS
  fs.SATUR= 'N'
  fs.BRIGHT = 'N'
  fs.EDGE = 'N'
  colorindex = 2
  
  flag_select, struct[s], fs, colorindex, s2, objc=1
  IF s2[0] EQ -1 THEN BEGIN
      ngood = 0
      return
  ENDIF 
  s = s[s2]

  ;; This gets all objects that aren't duplicates
  ss.SECONDARY = 'Y'
  status_select,struct[s],ss,s3
  IF s3[0] EQ -1 THEN BEGIN 
      ngood = 0
      return
  ENDIF 
  s=s[s3]
  
  petfix = struct[s].petrocounts[clr] - struct[s].reddening[clr]
  s4 = where(petfix GT minmag[clr] AND petfix LT maxmag[clr], ngood)
  IF ngood NE 0 THEN good = s[s4]

  return

END 

PRO make_mywthetacat, stripe, clr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: make_mywthetacat, stripe, clr'
      return
  ENDIF 

  sdss_goodstripe, stripe, runs, reruns

  maxmag = [21.0, 21.0, 21.0, 21.0, 21.0]
  maxmagstr = ['21','21','21','21','21']
  minmag = [15.5, 15.5, 15.5, 15.5, 15.5]

  colmin=1
  colmax=6
  ncol = (colmax-colmin+1)
;runs=runs[0:1]
;  print,runs
  nrun = n_elements(runs)
  strips = strarr(nrun)

  outdir = '/sdss6/data0/wtheta/'
  bigoutfile = outdir+'wtheta-stripe'+ntostr(stripe)+$
        '-gal-allmag-'+!colors[clr]+maxmagstr[clr]+'.fit'

  tags = ['run','rerun','camcol','field','id',$
          'lambda','eta','petrocounts','counts_model',$
          'reddening','objc_flags','status', 'objc_type']

  IF nrun NE n_elements( rem_dup(runs) ) THEN message,'We have more than one rerun for some of these runs!!'

  keepst = create_struct(name='keepst_struct',$
                         'run',0L,$
                         'lambda',0d,$
                         'eta',0d,$
                         'petrocounts', fltarr(5),$
                         'gr',0.0)
;GOTO,jump
  FOR irun=0L, nrun-1 DO BEGIN 

      run = runs[irun]
      rerun = reruns[irun]
      
      ii=0L
      FOR camcol = colmin, colmax DO BEGIN 
          
          fetch_dir, run, camcol, rerun, dir
          fetch_file_list, dir, files, fnums

          nfile = n_elements(files)

          IF camcol EQ colmin THEN BEGIN
              head=headfits(files[0])
              strips[irun] = ntostr( sxpar(head,'strip') )
              print
              print,'Run: ',ntostr(runs[irun]),'  strip: ',strips[irun]
              print
              totfiles = ncol*nfile
              ptrlist = ptrarr(totfiles)
              numlist = lonarr(totfiles)
          ENDIF 

          delvarx, tsobjstr
          FOR fi = 0L, nfile-1 DO BEGIN 
              
              read_tsobj, dir, tmp, start=fnums[fi], tsobjstr=tsobjstr, verbose=0, $
                taglist=tags
              
              make_mywthetacat_selectgal, tmp, clr, minmag, maxmag, good, ngood

              IF ngood NE 0 THEN BEGIN 
                  
                  keep = replicate(keepst, ngood)
                  
                  keep.lambda = tmp[good].lambda
                  keep.eta = tmp[good].eta
                  keep.petrocounts = tmp[good].petrocounts-$
                    tmp[good].reddening
                  model = tmp[good].counts_model - $
                    tmp[good].reddening
                  keep.gr = model[1]-model[2]
                  keep.run = run

                  ptrlist[ii] = ptr_new(keep, /no_copy)
                  numlist[ii] = ngood
                  
              ENDIF 
              
              ii = ii+1
              
              IF ((fi MOD 20) EQ 0) OR (fi EQ nfile-1) THEN print,'Run: '+ntostr(run)+'  Camcol: '+ntostr(camcol)+'  Field: '+ntostr(fnums[fi])+'/'+ntostr(fnums[nfile-1]) 
          ENDFOR 
      ENDFOR 
      ntotal = long(total(numlist))
      comb_ptrstruct, keepst, ptrlist, numlist, outstruct

      outfile = outdir+'wtheta-'+ntostr(runs[irun])+$
        '-gal-allmag-'+!colors[clr]+maxmagstr[clr]+'.fit'

      print
      print,'Outputting file: ',outfile
      print,'Number of objects: '+ntostr(ntotal)
      print
      mwrfits2, outstruct, outfile, /create, /destroy

  ENDFOR 
jump:
  print,'Combining all the files'
  ptrlist = ptrarr(nrun)
  numlist = lonarr(nrun)

  FOR irun = 0L, nrun-1 DO BEGIN 
          infile = outdir+'wtheta-'+ntostr(runs[irun])+$
            '-gal-allmag-'+!colors[clr]+maxmagstr[clr]+'.fit'
          tmp = mrdfits(infile, 1, structyp='keepst_struct')
          numlist[irun] = n_elements(tmp)
          ptrlist[irun] = ptr_new(tmp, /no_copy)

  ENDFOR 
  comb_ptrstruct, keepst, ptrlist, numlist, outstruct

  make_mywthetacat_combine, outstruct, runs, strips, stripe

  ;; cut off ends
  print,'Outputting to file ',bigoutfile
  mwrfits2, outstruct, bigoutfile, /create, /destroy

END 
