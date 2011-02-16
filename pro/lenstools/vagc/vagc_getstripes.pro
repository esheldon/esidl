PRO vagc_getstripes, stripes, nstripe, lss=lss, letter=letter, post=post, $
                     sample=sample

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: vagc_getstripes, stripes, nstripe, lss=lss, letter=letter, post=post, sample=sample'
      return
  ENDIF 

  read_runlist, runstruct, /silent
  
  catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)

  byrun_dir = sdssidl_config('spec_dir') + 'blanton/gal_collated/byrun/'
  files = findfile(byrun_dir, count=nf)
  w=where( strmatch(files,'*'+catname+'*'), nmatch)
  IF nmatch EQ 0 THEN BEGIN 
      message,'No '+catname+' files found'
  ENDIF 

  nf = nmatch
  files = files[w]
  runs = lonarr(nf)

  FOR i=0L, nf-1 DO BEGIN 
      t=( strsplit(files[i],'run',/extract) )[0]
      runstr = ( strsplit(t,'-',/extract) )[0]
      runs[i] = long(runstr)
  ENDFOR 

  match, runs, runstruct.run, mruns, mstruct

  rmd = rem_dup(runstruct[mstruct].stripe)
  stripes = runstruct[mstruct[rmd]].stripe
  w=where(stripes LE 86 AND $
          stripes NE 61 AND $
          stripes NE 62, nstripe)
  stripes = stripes[w]

END 
