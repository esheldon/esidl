PRO get_scat_old, run1, run2, clr, scat, indir=indir, typecut=typecut

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: get_scat_old, run1, run2, clr, scat, indir=indir, typecut=typecut'
      return
  ENDIF 

  colors = ['u','g','r','i','z']
  r1str = ntostr(run1)
  r2str = ntostr(run2)
  
  IF NOT keyword_set(typecut) THEN tc = '' ELSE tc = '_typecut'

  IF n_elements(indir) EQ 0 THEN BEGIN
      
      camcol=1
      fetch_rerun, run1, rerun
      fetch_dir, run1, camcol, rerun, corrdir=corrdir1
      base = ( str_sep(corrdir1, 'calibChunks/1/') )[0]
      corrdir1 = base + 'combined/'

      fetch_rerun, run2, rerun
      fetch_dir, run2, camcol, rerun, corrdir=corrdir2
      base = ( str_sep(corrdir2, 'calibChunks/1/') )[0]
      corrdir2 = base + 'combined/'
  ENDIF ELSE BEGIN
      corrdir1=indir
      corrdir2=indir
  ENDELSE 

  sname=corrdir1+'run'+r1str+'_'+r2str+'_srcgal_'+colors[clr]+tc+'_overlap.fit'
      
  IF NOT exist(sname) THEN BEGIN 
      sname=corrdir2+'run'+r1str+'_'+r2str+'_srcgal_'+colors[clr]+tc+'_overlap.fit'
      IF NOT exist(sname) THEN BEGIN 
          sname=corrdir1+'run'+r2str+'_'+r1str+'_srcgal_'+colors[clr]+tc+'_overlap.fit'
          IF NOT exist(sname) THEN BEGIN 
              sname=corrdir2+'run'+r2str+'_'+r1str+'_srcgal_'+colors[clr]+tc+'_overlap.fit'
              IF NOT exist(sname) THEN BEGIN
                  print,'No overlap file found for runs '+r1str+'/'+r2str
                  return
              ENDIF 
          ENDIF
      ENDIF 
  ENDIF 

  scat = mrdfits(sname, 1)

  return
END 
