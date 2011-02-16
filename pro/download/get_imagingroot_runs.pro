PRO get_imagingroot_runs

  ;; parse the dump from the imagingRoot page
  indir = sdssidl_config('data_dir')+'/dbm/imagingRoot/'
  infile = indir + 'imagingRoot.dat'
  outfile = indir + 'imagingRootRuns2Get.dat'

  format = 'L,L,L,L,L,L,L,L,L,L,L'

  readcol, infile, $
    runs, reruns, $
    nfpC, nfpBIN, nfpM, nfpObjc, nfpFieldStat, $
    nfpAtlas, npsField, ntsField, ntsObj, $
    format=format, /silent

  nruns = n_elements(runs)

  minrerun = 40
  maxrerun = 49

  ;; Now get the ones we want
  w = where(reruns GE minrerun AND reruns LE 50 AND $
            ntsObj GT 0, nw)
  
  w2 = where(reruns GE 40 AND reruns LT 50 AND $
             nfpAtlas GT 0 AND npsField GT 0 AND $
             ntsField GT 0 AND ntsObj GT 0, nw2)



  print,'Minrerun: '+ntostr(minrerun)+'  Maxrerun: '+ntostr(maxrerun)
  print,' ntsObj > 0: '+ntostr(nw)
  print,' ntsObj,nfpAtlas,npsField,ntsField,ntsObj > 0: '+ntostr(nw2)

  IF nw NE nw2 THEN message,'Differ, you should address this'

  hh = histogram(runs[w], rev=rev)
  nhh = n_elements(hh)
  
  keep = bytarr(nruns)

  FOR i=0L, nhh-1 DO BEGIN 

      IF rev[i] NE rev[i+1] THEN BEGIN 
          wrun = rev[ rev[i]:rev[i+1]-1 ]
          wrun = w[wrun]

          nwrun = n_elements(wrun)
          IF nwrun GT 1 THEN BEGIN 
              print,'Run: '+ntostr(runs[wrun[0]])
              print,'  Reruns: ',ntostr(reruns[wrun])

              wkeep = where( reruns[wrun] EQ max(reruns[wrun]) )
              wkeep = wrun[wkeep]
             
              keep[wkeep] = 1
              print,'  Keeping: '+ntostr(reruns[wkeep])
          ENDIF ELSE BEGIN 
              keep[wrun] = 1
          ENDELSE 
      ENDIF 

  ENDFOR 

  
  wkeep = where(keep)
  print
  print,'Keeping total of '+ntostr(n_elements(wkeep))+' reruns'

  print
  print,'Writing to file: ',outfile
  header = string('run','rerun','#fpC','#fpBIN','#fpM','#fpObjc','#fpFieldStat','#fpAtlas','#psField','#tsField','#tsObj', format='(11A15)')
  len=strlen(header)
  line=""
  FOR i=0L, len-1 DO line=line+"-"

  openw, lun, outfile, /get_lun
;  lun=-1
  printf, lun, header
  printf, lun, line
  colprint,runs,reruns,nfpC, nfpBIN, nfpM, nfpObjc, nfpFieldStat, $
    nfpAtlas, npsField, ntsField, ntsObj, format='(11I15)', lun=lun, $
    wuse=w
  free_lun, lun

END 
