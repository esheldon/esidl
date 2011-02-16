PRO scat_create_dbdfile, chunk, scat=scat

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: create_scat_dbdfile, chunk, scat=scat'
      return
  ENDIF 
  IF chunk LT 0 OR chunk GT 3 THEN BEGIN 
      print,'Chunks must be in [1,3]'
      return
  ENDIF 

  CASE chunk OF
      1: maxentries = '6000000'
      2: maxentries = '7000000'
      3: maxentries = '4000000'
  ENDCASE 
  clrstr = clrarr2string(clr)

  zdbase = getenv('ZDBASE') + '/'

  ;; Read in a source catalog to get tags
  clr = [1,2,3]
  IF n_elements(scat) EQ 0 THEN get_scat,15,clr,scat,/hirata
  
  scat = scat[0]

  db_struct2def, scat, names, typedefs

  file = zdbase + 'scat'+ntostr(chunk)+'.dbd'
  openw, lun, file, /get_lun
;  lun=-1

  printf,lun,'#title'
  printf,lun,'Source catalog database for '+clrstr
  printf,lun
  printf,lun,'#maxentries'
  printf,lun, maxentries
  printf,lun
  printf,lun,'#items'
  
  ;; First some extra items
  printf,lun,'photoid                U*8     SDSSid'
  printf,lun,'stripe                 I*2     SDSS stripe'
  printf,lun,'photoz_use             B*1     Passes photoz cuts?'
  maxlen = max(strlen(names))
  FOR i=0L, n_elements(names)-1 DO BEGIN 

      tn = names[i]
      tnl = strlen(tn)
      IF tnl LT maxlen THEN BEGIN 
          tn = tn + strjoin(replicate(' ', maxlen-tnl))
      ENDIF 
      printf,lun,tn,typedefs[i],format='(a,a7)'

  ENDFOR 
  printf,lun
  printf,lun,'#formats'
  printf,lun
  printf,lun,'#index'
  printf,lun,'photoid    sort'
  ;printf,lun,'spec_entry index'
  printf,lun,'leafid     sort'
  printf,lun,'stripe     sort'

  free_lun, lun

END 
