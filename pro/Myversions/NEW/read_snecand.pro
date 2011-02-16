PRO read_snecand_getstruct, nCand, struct

  struct = create_struct('cand','',$
                         'run',   0L, $
                         'rerun', 0L, $
                         'camcol',0,  $
                         'field', 0,  $
                         'id',    0L, $
                         'plate', 0L, $
                         'fiber', 0,  $
                         'mjd',   0L, $
                         'ra',    0d, $
                         'dec',   0d)

  struct = replicate(struct, nCand)

return
  struct = create_struct('cand','',$
                         'location','',$
                         'template','',$
                         'mjd',0L,$
                         'plate',0L,$
                         'fiber',0,$
                         'snage',0L,$
                         $
                         'snMinAge', 0L, $
                         'snMaxAge', 0L, $
                         'AbsMagSN', fltarr(nAbsMagSN), $
                         'AbsMagTmp', fltarr(nAbsMagTmp), $
                         $
                         'snType','',$

                         'ra',0d,$
                         'dec',0d,$
                         'fibermag',fltarr(nfib),$
                         'synthmag',fltarr(nsyn),$
                         'z',0.0,$
                         'zconf',0.0,$
                         'eclass',0.0,$
                         'galaxyfit',fltarr(ngfit),$
                         'galaxychisq',0.0,$
                         'galaxysnefit',fltarr(ngsnefit),$
                         'galsnechisq',0.0,$
                         'rank',0)

END 

PRO read_snecand_parseline, keyword, val, struct, icand

  COMMON read_snecand_block, $
    nfib, nsyn, ngfit, ngsnefit, $
    nAbsMagSN, nAbsMagTmp

  CASE strlowcase(keyword) OF
      'cand':BEGIN 
          icand = icand + 1
          struct[icand].cand = val
      END 
      'run':      struct[icand].run      = long(val)
      'rerun':    struct[icand].rerun    = long(val)
      'camcol':   struct[icand].camcol   = fix(val)
      'field':    struct[icand].field    = fix(val)
      'objid':    struct[icand].id       = long(val)
      'ra':       struct[icand].ra       = double(val)
      'dec':      struct[icand].dec      = double(val)
      'plate':    struct[icand].plate    = long(val)
      'fiber':    struct[icand].fiber    = fix(val)
      'mjd':      struct[icand].mjd      = long(val)
      ELSE: 
  ENDCASE 

return
  ;; old way where we tried to read everything
  CASE strlowcase(keyword) OF
      'cand':BEGIN 
          icand = icand + 1
          struct[icand].cand = val
      END 
      'location': struct[icand].location = val
      'template': struct[icand].template = val
      $
      'run':      struct[icand].run      = long(val)
      'rerun':    struct[icand].rerun    = long(val)
      'camcol':   struct[icand].camcol   = fix(val)
      'field':    struct[icand].field    = fix(val)
      'objid':    struct[icand].id       = long(val)
      'ra':       struct[icand].ra       = double(val)
      'dec':      struct[icand].dec      = double(val)
      'plate':    struct[icand].plate    = long(val)
      'fiber':    struct[icand].fiber    = fix(val)
      'mjd':      struct[icand].mjd      = long(val)
      $
      'sntype':   struct[icand].snType   = val
      'snage':    struct[icand].SNAge    = long(val)
      'snminage': struct[icand].snMinAge = long(val)
      'snmaxage': struct[icand].snMaxAge = long(val)
      'absmagsn': BEGIN 
          el = strsplit(strtrim(val,2), /extract)
          IF n_elements(el) NE nAbsMagSN THEN $
            message,'AbsMagSN must have '+ntostr(nAbsMagSN)+' elements'
          struct[icand].AbsMagSN = float(el)
      END 
      'absmagtmp': BEGIN 
          el = strsplit(strtrim(val,2), /extract)
          IF n_elements(el) NE nAbsMagTmp THEN $
            message,'AbsMagTmp must have '+ntostr(nAbsMagTmp)+' elements'
          struct[icand].AbsMagTmp = float(el)
      END 
      'fibermag': BEGIN 
          el = strsplit(strtrim(val,2), /extract)
          IF n_elements(el) NE nfib THEN $
            message,'fiber mag must have '+ntostr(nfib)+' elements'
          struct[icand].fibermag = float(el)
      END 
      'synthmag': BEGIN 
          el = strsplit(strtrim(val,2), /extract)
          IF n_elements(el) NE nsyn THEN $
            message,'synth mag must have '+ntostr(nsyn)+' elements'
          struct[icand].synthmag = float(el)
      END 
      'zvalue': struct[icand].z = float(val)
      'zconf': struct[icand].zconf = float(val)
      'eclass': struct[icand].eclass = float(val)
      'galaxyfit': BEGIN 
          el = strsplit(strtrim(val,2), /extract)
          IF n_elements(el) NE ngfit THEN $
            message,'galaxy fit must have '+ntostr(ngfit)+' elements'
          struct[icand].galaxyfit = float(el)
      END 
      'galaxychisq': struct[icand].galaxychisq = float(val)
      'galaxysnefit': BEGIN 
          el = strsplit(strtrim(val,2), /extract)
          IF n_elements(el) NE ngsnefit THEN $
            message,'galaxy+sne fit must have '+ntostr(ngsnefit)+' elements'
          struct[icand].galaxysnefit = float(el)
      END
      'galsnechisq': struct[icand].galsnechisq = float(val)
      'rank': struct[icand].rank = long(val)
      ELSE: ;;message,'Unknown keyword: '+ntostr(keyword)
  ENDCASE 

END 

PRO read_snecand, struct, snefile=file

  COMMON read_snecand_block, $
    nfib, nsyn, ngfit, ngsnefit, $
    nAbsMagSN, nAbsMagTmp

  ;; dir = '~garyk/sne/newcode/'
  ;;dir = "/net/cheops1/home/products/sncode/"
  dir = "/net/cheops1/home/www/html/SN/"
  IF n_elements(file) EQ 0 THEN BEGIN 
      file = dir + 'SNeCand.txt'
  ENDIF 

  spawn,'grep CAND '+file+' | wc -l',ans
  ncand = long(ans[0])

  nfib = 3
  nsyn = 3
  ngfit = 6
  ngsnefit = 7

  nAbsMagSN = 3
  nAbsMagTmp = 3

  read_snecand_getstruct, nCand, struct
 
  openr, lun, file, /get_lun

  line=' '
  icand = -1L
  WHILE NOT eof(lun) DO BEGIN 

      readf, lun, line
      ;print,line

      el = strsplit(line, '=', /extract)
      nel = n_elements(el)

      IF nel GE 2 THEN BEGIN 

          read_snecand_parseline, el[0], el[1], struct, icand          

      ENDIF 

  ENDWHILE 

  free_lun, lun

END 
