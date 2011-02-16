PRO readngc, struct

  dir='/sdss4/data1/esheldon/ngc/'
  file=dir+'ngc2000.dat'

  openr, lun, file, /get_lun

  nl = numlines(file)

  ;; not using l_size

  s=create_struct('name', '', $
                  'catalog', '',    $
                  'number', 0L,     $
                  'type', '',      $
                  'rah', 0,         $
                  'ram', 0.,        $
                  'decd', 0,        $
                  'decm', 0,        $
                  'ra', 0d,         $
                  'dec',0d,         $
                  'source', '',     $
                  'const', '',      $
                  'size', 0.,       $
                  'mag', 0.,        $
                  'n_mag', '',      $
                  'description','' )
  struct = replicate(s, nl)
  string=''

  beg =     [0, 1, 6, 10, 13, 19, 23, 26, 29, 33, 40, 44, 46]
  lengths = [1, 5, 3,  2,  4,  3,  2,  1,  3,  5,  4,  1, 50]
  FOR i=0L, nl-1 DO BEGIN 
      readf, lun, string, format='(A96)'

      nn=0
      catalog = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      IF catalog EQ 'I' THEN catalog = 'IC' ELSE catalog = 'NGC'
      number = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      type = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      rah = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      ram = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      decd = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      decm = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      source = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      const = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      size = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      mag = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      n_mag = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1
      description = strmid(string,beg[nn], lengths[nn] ) & nn=nn+1

      struct[i].name = catalog+ntostr(number)
      struct[i].catalog = catalog
      struct[i].number = fix(ntostr(number))
      struct[i].type = ntostr(type)
      struct[i].rah = fix(ntostr(rah))
      struct[i].ram = float(ntostr(ram))
      struct[i].decd = fix(ntostr(decd))
      struct[i].decm = fix(ntostr(decm))
      struct[i].ra = (struct[i].rah + struct[i].ram/60.)*15.
      struct[i].dec = ten( [struct[i].decd, struct[i].decm] )
      struct[i].source = ntostr(source)
      struct[i].const = ntostr(const)

      ;; sometimes size is not listed
      IF ntostr(size) EQ '' THEN size = -1. ELSE size = float(ntostr(size))
      struct[i].size = size
      ;; sometimes mag is not listed
      IF ntostr(mag) EQ '' THEN mag=0. ELSE mag = float(ntostr(mag))
      struct[i].mag = mag
      struct[i].n_mag = ntostr(n_mag)
      struct[i].description = ntostr(description)

;      IF i EQ 1 THEN BEGIN
;          help,struct[i],/str
;          return
;      ENDIF 
  ENDFOR 

  mwrfits, struct, dir+'ngc2000.fit',/create

  return
END 
