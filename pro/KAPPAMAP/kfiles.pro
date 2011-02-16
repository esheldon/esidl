PRO kfiles, dir, file, image

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: kfiles, dir, file, image'
      return
  ENDIF 

  IF NOT keyword_set(write) THEN write=0
  fmt='A10, I, A8, A6, I, A3, I, A5, A49'

  readcol, dir+file, a1,a2,a3,a4,a5,a6,a7,a8, names,format=fmt,/silent
  names = dir+names
    
  nf = n_elements(names)

  t =  mrdfits(names[0],/silent)

  ;; This assumes that they stack in y-direction
  s=size(t)
  sx = s[1]
  sy = s[2]

  i=0L
  image = fltarr(sx, nf*sy)
  image[*, i*sy:(i+1)*sy-1] = t
  FOR i=1L, nf-1 DO BEGIN
      
      t = mrdfits(names[i], /silent)
      
      image[*, i*sy:(i+1)*sy-1] = t
  ENDFOR 


return
END 
