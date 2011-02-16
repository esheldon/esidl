
;;;;;;;  Need to deal with more that just r band

PRO  mockcat, run

  if N_params() eq 0 then begin
     print,'-Syntax: mockcat, run, color'
     print,''
     print,'Use doc_library,"mockcat"  for more help.'  
     return
  endif


  p = create_struct('field', 0, $
                    'camcol',0, $
                    'run',0, $
                    'id',0, $
                    'parent',0, $
                    'nchild', 0, $
                    'objc_type', 0, $
                    'objc_flags', 0, $
                    'objc_rowc', 0., $
                    'objc_colc', 0., $
                    'rowc', dblarr(5), $
                    'colc', dblarr(5), $
                    'fibercounts', dblarr(5), $
                    'fibercountserr', dblarr(5), $
                    'petrorad', dblarr(5), $
                    'petroraderr', 0., $
                    'petrocounts', 0., $
                    'petrocountserr', 0., $
                    'reddening', 0., $
                    'flags', dblarr(5), $
                    'ra', 0., $
                    'dec', 0., $
                    'ixx',    dblarr(5), $
                    'iyy',    dblarr(5), $
                    'ixy',    dblarr(5), $
                    'momerr', dblarr(5), $
                    'rho4',   dblarr(5), $
                    'whyflag',dblarr(5), $
                    'momag',  dblarr(5), $
                    'e1',     dblarr(5), $
                    'e2',     dblarr(5), $
                    'r',      dblarr(5), $)

p.e1 = -10.0
p.e2 = -10.0
p.r = -10.0

dir=run_dir(run)

FOR i=1, 1 DO BEGIN 
  tsobj_name, run, i, name
  fname = dir+name
  
  ;; do it in chunks
  fits_info, fname, /silent, n_ext=n
  
  nchunks = n/100
  left = n MOD 100

  FOR j=0, nchunks-1 DO BEGIN
    str=0
    s=0
    read_photo_col, fname, str, start=i*100+1,nframes=100
    n=n_elements(str)
    s=replicate(p, n)
    
    copy_struct, str, s
    
    qu2e, str.q, str.u, e1, e2
    s.e1 = e1
    s.e2 = e2
    
    s.rho4 = 2.0
    s.ixx = s.petrorad[2]^2/2.0
    s.iyy = s.ixx

    IF j EQ 0 THEN struct = s ELSE struct=[struct,s]
  ENDFOR 
  
  str=0
  s=0
  read_photo_col, fname, str, start=nchunks*100+1,nframes=left
  n=n_elements(str)
  s=replicate(p, n)

  copy_struct, str, s

  qu2e, str.q, str.u, e1, e2
  s.e1 = e1
  s.e2 = e2

  s.rho4 = 2.0
  s.ixx = s.petrorad[2]^2/2.0
  s.iyy = s.ixx

ENDFOR

return
END















