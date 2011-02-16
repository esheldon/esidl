pro make_aligntags,taglist, default=default

  if n_params() eq 0 then begin
    print,'-syntax make_lenstags, taglist'
    return
  endif
  
  taglist=['ID','RUN','RERUN','CAMCOL','FIELD',$
           'PARENT','NCHILD',$
           'OBJC_TYPE','TYPE',$
           'FLAGS','FLAGS2',$
           'OBJC_FLAGS','OBJC_FLAGS2',$
           'OBJC_ROWC','OBJC_COLC',$
           'ROWC','COLC',$
           'FIBERCOUNTS','FIBERCOUNTSERR',$
           'PETROCOUNTS','PETROCOUNTSERR',$
           'PSFCOUNTS','PSFCOUNTSERR',$
           'PETRORAD','PETRORADERR',$
           'REDDENING',$
           'RA','DEC', $
           'LAMBDA','ETA','ROTATION',$
	   'PLATE','FIBERID','OBJID']

  IF NOT keyword_set(default) THEN BEGIN
      addtags = [$
                  'SEEING','STARFLAG', $
                  'IXX','IYY','IXY','RHO4',$
                  'WHYFLAG',$
                  'MOMMAG',$
                  'E1','E2','MOMERR',$
                  'R',$
                  'PHOTOZ']
      taglist = [taglist, addtags]
  ENDIF 

  return
end
