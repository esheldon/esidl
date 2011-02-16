pro run_smooth, run

;if n_params() eq 0 then begin
;  print,'usage: run_smooth,run'
;  return
;endif

  dirone='/sdss3/usrdevel/philf/'
  
  CASE run OF 
    752: rr = 'run'+ntostr(run)+'/'
    756: rr = 'run'+ntostr(run)+'/'
    ELSE : BEGIN
      print,'Cant do run',run
      return
    END 
  ENDCASE 

  dir=dirone + rr

  print,'-----------------------'
  print,' Working on '+rr
  print,'-----------------------'
  FOR i=3,6,1 DO BEGIN 

    istr=ntostr(i)
    print,'-----------------------'
    print,' Processing column: ',istr
    print,'-----------------------'
    fname = dir + 'adat' + istr + 'c.fit'
    oname = dir + 'adat'+ istr + 'c'
    smooth_col,fname,oname, groupn=100

  ENDFOR 


  return
end
