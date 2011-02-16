pro qcorr,str,binsize,sm,lag,corr_func,step,smq,rows,doplots=doplots,range=range

  if n_params() eq 0 then begin
      print,'-Syntax: qcorr,str,binsize,sm,lag,corr_func,step,smq,rows,doplots=doplots,'
      print,'               range=range   STR SHOULD BE ONLY STARS'
      print,'Rebins and smooths f.  Calls a_correlate.  lag should be in'
      print,'range [-(rbin-2),(rbin-2)]' 
      return
  endif	

  print,'----------------------------------------'

  rbin = n_elements(str)/binsize ;will round off
  print,'New array size: ',rbin
  newn = rbin*binsize
  print,'newn: ',newn

  IF n_elements(lag) EQ 0 THEN BEGIN
      print,'lag should have '+ntostr(rbin)+' elements and'+$
            'be in range [-'+ntostr(rbin-2)+', '+ntostr(rbin-2)+']'
      return
  ENDIF 

  step = (max(str.objc_rowc) - min(str.objc_rowc))/rbin
  print,'Step is: ',step

  s=sort(str.objc_rowc)
  f=str[s]

  q=f[0:newn-1].q[2]
  rows=f[0:newn-1].objc_rowc

  q=rebin(q,rbin)
  rows=rebin(rows,rbin)

  smq=smooth(q,sm)
  mn=mean(smq)
  corr_func = a_correlate(smq,lag)

  print,'----------------------------------------'
  if keyword_set(doplots) then begin

      print,'Plotting smoothed q'
      !p.multi=[0,1,2]
      if keyword_set(range) then begin
          print,'Applying user defined range'
      endif else begin
          range=[min(rows),max(rows)]
      endelse
      
      plot,rows,smq,psym=3,xrange=range,ytitle='q(2)',xtitle='objc_rowc',$
        title='Smoothed q'
      
      print,'Plotting auto correlation vs. lag in frames'
      print,'----------------------------------------'
      xpos=lag*step/1361.0
      plot,xpos,corr_func,psym=1,$
        title='Autocorrelation Function', xtitle='Lag(frames)       1frame=1361 rows'
      oplot,[0,max(xpos)],[0,0]
      !p.multi=[0,1,1]
      
  endif
  
  return
end




