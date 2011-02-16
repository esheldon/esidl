pro plot1,x,y,sig,a,b,xtitle,ytitle,xrange,yrange,title,$
          siga=siga,sigb=sigb

  myusersym, 'fill_circle'
  psym=8

  xx=[-1,1]

  srt=sort(x)

  fitlin,x,y,sig,a,siga,b,sigb,/silent
  aploterror,1,x[srt],y[srt],sig[srt],xtitle=xtitle,ytitle=ytitle,yrange=yrange,xrange=xrange,title=title,psym=psym
  yy=a+b*xx
  oplot,xx,yy
  sa=strmid(strtrim(string(a,format='(e8.1)'),2),0,8)
  ssiga=strmid(strtrim(string(abs(siga/a)),2),0,5)
  sb=strmid(strtrim(string(b,format='(e8.1)'),2),0,8)
  ssigb=strmid(strtrim(string(abs(sigb/b)),2),0,5)
  stta='a='+sa+' '+ssiga
  sttb='b='+sb+' '+ssigb
;  xyouts,xrange[0]+0.01,yrange[0]+0.8*(yrange[1]-yrange[0]),stt,font=-1
  legend, [stta,sttb], /right, /clear
;  print,'yaya',a,b,strtrim(string(a),2),strtrim(string(b,format='(e8.1)'),2)
  
return

end

