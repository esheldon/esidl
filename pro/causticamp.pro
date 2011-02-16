pro causticamp, thetacrit, maxplot=maxplot,noplot=noplot

  if N_params() eq 0 then begin
	print,'Syntax: causticamp, thetacrit, maxplot=maxplot,noplot=noplot'
	return
  endif

testing

R=6.96e5    ;radius of sun in km
Dos=1256.0  ;distance to source (star) in Mpc

psixxx = 1.0/thetacrit

;;;;;;;;;;;;   This should be changed to fit needs     ;;;;;;;;;;;;
psiyy = 0.5
psixxx = 1.0/thetacrit

n=40.0
d1 = R*findgen(n)/n
d2= d1 + R

fac=1.2e7      ;factor so we can use Mpc, km and ratio in radians while
               ;using thetacrit in arcseconds

mag1 =  fac*sqrt(2.0/psixxx)/(1.0-psiyy)*sqrt(d1)/R
mag2 =  fac*sqrt(2.0/psixxx)/(1.0-psiyy)*[ sqrt(d2) - sqrt(d2-R) ]/R
print,'Max: ',strtrim(string(max(mag1)),2)

d=[d1,d2]/R
mag = [mag1,mag2]

title = 'Thetacrit = '+strtrim(string(thetacrit),2)
xtitle='d/R'
ytitle='Magnification.  Solar object'
if (not keyword_set(noplot) ) then $
		plot,d,mag,xtitle=xtitle,ytitle=ytitle,title=title

if keyword_set(maxplot) then begin
   old=!p.multi
   !p.multi=[0,0,1]
   maxcrit=30.0 ;arcseconds
   thetacrit=(findgen(n) + .2)
   thetacrit=(maxcrit/max(thetacrit))*thetacrit
   psixxx = 1.0/thetacrit
   max = fac*sqrt(2.0/psixxx)/(1.0-psiyy)*sqrt(d1[n-1])/R
   ytitle='Maximum Magnification.  Solar object.'
   xtitle='Thetacrit (arcseconds)'
   if (not keyword_set(noplot)) then begin
      print,'Hit any key'
      s=get_kbrd(20)
   endif
   plot,thetacrit,max,xtitle=xtitle,ytitle=ytitle
   !p.multi=old
endif






return
end

