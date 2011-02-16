PRO bincheck

;;; use fwhm = 1.5 arcseconds => sigma = 1.61 pixels in sloan
;;;                                      2*1.61= in finer grid
;;; Use aratio=0.816 and theta=0.0 for e1=.2

;;; go out 3*sigma on either side => 2*3*sigma wide = 19 pixels wide

fwhm = 1.5
sigma = .43*fwhm ;; arcseconds
sigma = sigma/.2  ;; pixels in finer grid
size = [fix(2*3*sigma)+1,fix(2*3*sigma)+1]
print,'using size = ',fix(2*3*sigma+1)

n=8

pold = !p.multi
!p.multi = [0,4,4]
FOR i=0,n-1 DO BEGIN
  
  cen = [(size[0]-1.)/2.0 + (randomu(seed)-.5),$
         (size[1]-1.)/2.0 + (randomu(seed)-.5)]
  makegauss,g,size,sigma,aratio=.816,cen=cen
  tvim2,g

  gb=binarr(g,2)
  tvim2,gb

ENDFOR

x=!d.x_size/2.0 - !d.x_size/4.0
y=!d.y_size/2.0

xyouts, x, y ,/device, $
  'FWHM = 1.5 arcseconds   e1=.2   e2=0.0   Before and After binning'

return
end






