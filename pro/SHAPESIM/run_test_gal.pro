pro run_test_gal,galtype,smear,gal,tot

if (n_params() eq 0) then begin
	print,'syntax: run_test_gal,galtype,smear,gal,tot'
	return
endif

path = '/sdss4/data1/esheldon/'
numfr = 500
theta = 0.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aratio = 1.0
newtest_gal,numfr,galtype,aratio,theta,smear,galaxy,gal,tot

mwrfits,gal,path+'mids2nexpgal1.fits'
mwrfits,tot,path+'mids2nexptot1.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aratio = .95
newtest_gal,numfr,galtype,aratio,theta,smear,galaxy,gal,tot

mwrfits,gal,path+'mids2nexpgal95.fits'
mwrfits,tot,path+'mids2nexptot95.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aratio = .9
newtest_gal,numfr,galtype,aratio,theta,smear,galaxy,gal,tot

mwrfits,gal,path+'mids2nexpgal9.fits'
mwrfits,tot,path+'mids2nexptot9.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aratio = .8
newtest_gal,numfr,galtype,aratio,theta,smear,galaxy,gal,tot

mwrfits,gal,path+'mids2nexpgal8.fits'
mwrfits,tot,path+'mids2nexptot8.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aratio = .6
newtest_gal,numfr,galtype,aratio,theta,smear,galaxy,gal,tot

mwrfits,gal,path+'mids2nexpgal6.fits'
mwrfits,tot,path+'mids2nexptot6.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aratio = .4
newtest_gal,numfr,galtype,aratio,theta,smear,galaxy,gal,tot

mwrfits,gal,path+'mids2nexpgal4.fits'
mwrfits,tot,path+'mids2nexptot4.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aratio = .3
newtest_gal,numfr,galtype,aratio,theta,smear,galaxy,gal,tot

mwrfits,gal,path+'mids2nexpgal3.fits'
mwrfits,tot,path+'mids2nexptot3.fits'

return
end




