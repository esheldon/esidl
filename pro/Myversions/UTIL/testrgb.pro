PRO testrgb,g,r,i,bmap,gmap,rmap,jpegfile=jpegfile

type = 5

IF type EQ 1 THEN BEGIN 
    fg='/sdss3/usrdevel/esheldon/JPG/n2903_pg.fits'
    fr='/sdss3/usrdevel/esheldon/JPG/n2903_pr.fits'
    fi='/sdss3/usrdevel/esheldon/JPG/n2903_pi.fits'
    gunn=1
    ;jpegfile = 'n2903_2.jpeg'
    ;giffile = 'n2903_2.gif'
ENDIF 

IF type EQ 2 THEN BEGIN 
    fg='/sdss3/usrdevel/esheldon/JPG/n3031_pg.fits'
    fr='/sdss3/usrdevel/esheldon/JPG/n3031_pr.fits'
    fi='/sdss3/usrdevel/esheldon/JPG/n3031_pi.fits'
    gunn=1
    low_cut=2.
    ;sat=1.9
    ;sat=0.
    ;jpegfile = 'n3031_2.jpeg'
    ;giffile = 'n3031_2.gif'
ENDIF 

IF type EQ 3 THEN BEGIN 
    fg='/sdss3/usrdevel/esheldon/JPG/n4559_pg.fits'
    fr='/sdss3/usrdevel/esheldon/JPG/n4559_pr.fits'
    fi='/sdss3/usrdevel/esheldon/JPG/n4559_pi.fits'
    gunn=1
    ;jpegfile = 'n4559_1.jpeg'
    ;giffile = 'n4559_1.gif'
ENDIF 

IF type EQ 4 THEN BEGIN 
    fg='/sdss3/usrdevel/esheldon/JPG/n4303_pg.fits'
    fr='/sdss3/usrdevel/esheldon/JPG/n4303_pr.fits'
    fi='/sdss3/usrdevel/esheldon/JPG/n4303_pi.fits'
    gunn=1
    ;low_cut=2.
    ;sat = .5
    jpegfile = 'n4303_1.jpeg'
    ;giffile = 'n4303_1.gif'
ENDIF 

IF type EQ 5 THEN BEGIN 
    fg='/sdss3/usrdevel/esheldon/JPG/NGC3521_pg.fits'
    fr='/sdss3/usrdevel/esheldon/JPG/NGC3521_pr.fits'
    fi='/sdss3/usrdevel/esheldon/JPG/NGC3521_pi.fits'
    sdss=1
    ;jpegfile = 'NGC3521_2.jpeg'
    ;giffile = 'NGC3521_2.gif'
ENDIF 

g=mrdfits(fg,0,ghdr)
r=mrdfits(fr,0,rhdr)
i=mrdfits(fi,0,ihdr)

gsky = sxpar(ghdr,'SKY')
rsky = sxpar(rhdr,'SKY')
isky = sxpar(ihdr,'SKY')

gsig = sxpar(ghdr,'SKYSIG')
rsig = sxpar(rhdr,'SKYSIG')
isig = sxpar(ihdr,'SKYSIG')


;rgbview2, i, r, g, /sqrt, /log
;return

;rgbview4, i, r, g, jpegfile=jpegfile
;return

rgbview4, i, r, g, rsig=isig, gsig=rsig, bsig=gsig, $
  rsky=isky, gsky=rsky, bsky=gsky, $
  low_cut=low_cut, contrast=50., sat=sat, gamma=gamma, sdss=sdss, gunn=gunn, $
  jpegfile=jpegfile, giffile=giffile, $
  rmap=rmap,gmap=gmap,bmap=bmap

return
END 
