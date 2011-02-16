pro build_system, rad, theta, sigcrit, $
                  e1, e2, uncert, $
                  rmin, rmax, binsize, nbins, $
                  ata, btb, av, bv, npair, rsum, lens_npair, lens_rad

;sets up the normal equations to solve for the shear profile

if n_params() LT 18 then begin
	print,'-syntax build_system, rad, theta, sigcrit, '
        print,'  e1, e2, uncert, '
        print,'  rmin, rmax, binsize, nbins, '
        print,'  ata, btb, av, bv, npair, rsum'
	return
endif

costh=-1.0*cos(2.0*theta)
sinth=-1.0*sin(2.0*theta)

weight=2.0/uncert               ;convert uncert(e) to uncert(gamma)

a=fltarr(nbins)
b=a

hist=histogram(rad, binsize=binsize, min=rmin, $
               max=rmax,rever=rev_ind)
numbin=n_elements(hist)
whist = where(hist NE 0, nhist)

;; Check if there are any in this annulus rmin-rmax
IF nhist NE 0 THEN BEGIN 
         
    FOR i=0L, nhist-1 DO BEGIN 
        binnum = whist[i]    
        w=rev_ind( rev_ind(binnum):rev_ind(binnum+1)-1 )

        lens_npair[w, binnum] = 1
        lens_rad[w,binnum] = rad[w]

        np=n_elements(w)
        npair[binnum]=npair[binnum]+np
        rsum[binnum]=rsum[binnum]+total(rad[w])
        a(binnum)=total(costh[w]/sigcrit[w])*weight
        b(binnum)=total(sinth[w]/sigcrit[w])*weight
        av(binnum)=av(binnum)+total(costh[w])*e1*(weight^2)*.5
        bv(binnum)=bv(binnum)+total(sinth[w])*e2*(weight^2)*.5
    ENDFOR 
    FOR i=0L, nhist-1 DO BEGIN 
        bini = whist[i] 
        FOR j=0L, nhist-1 DO BEGIN 
            binj = whist[j]
            ata(binj,bini)=ata(binj,bini)+a(bini)*a(binj)
            btb(binj,bini)=btb(binj,bini)+b(bini)*b(binj)
        ENDFOR
    ENDFOR     
    
ENDIF 

return
END 

