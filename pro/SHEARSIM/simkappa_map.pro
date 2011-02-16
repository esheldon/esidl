PRO simkappa_map, sigma, zlens, zsource, $
                  nstack=nstack, $
                  slength=slength, $
                  rfac=rfac, $
                  gridsize=gridsize, $
                  stepfac=stefpac, $
                  verbose=verbose, $
                  write=write, $
                  error=error,$
                  ntrial=ntrial, $
                  random=random, $
                  rotate=rotate, $
                  surface=surface, $
                  noprompt=noprompt, $
                  clust=clust, $
                  scat=scat, $
                  no_ps=no_ps, $
                  sra=sra, sdec=sdec, $
                  nonoise=nonoise, $
                  _extra=extra

IF n_params() LT 1 THEN BEGIN 
    print,'-Syntax: simkappa_map, sigma, zlens, zsource, '
    print,'               nstack=nstack, '
    print,'               slength=slength, '
    print,'               rfac=rfac, '
    print,'               gridsize=gridsize, '
    print,'               stepfac=stefpac, '
    print,'               verbose=verbose, '
    print,'               write=write, '
    print,'               error=error,'
    print,'               ntrial=ntrial, '
    print,'               random=random, '
    print,'               rotate=rotate, '
    print,'               surface=surface, '
    print,'               noprompt=noprompt, '
    print,'               clust=clust, '
    print,'               scat=scat, '
    print,'               no_ps=no_ps, '
    print,'               nonoise=nonoise, '
    print,'               _extra=extra'
    return
ENDIF 
time=systime(1)

IF NOT keyword_set(nstack) THEN nstack=1.
density = 7200.                 ;per square degree
density = density*nstack

IF NOT keyword_set(nonoise) THEN nonoise=0
IF NOT keyword_set(write) THEN write=0
IF NOT keyword_set(random) THEN random=0

IF n_elements(zlens) EQ 0 THEN zlens = .15
IF n_elements(zsource) EQ 0 THEN zsource = .4
IF n_elements(ntrial) EQ 0 THEN ntrial=1

;dummy variables
run1=1
run2=1
clr=1

IF n_elements(sra) EQ 0 THEN sra = 2.7 ;degrees
IF n_elements(sdec) EQ 0 THEN sdec = 2.7 ;degrees


;; use the structure from annis
IF n_elements(clust) EQ 0 THEN BEGIN 
    rdannis, tmp
    clust = tmp[0]
ENDIF 
clust.z = zlens

; Create the source catalog.

outdir='/sdss4/data1/esheldon/CLUSTER/SIM/'

IF random THEN BEGIN 

    print,'Sigma = ',sigma
    sis_shear, zlens, zsource, 0., sdec, sra, tmp, cen, kappa, $
               density=density, error=error

    sss = sort(tmp.ra)

    gridsize = 2000.
    nx = long(gridsize/step) + 1
    ny=nx

    clust.dec = cen[0]
    clust.ra  = cen[1]
    clust.name = 'sim_rand_sig'+ntostr(long(sigma))
    kappa_map, run1, run2, clr, clust, $
      slength=slength, $
      rfac=rfac, $
      gridsize=gridsize, $
      stepfac=stepfac, $
      scat=tmp[sss], $
      write=write, $
      /abs, $
      rotate=rotate, $
      surface=surface, $
      noprompt=noprompt, $
      verbose=verbose, $
      outdir=outdir, $
      no_ps=no_ps, $
      _extra=extra
      
ENDIF ELSE nse = .053           ;From actual data was .0455


IF ntrial NE 0 THEN BEGIN 
    FOR i=0, ntrial-1 DO BEGIN 

        print,'Sigma = ',sigma
        IF ntrial GT 1 THEN print,'Trial :',ntostr(i+1),'/',ntostr(ntrial)
        
        sis_shear, zlens, zsource, sigma, sdec, sra, scat, cen, $
                   density=density, error=error

        ss = sort(scat.ra)

        clust.dec = cen[0]
        clust.ra  = cen[1]
        clust.name = 'sim_sig'+ntostr(long(sigma))

        kappa_map, run1, run2, clr, clust, $
          slength=slength, $
          rfac=rfac, $
          gridsize=gridsize, $
          stepfac=stepfac, $
          scat=scat, $
          write=write, $
          /abs, $
          rotate=rotate, $
          surface=surface, $
          noprompt=noprompt, $
          verbose=verbose, $
          no_ps=no_ps, $
          outdir=outdir, $
          _extra=extra

        IF nonoise THEN BEGIN
    
            print,'Now doing no noise'
            sis_shear, zlens, zsource, sigma, sdec, sra, scat, cen, shear, $
              density = 25000.

            clust.name = 'sim_nonoise'

            s = sort(shear.ra)
            shear.e1 = 2*shear.e1 ;convert shear to polarization
            shear.e2 = 2*shear.e2
            kappa_map, run1, run2, clr, clust, $
              slength=slength, $
              rfac=rfac, $
              gridsize=gridsize, $
              stepfac=stepfac, $
              scat=shear, $
              write=write, $
              /abs, $
              rotate=rotate, $
              surface=surface, $
              noprompt=noprompt, $
              verbose=verbose, $
              no_ps=no_ps, $
              outdir=outdir, $
              _extra=extra


        ENDIF 
    ENDFOR 
ENDIF 

print
ptime, systime(1)-time

return
END 
