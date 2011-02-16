PRO run_sim, slength=slength, rfac=rfac

ns = n_elements(slength)
IF ns EQ 0 THEN slength = 120.
IF n_elements(rfac) EQ 0 THEN rfac=10.

time = systime(1)

zlens = .15
sigma = [0., 900., 1000., 1100., 1200.]
nsig = n_elements(sigma)
ntrial = 130L
gridsize = 1000.

FOR j=0, n_elements(slength)-1 DO BEGIN 
    sra = 1.05*2.0*(gridsize/2.0 + 20.*slength[j] )/3600.
    sdec = sra

    FOR i=0L, nsig-1 DO BEGIN 

        simkappa_map, sigma[i], zlens, $
                      slength=slength[j], $
                      ntrial=ntrial, $
                      /error, $
                      clust=clust, $
                      rfac=rfac, $
                      verbose=0, $
                      sra=sra, sdec=sdec, $
                      /no_ps, $
                      /write

    ENDFOR 
    ptime,systime(1)-time
ENDFOR 
ptime,systime(1)-time

return
end
