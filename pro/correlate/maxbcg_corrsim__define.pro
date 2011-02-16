; I separated this out from maxbcg_correlate__define just 
; because there were so many specific routines.  But you can
; read this sample there as well as sample 5
function maxbcg_corrsim::init, sample

    if n_elements(sample) eq 0 then begin 
        message,'You must send the sample number',/inf
        message,"ml = obj_new('maxbcg_corrsim', sample)"
    endif 

    par = correlate_pardef()
    case sample of
        1: begin 
            ;; see correlate_pardef for other defaults
            ;; Sample 5 in maxbcg_correlate
            maxbcg_catalog = 'hv1'
            par.psample = 'p03'   
            par.DPsample = 'hv1bcg'
            par.RPsample = 'rhv1bcg'
            par.DSsample = 'hv1gals'
            par.RSsample = 'rhv1gals'

            par.omega_m = 0.3
        end 
        2: begin 
            ; we only need to re-do the measurements for DD and DR
            ; Sample 6 in maxbcg_correlate
            maxbcg_catalog = 'hv1halo'
            par.psample = 'p03'   
            par.DPsample = 'hv1halo'
            ; We can re-use all these since we re-sample the redshifts
            ; anyway and the secondaries are the same.
            par.RPsample = 'rhv1bcg'
            par.DSsample = 'hv1gals'
            par.RSsample = 'rhv1gals'

            par.omega_m = 0.3
        end

      else: message,'Unsupported sample number: '+ntostr(sample)
  endcase 

  cinit = self->correlate::init(par)
  mbinit = self->maxbcg::init(maxbcg_catalog)
  return, (cinit+mbinit) eq 2

end

function maxbcg_corrsim::zrange
  return, [0.1, 0.3]
end 

function maxbcg_corrsim::radec_cuts, ra, dec, nw
    w=where(ra ge 0d and ra le 90d and dec ge 0d and dec le 90d, nw)
    return, w
end

; Ben's crazy 1Mpc cut plus a redshift cut
function maxbcg_corrsim::cuts, ra, dec, z, nkeep, comp=comp

    par = self->par_struct()
    n = n_elements(ra)
    ind = lindgen(n)

    ; Redshift cuts
    zrange = self->zrange()
    print
    print,'applying redshift cuts: '+strjoin(ntostr(zrange),', ')
    wz = where(z[ind] ge zrange[0] and z[ind] le zrange[1], nz)

    if nz eq 0 then message,'None passed redshift cuts'
    print,'  kept: '+ntostr(nz)
    ind = ind[wz]

    print
    print,'Getting max angle'
    rmpc = 1.0
    max_angle_rad = rmpc/angdist(0.0, z[ind], omega_m=par.omega_m)
    max_angle_deg = max_angle_rad*180.0/!pi
 
    ; A cut 1Mpc above the dec=0 line
    print
    print,'applying 1Mpc dec=0 cuts'
   
    w0 = where( (dec[ind] - 0.0) gt max_angle_deg, n0)
    if n0 eq 0 then message,'None passed dec=0 edge cut'
    print,'  kept: '+ntostr(n0)
    ind = ind[w0]

    ; Kludgy edge cuts. First split into 2.5 degree stripes
    print
    print,'applying 1Mpc ra cuts'
    print,'  binning...', format='($,a)'
    bs=binner(dec[ind], min=0.0, max=90.0, binsize=2.5, rev=rev)

    ; Now for each stripe do the 1Mpc cut
    print,'checking...'
    nbin = n_elements(bs.hist)
    keepbool = replicate(1, n_elements(ind))
    for i=0L, nbin-1 do begin
        if rev[i] ne rev[i+1] then begin
            w=rev[ rev[i]:rev[i+1]-1 ]
            w=ind[w]
         
            ; 1Mpc from edge at ra=0
            decmean = mean(dec[w])   
            mygcirc, 0d, decmean, ra[w], dec[w], disrad0, /radians_out
            mygcirc, 90d, decmean, ra[w], dec[w], disradmax, /radians_out
            wkeep = where( disrad0 gt max_angle_rad[w] $
                            and disradmax gt max_angle_rad[w], $
                            nkeep, comp=comp, ncomp=ncomp)

            if ncomp ne 0 then begin
                keepbool[w[comp]] = 0
            endif
        endif
    endfor
    keep = where( keepbool eq 1, nkeep)
    if nkeep eq 0 then message,'No objects passed crazy edge cuts'
    print,'  kept: '+ntostr(nkeep)
    ind = ind[keep]
    if arg_present(comp) then begin
        comp = lindgen(n_elements(ra))
        remove, ind, comp
    endif
    return, ind

end


function maxbcg_corrsim::nperfile
    return, 1000000
end
pro maxbcg_corrsim::genrand_radec, num, rra, rdec, nperfile=nperfile

    ; For the sims this is just a quadrant
    ;   0 < ra < 90
    ;   0 < dec < 90

    minra = 0d
    maxra = 90d
    mindec = 0d
    maxdec = 90d

    ; generate uniformly in sin(dec)
    sinfdec = sin( minra*!d2r )
    sinldec = sin( maxra*!d2r )

    rsindec = arrscl( randomu(seed, num, /double), $
        sinfdec, sinldec, arrmin=0d, arrmax=1d )
    rdec = asin(rsindec)*!r2d
    rra = arrscl( randomu(seed, num, /double), $
        minra, maxra, arrmin=0d, arrmax=1d )

end
function maxbcg_corrsim::genrand, type

    par = self->par_struct()
    nrand = self->nperfile()
    if type eq 'primary' then begin

        self->genrand_radec, nrand, rra, rdec

        print
        print,'Generating redshifts'
        cosmo = obj_new('cosmology')

        zrange = self->zrange()
        zrand = cosmo->genrandz(nrand, $
            zrange[0], zrange[1], $
            omega_m=par.omega_m)
        obj_destroy, cosmo

        print
        print,'Making edge cut'
        rkeep = self->cuts(rra, rdec, zrand, nrkeep)

        print,'Kept: ',nrkeep
        struct = self->primary_struct(nrkeep)

        rra = rra[rkeep]
        rdec = rdec[rkeep]
        zrand = zrand[rkeep]
        bcg_id = lindgen(nrkeep)

        struct.bcg_id = bcg_id
        struct.ra = rra
        struct.dec = rdec
        struct.z = zrand

            
    endif else begin

        self->genrand_radec, nrand, rra, rdec

        struct = self->secondary_struct(nrand, /random)
        struct.ra = rra
        struct.dec = rdec
        struct.htm_index = htm_index(rra, rdec, par.depth)

    endelse

    return, struct

end


function maxbcg_corrsim::getz, struct
    if tag_exist(struct, 'photoz_cts') then return, struct.photoz_cts
    if tag_exist(struct, 'z') then return, struct.z
    message,'no z tag found'
end
pro maxbcg_corrsim::primary_input, randnum=randnum

    nrand=n_elements(randnum)
    if nrand eq 0 then begin

        bsim=self->get()
        z = self->getz(bsim)

        ntot=n_elements(bsim)
        ngals_min = 3
        print
        print,'Selecting ngals >= '+ntostr(ngals_min)
        keep=where(bsim.ngals ge ngals_min, nkeep)

        print
        print,'Kept '+ntostr(nkeep)+'/'+ntostr(ntot)

        keep2 = self->cuts($
            bsim[keep].ra, bsim[keep].dec, z[keep], nkeep)
        print,'Kept '+ntostr(nkeep)+'/'+ntostr(ntot)
        keep = keep[keep2]

        outst = self->primary_struct(nkeep)
        outst.bcg_id = keep
        outst.ra = bsim[keep].ra
        outst.dec = bsim[keep].dec
        outst.z = z[keep]

        outfile = self->corrfile('primary', 'input', /createdir)
        print
        print,'Writing to file: ',outfile
        write_idlstruct, outst, outfile

    endif else begin

        for i=0l, nrand-1 do begin

            print,'------------------------------------------------------------'
            rstruct = self->genrand('primary')
            routfile = $
                self->corrfile('primary_random','input',$
                                primary_randnum=randnum[i],/createdir)
            print
            print,'Writing random file: ',routfile
            write_idlstruct, rstruct, routfile

        endfor 

    endelse
end

pro maxbcg_corrsim::secondary_input, randnum=randnum

    nrand = n_elements(randnum)
    par = self->par_struct()

    if nrand eq 0 then begin
        outfile = self->corrfile('secondary','input',/createdir)
        print
        print,'Will write to file: ',outfile

        gals=self->get(/gals)
        ntot=n_elements(gals)
        print
        print,'Making basic ra/dec cuts'
        keep = self->radec_cuts(gals.ra, gals.dec, nkeep)
        print,'Kept: '+ntostr(nkeep)+'/'+ntostr(ntot)

        outst = self->secondary_struct(nkeep)

        print,'  copying...',f='($,a)'
        outst.ra = gals[keep].ra
        outst.dec = gals[keep].dec

        nmgy=lups2nmgy(gals[keep].omag)
        outst.gflux = reform(nmgy[1,*])
        outst.rflux = reform(nmgy[2,*])
        outst.iflux = reform(nmgy[3,*])

        print,'htm_index...'
        outst.htm_index = htm_index(outst.ra, outst.dec, par.depth)

        print
        print,'Writing to file: ',outfile
        stop
        write_idlstruct, outst, outfile
    endif else begin

        for i=0l, nrand-1 do begin

            rstruct = self->genrand('secondary')
            routfile = $
                self->corrfile('secondary_random','input',$
                                secondary_randnum=randnum[i],/createdir)
            print
            print,'Writing random file: ',routfile
            write_idlstruct, rstruct, routfile

        endfor 
 
    endelse
end


pro maxbcg_corrsim__define

  struct = { $
             maxbcg_corrsim,     $
             sample: 0,          $
             inherits correlate, $
             inherits maxbcg     $
           }

end 
