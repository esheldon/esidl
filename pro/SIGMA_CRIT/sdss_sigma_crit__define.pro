; pzsample is the sample for the zweight class
function sdss_sigma_crit::init, zwsample
    if n_elements(zwsample) eq 0 then begin
        message,'zwsample not entered, using 2',/inf
        zwsample=2
    endif

    retval = self->zweight::init(zwsample)
    return, retval

end




; return inverse critical density in units of (msolar/pc^2)^(-1)
; for a range of zl given zs,pofzs
function sdss_sigma_crit::calc_scinv, zl, zs, pzs, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h, $
        npts=npts

    nzl = n_elements(zl) & nzs=n_elements(zs) & npzs=n_elements(pzs)
    if nzl eq 0 or nzs eq 0 or npzs eq 0 then begin
        print,'-Syntax: meanscinv=sc->calc_scinv(zl, zs, pzs, omega_m=, omega_l=, omega_k=, flat=, h=)'
        on_error, 2
        message,'Halting'
    endif

    if n_elements(npts) eq 0 then npts = 200
    gauleg, min(zs), max(zs), npts, ZZi, WWi

    scinv = fltarr(nzl)

    norm = total( WWi*pzs )
    pofzs = pzs/norm

    for i=0L, nzl-1 do begin
        tscinv = sigmacritinv(zl[i], zs, $
            omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h)
        
        scinv[i] = total( WWi*tscinv*pofzs )

    endfor

    return, scinv

end



;; generate sigma crit inv at a bunch of grid points
FUNCTION sdss_sigma_crit::sigmacritinv_func, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h, $
        photoz=photoz

  sdssidl_setup, /silent
  esheldon_setup,/silent

  IF n_elements(npts) EQ 0 THEN npts = 100

  ;; Read the source redshift distribution
  ;; defaults to Hybrid

  pz = obj_new('photoz_uchicago',type='nn')
  nzstruct = pz->deconv_read(photoz=photoz)

  ;; Get normalized P(z)
  z1 = min(nzstruct.z, max=z2)

  gauleg, z1, z2, npts, ZZi, WWi
  nofz = interpol( nzstruct.dndz, nzstruct.z, ZZi )

  norm = total( WWi*nofz )
  pofz = nofz/norm

  obj_destroy,pz

  ;; output lens redshifts

  nlens = 1000L
  zlens = dblarr(nlens)
  scritinv = dblarr(nlens)

  FOR i=0, nlens-1 DO BEGIN 
      zlens[i] = .005 + i*.001
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate mean inverse critical density at each zlens
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
  FOR i=0L, nlens-1 DO BEGIN 
      sig_inv = sigmacritinv(zlens[i], ZZi, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h)

      scritinv[i] = total(WWi*sig_inv*pofz)
  ENDFOR 

  ;; output structure
  s=create_struct('zlens', 0., 'sigcritinv', 0.)
  struct = replicate(s, nlens)
  struct.zlens = zlens
  struct.sigcritinv = scritinv

  return,struct

END 

;; generate sigma crit inv at a bunch of grid points
;; In this one, we assume the user chooses only sources with
;; z > zL so the inverse critical density will be larger
FUNCTION sdss_sigma_crit::sigmacritinv_func_limit, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h, $
        photoz=photoz

  sdssidl_setup, /silent
  esheldon_setup,/silent

  IF n_elements(npts) EQ 0 THEN npts = 100

  ;; Read the source redshift distribution
  ;; defaults to Hybrid

  pz = obj_new('photoz_uchicago',type='nn')
  nzstruct = pz->deconv_read(photoz=photoz)

  ;; output lens redshifts and scritinv
  nlens = 1000L
  zlens = dblarr(nlens)
  scritinv = dblarr(nlens)

  FOR i=0, nlens-1 DO BEGIN 
      zlens[i] = .005 + i*.001
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate mean inverse critical density at each zlens
  ;; Because we limit each time, must always re-calculate the weights
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
  FOR i=0L, nlens-1 DO BEGIN 

      z1 = zlens[i]
      z2 = max(nzstruct.z)

      gauleg, z1, z2, npts, ZZi, WWi
      nofz = interpol( nzstruct.dndz, nzstruct.z, ZZi )

      sig_inv = sigmacritinv(zlens[i], ZZi, $
            omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h)
      scritinv[i] = total(WWi*sig_inv*nofz)/total(WWi*nofz)
  ENDFOR 

  obj_destroy,pz


  ;; output structure
  s=create_struct('zlens', 0., 'sigcritinv', 0.)
  struct = replicate(s, nlens)
  struct.zlens = zlens
  struct.sigcritinv = scritinv

  return,struct

END 


PRO sdss_sigma_crit::plot, dops=dops

  psfile = '~/plots/photoz/uchicago/sigmacrit_compare_hybrid.eps'
  botytitle = 'deconvolved/photoz'
  leg = 'Deconvolved'


  IF keyword_set(dops) THEN BEGIN 
      begplot,name=psfile,/encapsulated
  ENDIF 
  oline = 2

  omega_m = 0.27
  sc_calc   = self->sigmacritinv_func(omega_m)
  sc_photoz = self->sigmacritinv_func(omega_m, /photoz)

  ratio = sc_calc.sigcritinv/sc_photoz.sigcritinv

  siginv = !csym.sigma_cap+ '!S!U'+!csym.minus+'1!N!R!Dcrit!N'
  ytitle = siginv + ' [10!U-4!N pc!U2!N / M'+sunsymbol()+' ]'
  xtitle = 'z!DLens!N'

  xrange = [0.0, 0.4]
  topyrange = [0,1.05]
  botyrange = [0.95, 1.05]
  plottwo,sc_calc.zlens,sc_calc.sigcritinv/1.e-4, sc_calc.zlens, ratio,$
          xoplot=sc_photoz.zlens,yoplot=sc_photoz.sigcritinv/1.e-4,oplotline=oline, $
          frac1=0.75, topaspect=1, $
          /ynozero,$
          xrange=xrange, xstyle=1, $
          topyrange=topyrange, ystyle=1+2,$
          botyrange=botyrange, $
          xtitle=xtitle, $
          topytitle=ytitle, $
          botytitle = botytitle, $
          xwindow1=xwindow1, ywindow1=ywindow1

  oplot,[0,10],[1,1]

  xwold = !x.window & !x.window=xwindow1
  ywold = !y.window & !y.window=ywindow1
          
  mess = [leg,'Photoz']
  line = [0,oline]

  legend, mess, /right, box=0, charsize=1, line=line

  !x.window=xwold & !y.window=ywold

  IF keyword_set(dops) THEN endplot,/trim_bbox

END 



; This is an integral over all sources
FUNCTION sdss_sigma_crit::sigmacritinv, omegamat, zlens, $
                        h=h, $
                        scinvstruct=scinvstruct, $
                        npts=npts, $
                        wgood=wgood, cgs=cgs, $
                        photoz=photoz, $
                        limit=limit

  on_error,2

  IF N_params() LT 2 THEN BEGIN 
     print,'-Syntax: result = sc->sigmacritinv(omegmat, zlens, h=, npts=, scinvstruct=, wgood=, /cgs, /photoz)'
     return,-1
  ENDIF 

  ntot = n_elements(zlens)

  ;; parameters
  ;; to convert to gm^-1 cm^2 from Msolar^-1 pc^2, mult. by fac2/fac1
  fac1 =   .3474                ;g/cm^2
  fac2 = 1.663e3                ;Msolar/parsec^2
  tol = 1.e-8  
  IF keyword_set(cgs) THEN fac = fac2/fac1 ELSE fac = 1.0

  ;; calculate for this cosmology
  IF keyword_set(limit) THEN BEGIN 
      scinvstruct = self->sigmacritinv_func_limit(omegamat, h=h, npts=npts, photoz=photoz)
  ENDIF ELSE BEGIN 
      scinvstruct = self->sigmacritinv_func(omega_m=omegamat, h=h, photoz=photoz)
  ENDELSE 
  ;; range checking
  maxz = max(scinvstruct.zlens)
  minz = min(scinvstruct.zlens)

  wgood = where(zlens LE maxz, ngood)
  IF ngood EQ 0 THEN BEGIN
      print,'No lenses with z < maxz.'
      return,-1
  ENDIF 
  wgood2 = where(zlens[wgood] GE minz, ngood)
  IF ngood EQ 0 THEN BEGIN
      print,'No lenses with z > minz.'
      return,-1
  ENDIF 
  wgood = wgood[wgood2]

  ;; interpolate
  IF ntot EQ 1 THEN sigcritinv = -1000d $
  ELSE sigcritinv = replicate(-1000d, ntot)

  sigcritinv[wgood] = $
    interpol( scinvstruct.sigcritinv, scinvstruct.zlens, zlens[wgood] )*fac

  ;; check versus tolerance
  wgood3 = where(sigcritinv[wgood] GT tol, ngood)
  IF ngood EQ 0 THEN BEGIN 
      print,'No lenses with 1/sig_crit > ',ntostr(tol)
      return,-1
  ENDIF 

  wgood = wgood[wgood3]
  return, sigcritinv


END 




; Generate the mean inverse critical density 
;    <scinv>(zL, zphotz0
; given the input ztrue distributions as a function of zphot.
;
; We have characterized ztrue in various photoz bins, which are irregularly
; gridded.  We will first calculae the mean scinv in each of these bins for
; a grid of zlens.  Then we will go back and interpolate to a regular grid
; of zsource for each zlens

pro sdss_sigma_crit::generate_meanscinv, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h

    file = self->file('output','meanscinv', project='scinv', /createdir)
    print
    print,'Will write to file: ', file

    zpbins = self->zweight::read('output','zbinhist')
    nzp = n_elements(zpbins)

    npts = 200

    nzl = 40
    nzs = 40

    zlmin = 0.044
    zlmax = 0.3
    zlstep = (zlmax-zlmin)/(nzl-1)

    zsmin = 0.1
    zsmax = 1.0
    zsstep = (zsmax-zsmin)/(nzl-1)

    zl = arrscl( findgen(nzl), zlmin, zlmax )
    zs = arrscl( findgen(nzs), zsmin, zsmax )


    ; First pass over irregular grid of photoz bins
    print
    print,'Doing irregular grid on zphot bins'
    tscinv = dblarr(nzl, nzp)
    for is=0L, nzp-1 do begin
        tzs = zpbins[is].xcenter
        tpzs = zpbins[is].whist
        tscinv[*, is] = self->calc_scinv(zl, tzs, tpzs , $
            omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h, $
            npts=npts)
    endfor

    ; Now go through and re-interpolate to the grid
    print,'Interpolating to a regular grid'
    scinv = dblarr(nzl, nzs)
    for il=0L, nzl-1 do begin
        scinv[il,*] = interpol( tscinv[il,*], zpbins.zpmean, zs )
    endfor

    ; Now write them all out
    print,'Writing to file: ',file
    openw, lun, file, /get_lun

    writeu, lun, long(npts)

    writeu, lun, long(nzl)
    writeu, lun, long(nzs)

    writeu, lun, zlstep
    writeu, lun, zlmin
    writeu, lun, zlmax
    writeu, lun, zl
    writeu, lun, findgen(nzl)

    writeu, lun, zsstep
    writeu, lun, zsmin
    writeu, lun, zsmax
    writeu, lun, zs
    writeu, lun, findgen(nzs)

    for i=0L, nzl-1 do begin
        for j=0L, nzs-1 do begin
            writeu, lun, scinv[i,j]
        endfor
    endfor

    free_lun, lun

end

function sdss_sigma_crit::meanscinv_struct, nzl, nzs
    st = {  $
        npts: 0L, $
        nzl: long(nzl), $
        nzs: long(nzs), $
        zlstep: 0.0,    $
        zlmin:  0.0,    $
        zlmax:  0.0,    $
        zl: fltarr(nzl), $
        zli: fltarr(nzl), $
        $
        zsstep: 0.0,    $
        zsmin:  0.0,    $
        zsmax:  0.0,    $
        zs: fltarr(nzs), $
        zsi: fltarr(nzs), $
        $
        scinv: dblarr(nzl, nzs) $
    }

    return, st
 
end
function sdss_sigma_crit::read_meanscinv
     
    file = self->file('output','meanscinv', project='scinv')
    print,'Reading from file: ', file
    openr, lun, file, /get_lun

    npts = 0L
    nzl = 0L
    nzs = 0L

    readu, lun, npts
    readu, lun, nzl
    readu, lun, nzs

    print,'nzl = ', nzl
    print,'nzs = ', nzs
    st = self->meanscinv_struct(nzl, nzs)
    st.npts = npts

    tflt = 0.0
    readu, lun, tflt
    st.zlstep = tflt
    readu, lun, tflt
    st.zlmin = tflt
    readu, lun, tflt
    st.zlmax = tflt
    for i=0L, nzl-1 do begin
        readu, lun, tflt
        st.zl[i] = tflt
    endfor
    for i=0L, nzl-1 do begin
        readu, lun, tflt
        st.zli[i] = tflt
    endfor

    readu, lun, tflt
    st.zsstep = tflt
    readu, lun, tflt
    st.zsmin = tflt
    readu, lun, tflt
    st.zsmax = tflt
    for i=0L, nzl-1 do begin
        readu, lun, tflt
        st.zs[i] = tflt
    endfor
    for i=0L, nzl-1 do begin
        readu, lun, tflt
        st.zsi[i] = tflt
    endfor

    tdbl = 0d
    for i=0L, nzl-1 do begin
        for j=0L, nzs-1 do begin
            readu, lun, tdbl
            st.scinv[i,j] = tdbl
        endfor
    endfor


    free_lun, lun

    return, st
    
end


function sdss_sigma_crit::interp, scstruct, zl, zs

    nst=n_elements(scstruct) & nl=n_elements(zl) & ns=n_elements(zs)
    if nst eq 0 or nl eq 0 or ns eq 0 then begin
        print,'-Syntax: scinv=sc->interp(scstruct, zl, zs)'
        on_error, 2
        message,'Halting'
    endif
    if nl ne 1 then message,'Only scalar zl allowed for now'
    if ns ne 1 then message,'Only scalar zs allowed for now'
    
    wl=where(zl lt scstruct.zlmin or zl gt scstruct.zlmax, nwl)
    if nwl ne 0 then begin
        message,'zl out of range values encountered', /inf
    endif
    ws=where(zs lt scstruct.zsmin or zs gt scstruct.zsmax, nws)
    if nws ne 0 then begin
        message,'zs out of range values encountered', /inf
    endif
    zli = interpol( scstruct.zli, scstruct.zl, zl)
    zsi = interpol( scstruct.zsi, scstruct.zs, zs)

    return, interpolate( scstruct.scinv, zli, zsi, missing=-1d)

end
PRO sdss_sigma_crit__define

  struct = {$
             sdss_sigma_crit,   $
             ssc_dummy: 0,          $
             inherits zweight   $
           }

END 
