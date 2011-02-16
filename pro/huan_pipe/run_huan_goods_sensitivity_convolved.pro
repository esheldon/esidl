pro run_huan_goods_sensitivity_convolved, dops=dops, type=type

    if n_elements(type) eq 0 then type='des5yr'

    outdir = esheldon_config('des_sensitivity_dir')
    outdir = path_join(outdir, 'sensitivity.new')
    plotdir = path_join(outdir, 'plots')

    if not fexist(outdir) then file_mkdir, outdir
    if not fexist(plotdir) then file_mkdir, plotdir

    stfile = path_join(outdir, type+'_convolved.st')
    fitsfile = path_join(outdir, type+'_convolved.fits')

;  seeing = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, $
;            0.50, 0.55, 0.60, 0.65, 0.70, 0.76, $
;            0.80, 0.85, 0.90, 0.95, 1.00, 1.05, $
;            1.10, 1.15, 1.20]

    ;; for 0.7 use 0.675
    ;; for 0.8 use 0.78
    ;; for 0.9 use 0.87

    seeing = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, $
            0.50, 0.55, 0.60, 0.65, $
            0.675, $
            0.70, 0.76, $
            0.78, $
            0.80, 0.85, $
            0.87, $
            0.90, 0.95, 1.00, 1.05, $
            1.10, 1.15, 1.20]
    psfiles = path_join(plotdir, $
        type+'_convolved'+'_'+ntostr(seeing, f='(f0.3)')+'.ps')

    nseeing = n_elements(seeing)
    for i=0l, nseeing-1 do begin 

        if keyword_set(dops) then begplot,psfiles[i],/color
        tsens_struct = $
            huan_goods_sensitivity_convolved(seeing[i], type=type)
        if keyword_set(dops) then endplot

        if i eq 0 then begin 
            sens_struct = replicate(tsens_struct, nseeing)
        endif else begin 
            tmp = sens_struct[i]
            struct_assign, tsens_struct, tmp
            sens_struct[i] = tmp
        endelse 

    endfor 

    print
    print,'Writing to idlstruct file: ',stfile
    write_idlstruct, sens_struct, stfile
    print
    print,'Writing to fits file: ', fitsfile
    mwrfits, sens_struct, fitsfile, /create

end 
