function _maxbcg_lensing_m2l_table_getsubstr, subtype, struct
    case subtype of
        'ngals200_12': begin
            tagname = 'mean_ngals200'
            format='(I0)'
        end
        'ilum200_16': begin
            tagname = 'mean_ilum200'
            format='(F5.1)'
        end
        else: message,'unknown type '+subtype
    endcase
    if not tag_exist(struct, tagname, ind=ind) then message,'tag not found'

    str = string(struct.(ind), format=format)
    return,str
    
end

pro _maxbcg_lensing_m2l_table_name_label, subtype, name, label
    case subtype of
        'ngals200_12': begin
            name = '\nvir'
            label = 'm2lngals'
        end
        'ilum200_16': begin
            name = '\lvir'
            label = 'm2llum'
        end
        else: message,'unknown type '+subtype
    endcase
end

pro maxbcg_lensing_m2l_table, subtype, m2lstruct=m2lstruct

    if n_elements(subtype) eq 0 then begin
        print,'maxbcg_lensing_plot_m2l, subtype'
        return
    endif

    ; hard wiring the samples for simplicity
    mass_sample = [21,22]
    lsample = 4

    ; build the basic dir and file names
    dir = '~/plots/m2l/maxbcg/'+subtype
    if not fexist(dir) then file_mkdir, dir

    mstr = strjoin( ntostr(mass_sample), '-')
    lstr = ntostr(lsample)

    front = 'm2l-'+subtype+'-m'+mstr+'-l'+lstr


    output_name = path_join(dir, front+'-m2l.fits')
    print,'Reading ',output_name
    m2lstruct = mrdfits(output_name,1)


    mobj = obj_new('maxbcg_lensing', mass_sample[0]) ; vers doesn't matter here
    
    m=mobj->lensread('jackknife_invert', sub=subtype, samp=mass_sample)
    n=n_elements(m2lstruct)

    ; asymptotic had some instability
    wbad=where(m2lstruct.m2lasym gt 700,nbad,comp=comp)

    m2la = m2lstruct.m2lasym
    m2la_err = m2lstruct.m2lasym_err

    if nbad ne 0 then begin
        ;; get average conversion of regular error to mcmc error
        errfac = mean( m2la_err[comp]/m2lstruct[comp].m2lasym_err_0 )
        m2la[wbad] = m2lstruct[wbad].m2lasym_0
        m2la_err[wbad] = m2lstruct[wbad].m2lasym_err_0*errfac
    endif

    _maxbcg_lensing_m2l_table_name_label, subtype, name, label
    print
    print
    print,'\begin{deluxetable}{cccccc}'
    print,'\tabletypesize{\small}'
    print,'\tablecaption{M/L Statistics for '+name+'\ Bins \label{tab:'+label+'}}'
    print,'\tablewidth{0pt}'
    print,'\tablehead{'
    print,'    \colhead{Mean '+name+'}       &'
    print,'    \colhead{$L_{200}^{tot}$}         &'
    print,'    \colhead{$M_{200}$}  &'
    print,'    \colhead{$M/L_{200}$}  &'
    print,'    \colhead{$M/L_{22}$}  &'
    print,'    \colhead{$M/L_{asym}$}  '
    print,'}'
    print,'\'
    print,'\startdata'

    ; substr & l200+/-err & m200+/-err & M/L200+/-err & M/Llast+/-err & M/Lasymp+/-err
    for i=0L, n-1 do begin
        substr = _maxbcg_lensing_m2l_table_getsubstr(subtype, m[i])

        l200=string(m2lstruct[i].l200/1.e10, format='(F6.2)')
        l200err=string(m2lstruct[i].l200_err/1.e10, format='(F6.2)')

        m200=string(m2lstruct[i].m200/1.e12, format='(F6.2)')
        m200err=string(m2lstruct[i].m200_err/1.e12, format='(F6.2)')

        m2l200=string(m2lstruct[i].m2l200, format='(I0)')
        m2l200err=string(m2lstruct[i].m2l200_err, format='(I0)')

        m2lmax=string(m2lstruct[i].m2lmax, format='(I0)')
        m2lmaxerr=string(m2lstruct[i].m2lmax_err, format='(I0)')

        m2lasym=string(m2la[i], format='(I0)')
        m2lasymerr=string(m2la_err[i], format='(I0)')

        print,$
            substr+' & '+$
            '$'+l200+' \pm '+l200err+'$ & '+$
            '$'+m200+' \pm '+m200err+'$ & '+$
            '$'+m2l200+' \pm '+m2l200err+'$ & '+$
            '$'+m2lmax+' \pm '+m2lmaxerr+'$ & '+$
            '$'+m2lasym+' \pm '+m2lasymerr+'$  ', format='($,a)'

        if i eq (n-1) then print,'' else print,' \\'

    endfor
    print,'\enddata'
    print,'\end{deluxetable}'

end
