;+
; NAME:
;   icl_gonzales_read
;
; PURPOSE:
;   read in parameters for ICL as measured by Gonzales et al. 2005
;
; CALLING SEQUENCE:
;   st=icl_gonzales_read(type)
;
; INPUTS:
;   type: either 'sersic' or '2dev'
;
; OUTPUTS:
;   a structure representing the table in that paper.
;
; MODIFICATION HISTORY:
;   Documented: 2007-04-15 Erin Sheldon, NYU
;-

function icl_gonzales_sersic_read

    ;; note they used h=0.7, so we need to convert
    hratio = (0.7/1.0)
    magdiff = 5.0*alog10(hratio)

    ;dir = '/global/early2/esheldon/icl/gonzales'
    ;file = concat_dir(dir, 'sersic.dat')

    icl=obj_new('icl')
    file = icl->gonzales_file('sersic')
    obj_destroy, icl

    sdef =                  $
        {                   $
            name:mkstr(11), $
            absmag: 0.0,       $
            absmag_err: 0.0,   $
            r50: 0.0,        $
            r50_err: 0.0,    $
            r0: 0.0,        $
            n: 0.0,         $
            n_err: 0.0,     $
            ab: 0.0,        $
            ab_err: 0.0,    $
            dchi2: 0L       $
        }

    readcol, file, format='a11,f,f,f,f,f,f,f,f,I', name, absmag, absmag_err, r50, r50_err, n, n_err, ab, ab_err, dchi2

    ; convert sizes and mags
    r50 = r50*hratio
    r50_err = r50_err*hratio

    r0 = sersic_r0(r50, n) 

    absmag = absmag+magdiff

    nstruct=n_elements(name)
    struct = replicate(sdef, nstruct)
    struct.name=name
    struct.absmag=absmag
    struct.absmag_err=absmag_err
    struct.r50=r50
    struct.r50_err=r50_err
    struct.r0=r0
    struct.n=n
    struct.n_err=n_err
    struct.ab=ab
    struct.ab_err=ab_err
    struct.dchi2=dchi2
    return, struct

end

function icl_gonzales_2dev_read

    ;; note they used h=0.7, so we need to convert
    hratio = (0.7/1.0)
    magdiff = 5.0*alog10(hratio)

    ;dir = '/global/early2/esheldon/icl/gonzales'
    ;file = concat_dir(dir, '2dev.dat')

    icl=obj_new('icl')
    file = icl->gonzales_file('2dev')
    obj_destroy, icl


    sdef =                   $
        {                    $
            name:mkstr(11),  $
            absmag: 0.0,     $
            dphi: 0.0,       $
            dphi_err:0.0,    $
            absmag1: 0.0,    $
            absmag1_err:0.0, $
            r50_1: 0.0,      $
            r50_1_err:0.0,   $
            ab1: 0.0,        $
            ab1_err: 0.0,    $
            dchi2d: 0L,      $
            dchi2s: 0L,      $
            absmag2: 0.0,    $
            absmag2_err:0.0, $
            r50_2: 0.0,      $
            r50_2_err: 0.0,  $
            ab2: 0.0,        $
            ab2_err: 0.0     $
       }

    format='a11, f,f,f, f,f,f,f,f,f, I,I, f,f,f,f,f,f'
    readcol, file, format=format, $
        name, $
        absmag, dphi, dphi_err, $
        absmag1, absmag1_err, r50_1, r50_1_err, ab1, ab1_err, $
        dchi2d, dchi2s, $
        absmag2, absmag2_err, r50_2, r50_2_err, ab2, ab2_err
    

    ; convert sizes and mags
    r50_1 = r50_1*hratio
    r50_1_err = r50_1_err*hratio
    r50_2 = r50_2*hratio
    r50_2_err = r50_2_err*hratio
    absmag = absmag+magdiff
    absmag1 = absmag1+magdiff
    absmag2 = absmag2+magdiff

    nstruct=n_elements(name)
    struct = replicate(sdef, nstruct)
    struct.name=name
    struct.absmag=absmag
    struct.dphi = dphi
    struct.dphi_err = dphi_err

    struct.absmag1 = absmag1
    struct.absmag1_err = absmag1_err
    struct.r50_1=r50_1
    struct.r50_1_err=r50_1_err
    struct.ab1 = ab1
    struct.ab1_err = ab1_err

    struct.dchi2d = dchi2d
    struct.dchi2s = dchi2s

    struct.absmag2 = absmag2
    struct.absmag2_err = absmag2_err
    struct.r50_2=r50_2
    struct.r50_2_err=r50_2_err
    struct.ab2 = ab2
    struct.ab2_err = ab2_err

    return, struct

end


function icl_gonzales_read, type
    if n_elements(type) ne 1 then begin
        print,'-Syntax: st=icl_gonzales_read(type)'
        print,'type = "sersic" or "2dev"'
        on_error, 2
        message,'Halting'
    endif
    case strlowcase(type) of
        'sersic': return, icl_gonzales_sersic_read()
        '2dev': return, icl_gonzales_2dev_read()
        else: message,'Bad type: '+ntostr(type)
    endcase
    
end
