;+
; NAME:
;   sersic_r0
;
; PURPOSE:
;   Convert half-light radius r50 to sersic r0 parameters.
;      sersic = exp( -(r/r0)^(1/n) )
;   Thus the r50 is related to r0 through the roots of the incomplete gamma function. 
;
; CALLING SEQUENCE:
;   r0 = sersic_r0(r50, n, /inverse, status=)
;
; INPUTS:
;   r50: the half light radius. 
;   n: the index. Range of tested values [0.05, 20]
;
; KEYWORD PARAMETERS:
;   /inverse: The input is actually r0 and r50 should be computed.
;
; OUTPUTS:
;   r0: r0 in the sersic profile.
; OPTIONAL OUTPUTS:
;   status: 0 if success, 1 if failed to converge.
;
; COMMENTS:
;   Has been tested to work for n ranging [0.05,20], including inverse.
;
; MODIFICATION HISTORY:
;   Created: 2007-04-11 Erin Sheldon, NYU 
;-

function sersic_igamma, x
    common sersic_igamma_block, a
    tol = 1.e-4
    return, igamma(a, x > tol, /double) - 0.5
end


function sersic_r0guess, n

    ; use these known values to get good guesses from n=0.1 to 20
    ng   = [0.05,   0.1,  0.5, 2.0, 3.0, 4.0, 10.0, 15.0, 20.0]
    glow = [0.0002, 0.01, 0.4, 3.0, 5.0, 6.0, 15.0, 25.0, 35.0]
    g    = [0.0006, 0.02, 0.7, 3.7, 5.7, 7.7, 19.7, 29.7, 39.7]
    ghigh= [0.001,  0.03, 1.0, 4.4, 6.4, 8.0, 25.0, 35.0, 45.0]

    guessl = interpol(glow,ng,n)
    guess = interpol(g,ng,n)
    guessh = interpol(ghigh,ng,n)

    guess = [guessl, guess, guessh]

    return, guess
end

function sersic_r0, r50_in, n_in, inverse=inverse, root=root, guess=guess_in, status=status
    
    common sersic_igamma_block, a


    status=1
    n50=n_elements(r50_in) & nn=n_elements(n_in)
    if n50 eq 0 or nn eq 0 then begin
        print,'-Syntax: r0 = sersic_r0(r50, n, /inverse, status=)'
        on_error, 2
        message,'Halting'
    endif

    if n50 ne nn then message,'r50 and n must be same length'

	r50 = double(r50_in)
	n = double(n_in)

    rout = r50
    for i=0L, n50-1 do begin

        a = 2.0*n[i]

        if n_elements(guess_in) eq 0 then guess=sersic_r0guess(n[i]) else guess=guess_in

        ; we need to catch the error
        command = "root = fx_root(guess,'sersic_igamma', /double)"
        if not execute(command) then begin
            print,'Could not calculate r0 from (r50,n) = ',r50[i],n[i]
            rout[i] = -1
        endif

        ; This root is in the (r/r0)^(1/n) space
        ;   i.e.  root = (r50/r0)^(1/n)

        status = 0
        if not keyword_set(inverse) then begin
            ; r50 was input
            tmp_r0 = r50[i]/root^n[i]
            rout[i] = tmp_r0
        endif else begin
            ; r0 was input
            tmp_r50 = r50[i]*root^n[i]
            rout[i] = tmp_r50
        endelse

    endfor  

    return, rout

end
