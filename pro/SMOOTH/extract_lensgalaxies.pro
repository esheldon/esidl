pro extract_lensgalaxies, pstruct, clr, indices, min_mag=min_mag, max_mag=max_mag, plot=plot, silent=silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Extracts a clean set of galaxies
; 
; Inputs:  pstruct: input photo structure 
;	   clr: bandpass to select on:
;          ostruct: output photo structure containing galaxies
;          max_mag: maximum magnitude to use (default=24.0)
;	   silent: don't make plots
;
; Outputs: Plots flags for these objects....
;	   indices: returns the indices of the galaxies in the original
;			struct
;
; Author:  Phil Fischer (extract_stars)
; Date: 1/14/99
; Altered to get galaxies: Erin Scott Sheldon
; Date: 2/19/99
; Added silent and indices options: Tim McKay
; Date: 6/10/99
; Added altered to make cuts based on Phil's r parameter: Tim McKay
; Date: 6/23/99
; No longer creates a new struct.  Instead just returns indices.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Help message
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if n_params() eq 0 then begin
    print,'-Syntax extract_galaxies, pstruct, clr, indices, min_mag=min_mag, max_mag=max_mag, plot=plot, silent=silent'
    return
  ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(min_mag) THEN min_mag = 16.0
  IF NOT keyword_set(max_mag) THEN max_mag = 22.0
  IF NOT keyword_set(silent) THEN silent = 0 ELSE silent = 1
  IF keyword_set(plot) THEN plot = 1 ELSE plot = 0
  colors=['u','g','r','i','z']

  IF NOT silent THEN print,'Number of entries: ',n_elements(pstruct)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Objects must have well measured shapes. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(pstruct.e1(clr) NE 1.e10 AND pstruct.e2(clr) NE 1.e10, nw)
  indices=w
  IF NOT silent THEN print,'Objects with well measured shapes: ',ntostr(nw)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Cut on polarizability
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(pstruct[indices].r(clr) lt 0.8 AND $
            pstruct[indices].r(clr) NE -1, nw)
  indices=indices(w)
  IF NOT silent THEN print,'After polarizability cuts: ',ntostr(nw)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; cut on maximum and minimum magnitudes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT silent THEN print,'Selecting gal between '+$
                           colors[clr]+' = '+ntostr(min_mag)+' and '+$
                           ntostr(max_mag)
  w=where(pstruct[indices].petrocounts(clr) lt max_mag AND $
            pstruct[indices].petrocounts(clr) GT min_mag, nw)
  indices=indices[w]

  IF NOT silent THEN print,'After magnitude cuts: ', ntostr(nw) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if plot then begin
    title='Extracted Galaxies'
    xtitle='petrocounts'
    ytitle='petrorad'

    plot,pstruct[indices].petrocounts(clr), $
         pstruct[indices].petrorad(clr), $
;         pstruct[indices].r(clr), $
         yrange=[0,10],$
         psym=3,$
         title=title,xtitle=xtitle,ytitle=ytitle
  endif

  return 
END 










