pro phextract_stars, pstruct, color_index, indices, max_mag=max_mag, sig_clip=sig_clip, plot1=plot1, silent=silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Extracts a clean set of stars suitable for PSF determination
; 
; Inputs:  pstruct: photo structure 
;	   color_index: bandpass to select on:
;          ostruct: output photo structure containing stars
;          max_mag: maximum magnitude to use (default=20)
;	   sig_clip: Number of sigma for clipping large radius objects (default=3)
; Outputs: ostruct: structure containing stars
;
; Author:  Phil Fischer
; Date: 1/14/99
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Help message
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  if n_params() eq 0 then begin
      print,'-syntax extract_stars, pstruct, color_index, ostruct, max_mag=max_mag'
      return
  endif
  
;extract objects classified as stars
  
  IF NOT keyword_set(silent) THEN sil=0 ELSE sil=1
  if not keyword_set(plot1) then plot1 = 0 else plot1 = 1
  
  indices=where( (pstruct.type[1] eq 6) and (pstruct.type[2] EQ 6) AND $
                 (pstruct.ixx[color_index] gt 0.) , nss)
  
  IF nss EQ 0 then begin
      print,'Not enough stars 1'
      return
  endif 
  
  ;;;;;;;;;;  extract based on object1 flags in selected bandpass
  
  make_flag_struct, fs
  fs.canonical_center='N' 
  fs.edge='N' 
  fs.blended='N' 
;fs.child='N' 
  fs.peakcenter='N' 
  fs.nodeblend='N' 
  fs.nopetro='N'   
;fs.manypetro='N'   
  fs.manyr50='N'   
  fs.manyr90='N'   
  fs.incomplete_profile='N'   
  fs.interp='N'   
  fs.notchecked='N'   
  fs.subtracted='N'   
  fs.nostokes='N'   
  fs.badsky='N'   
;fs.petrofaint='N'   
;fs.too_large='N'   
  fs.deblended_as_psf='N'
  fs.deblend_pruned='N'   
  fs.ellipfaint='N'   
  fs.moved='N'   

  flag_select,pstruct[indices],fs,color_index,ind
  if(ind[0] EQ -1) then begin
      print,'Not enough stars 2'
      indices=ind
      return
  endif
  indices = indices[ind]
  
  ;;;;;;  objc flags  
  fs.bright='N'
  fs.blended='N'
;fs.child='N'
  make_flag_struct, fs
  flag_select,pstruct[indices],fs,color_index,ind,objc=1
  make_flag_struct, fs
  
  if(n_elements(indices) eq 1) then begin
      print,'Not enough stars 3'
      indices=ind
      return
  ENDIF
  indices = indices[ind]
  
  ;;;;; Magnitude cuts
  maxm=20.
  if keyword_set(max_mag) then begin
      maxm=max_mag
  endif
  ind=where(pstruct[indices].petrocounts[color_index] lt maxm, nind)
  if(nind LT 2) then begin
      print,'Not enough stars 4'
      indices=ind
      return
  endif
  indices = indices[ind]
  
  sigc=3
  if keyword_set(sig_clip) then sigc=sig_clip
  
  ;;;;; Cut on size
  IF NOT sil THEN print,' pre-clip ',n_elements(indices)
  size=pstruct[indices].ixx[color_index]+pstruct[indices].iyy[color_index]
  ind=where(size gt 0, nind)
  IF nind LT 2 THEN BEGIN
      print,'Not enough stars 5'
      indices=ind
  ENDIF 
  indices = indices[ind]
  
  size=pstruct[indices].ixx[color_index]+pstruct[indices].iyy[color_index]
  npt=nind
  fsize=sqrt(size/2.)*0.4*2.35  ;???
  
  IF (plot1 eq 1) THEN BEGIN
      !p.multi=[0,1,2]
      plot,pstruct[indices].petrocounts[color_index],fsize,yrange=[0,4],psym=3
  ENDIF 
  
  _dd=moment(size)
  ind=where(size lt _dd(0)+sigc*sqrt(_dd(1)), nind)
  IF nind LT 2 then begin
      print,'Not enough stars 6'
      indices=ind
      return
  ENDIF
  indices = indices[ind]
  
  size=pstruct[indices].ixx[color_index]+pstruct[indices].iyy[color_index]
  ind=where(size gt _dd(0)-sigc*sqrt(_dd(1)), nind)
  if(nind LT 2) then begin
      print,'Not enough stars 7'
      indices=ind
      return
  endif
  indices = indices[ind]
  
  i=0
  npt=n_elements(indices)+1
  while (n_elements(indices) lt npt) do begin
      i=1+1
      npt=n_elements(indices)
      size=pstruct[indices].ixx[color_index]+pstruct[indices].iyy[color_index]
      _dd=moment(size)
      ind=where(size lt _dd[0]+sigc*sqrt(_dd[1]), nind)
      IF nind LT 2 then begin
          print,'Not enough stars'
          indices=ind
          return
      ENDIF
      indices=indices[ind]
      
      size=pstruct[indices].ixx[color_index]+pstruct[indices].iyy[color_index]
      ind=where(size gt _dd[0]-sigc*sqrt(_dd[1]), nind)
      IF nind LT 2 then begin
          print,'Not enough stars'
          indices=ind
          return
      endif 
      indices=indices[ind]
  endwhile
  
  IF NOT sil THEN print,' after clip ',n_elements(indices),' iterations ',i,' size= ',_dd[0],sqrt(_dd[1])
  
  
  size=pstruct[indices].ixx[color_index]+pstruct[indices].iyy[color_index]
  
  fsize=sqrt(size/2.)*0.4*2.35
  IF (plot1 eq 1) THEN BEGIN
      plot,pstruct[indices].petrocounts[color_index],fsize,yrange=[0,4],psym=3
      
      !p.multi=0
  ENDIF 
end
