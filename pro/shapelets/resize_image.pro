PRO resize_image, image, x0, sz

; Author: Brandon Kelly June 2002
; Routine to resize and image, centered around a point x0, with a given
; size (sz)

if n_params() lt 3 then begin
    print, 'Syntax- resize_image, image, x0, sz'
    return
endif

sz1=sz[0]
sz2=sz[1]
nx=size(image)

IF x0[0]-sz1/2 LE 0 THEN BEGIN
    newx=fix(sz1/2-x0[0])
    image1=fltarr(nx[1]+newx, nx[2])
    if n_elements(image[*,0]) ne nx[1] or $
      n_elements(image[0,*]) ne nx[2] then begin
        print, 'Size of Image does not match induces for Image 1, returning...'
        return
    endif else image1[newx : nx[1]+newx-1, 0 : nx[2]-1]=image
ENDIF
IF x0[0]-sz1/2 GT 0 THEN BEGIN
    newx=fix(x0[0]-sz1/2)
    if n_elements(image[*,0]) lt nx[1] or $
      n_elements(image[0,*]) lt nx[2] then begin
        print, 'Size of Image is too small compared to Image 1, returning...'
        return
    endif else image1=image[newx : nx[1]-1, 0:nx[2]-1]
ENDIF
IF nx[1]-1-x0[0]-sz1/2 LE 0 THEN BEGIN
    nx1=size(image1)
    newx=fix(sz1/2-(nx[1]-1-x0[0]))
    image2=fltarr(nx1[1]+newx,nx[2])
    if n_elements(image1[*,0]) ne nx1[1] or $
      n_elements(image1[0,*]) ne nx1[2] then begin
        print, 'Size of Image 1 does not match induces for Image 2, returning...'
        return
    endif else image2[0:nx1[1]-1, 0:nx1[2]-1]=image1
ENDIF
IF nx[1]-1-x0[0]-sz1/2 GT 0 THEN BEGIN
    nx1=size(image1)
    newx=fix(nx[1]-1-x0[0]-sz1/2)
    if n_elements(image1[*,0]) lt nx1[1]-newx or $
      n_elements(image1[0,*]) lt nx1[2] then begin
        print, 'Size of Image 1 is too small compared to Image 2, returning...'
        return
    endif else image2=image1[0:nx1[1]-newx-1,0:nx1[2]-1]
ENDIF
IF x0[1]-sz2/2 LE 0 THEN BEGIN
    nx2=size(image2)
    newy=fix(sz2/2-x0[1])
    image3=fltarr(nx2[1], nx2[2]+newy)
    if n_elements(image2[*,0]) ne nx2[1] or $
      n_elements(image2[0,*]) ne nx2[2] then begin
        print, 'Size of Image 2 does not match induces for Image 3, returning...'
        return
    endif else image3[0 : nx2[1]-1, newy : newy+nx2[2]-1]=image2
ENDIF
IF x0[1]-sz2/2 GT 0 THEN BEGIN
    nx2=size(image2)
    newy=fix(x0[1]-sz2/2)
    if n_elements(image2[*,0]) lt nx2[1] or $
      n_elements(image2[0,*]) lt nx2[2] then begin
        print, 'Size of Image 2 is too small compared to Image 3, returning...'
        return
    endif else image3=image2[0 : nx2[1]-1, newy : nx2[2]-1]
ENDIF
IF nx[2]-1-x0[1]-sz2/2 LE 0 THEN BEGIN
    nx3=size(image3)
    newy=fix(sz2/2-(nx[2]-1-x0[1]))
    image4=fltarr(nx3[1],nx3[2]+newy)
    if n_elements(image3[*,0]) ne nx3[1] or $
      n_elements(image3[0,*]) ne nx3[2] then begin
        print, 'Size of Image 3 does not match induces for Image 4, returning...'
        return
    endif else image4[0:nx3[1]-1, 0:nx3[2]-1]=image3
ENDIF
IF nx[2]-1-x0[1]-sz2/2 GT 0 THEN BEGIN
    nx3=size(image3)
    newy=fix(nx[2]-1-x0[1]-sz2/2)
    if n_elements(image3[*,0]) lt nx3[1]-1 or $
      n_elements(image3[0,*]) lt nx3[2]-newy-1 then begin
        print, 'Size of Image 3 is too small compared to Image 4, returning...'
        return
    endif else image4=image3[0:nx3[1]-1, 0:nx3[2]-newy-1]
ENDIF
image=image4

return
END



