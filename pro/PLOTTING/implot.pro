;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:   
;    tvim2
;
; USEAGE:   tvim2,a
;
;           tvim2,a,title=,xtitle=,ytitle=,xrange=,yrange=,subtitle=,$
;                   scale=,range=,/noframe,aspect=
;
; PURPOSE:  
;    Display an image.
;            1. numbered color scale 
;            2. plot title
;            3. annotated x and y axis 
;            4. simplified OPLOT capability
;
; INPUT    a           image quantity
;
; Optional keyword input:
;
;	   **NOTE**: This routine now uses the _extra mechanism to pass
;		     keywords for put_clr_scl2.
;
;          title       plot title
;
;          xtitle      x axis title
; 
;          ytitle      y axis title
;
;	   subtitle    x axis subtitle
;
;          xrange      array spanning entire x axis range.  
;                      NOTE:  TVIM2 only uses XRANGE(0) and
;				XRANGE(N_ELEMENTS(XRANGE)-1).
;
;          yrange      array spanning entire y axis range.  
;                      NOTE:  TVIM2 only uses YRANGE(0) and
;				YRANGE(N_ELEMENTS(YRANGE)-1).
;
;          scale       if set draw color scale.  SCALE=2 causes steped
;                      color scale
;
;          range       two or three element vector indicating physical
;                      range over which to map the color scale.  The
;                      third element of RANGE, if specified, sets the
;                      step interval of the displayed color scale.  It
;                      has no effect when SCALE is not set. E.g.,
;                      RANGE=[0., 1., 0.1] will cause the entire color
;                      scale to be mapped between the physical values
;                      of zero and one; the step size of the displayed 
;                      color scale will be set to 0.1.
;
;          aspect      the x over y aspect ratio of the output image
;                      if not set aspect is set to (size_x/size_y) of the
;                      input image. If set it -1, no aspect ratio is forced. 
;
;          /noCenter:  if set, don't center the display
;
;          noframe     if set do not draw axis box around image
;
;	   nolabels    If set, inhibits labels on plot and scale.
;	  
;	   invbw       To invert the color table ie. array=!d.n_colors-1-array
;		       before calling tv. The deault is invbw=1 for the
;		       postscript device and invbw=0 else.		
;
;
; SIDE EFFECTS:        Setting SCALE=2 changes the color scale using the
;                      STEP_CT procedure.  The color scale may be returned
;                      to its original state by the command, STEP_CT,/OFF
;
; PROCEDURE            TVIM first determins the size of the plot data window
;                      with a dummy call to PLOT.  When the output device is
;                      "X", CONGRID is used to scale the image to the size of
;                      the available data window.  Otherwise, if the output
;                      device is Postscript, scaleable pixels are used in the
;                      call to TV.  PUT_COLOR_SCALE draws the color scale and
;                      PLOT draws the x and y axis and titles.
;
; DEPENDENCIES:        PUT_COLOR_SCALE, STEP_CT
;
; AUTHOR:              Paul Ricchiazzi    oct92 
;                      Earth Space Research Group, UCSB
;
;                      Modified version of TVIM by Jeff Bloch, SST-9, LANL
;
;			           1.12	10/13/95
;		               Added invbw keyword to invert the color table. 
;		               The default is invbw=1 for postscript and invbw=0 else.
;		               Changed a line to use bytscl function.
;				              David Johnston UChicago July 1999
;                      Centered up the display. Comments, Cleaned up. 
;                             Erin Scott Sheldon UMich 11/24/99
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function implot_asinh_scale, image, alpha=alpha, nonlinearity=nonlinearity
    if n_elements(nonlinearity) eq 0 then nonlinearity = 8
    if n_elements(alpha) eq 0 then alpha = 0.02
 
    image_out = asinh( alpha*nonlinearity*image )/nonlinearity
    image_out = bytscl( temporary(image_out) )
    return, image_out

end 

function implot_byte_scale, image, range=range, max_color=max_color
    if n_elements(max_color) eq 0 then max_color=!d.n_colors-1

    amin = min(image, max=amax)
    amin = float(amin) & amax = float(amax)
    if amin eq amax then amax = amin+1
    if n_elements(range) eq 0 then range=[amin, amax]
    if range(0) eq range(1) then range(1)=range(1)+1

    image_scaled = bytscl(image,min=range(0),max=range(1),top=max_color)
  
    return, image_scaled 

end

function implot_scale, image, type=type, $
    alpha=alpha, nonlinearity=nonlinearity,$ ; asinh parameters
    range=range, max_color=max_color ; regular parameters

    if n_elements(type) eq 0 then type='byte'
    case strlowcase(type) of
        'asinh':begin
            imout=implot_asinh_scale(image, alpha=alpha, $
                                     nonlinearity=nonlinearity)
        end
        'byte': begin
            imout=implot_byte_scale(image, range=range, max_color=max_color)
        end
        else: message,'Unknown type: '+string(type)
    endcase
    return, imout
end

PRO implot, image, $
    type=type, $
    alpha=alpha, nonlinearity=nonlinearity, $
    scale=scale, range=range, xrange=xrange, yrange=yrange, $
    aspect=aspect,$
    image_scaled=image_scaled, $
    title=title, xtitle=xtitle, ytitle=ytitle, subtitle=subtitle, $
    noframe=noframe, nolabels=nolabels,$
    invbw=invbw, max_color=max_color, position=pos, $
    noCenter=noCenter, $
    _extra=extra_key


    if n_params() eq 0 then begin
      print,'-Syntax: implot, image, '
      print,'   scale=scale, range=range, xrange=xrange, yrange=yrange, '
      print,'   aspect=aspect,'
      print,'   title=title, xtitle=xtitle, ytitle=ytitle, subtitle=subtitle, '
      print,'   noframe=noframe, nolabels=nolabels,'
      print,'   invbw=invbw, '
      print,'   _extra=extra_key'
      print
      print,'For more help use doc_library,"tvim2"'
      return
    endif 

    if n_elements(invbw) eq 0 then begin 
        ; if device is postscript, reverse color by default
        if !d.name eq 'PS' then invbw = 1 else invbw = 0
    endif 

    if n_elements(max_color) eq 0 then max_color=!d.n_colors-1
    if n_elements(title) eq 0 then title = ''
    if n_elements(xtitle) eq 0 then xtitle=''
    if n_elements(ytitle) eq 0 then ytitle=''
    if n_elements(subtitle) eq 0 then subtitle=''
    if keyword_set(scale) then doscale = 1 else doscale=0

    sz=size(image)
    nx=sz(1)
    ny=sz(2)

    nxm=nx-1
    nym=ny-1

    if n_elements(xrange) eq 0 then begin
        xrng=[ -0.5, nxm+0.5]
    endif else begin
        xrng=[xrange(0), xrange(n_elements(xrange)-1)]
    endelse 

    if n_elements(yrange) eq 0 then begin
        yrng = [-0.5,nym+0.5]
    endif else begin
        yrng = [yrange(0), yrange(n_elements(yrange)-1)]
    endelse 

    ; normalized position
    pos = plotposition()

    ; Size of region in normalized coords
    xsize = pos[2] - pos[0]
    ysize = pos[3] - pos[1]

    image_scaled = $
        implot_scale(image, type=type, $
                     range=range, max_color=max_color, $
                     alpha=alpha, nonlinearity=nonlinearity)
    if keyword_set(invbw) then image_scaled = max_color-image_scaled

    if !d.name eq 'X' then begin ;x window
        ; number of pixels in each direction
        xwinsize = !d.x_vsize
        ywinsize = !d.y_vsize
        ; now times size in normalized coords gives number of pixels
        ; used on the device
        new_nx = xsize*xwinsize
        new_ny = ysize*ywinsize

        tv, congrid(image_scaled,new_nx,new_ny), pos[0], pos[1], /norm
        ;pos = [px(0), py(0), px(1), py(1)]
    endif else begin 
        ;pos = [px(0), py(0), px(1), py(1)] ;postscript
        tv, image_scaled, pos[0], pos[1], xsize=xsize, ysize=ysize, /norm
    endelse 

    if keyword_set(noframe) or keyword_set(nolabels) then begin 
        plot, [0,0], [0,0], xstyle=5, ystyle=5, $
            title=title,xtitle=xtitle,ytitle=ytitle, subtitle=subtitle, $
            xrange=xrng, yrange=yrng, position=pos, /noerase, /nodata,$
            _extra=extra_key
    endif else begin 
        plot, [0,0], [0,0], xstyle=1, ystyle=1, $
            title=title, xtitle=xtitle, ytitle=ytitle, subtitle=subtitle, $
            xrange=xrng, yrange=yrng, position=pos, /noerase, /nodata,$
            _extra=extra_key
    endelse 
  
    if (not keyword_set(noframe)) and keyword_set(nolabels) then begin 
        axis,xaxis=1,xtickname=strarr(10)+" ",$
            _extra=extra_key
        axis,xaxis=0,xtickname=strarr(10)+" ",$
            _extra=extra_key
        axis,yaxis=1,ytickname=strarr(10)+" ",$
            _extra=extra_key
        axis,yaxis=0,ytickname=strarr(10)+" ",$
            _extra=extra_key
    endif 

end 





