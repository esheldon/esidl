PRO charsets

; This MAIN_PROGRAM creates a printed Character Set list for a Postscript
; printer.  This program requires the 'setps' and 'psclose' programs,
; although it can be modified to just use set_plot,'ps' and then
; appropriate scaling commands (which is what setps does.)
; Send comments or changes to: deutsch@astro.washington.edu (Eric Deutsch)
;
; Simply run with:
;    IDL> .run charsets


  setps,8,10.5,.25,.25          ; Set 8"x10" Portrait plot window
  !p.font=0                     ; harware fonts
  
  fontnum=2                     ; select a font to print
; below are some examples of 4 interesting fonts.  Change or add the ones
; you want to print.
  
  if (fontnum eq 1) then begin 
      device,/helv,/isolatin1
      title='Helvetica ISOLatin1'
  endif
    
  if (fontnum eq 2) then begin 
      device,/symbol 
      title='Symbol' 
  endif
    
  if (fontnum eq 3) then BEGIN
      device,/times,/isolatin1
      title='Times Roman ISOLatin1' 
  endif
    
  if (fontnum eq 4) then begin 
      device,/palatino,/italic,/isolatin1
      title='Palatino Italic ISOLatin1'
  endif
    
  plot,[0,60],[0,60],/nodata,xsty=7,ysty=7,$
    position=[0,0,1,1]          ; set up coordinates
  xyouts,30,58,/data,align=.5,title ; print title
    
  col=0
  row=0
  for i=0,255 do begin          ; loop over all characters
      xyouts,col+5,55-row,strn(i)+'= '+string(byte(i)),/data
      row=row+1
      if (row eq 52) then begin
          col=col+11
          row=0
      endif
  endfor
  
  psclose                       ; close Postscript channel
    
end








