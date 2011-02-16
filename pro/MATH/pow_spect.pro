pro pow_spect,f,step,direction,spect,freq,plot=plot

  if n_params() eq 0 then begin
    print,'-Syntax: pow_spect,f,step,direction,spect,freq,yesplot=yesplot'
    print,' direction = -1 for transform, +1 for inverse'
    print,' Computes power spectrum (spect) for function f from its fft.  Plotted with most negative first.'  
    return
  endif

  n = n_elements(f)
  n21 = n/2 + 1

  trans = fft(f,direction)
  spect = abs(trans)

  freq = indgen(n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;insert negative frequencies from n/2+1....n-1 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  freq[n21] = n21 - n + findgen(n21 -2)
  
  freq = freq/(step*n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;using shift function to get most negative frequencies first in case
;function is not nice(should wrap around automatically, thats the point
;of doing above)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if keyword_set(plot) then begin
    mn = min(shift(freq,-n21))
    mx = -mn
    plot,/ylog, shift(freq,-n21), shift(spect, -n21), psym=3,$
      title='Power spectrum',xtitle='freq  row^-1',$
      xrange=[mn,mx]
  endif

  return
end


