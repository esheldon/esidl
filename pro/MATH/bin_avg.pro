pro  bin_avg, x, y, nbins, mean_x, mean_y, median_y, bin_sig

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:
;       
; PURPOSE:
;	
;
; CALLING SEQUENCE:
;      
;                 
;
; INPUTS: 
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUT ARRAYS:
;
; INPUT KEYWORD PARAMETERS:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
	print,'Syntax: bin_avg, x, y, nbins, mean_x, mean_y, median_y, bin_sig'
	return
  endif

;;;;;;; cast nbins into integer.  fix will round down  ;;;;;;;;;
nbins = fix(nbins)

xmin = min(x)
xmax = max(x)
binsize = (xmax - xmin)/nbins

;;;;; find average and standard deviation in each bin ;;;;;;;;


mean_x = fltarr(nbins)
mean_y = fltarr(nbins)
median_y = fltarr(nbins)
bin_sig = fltarr(nbins)

for i=0,nbins-1 do begin
  min = i*binsize + xmin
  max = min + binsize
  w = where( x ge min and x le max)
  x_mom = moment(x[w])
  y_mom = moment(y[w])
  even = 1
  if ( (n_elements(w) mod 2) ) then begin
	even = 0
  endif
  y_med = median(y[w], even=even)

  mean_x(i) = x_mom(0)
  mean_y(i) = y_mom(0)
  median_y(i) = y_med
  bin_sig(i) = sqrt( y_mom(1) )
endfor

plot,x,y,psym=3
oplot,mean_x,mean_y,psym=1
errplot, mean_x, mean_y-bin_sig, mean_y+bin_sig

oplot, mean_x, median_y,psym=2
legend,['Mean', 'Median'],psym=[1,2]
return
end


























