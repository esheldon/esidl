pro fukugita_k_corrections,gmr,z,kcorr,type=type

   ;First set up the information
   type=findgen(2,6)
   color=findgen(2,6)
   color(0,*)=[0.77,0.68,0.66,0.52,0.48,0.20]	
   color(1,*)=[1.31,1.13,1.02,0.71,0.62,0.32]
   slope_array=(color(1,*)-color(0,*))/0.2

   ;First use color, and redshift to pick type
   nobj=n_elements(gmr)
   type=fltarr(nobj)
   for i=0,nobj-1 do begin
	;Find color of each type at this redshift
	carr=reform(color(0,*)+z(i)*slope_array)
	type(i)=interpol(findgen(6),carr,gmr(i))
	;if (type(i) lt 0.) then type(i)=0.
	;if (type(i) gt 5.) then type(i)=5.
   endfor

   ;Make an array in type redshift space which describes the k corrections
   karray=fltarr(6,40)
   za=findgen(40)*0.01
   karray(0,*)=(0.6/0.4)*za
   karray(1,*)=(0.55/0.4)*za
   karray(2,*)=(0.5/0.4)*za
   karray(3,*)=(0.3/0.4)*za
   karray(4,*)=(0.15/0.4)*za
   karray(5,*)=(-0.15/0.4)*za

   kcorr=interpolate(karray,type,z/0.01)   
   kcorr=-1.0*kcorr

   return
   end
