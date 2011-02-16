pro fukugita_k_corrections_array,gmr,z,kcorr,type=type,plot=plot

   if (n_params(0) eq 0) then begin
	print,'Syntax: fukugita_k_corrections_array,gmr,z,kcorr,type=type,plot=plot'
	return
   endif

   ;Make an array in type redshift space which describes the k corrections
   karray=mrdfits('/sdss3/data5/mckay/fukugita_k_corrections_array.fit',0)
   ;This array has three indices:
   ;	First: Color u,g,r,i,z
   ;	Second: Redshift bin
   ;	Third: z, then six types (E,SO,Sab,Sbc,Scd,Im)
   ;
   ;	For example karray(2,*,0) is the redshift values for the r data
   ;		    karray(2,*,2) is the r' k corrections for SO galaxies


   ;First set up the color information required to type the galaxies
   gmr_color=karray(1,*,1:6)-karray(2,*,1:6)
   zval=karray(1,*,0)
   gmr_color=reform(gmr_color)
   ;Now correct for zero redshift color
   gmr_color(*,0)=gmr_color(*,0)+0.77
   gmr_color(*,1)=gmr_color(*,1)+0.68
   gmr_color(*,2)=gmr_color(*,2)+0.66
   gmr_color(*,3)=gmr_color(*,3)+0.52
   gmr_color(*,4)=gmr_color(*,4)+0.48
   gmr_color(*,5)=gmr_color(*,5)+0.20

   ;First use color, and redshift to pick type
   nobj=n_elements(gmr)
   type=fltarr(nobj)
   for i=0,nobj-1 do begin
	;Find color of each type at this redshift
	ce=interpol(gmr_color(*,0),zval,z(i))
	cso=interpol(gmr_color(*,1),zval,z(i))
	csab=interpol(gmr_color(*,2),zval,z(i))
	csbc=interpol(gmr_color(*,3),zval,z(i))
	cscd=interpol(gmr_color(*,4),zval,z(i))
	cim=interpol(gmr_color(*,5),zval,z(i))
	carr=[ce,cso,csab,csbc,cscd,cim]
	if (gmr(i) ne 0.0) then begin
	  type(i)=interpol(findgen(6),carr,gmr(i))
	endif else begin
	  type(i)=2.0
        endelse
	if (type(i) lt 0.) then type(i)=0.
	if (type(i) gt 5.) then type(i)=5.
   endfor

   if keyword_set(plot) then begin
     plot,zval,gmr_color(*,0),linestyle=1
     oplot,zval,gmr_color(*,1),linestyle=1
     oplot,zval,gmr_color(*,2),linestyle=1
     oplot,zval,gmr_color(*,3),linestyle=1
     oplot,zval,gmr_color(*,4),linestyle=1
     oplot,zval,gmr_color(*,5),linestyle=1
     oplot,z,gmr,psym=3
     res=get_kbrd(10)
   endif

   kcorr=fltarr(5,nobj)
   kcorr(0,*)=interpolate(reform(karray(0,*,1:6)),z/0.05,type)
   kcorr(1,*)=interpolate(reform(karray(1,*,1:6)),z/0.05,type)
   kcorr(2,*)=interpolate(reform(karray(2,*,1:6)),z/0.05,type)
   kcorr(3,*)=interpolate(reform(karray(3,*,1:6)),z/0.05,type)
   kcorr(4,*)=interpolate(reform(karray(4,*,1:6)),z/0.05,type)
   
   if keyword_set(plot) then begin
     ;Now make a plot to check the results
     plot,karray(2,*,0),karray(0,*,1),linestyle=1,yrange=[-0.0,1.5],$
	xrange=[0,0.5],ystyle=1,xstyle=1,$
	title='K corrections for u band'
     oplot,karray(2,*,0),karray(0,*,2),linestyle=1
     oplot,karray(2,*,0),karray(0,*,3),linestyle=1
     oplot,karray(2,*,0),karray(0,*,4),linestyle=1
     oplot,karray(2,*,0),karray(0,*,5),linestyle=1
     oplot,karray(2,*,0),karray(0,*,6),linestyle=1
     oplot,z,kcorr(0,*),psym=3
     res=get_kbrd(10)

     plot,karray(2,*,0),karray(1,*,1),linestyle=1,yrange=[-0.2,1.0],$
	xrange=[0,0.5],ystyle=1,xstyle=1,$
	title='K corrections for g band'
     oplot,karray(2,*,0),karray(1,*,2),linestyle=1
     oplot,karray(2,*,0),karray(1,*,3),linestyle=1
     oplot,karray(2,*,0),karray(1,*,4),linestyle=1
     oplot,karray(2,*,0),karray(1,*,5),linestyle=1
     oplot,karray(2,*,0),karray(1,*,6),linestyle=1
     oplot,z,kcorr(1,*),psym=3
     res=get_kbrd(10)

     plot,karray(2,*,0),karray(2,*,1),linestyle=1,yrange=[-0.2,1.0],$
	xrange=[0,0.5],ystyle=1,xstyle=1,$
	title='K corrections for r band'
     oplot,karray(2,*,0),karray(2,*,2),linestyle=1
     oplot,karray(2,*,0),karray(2,*,3),linestyle=1
     oplot,karray(2,*,0),karray(2,*,4),linestyle=1
     oplot,karray(2,*,0),karray(2,*,5),linestyle=1
     oplot,karray(2,*,0),karray(2,*,6),linestyle=1
     oplot,z,kcorr(2,*),psym=3
     res=get_kbrd(10)

     plot,karray(2,*,0),karray(3,*,1),linestyle=1,yrange=[-0.2,0.6],$
	xrange=[0,0.5],ystyle=1,xstyle=1,$
	title='K corrections for i band'
     oplot,karray(2,*,0),karray(3,*,2),linestyle=1
     oplot,karray(2,*,0),karray(3,*,3),linestyle=1
     oplot,karray(2,*,0),karray(3,*,4),linestyle=1
     oplot,karray(2,*,0),karray(3,*,5),linestyle=1
     oplot,karray(2,*,0),karray(3,*,6),linestyle=1
     oplot,z,kcorr(3,*),psym=3
     res=get_kbrd(10)

     plot,karray(2,*,0),karray(4,*,1),linestyle=1,yrange=[-0.2,0.6],$
	xrange=[0,0.5],ystyle=1,xstyle=1,$
	title='K corrections for z band'
     oplot,karray(2,*,0),karray(4,*,2),linestyle=1
     oplot,karray(2,*,0),karray(4,*,3),linestyle=1
     oplot,karray(2,*,0),karray(4,*,4),linestyle=1
     oplot,karray(2,*,0),karray(4,*,5),linestyle=1
     oplot,karray(2,*,0),karray(4,*,6),linestyle=1
     oplot,z,kcorr(4,*),psym=3
     res=get_kbrd(10)

   endif
   
   kcorr=-1.0*kcorr

   return
   end