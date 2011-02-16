PRO plot_fit_wtheta


  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
  ytitle='Overdensity (# Mpc!U'+!tsym.minus+'2!N)'

  fit_wtheta_tot, tmeanr, tdiffdensecomb, terrorcomb, tcumul, txx, tyycomb
  fit_wtheta_high, hmeanr, hdiffdensecomb, herrorcomb, hcumul, hxx, hyycomb
  fit_wtheta_low, lmeanr, ldiffdensecomb, lerrorcomb, lcumul, lxx, lyycomb

  outdir='/sdss5/data0/lensout/wtheta_conv/'
  psfile=outdir+'wtheta_all.ps'
  begplot,name=psfile

;  yrange=prange(hdiffdensecomb,ldiffdensecomb,$
;                herrorcomb,lerrorcomb)
  erase & multiplot,[1,3]

  hxx=hxx*1000. & txx=txx*1000.

  ploterror, tmeanr, tdiffdensecomb, terrorcomb, psym=4;,yrange=yrange
  oplot,txx,tyycomb
  oplot,tmeanr,tcumul,line=2
  oplot,[0,10000],[0,0]
  legend,'All Galaxies',/right,box=0
  multiplot

  ploterror, hmeanr, hdiffdensecomb, herrorcomb, psym=4;,yrange=yrange
  oplot,hxx,hyycomb
  oplot,hmeanr,hcumul,line=2
  oplot,[0,10000],[0,0]
  legend,'High Density Regions',/right,box=0
  multiplot

  ploterror, lmeanr, ldiffdensecomb, lerrorcomb, psym=4,$
    ytitle=ytitle,xtitle=xtitle;, yrange=yrange
  oplot,lxx,lyycomb
  oplot,[0,10000],[0,0]
  legend,'Low Density Regions',/left,box=0
  multiplot,/reset

  key=get_kbrd(1)



  aploterror, !gratio, hmeanr, hdiffdensecomb, herrorcomb, psym=4,yrange=yrange,$
    ytitle=ytitle,xtitle=xtitle
;  oplot,hmeanr,hcumul,line=2
  oplot,hxx,hyycomb
  oploterror, tmeanr, tdiffdensecomb, terrorcomb, psym=5
;  oplot,tmeanr,tcumul,line=2
  oplot,txx,tyycomb
  oploterror, lmeanr, ldiffdensecomb, lerrorcomb, psym=1
  oplot,lxx,lyycomb
  oplot,[0,10000],[0,0]

  legend,['All Galaxies','High Density Regions',$
          'Low Density Regions'],/right,psym=[5,4,1],thick=replicate(!p.thick,3)

  endplot

END 
