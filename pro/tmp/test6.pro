PRO test6, admomstruct, psf, get=get


  ;; compare psf measurement using old way and 
  ;; way where we measure moments using weight of galaxy

  indir='/sdss5/data2/corrected/corr1033/objcs/1/'


  IF keyword_set(get) THEN BEGIN 

      setzero,admomstruct,psf

      ;; read in admom info
      fetch_dir,1033,1,1,dir,corrdir=corrdir
      admomatlas_infile, 1033, 1, 1, fname_in, fname_out
      fname_out = repstr(fname_out, '.dat', '.fit')
      admomstruct = mrdfits(corrdir + fname_out, 1)

      FOR i=65L,131L DO BEGIN 
          
          istr=field2string(i)
          t2=mrdfits(indir+'psffit-001033-1-1-'+istr+'.fit_test',1,/silent)
          
          IF i EQ 65 THEN BEGIN 
              psf = t2
          ENDIF ELSE BEGIN 
              concat_structs, temporary(psf), t2, psf
          ENDELSE 
          
          setzero,t2
      ENDFOR 
  ENDIF 
  sphoto_match, admomstruct, psf, mad, mnew

  clr=2
  admomsize=admomstruct[mad].ixx[clr]+admomstruct[mad].iyy[clr]
  e1admom = (admomstruct[mad].ixx[clr] - admomstruct[mad].iyy[clr])/admomsize
  e2admom = 2.*admomstruct[mad].ixy[clr]/admomsize
  
  newsize=psf[mnew].psfixx[clr]+psf[mnew].psfiyy[clr]
  e1new = (psf[mnew].psfixx[clr] - psf[mnew].psfiyy[clr])/newsize
  e2new = 2.*psf[mnew].psfixy[clr]/newsize

  xx=arrscl(findgen(100), -5., 10.)
  setupplot, dtype, /true
;  !p.multi=[0,0,2]
;  myusersym,'fill_circle'
;  aplot,1,e1admom,e1new,psym=3,symsize=0.5,$
;    xtitle='!7psf_e1  admom',ytitle='!7psf_e1  new',$
;    yrange=[-0.2,0.1],ystyle=1,title='!7Round weight',/center
;  oplot,xx,xx,color=!red
;  aplot,1,e2admom,e2new,psym=3,symsize=0.5,$
;    xtitle='psf_e2  admom',ytitle='psf_e2  new',$
;    yrange=[-0.04,0.1],ystyle=1,/center
;  oplot,xx,xx,color=!red

  range=[-1,1]
  histogram_2d,e1admom,e1new,hist,range,range
  IF dtype EQ 'X' THEN BEGIN
      padmom = !p.charsize & !p.charsize=2.0
  ENDIF 
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='!7e!D1!N  galaxy',ytitle=,title='!7Galaxy Weighted',/invbw,$
    range=[0,11.5]
  legend,['Round Input PSF'],/left,box=0
  ;oplot,xx,xx
  key=get_kbrd(1)
  histogram_2d,e2admom,e2new,hist,range,range
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='!7e!D2!N  galaxy',ytitle='!7e!D2!N PSF measured',title='Galaxy Weighted',/invbw,$
    range=[0,28.8]
  legend,['Round Input PSF'],/left,box=0
  ;oplot,xx,xx
 
  !p.multi=0
return
  key=get_kbrd(1)

;  plot,admomsize, newsize, psym=3,symsize=0.5,$
;    yrange=[2.0,6.0],xtitle='(ixx+iyy)_admom',$
;    ytitle='(ixx+iyy)_new',title='Round Weight'
;  oplot,xx,xx,color=!red

  histogram_2d, admomsize, newsize, hist, [0.0,8.0],[0.0,8.0]
  rdis, hist.map, xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='(ixx+iyy)  galaxy',$
    ytitle='(ixx+iyy) PSF  measured',title='Galaxy Weight',/invbw,$
    range=[0,28.8]
  legend,['Input PSF (ixx+iyy)=3.0'],/left,box=0
  ;oplot,xx,xx
  IF dtype EQ 'X' THEN !p.charsize=padmom

END 
