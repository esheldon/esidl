pro binner, xi,yi,xo,yo,sig,bin
  
  histo=histogram(xi,bin=0.01,reverse_indices=r)
  nn=n_elements(histo)
  xo=fltarr(nn)
  yo=fltarr(nn)
  sig=fltarr(nn)
  kk=0
  for i=0,nn-1 do begin
    if(r(i+1)-r(i) gt 2) then begin
      result=moment(yi(r(r(i):r(i+1)-1)))
      yo[kk]=result[0]
      sig[kk]=sqrt(result[1]/(r(i+1)-r(i)))
      result=moment(xi(r(r(i):r(i+1)-1)))
      xo[kk]=result[0]
;    print,xo[kk],yo[kk],sig[kk],r(i+1)-r(i)
      kk=kk+1
    endif
  endfor
  yo=yo[0:kk-1]
  xo=xo[0:kk-1]
  sig=sig[0:kk-1]
  
  return
  
end


PRO compare_testbed, clr, camcol, orig, new, matcho, matchn, notot, matchtot, nomatch

  ;; color for testing
  colors=['u','g','r','i','z']
  ;clr=2

  run=745
  rerun=1
  field0=395L
  fieldL=514L
  nf=fieldL-field0+1
  ;camcol=3
  fetch_dir,run,camcol,rerun,dir,corrdir=corrdir

  indir='/sdss3/data6/745/testbed/'+ntostr(camcol)+'/'

  arrval=fltarr(5)
  newstruct=create_struct('field',0L,$
                          'id',0L,$
                          'colc',arrval,$
                          'rowc',arrval,$
                          'mcc',arrval,$
                          'mccerr',arrval,$
                          'mrr',arrval,$
                          'mrrerr',arrval,$
                          'mcr',arrval,$
                          'mcrerr',arrval,$
                          'mcc_psf',arrval,$
                          'mrr_psf',arrval,$
                          'mcr_psf',arrval,$
                          'mcr4',arrval,$
                          'mcr4_psf',arrval,$
                          'flags2', lonarr(5),$
                          'type',lonarr(5))

  nomatchstruct = create_struct('run', 0L, $
                                'camcol',0,$
                                'field',0L,$
                                'rowc', fltarr(5), $
                                'colc', fltarr(5))
  
  ptrlist = ptrarr(nf)
  numlist = lonarr(nf)
  ntotal = 0L

  IF n_elements(orig) EQ 0 THEN BEGIN 
      tag=['field','id','rowc','colc','r','ixx','iyy',$
           'ixy','e1','e2','rho4','momerr','petrocounts']
      read_tsobj,corrdir,orig,start=field0,nframes=nf,tag=tag,front='adatc'
      
      ff=0
      FOR i=field0,fieldL DO BEGIN 
          
          file=indir+'fpObjc-000745-'+ntostr(camcol)+'-'+field2string(i)+'.fit'
          openr,lun, file, /get_lun
          IF i EQ field0 THEN BEGIN 
              tmp=mrdfits3(lun, 1, 0, /silent)
          ENDIF ELSE BEGIN 
              tmp=mrdfits3(lun, 1, 0, /silent, /deja_vu)
          ENDELSE 
          free_lun,lun
          print,'.',format='(a,$)'

          ntmp = n_elements(tmp)
          newtmp = replicate(newstruct, ntmp )
          copy_struct, tmp, newtmp
          newtmp.field =i

          ptrlist[ff] = ptr_new(newtmp)
          numlist[ff] = ntmp
          ntotal=ntotal+ntmp
          ff=ff+1
          setzero, tmp, newtmp

      ENDFOR 

      new = replicate(newstruct, ntotal)
      beg=0L
      FOR fi=0L, nf-1 DO BEGIN 
          IF (numlist[fi] NE 0) THEN BEGIN 
              new[beg:beg+numlist[fi]-1] = *ptrlist[fi]
          ENDIF 
          ptr_free, ptrlist[fi]
          beg=beg+numlist[fi]
      ENDFOR 

  ENDIF 

  

  IF n_elements(matcho) EQ 0 THEN BEGIN 

      notot=0L
      matchtot=0L
      print
      print,'matching'
      print
      
      AMOMENT_FAINT = '200000'X
      AMOMENT_SHIFT = '400000'X
      AMOMENT_MAXITER = '800000'X

      delvarx,matcho,matchn,nomatch
      FOR i=field0,fieldL DO BEGIN 
          

          wo=where( (orig.field EQ i) AND (orig.ixx[clr] GT 0.0) AND $
                    (orig.petrocounts[clr] LT 22.),no)
          wn2 = where(new.field EQ i,nn2)
          wn=where( (new[wn2].field EQ i) AND $
                    ( (new[wn2].flags2[clr] AND AMOMENT_FAINT) EQ 0) AND $
                    ( (new[wn2].flags2[clr] AND AMOMENT_SHIFT) EQ 0) AND $
                    ( (new[wn2].flags2[clr] AND AMOMENT_MAXITER) EQ 0) ,nn1)
          wn2_keep = wn2
          remove, wn, wn2_keep
          wn = wn2[wn]
          wn2 = wn2_keep

          close_match, orig[wo].colc[clr], orig[wo].rowc[clr],$
            new[wn].colc[clr],  new[wn].rowc[clr],$
            mo1, mn1, 1.5, 1, /silent
          
          nmatch = n_elements(mn1)
          IF mo1[0] NE -1 THEN BEGIN 
              mo = wo[mo1]
              mn = wn[mn1]
              add_arrval, mo, matcho
              add_arrval, mn, matchn
          ENDIF ELSE nmatch=0

          ;; See what else matches
          IF (nmatch NE no) AND (nmatch NE 0) THEN remove, mo1, wo

          close_match, orig[wo].colc[clr], orig[wo].rowc[clr],$
            new[wn2].colc[clr],  new[wn2].rowc[clr],$
            mo2, mn2, 1.5, 1, /silent

          IF mn2[0] NE -1 THEN BEGIN 
              mn2 = wn2[mn2]
              nmatch2 = n_elements(mn2)
              nomatchtmp=replicate(nomatchstruct, nmatch2)

              copy_struct, new[mn2], nomatchtmp
              nomatchtmp.camcol = camcol
              nomatchtmp.run = run

              IF n_elements(nomatch) EQ 0 THEN nomatch=temporary(nomatchtmp) $
              ELSE concat_structs, temporary(nomatch), temporary(nomatchtmp), $
                nomatch
              
          ENDIF 
          
          print,'-'+ntostr(nmatch)+'/'+ntostr(no),format='(a,$)'
          notot=notot+no
          matchtot = matchtot+nmatch

      ENDFOR 
  ENDIF
 
  print
  print
  print,'Total in orig = ',notot
  print,'Total matched = ',matchtot
  print,'Fraction matched = ',matchtot/float(notot)

  ;; compare shape measurements

  ;; generate e1,e2,momerr for testbed

  dtype = !d.name


  tsize = new[matchn].mcc[clr] + new[matchn].mrr[clr]
  psfsize = new[matchn].mcc_psf[clr] + new[matchn].mrr_psf[clr]
  psfrho4 = new[matchn].mcr4_psf[clr]
  corr = (psfsize/tsize)*(4./psfrho4-1.)/(4./new[matchn].mcr4[clr]-1.)
  gal = where(corr GT 0.0 AND corr LT 0.8,ngal)
  galn = matchn[gal]
  galo = matcho[gal]


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get some stars. Compare shape of psf reconstruction to
  ;; shape of star, check error formula
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  tmpmn=matchn
  tmpmo = matcho
  remove,gal,tmpmn
  remove,gal,tmpmo
  star=where(new[tmpmn].type[2] EQ 6 AND new[tmpmn].type[1] EQ 6 AND $
             orig[tmpmo].petrocounts[clr] LT 19,nstar)
  star = tmpmn[star]
  starsize = new[star].mcc[clr] + new[star].mrr[clr]
  startop = (new[star].mcc[clr] - new[star].mrr[clr])
  starpsfsize = new[star].mcc_psf[clr] + new[star].mrr_psf[clr]
  stare1 = startop/starsize
  stare2 = 2.*new[star].mcr[clr]/starsize
  starpsfe1 = (new[star].mcc_psf[clr] - new[star].mrr_psf[clr])/starsize
  starpsfe2 = 2.*new[star].mcr_psf[clr]/starpsfsize
  stare1err = abs(stare1)*sqrt( new[star].mccerr[clr]^2 + new[star].mrrerr[clr]^2)*$
    sqrt( 1./starsize^2 + 1./startop^2 )

  stare2err = abs(stare2)*sqrt((new[star].mcrerr[clr]/new[star].mcr[clr])^2 + $
    (new[star].mccerr[clr]^2 + new[star].mrrerr[clr]^2)/tsize^2)


  ;; galaxies
  tsize = tsize[gal]
  psfsize = psfsize[gal]
  psfrho4=psfrho4[gal]
  corr = corr[gal]

  ttop = (new[galn].mcc[clr] - new[galn].mrr[clr])

  te1 = ttop/tsize
  te2 = 2.*new[galn].mcr[clr]/tsize

  te1err = abs(te1)*sqrt( new[galn].mccerr[clr]^2 + new[galn].mrrerr[clr]^2)*$
               sqrt( 1./tsize^2 + 1./ttop^2 )

  te2err = abs(te2)*sqrt((new[galn].mcrerr[clr]/new[galn].mcr[clr])^2 + $
                    (new[galn].mccerr[clr]^2 + new[galn].mrrerr[clr]^2)/tsize^2)

  ;; correct e1,e2


  psfe1 = (new[galn].mcc_psf[clr] - new[galn].mrr_psf[clr])/psfsize
  psfe2 = 2.*new[galn].mcr_psf[clr]/psfsize

  te1_corr = te1 - corr*psfe1
  te2_corr = te2 - corr*psfe2

  te1_corr2 = te1_corr/(1.-corr)
  te2_corr2 = te2_corr/(1.-corr)

  xx=arrscl( findgen(100), -100, 100 )

  !p.multi=[0,0,2]
  mapfac=25.

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; stare1-psfe1
  ;;;;;;;;;;;;;;;;;;;;;;;;

  plothist,stare1err,bin=.001

  plothist,stare1-starpsfe1,bin=.01
  key=get_kbrd(1)
!p.multi=0
  histogram_2d,stare1err,stare1-starpsfe1,hist,[0,.015],[-0.2,0.2]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='erre1',ytitle='e1star-e1psf',/invbw,range=range
  oplot, xx, xx
  key=get_kbrd(1)
  rotate_plot,hist.map
;return
  ;; ixx, iyy, ixy
  histogram_2d,orig[galo].ixx[clr],new[galn].mcc[clr],hist,[0,10],[0,10]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='I!Dxx!N  orig',ytitle='M!Dcc!N',/invbw,range=range
  oplot, xx, xx

  histogram_2d,orig[galo].iyy[clr],new[galn].mrr[clr],hist,[0,10],[0,10]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='I!Dyy!N  orig',ytitle='M!Drr!N',/invbw,range=range
  oplot, xx, xx

  if dtype eq 'X' then key=get_kbrd(1)
  histogram_2d,orig[galo].ixy[clr],new[galn].mcr[clr],hist,[-5,5],[-5,5]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='I!Dxy!N  orig',ytitle='M!Dcr!N',/invbw,range=range
  oplot, xx, xx

  ;; rho4
  histogram_2d,orig[galo].rho4[clr],new[galn].mcr4[clr],hist,[0,5],[0,5]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle=!tsym.rho+'!U4!N  orig!N',ytitle='M!Dcr!U4',/invbw,range=range
  oplot, xx, xx

  ;; e1/e2
  if dtype eq 'X' then key=get_kbrd(1)
  !p.multi=[0,0,2]
  histogram_2d,orig[galo].e1[clr],te1_corr,hist,[-1.0,1.0],[-1.0,1.0]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='e!D1!N orig',ytitle='e!D1!N new',/invbw,range=range
  oplot, xx, xx
  histogram_2d,orig[galo].e2[clr],te2_corr,hist,[-1.0,1.0],[-1.0,1.0]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='e!D2!N orig',ytitle='e!D2!N new',/invbw,range=range
  oplot, xx, xx

  ;; error
  if dtype eq 'X' then key=get_kbrd(1)
  histogram_2d,orig[galo].momerr[clr],te1err,hist,[0,0.4],[0,0.4]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='err old',ytitle='e!D1!N err new',/invbw,range=range
  oplot, xx, xx
  histogram_2d,orig[galo].momerr[clr],te2err,hist,[0,0.4],[0,0.4]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='err old',ytitle='e!D2!N err new',/invbw,range=range
  oplot, xx, xx

  ;;error*2
  if dtype eq 'X' then key=get_kbrd(1)
jump:
  histogram_2d,orig[galo].momerr[clr],2./sqrt(1.+2.*te1^2)*te1err,hist,[0,0.4],[0,0.4]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='err old',ytitle='2*e!D1!N err new',/invbw,range=range
  oplot, xx, xx
  histogram_2d,orig[galo].momerr[clr],2./sqrt(1.+2.*te2^2)*te2err,hist,[0,0.4],[0,0.4]
  range=[0,max(hist.map)/mapfac]
  rdis,hist.map,xrange=hist.xrange,yrange=hist.yrange,$
    xtitle='err old',ytitle='2*e!D2!N err new',/invbw,range=range
  oplot, xx, xx

  !p.multi=0
  !p.charsize=1.0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now e_corr vs psf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if dtype eq 'X' then key=get_kbrd(1)
  xrange=[-0.4,0.4]
  yrange=[-0.1,0.1]

;  erase & multiplot,[0,2,2,0,0],/square
;  ;; e1 vs psf_e1
;  binner, psfe1, te1_corr2, xo, yo, sig
;  plot2,xo,yo,sig,aa,bb,'','e!D1!N',xrange,yrange,'After Correction smear'
;  if dtype eq 'X' then key=get_kbrd(1)
;  ;; e1 vs psf_e2
;  binner, psfe2, te1_corr2, xo, yo, sig
;  multiplot & plot2,xo,yo,sig,aa,bb,'','',xrange,yrange
;  if dtype eq 'X' then key=get_kbrd(1)
;  ;; e2 vs psf_e1
;  binner, psfe1, te2_corr2, xo, yo, sig
;  multiplot & plot2,xo,yo,sig,aa,bb,'e!D1!N!UPSF!N','e!D2!N',xrange,yrange
;  if dtype eq 'X' then key=get_kbrd(1)
;  ;; e2 vs psf_e2
;  binner, psfe2, te2_corr2, xo, yo, sig
;  multiplot & plot2,xo,yo,sig,aa,bb,'e!D2!N!UPSF!N','',xrange,yrange
;  multiplot,/reset
;
;  if dtype eq 'X' then key=get_kbrd(1)

  erase & multiplot,[0,2,2,0,0],/square
  ;; e1 vs psf_e1
  binner, psfe1, te1_corr, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'','e!D1!N',xrange,yrange,'After Correction'
  ;; e1 vs psf_e2
  multiplot
  binner, psfe2, te1_corr, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'','',xrange,yrange
  ;; e2 vs psf_e1
  multiplot
  binner, psfe1, te2_corr, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'e!D1!N!UPSF!N','e!D2!N',xrange,yrange
  ;; e2 vs psf_e2
  multiplot
  binner, psfe2, te2_corr, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'e!D2!N!UPSF!N','',xrange,yrange
  multiplot,/reset


  if dtype eq 'X' then key=get_kbrd(1)
  xrange=[-0.5,0.5]
  yrange=[-0.5,0.5]

  erase & multiplot,[0,2,2,0,0],/square
  ;; e1 vs psf_e1
  binner, psfe1, te1, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'','e!D1!N',xrange,yrange,'Before Correction'
  ;; e1 vs psf_e2
  multiplot
  binner, psfe2, te1, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'','',xrange,yrange
  ;; e2 vs psf_e1
  multiplot
  binner, psfe1, te2, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'e!D1!N!UPSF!N','e!D2!N',xrange,yrange
  ;; e2 vs psf_e2
  multiplot
  binner, psfe2, te2, xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'e!D2!N!UPSF!N','',xrange,yrange
  multiplot,/reset

  ;; now e_corr original vs psf
  if dtype eq 'X' then key=get_kbrd(1)
  xrange=[-0.4,0.4]
  yrange=[-0.1,0.1]


  corr_o=orig[galo].r[clr]

;  erase & multiplot,[0,2,2,0,0],/square
;  ;; e1 vs psf_e1
;  binner, psfe1, orig[galo].e1[clr]/(1-corr_o), xo, yo, sig
;  plot2,xo,yo,sig,aa,bb,'','e!D1!N',xrange,yrange,'After Correction smear orig'
;  if dtype eq 'X' then key=get_kbrd(1)
;  ;; e1 vs psf_e2
;  binner, psfe2, orig[galo].e1[clr]/(1-corr_o), xo, yo, sig
;  multiplot & plot2,xo,yo,sig,aa,bb,'','',xrange,yrange
;  if dtype eq 'X' then key=get_kbrd(1)
;  ;; e2 vs psf_e1
;  binner, psfe1, orig[galo].e2[clr]/(1-corr_o), xo, yo, sig
;  multiplot & plot2,xo,yo,sig,aa,bb,'e!D1!N!UPSF!N','e!D2!N',xrange,yrange
;  if dtype eq 'X' then key=get_kbrd(1)
;  ;; e2 vs psf_e2
;  binner, psfe2, orig[galo].e2[clr]/(1-corr_o), xo, yo, sig
;  multiplot & plot2,xo,yo,sig,aa,bb,'e!D2!N!UPSF!N','',xrange,yrange
;  multiplot,/reset
;
;  if dtype eq 'X' then key=get_kbrd(1)

  ;; e1 vs psf_e1
  erase & multiplot,[0,2,2,0,0],/square
  binner, psfe1, orig[galo].e1[clr], xo, yo, sig
  plot2,xo,yo,sig,aa,bb,'','e!D1!N',xrange,yrange,'After Correction orig'
  ;; e1 vs psf_e2
  binner, psfe2, orig[galo].e1[clr], xo, yo, sig
  multiplot & plot2,xo,yo,sig,aa,bb,'','',xrange,yrange
  ;; e2 vs psf_e1
  binner, psfe1, orig[galo].e2[clr], xo, yo, sig
  multiplot & plot2,xo,yo,sig,aa,bb,'e!D1!N!UPSF!N','e!D2!N',xrange,yrange
  ;; e2 vs psf_e2
  binner, psfe2, orig[galo].e2[clr], xo, yo, sig
  multiplot & plot2,xo,yo,sig,aa,bb,'e!D2!N!UPSF!N','',xrange,yrange
  multiplot,/reset

  !p.multi=0

  outdir='/sdss2/data0/sdss/tmp/'
  outname = outdir+'matchad-'+colors[clr]+'-'+run2string(run)+$
    '-'+ntostr(camcol)+'.fit'
  print,'output fits file: outname'
  mwrfits, nomatch, outname, /create

END 
