PRO rsmear_cuts_get_rsmear, struct, clr, rsmear, hirata=hirata

  nclr = n_elements(clr)

  tt = tag_names(struct[0])
  IF keyword_set(hirata) THEN BEGIN 
      rind = where(tt EQ 'M_R_H')
  ENDIF ELSE BEGIN 
      rind = where(tt EQ 'M_R')
  ENDELSE 

  IF nclr EQ 1 THEN BEGIN
      rsmear = struct.(rind)[clr]
  ENDIF ELSE IF nclr EQ 3 THEN BEGIN 

      boloff= 0.0
      
      off=.15
      power=.2
      
      anone = 0

      ;; check for binned1 (detected) in 
      ;; each bandpass, as well as bad
      ;; measurements (-9999.), etc.

      make_flag_struct,fs
      fs.binned1='Y'

      det1=fltarr(n_elements(struct))
      det2=det1
      det3=det1

      ;; min mag val
      minmagval = 0.0
      ;; min moment val
      minmomval = 0.0

      ;; smear polarizabilities
      grsmear = struct.(rind)[1]
      rrsmear = struct.(rind)[2]
      irsmear = struct.(rind)[3]

      clr = 1
      flag_select,struct,fs,clr,index
      IF index[0] NE -1 THEN BEGIN 
          w=where(struct[index].psfcounts[clr] GT minmagval AND $
                  struct[index].counts_exp[clr] GT minmagval AND $
                  struct[index].m_rr_cc_psf[clr] GT minmomval AND $
                  grsmear[index] GT 0.0, nw)
      ENDIF ELSE nw=0
      IF nw NE 0 THEN BEGIN
          index = index[w]
          det1[index]=1.0
      ENDIF 

      clr = 2
      flag_select,struct,fs,clr,index
      IF index[0] NE -1 THEN BEGIN 
          w=where(struct[index].psfcounts[clr] GT minmagval AND $
                  struct[index].counts_exp[clr] GT minmagval AND $
                  struct[index].m_rr_cc_psf[clr] GT minmomval AND $
                  rrsmear[index] GT 0.0, nw)
      ENDIF ELSE nw=0
      IF nw NE 0 THEN BEGIN
          index = index[w]
          det2[index]=1.0
      ENDIF 

      clr = 3
      flag_select,struct,fs,clr,index
      IF index[0] NE -1 THEN BEGIN 
          w=where(struct[index].psfcounts[clr] GT minmagval AND $
                  struct[index].counts_exp[clr] GT minmagval AND $
                  struct[index].m_rr_cc_psf[clr] GT minmomval AND $
                  irsmear[index] GT 0.0, nw)
      ENDIF ELSE nw=0
      IF nw NE 0 THEN BEGIN
          index = index[w]
          det3[index]=1.0
      ENDIF 

      detnum=det1+det2+det3
      wnone=where(detnum EQ 0,anone)

      gpsf=struct.psfcounts(1)
      rpsf=struct.psfcounts(2)
      ipsf=struct.psfcounts(3)
      
      gexp=struct.counts_exp(1)
      rexp=struct.counts_exp(2)
      iexp=struct.counts_exp(3)

      gmom=struct.m_rr_cc_psf(1)
      rmom=struct.m_rr_cc_psf(2)
      imom=struct.m_rr_cc_psf(3)

      ;; just sanity checks
      ggpsf=gpsf > 1.0 < 28.0
      iipsf=ipsf > 1.0 < 28.0
      rrpsf=rpsf > 1.0 < 28.0
      
      ggexp=gexp > 1.0 < 28.0
      iiexp=iexp > 1.0 < 28.0
      rrexp=rexp > 1.0 < 28.0
      
      lup2flux,ggpsf,fgpsf,b=1
      lup2flux,rrpsf,frpsf,b=2
      lup2flux,iipsf,fipsf,b=3
      
      lup2flux,ggexp,fgexp,b=1
      lup2flux,rrexp,frexp,b=2
      lup2flux,iiexp,fiexp,b=3
  
      f_psf=(  det1*(fgpsf>0.0) + det2*(frpsf>0.0) + det3*(fipsf>0.0) )/3.0
      f_exp=(  det1*(fgexp>0.0) + det2*(frexp>0.0) + det3*(fiexp>0.0) )/3.0
      
      flux2lup,f_psf,mag_psf,b=2
      flux2lup,f_exp,mag_exp,b=2
      
      mag_psf=mag_psf+boloff
      mag_exp=mag_exp+boloff
      
      con=mag_psf-mag_exp
      
      c=(((con+off) > 0)^power)-off^power
      c=c<1.0

      mom=(gmom*fgpsf+rmom*frpsf+imom*fipsf)/(fgpsf+frpsf+fipsf)
      mom2seeing,mom,see
      mag=mag_exp

      norm = ( (fgpsf>0.0)+(frpsf>0.0)+(fipsf>0.0))
      rsmear = (grsmear*(fgpsf>0.0)+rrsmear*(frpsf>0.0)+irsmear*(fipsf>0.0))/norm

;      norm = (det1*fgpsf+det2*frpsf+det3*fipsf)
;      rsmear = (det1*grsmear*(fgpsf>0.0)+det2*rrsmear*(frpsf>0.0)+det3*irsmear*(fipsf>0.0))/norm

;      norm = (det1 + det2 + det3)
;      rsmear = (det1*grsmear + det2*rrsmear + det3*irsmear)/(det1+det2+det3)

      IF anone GT 0 THEN BEGIN
          c[wnone]=-10.0
          mag[wnone]=25.0
          rsmear[wnone] = -9999.
      ENDIF


  ENDIF ELSE message,'nclr = 1 or 3 please'

END 

PRO rsmear_cuts_percentage_cut, probgal, rsmear, purity, clr, rcut, ngal

  ;; new way with probabilities
  nprob = n_elements(probgal)
  s = sort(rsmear)
      
  cumprob = total(probgal[s], /cumulative)
  num = findgen(nprob) + 1.

  frac = cumprob/num

  w=where(frac GE purity, nw)
  IF nw EQ 0 THEN BEGIN
      rcut = 0.0
      ngal = 0L
  ENDIF ELSE BEGIN 
      w = (max(w))[0]
      rcut = rsmear[s[w]]
      ngal = num[w]
  ENDELSE 

  ;; plot only a few points
  nplot = 5000

  ;; randomly sample nplot points
;  xx = arrscl(findgen(nplot), 0.0, 2.5)
;  yy = interpol(frac, rsmear[s], xx)

  randind = long( arrscl(randomu(seed, nplot), 0.0, nprob, $
                         arrmin=0.0, arrmax=1.0) )
  xx = rsmear[s[randind]]
  yy = frac[randind]

  sx = sort(xx)
  xx = xx[sx]
  yy = yy[sx]

  plot, xx, yy, $
        ytitle='Cumulative Galaxy Fraction', $
        xtitle='R!Dsmear!N['+!colors[clr]+']', charsize=1.25

  oplot, [rcut, rcut], [0, 10000.], color=!red
  mess = ['Purity Cut = '+ntostr(purity,4,/round), $
          'R!Dsmear!N cut = '+ntostr(rcut, 5, /round)]
  legend, mess, /right,charsize=0.7

END 

PRO rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts, $
                 camcol=camcol, overwrite=overwrite, nops=nops, $
                 hirata=hirata, outdir=outdir

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts,'
      print,' camcol=camcol, overwrite=overwrite, nops=nops, hirata=hirata'
      return
  ENDIF 

  IF n_elements(struct) EQ 0 THEN BEGIN 
      
      tags = ['run','rerun','camcol','field','id', $
              'm_e1', 'm_e2', 'm_e1e1err', 'm_e1e2err', 'm_e2e2err', $
              'm_e1_corr', 'm_e2_corr', $
              'm_r', 'm_r_h', 'petrocounts','reddening', 'objc_prob_psf', $
              'PSFCOUNTS', 'COUNTS_EXP','M_RR_CC_PSF', 'FLAGS', 'FLAGS2']

      IF n_elements(camcol) EQ 0 THEN BEGIN 
          colmin = 1
          colmax = 6
      ENDIF ELSE BEGIN 
          colmin = camcol
          colmax = camcol
      ENDELSE 
      FOR col=colmin,colmax DO BEGIN 
          read_tsobj, [run,rerun,col], tstruct, /all, /corr, $
                      taglist=tags, ex_struct={m_r_avg:0.0}

          IF n_elements(struct) EQ 0 THEN BEGIN
              struct=temporary(tstruct) 
          ENDIF ELSE BEGIN
              concat_structs, temporary(struct),temporary(tstruct),tmp
              struct=temporary(tmp)
          ENDELSE 
          
      ENDFOR 

;      add_tag, struct, 'm_r_avg', 0.0, new_struct
;      delvarx, struct
;      struct = temporary(new_struct)

      combine_ellip_cove1e2, struct.m_e1[1:3], struct.m_e2[1:3], $
                             struct.m_e1e1err[1:3], struct.m_e1e2err[1:3], $
                             struct.m_e2e2err[1:3], struct.m_r[1:3], $
                             new_e1, new_e2, $
                             new_e1e1err, new_e1e2err, new_e2e2err,$
                             new_smear, $
                             combine_flag, good
      
      struct.m_r_avg = new_smear

  ENDIF 

  IF n_elements(camcol) EQ 0 THEN BEGIN
      clstr = '' 
  ENDIF ELSE BEGIN
      clstr = ' Camcol: '+ntostr(camcol)
  ENDELSE 

  ;; output file names
  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, $
                     camcol=camcol, hirata=hirata, outdir=outdir

  IF fexist(fitfile) AND (NOT keyword_set(overwrite)) THEN BEGIN
      print,'Fits file: ',fitfile,' already exists. Returning'
      return
  ENDIF 

  tt = tag_names(struct[0])
  IF keyword_set(hirata) THEN BEGIN 
      rind = where(tt EQ 'M_R_H')
  ENDIF ELSE BEGIN 
      rind = where(tt EQ 'M_R')
  ENDELSE 

  IF NOT keyword_set(nops) THEN begplot, name=psfile, /color

  rmag = struct.petrocounts[2] - struct.reddening[2]
;  rmag = struct.counts_exp[2] - struct.reddening[2]

  step = 0.25
  minmag = 20.0
  nmag = 8
  
  maglow = fltarr(nmag)
  maghigh = fltarr(nmag)

  ;; bigger steps at low mag
  ;; smaller steps at high mag
  maglow[0] = minmag
  maghigh[0] = minmag + step
  FOR i=1L, nmag-1 DO BEGIN 
      maglow[i] = maglow[i-1] + step
      maghigh[i] = maghigh[i-1] + step
  ENDFOR 

  maglow = [18.0, 19.0, 19.5, maglow]
  maghigh = [19.0, 19.5, 20.0, maghigh]
  nmag = n_elements(maglow)

;forprint,maglow,maghigh
;stop
  meanmag = fltarr(nmag)
  ngal = fltarr(nmag)
  ngalcumul = fltarr(nmag)
  rcuts = fltarr(nmag)

  minrsmear = 0.0

;  rsmear_cuts_get_rsmear, struct, clr, rsmear, hirata=hirata
  rsmear = struct.m_r_avg

  FOR i=0L, nmag-1 DO BEGIN 
      ;;struct.objc_prob_psf LT 0.5 AND $
      IF n_elements(camcol) NE 0 THEN BEGIN  
          w=where(rmag LT maghigh[i] AND $
                  rmag GT maglow[i] AND $
                  rsmear GT minrsmear AND $
                  rsmear LT 2.0 AND $
                  struct.camcol EQ camcol AND $
                  struct.objc_prob_psf GE 0.0, ngood)
      ENDIF ELSE BEGIN 
          w=where(rmag LT maghigh[i] AND $
                  rmag GT maglow[i] AND $
                  rsmear GT minrsmear AND $
                  rsmear LT 2.0 AND $
                  struct.objc_prob_psf GE 0.0, ngood)
      ENDELSE 

      meanmag[i] = mean( rmag[w] )

      probgal = 1.-struct.objc_prob_psf

      min=0
      max=2.0
      bin=0.01
      
      !p.multi = [0,0,2]
      
      objc_prob_psf_cut = 0.5
      gal  = where( struct[w].objc_prob_psf LT objc_prob_psf_cut, ng)
      star = where( struct[w].objc_prob_psf GT objc_prob_psf_cut, ns)
      mess = ['All', 'Bayes Gal', 'Bayes Star']
      mess2 = ['Bayesian','Run: '+ntostr(run)+' Rerun: '+ntostr(rerun)+clstr, $
               ntostr(maglow[i],5,/round)+' < r < '+$
               ntostr(maghigh[i],5,/round)]
 
     
      ;; plot everything
      plothist,rsmear[w], xhist, yhist, min=min,max=max,bin=bin, $
               xtitle='R!Dsmear!N['+!colors[clr]+']',$
               title=title, charsize=1.25

      ;; "galaxies"
      plothist,rsmear[w[gal]],gal_xhist, gal_yhist, $
               min=min,max=max,bin=bin, $
               color=!red, /overplot

      ;; "stars"
      plothist,rsmear[w[star]],min=min,max=max,bin=bin,$
               color=!blue,/overplot
      legend, mess,$
              colors=[!p.color,!red,!blue],line=[0,0,0], $
              thick=replicate(!p.thick,3),charsize=0.7
      legend, mess2, /right, charsize=0.7

      rsmear_cuts_percentage_cut, probgal[w], rsmear[w], purity, clr, rcut, tngal


      ngal[i] = tngal
      IF i GT 0 THEN ngalcumul[i] = ngalcumul[i-1] + tngal $
      ELSE ngalcumul[i] = tngal

      rcuts[i] = rcut
      IF (!d.name EQ 'X') THEN key=get_kbrd(1)

      !p.multi=0

  ENDFOR 

  xtit = 'r petrosian'
  ytit = 'R!Dsmear!N['+!colors[clr]+'] Cut'
  myusersym, 'fill_circle'
  aplot, 1, meanmag, rcuts, psym=8, xtit=xtit,ytit=ytit
  oplot, meanmag, rcuts

  legend,['Run: '+ntostr(run)+' Rerun: '+ntostr(rerun)+clstr,$
          'Purity = '+ntostr(purity,4,/round), $
          'Ngal = '+ntostr(long(total(ngal)))],/right,charsize=0.7

  IF NOT keyword_set(nops) THEN endplot

  t=create_struct('meanmag', meanmag, $
                  'rsmear_cuts', rcuts,$
                  'ngal', ngal, $
                  'ngalcumul', ngalcumul, $
                  'ngaltot', long(total(ngal)))

  hdr = ['END']
  sxaddpar, hdr, 'PURITY', purity
  mwrfits2, t, fitfile, /create, hdr0=hdr

END 
