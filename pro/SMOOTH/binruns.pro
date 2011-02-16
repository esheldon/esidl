PRO binruns, run1, run2, clr, struct, $
             binsize0=binsize0, step=step, nstep=nstep,$
             map=map, $
             get=get, $
             star=star,$
             hist=hist, $
             dir=dir, $
             outdir=outdir 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
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
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 3 THEN BEGIN 
      print,'-Syntax: binruns, run1, run2, clrindex [, struct, binsize0=binsize0, step=step, nstep=nstep, map=map, combine=combine, get=get, stars=stars, hist=hist]'
      print,''
      print,'Use doc_library,"binradec"  for more help.'  
      return
  ENDIF 
  
  pold=!p.multi
  !p.multi=[0,1,2]

  tt = systime(1)

  COMMON par, minra,maxra,mindec,maxdec,radiff,decdiff

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Defaults
  ;; side of the initial binsize in degrees. 
  ;; .2276 = side of frame in column directon.  But: overlap reduces
  ;; what we can fit.  For 12 across use .21
  IF n_elements(binsize0) EQ 0 THEN binsize0 =  0.21     
  IF n_elements(nstep) EQ 0 THEN nstep = 11
  IF n_elements(step) EQ 0 THEN step = binsize0

  ;; take .32 as intrinsic spread in weight
  spread = .32

;  binsize0 = .1  ; 1 square degree on 10th bin
;  step = binsize0  ; how much to increase binsize by each time
;  nstep = 26

  IF keyword_set(map) THEN map=1 ELSE map=0
  IF keyword_set(get) THEN get=1 ELSE get=0
  IF keyword_set(hist) THEN hist = 1 ELSE hist = 0
  IF keyword_set(star) THEN star=1 ELSE star=0
  IF NOT keyword_set(dir) THEN dir = '/sdss4/data1/esheldon/CORRECTED/'
  IF NOT keyword_set(outdir) THEN dir = '/sdss4/data1/esheldon/SHAPEVAR/'
  IF star THEN nnm='star' ELSE nnm='srcgal'
  imdir = outdir+'SHAPEIM/'+nnm+'/'

  colors=['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get the 2 interleaved runs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  r1str = ntostr(run1)
  r2str = ntostr(run2)

  IF get THEN BEGIN
      sname=dir+'run'+r1str+'_'+r2str+'_'+nnm+'_'+colors[clr]+'_overlap.fit'
      IF NOT exist(sname) THEN $
        sname=dir+'run'+r2str+'_'+r1str+'_'+nnm+'_'+colors[clr]+'_overlap.fit'
      IF NOT exist(sname) THEN BEGIN
          print,'No overlap file exists for these two runs'
          return
      ENDIF 
      struct = mrdfits(sname, 1, hdr)
  ENDIF 
  minra = min(struct.ra)
  maxra = max(struct.ra)
  mindec = min(struct.dec)
  maxdec = max(struct.dec)
  radiff  = maxra - minra
  decdiff = maxdec - mindec

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Open the files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = 'run'+r1str+'_'+r2str+'_'+nnm+'_'+colors[clr]+'_bin'

  psname = outdir+prename+'_N1.ps'
  dataname = outdir+prename+'_N1.dat'
  WHILE exist(psname) OR exist(dataname) DO BEGIN 
      psname = newname(psname)
      dataname = newname(dataname)
      i = i+1
  ENDWHILE 
      
  print, 'postscript file: ',psname
  print, 'data file:       ',dataname

  logname = outdir+'log.txt'
  openw, lun1, logname, /get_lun


  printf, lun1, '  #bins  sq.deg/bin  sq.deg.  e1err   e2err  nobj/bin'
  printf, lun1, '-------------------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; save time by declaring arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  binsize = dblarr(nstep)
  tmpmaxra= binsize
  tmpmaxdec = binsize
  ta      = binsize
  nra     = lonarr(nstep)
  ndec    = nra
  nb      = nra

  bina    = binsize
  nobjavg = binsize
  nobjavgerr = binsize
  nobjtot = binsize
  
  e1       = binsize            ; -average e1 over all bins
  e2       = binsize
  e1uncert = binsize            ; -uncertainty in this average
  e2uncert = binsize            ;  Note this is VERY different from e1err
  e1sig    = binsize            ;  is average uncertainty for a given bin
  e2sig    = binsize            ; -spread around average
  e1sigerr = binsize            ; -Error in sigma.
  e2sigerr = binsize

  e1err    = binsize            ; -Average uncert over this sized bins
  e2err    = binsize            ;  Found uncert in each bin and averaged
  e1errerr = binsize
  e2errerr = binsize

  ae1err   = binsize
  ae2err   = binsize
          

  FOR ci=0, nstep-1 DO BEGIN
      binsize[ci] = binsize0 + ci*step

      ;; Number of bins in each direction
      nra[ci] = long( radiff/binsize[ci] ) > 1
      ndec[ci] = long( decdiff/binsize[ci] ) > 1
      tmpmaxra[ci] = minra + nra[ci]*binsize[ci]
      tmpmaxdec[ci] = mindec + ndec[ci]*binsize[ci]

      ;; bin area in sq. deg.
      bina[ci] = binsize[ci]^2

      ;; total no. of bins.
      nb[ci] = nra[ci]*ndec[ci]

      ;; total area used
      ta[ci] = nb[ci]*bina[ci]
  ENDFOR 

  ;; images
  IF map THEN BEGIN 
      ime1 = dblarr(ndec[0], nra[0])
      ime2 = ime1
      ime  = ime1
      imra = ime1
      imdec = ime1
      imnum = lonarr(ndec[0], nra[0])
  ENDIF 

  ;; temporary arrays
  be1 = dblarr(ndec[0], nra[0])
  be1 = be1
  be2 = be1
  be1sig = be1
  be2sig = be1
  be1err = be1
  be2err = be1
  be1sigerr = be1
  be2sigerr = be1
  be1errerr = be1
  be2errerr = be1
          
  nobj = lonarr(ndec[0], nra[0])
  index=lindgen(ndec[0], nra[0]) ;Used to subscript arrays

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bin by various sizes. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  FOR ci = 0, nstep-1 DO BEGIN 

      print
      print,ntostr(ci+1)+'/'+ntostr(nstep)
      print,ntostr(nra[ci])+' dots'

      ; Get averages in each bin in both directions

      binarr, struct.ra, binsize[ci], minra, tmpmaxra[ci], rev_ind1
      FOR ira = 0L, nra[ci]-1 DO BEGIN 
          print, format='($,A)', '.'
          w1 = rev_ind1( rev_ind1[ira]:rev_ind1[ira+1]-1 ) 

          binarr, struct[w1].dec, binsize[ci], mindec, tmpmaxdec[ci], rev_ind2
          FOR idec = 0L, ndec[ci]-1 DO BEGIN 

              w2 = rev_ind2( rev_ind2[idec]:rev_ind2[idec+1]-1 ) 
              w = w1[w2]
              tnobj=n_elements(w)
              ra1 = min(struct[w].ra)
              ra2 = max(struct[w].ra)
              dec1 = min(struct[w].dec)
              dec2 = max(struct[w].dec)
        
              wmom, struct[w].e1, sqrt( spread^2 + struct[w].uncert^2 ), $
                te1, te1sig, te1err, te1sigerr, te1errerr
              wmom, struct[w].e2, sqrt( spread^2 + struct[w].uncert^2 ), $
                te2, te2sig, te2err, te2sigerr, te2errerr

              be1[idec,ira] = te1
              be2[idec,ira] = te2
              be1sig[idec,ira] = te1sig
              be2sig[idec,ira] = te2sig
              be1err[idec,ira] = te1err
              be2err[idec,ira] = te2err
              be1sigerr[idec,ira] = te1sigerr
              be2sigerr[idec,ira] = te2sigerr
              be1errerr[idec,ira] = te1errerr
              be2errerr[idec,ira] = te2errerr
          
              nobj[idec,ira] = tnobj

              IF map THEN BEGIN 
                  ime1[idec, ira] = te1
                  ime2[idec, ira] = te2
                  ime[idec, ira] = sqrt( te1^2 + te2^2 )
                  imnum[idec, ira] = tnobj
                  imra[idec, ira]=(ra1+ra2)/2.
                  imdec[idec, ira]=(dec1+dec2)/2.
              ENDIF 
              
          ENDFOR 
      ENDFOR 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Output the images
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Indices to use
      jj = index[0:idec-1, 0:ira-1] ;; -1 because of exit status

      IF map THEN BEGIN 
          nrastr = ntostr(nra[ci])
          ndecstr = ntostr(ndec[ci])
          fe1=imdir+'e1im_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fe2=imdir+'e2im_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fe = imdir+'eim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fnum=imdir+'nim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fra=imdir+'raim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fdec=imdir+'decim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          print
          print,'Creating images: ',imdir+'*_'+$
                nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          writefits, fe1,   ime1[jj]
          writefits, fe2,   ime2[jj]
          writefits, fe,     ime[jj]
          writefits, fnum, imnum[jj]
          writefits, fra,   imra[jj]
          writefits, fdec, imdec[jj]
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Get averages from the individual bins
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      vobj = variance(nobj[jj])
      tnobjavgerr = sqrt(vobj/nb[ci])
      tnobjavg=median(nobj[jj])
      tnobjtot=total(nobj[jj])
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; find bin to bin spread
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ve1 = variance(be1[jj])
      ve2 = variance(be2[jj])
      s1=sqrt(ve1)
      s2=sqrt(ve2)

      var1 = be1err[jj]^2
      var2 = be2err[jj]^2
      
      v1int = (ve1 - median(var1) ) > 0.
      v2int = (ve2 - median(var2) ) > 0.
      IF v1int EQ 0 THEN print,'Dohh1'
      IF v2int EQ 0 THEN print,'Dohh2'
      
      wmom, be1[jj], sqrt(v1int + var1), tbe1, tbe1sig, tbe1uncert, tbe1sigerr
      wmom, be2[jj], sqrt(v2int + var2), tbe2, tbe2sig, tbe2uncert, tbe2sigerr

      ;; Using standard Unweighted formula
      tbe1sigerr = tbe1sig/sqrt(2*nb[ci])
      tbe2sigerr = tbe2sig/sqrt(2*nb[ci])

      IF hist THEN BEGIN 
          plothist, be1[jj], bin=.16*s1,xstyle=1, xrange=[-.2,.2],$
            xtitle='e1',title=title
          plothist, be2[jj], bin=.16*s2,xstyle=1,  xrange=[-.2,.2],$
            xtitle='e2',title=title
          key=get_kbrd(20)
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Find mean uncertainty for each binning
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      v1 = variance(be1err[jj])
      v2 = variance(be2err[jj])

      s1 = sqrt(v1)
      s2 = sqrt(v2) 

      ;; Weighted 
      ; Estimate intrinsic spread for weighting
      var1 = be1errerr[jj]^2
      var2 = be2errerr[jj]^2

      v1int = (v1 - median(var1) ) > 0.
      v2int = (v2 - median(var2) ) > 0.

      IF v1int EQ 0 THEN print,'Dohh2'
      IF v2int EQ 0 THEN print,'Dohh3'

      wmom, be1err[jj], sqrt( v1int + var1 ), tbe1err, tbe1errsig, tbe1errerr
      wmom, be2err[jj], sqrt( v2int + var2 ), tbe2err, tbe2errsig, tbe2errerr

      ;; Unweighted formulae

;      tbe1err = median(be1err)
;      tbe2err = median(be2err)
;      tbe1errsig = s1
;      tbe2errsig = s2
;      tbe1errerr = tbe1errsig/sqrt(nb[ci])
;      tbe2errerr = tbe2errsig/sqrt(nb[ci])

      IF hist THEN BEGIN 
          plothist, be1err[jj], bin=.16*s1,xrange=[0,.07],xstyle=1, $
            xtitle='e1err',title=title
          plothist, be2err[jj], bin=.16*s2,xrange=[0,.07],xstyle=1, $
            xtitle='e2err',title=title
          key=get_kbrd(20)
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Approximate values (a for approximate)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tae1err = 0.32/sqrt(tnobjavg)
      tae2err = 0.32/sqrt(tnobjavg)

      out='   '+ntostr(nb[ci],6)+'     '+ntostr(bina[ci],6)+'   '
      out=out+ntostr(ta[ci],6)+'   '
      out=out+ntostr(tbe1err,7)+'  '+ntostr(tbe2err,7)+'   '
      out=out+ntostr(tnobjavg,6) 
      printf, lun1, out

      nobjavg[ci] = tnobjavg
      nobjavgerr[ci] = tnobjavgerr
      nobjtot[ci] = tnobjtot

      e1[ci] = tbe1             ; -average e1 over all bins
      e2[ci] = tbe2
      e1uncert[ci] = tbe1uncert ; -uncertainty in this average
      e2uncert[ci] = tbe2uncert ;  Note this is VERY different from e1err which
      e1sig[ci] = tbe1sig       ;  is average uncertainty for a given bin
      e2sig[ci] = tbe2sig       ; -spread around average
      e1sigerr[ci] = tbe1sigerr ; -Error in sigma.
      e2sigerr[ci] = tbe2sigerr

      e1err[ci] = tbe1err       ; -Average uncertainty over this sized bins
      e2err[ci] = tbe2err       ;  Found uncertainty in each bin and averaged
      e1errerr[ci] = tbe1errerr
      e2errerr[ci] = tbe2errerr

      ae1err[ci] = tae1err
      ae2err[ci] = tae2err

  ENDFOR 
      

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Sort them by bin area
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  s=sort(bina)
  bina = bina(s)

  nobjavg = nobjavg(s)
  nobjavgerr = nobjavgerr(s)
  nobjtot = nobjtot(s)
  e1 = e1[s]
  e2 = e2[s]
  e1uncert = e1uncert[s]
  e2uncert = e2uncert[s]
  e1sig = e1sig[s]
  e2sig = e2sig[s]
  e1sigerr = e1sigerr[s]
  e2sigerr = e2sigerr[s]

  e1err = e1err(s)
  e2err = e2err(s)
  e1errerr = e1errerr(s)
  e2errerr = e2errerr(s)

  ae1err = ae1err(s)
  ae2err = ae2err(s)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Plot and print the results
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mes=    '      binarea'
  mes1=mes+'         e1err        e1errerr         e2err          e2errerr'
  mes2=mes+'         e1sig        e1sigerr         e2sig          e2sigerr'
  mes3=mes+'         <e1>         e1uncert         <e2>           e2uncert'
  mes4=mes+'        <nobj>      .32/sqrt(nobj)    nobjtot'
  
  openw, lun2, dataname, /get_lun
  !textunit = lun2

  printf, lun2, mes1
  forprint, bina, $
            e1err, e1errerr, e2err, e2errerr, $
            TEXT=5, $  ;; prints to datafile ( I set !textunit = lun2)
            /silent
  printf, lun2, mes2
  forprint, bina, $
            e1sig, e1sigerr, e2sig, e2sigerr, $
            TEXT=5, $
            /silent
  printf, lun2, mes3
  forprint, bina, $
            e1, e1uncert, e2, e2uncert, $
            TEXT=5, $
            /silent
  printf, lun2, mes4
  forprint, bina, $
            nobjavg, ae1err, nobjtot, $
            TEXT=5, $
            /silent

  title = 'Runs '+r1str+' and '+r2str+'  '+nnm+' in '+colors[clr]
  xt='square degrees'

  m = max([ 1.1*max(e1err), 1.1*max(e2err),1.1*max(e1sig), 1.1*max(e2sig) ])

  ; create postscript file for output
  makeps, psname, /noland

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Log plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  yt=''
  ploterr, bina, e1err, e1errerr, title=title, xtitle=xt, ytitle=yt,psym=1, $
           /xlog, /ylog
  oploterr, bina, e1sig, e1sigerr, psym=3
  oplot, bina, ae1err
  legend, ['<e1err>', 'bin2bin spread','.32/sqrt(N)'],psym=[1,3,0],$
          position=[.6,.06]

  yt=''
  ploterr, bina, e2err, e2errerr, title=title, xtitle=xt, ytitle=yt,psym=1, $
           /xlog, /ylog
  oploterr, bina, e2sig, e2sigerr, psym=3
  oplot, bina, ae2err
  legend, ['<e2err>', 'bin2bin spread','.32/sqrt(N)'],psym=[1,3,0],$
          position=[.6,.06]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Linear plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  yt=''
  ploterr, bina, e1err, e1errerr, title=title, xtitle=xt, ytitle=yt,psym=1, $
           yrange=[0,m]
  oploterr, bina, e1sig, e1sigerr, psym=3
  oplot, bina, ae1err
  legend, ['<e1err>', 'bin2bin spread','.32/sqrt(N)'],psym=[1,3,0],$
          position=[3.,.035]

  yt=''
  ploterr, bina, e2err, e2errerr, title=title, xtitle=xt, ytitle=yt,psym=1, $
           yrange=[0,m]
  oploterr, bina, e2sig, e2sigerr, psym=3
  oplot, bina, ae2err
  legend, ['<e2err>', 'bin2bin spread','.32/sqrt(N)'],psym=[1,3,0],$
          position=[3.,.035]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Plot of <e1> and <e2> vs binsize
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  yt='<e1>'
  ploterr, bina, e1, e1uncert, title=title, xtitle=xt, ytitle=yt,psym=1

  yt='<e2>'
  ploterr, bina, e2, e2uncert, title=title, xtitle=xt, ytitle=yt,psym=1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Plots of <nobjavg> and nobjtot vs binsize
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  yt='<nobj>'
  ploterr, bina, nobjavg, nobjavgerr, title=title, ytitle=yt, xtitle=xt,psym=1

  yt = 'Total Objects Used'
  plot, bina, nobjtot, title=title, ytitle=yt, xtitle=xt,psym=1

  ;close plot
  ep

  !p.multi=pold


  close, lun1
  free_lun, lun1
  close, lun2
  free_lun, lun2

  ptime, systime(1)-tt
  return 
END 































