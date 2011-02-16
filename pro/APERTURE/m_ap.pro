PRO m_ap, run1, run2, clr, struct, $
          binsize0=binsize0, $
          binmax = binmax, nstep=nstep,$
          domap=domap, $
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
      print,'-Syntax: m_ap, run1, run2, clrindex [, struct, binsize0=binsize0, binmax=binmax, nstep=nstep, domap=domap, combine=combine, get=get, stars=stars, hist=hist]'
      print,''
      print,'NEED TO FIX OUTPUT FROM BINSIZE TO RADIUS (including the logfile.)  FIX PLOTS'
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


  IF n_elements(binsize0) EQ 0 THEN binsize0 =  0.017   
  IF n_elements(binmax) EQ 0 THEN binmax = 1.67
  IF n_elements(nstep) EQ 0 THEN nstep = 10 


  ;; Make lagarithmically spaced bins.
  iii = lindgen(nstep)
  alpha = alog( binmax/binsize0 )/(nstep - 1)
  binsize = binsize0*exp(alpha*iii)


  ;; take .32 as intrinsic spread in weight
  spread = .32

;  binsize0 = .1  ; 1 square degree on 10th bin
;  step = binsize0  ; how much to increase binsize by each time
;  nstep = 26

  IF keyword_set(domap) THEN domap=1 ELSE domap=0
  IF keyword_set(get) THEN get=1 ELSE get=0
  IF keyword_set(hist) THEN hist = 1 ELSE hist = 0
  IF keyword_set(star) THEN star=1 ELSE star=0
  IF NOT keyword_set(dir) THEN dir = '/sdss4/data1/esheldon/CORRECTED/'
  IF NOT keyword_set(outdir) THEN outdir = '/sdss4/data1/esheldon/M_ap/'
  IF star THEN nnm='star' ELSE nnm='srcgal'
  imdir = dir+'SHAPEIM/'+nnm+'/'

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

  i=1
  istr = '_'+ntostr(i)
  psname = outdir+prename+istr+'_N1.ps'
  dataname = outdir+prename+istr+'_N1.dat'
  WHILE exist(psname) OR exist(dataname) DO BEGIN 
      psname = newname(psname)
      dataname = newname(dataname)
  ENDWHILE 
      
  print, 'postscript file: ',psname
  print, 'data file:       ',dataname

  logname = outdir+'log.txt'
  openw, lun1, logname, /get_lun


  printf, lun1, '  #bins  angsize/bin  sq.deg.  Map   Maperr   rad   raderr  nobj/bin'
  printf, lun1, '--------------------------------------------------------------------'
  flush,lun1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; save time by declaring arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nra     = lonarr(nstep)
  ndec    = nra
  nb      = nra
  ta      = dblarr(nstep)
  tmpmaxra= ta
  tmpmaxdec = ta


  bina    = ta
  nobjavg = ta
  nobjavgerr = ta
  nobjtot = ta
  
  Map       = ta           ; -average aperture mass over all bins
  rad       = ta           ; -average radial
  Maperr    = ta           ; -uncertainty in this average
  raderr    = ta           ;  Note this is VERY different from Maperr

  aMaperr   = ta
  araderr   = ta
          


  ;; Number of bins in each direction
  nra = long( radiff/binsize ) > 1
  ndec = long( decdiff/binsize ) > 1
  tmpmaxra = minra + nra*binsize
  tmpmaxdec = mindec + ndec*binsize

  ;; bin area in sq. deg.
  bina = binsize^2
  
      ;; total no. of bins.
  nb = nra*ndec
  
      ;; total area used
  ta = nb*bina

  ;; images
  IF domap THEN BEGIN 
      imMap = dblarr(ndec[0], nra[0])
      imrad = imMap
      imra = imMap
      imdec = imMap
      imnum = lonarr(ndec[0], nra[0])
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; temporary arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;

  binMap = dblarr(ndec[0], nra[0])
  binMap = binMap
  binrad = binMap
  binMaperr = binMap
  binraderr = binMap
          
  nobj = lonarr(ndec[0], nra[0])
  index=lindgen(ndec[0], nra[0]) ;Used to subscript arrays


  print,'Using binsizes:  ',binsize

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bin by various sizes. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  FOR ci = 0, nstep-1 DO BEGIN 

      print
      print,ntostr(ci+1)+'/'+ntostr(nstep),'  Binsize = ',ntostr(binsize[ci])
      print,ntostr(nra[ci])+' dots'

      ; Get averages in each bin in both directions
      
      Rmax = binsize[ci]/2.
      binarr, struct.ra, binsize[ci], minra, tmpmaxra[ci], rev_ind1
      nobj[*,*] = -1
      FOR ira = 0L, nra[ci]-1 DO BEGIN 
          print, format='($,A)', '.' ;+ntostr(ira)
          IF rev_ind1[ira] NE rev_ind1[ira+1] THEN BEGIN 
              w1 = rev_ind1( rev_ind1[ira]:rev_ind1[ira+1]-1 ) 

              binarr, struct[w1].dec, binsize[ci], mindec, tmpmaxdec[ci], $
                      rev_ind2
              FOR idec = 0L, ndec[ci]-1 DO BEGIN 

                  IF rev_ind2[idec] NE rev_ind2[idec+1] THEN BEGIN 
                      w2 = rev_ind2( rev_ind2[idec]:rev_ind2[idec+1]-1 ) 
                      w = w1[w2]
              
                      ra1 = min(struct[w].ra)
                      ra2 = max(struct[w].ra)
                      dec1 = min(struct[w].dec)
                      dec2 = max(struct[w].dec)

                      ra = (ra1 + ra2)/2.
                      dec = (dec1 + dec2)/2.

                      qmom_test, struct[w].e1, struct[w].e2, $
                        struct[w].dec, struct[w].ra, $
                        dec, ra, Rmax, sqrt( spread^2 + struct[w].uncert^2 ), $
                        tempMap, temprad, tempMaperr, tempraderr, nobj=tempnobj

                      IF tempnobj EQ 0 THEN tempnobj = -1
                  ENDIF ELSE BEGIN
                      tempMap = 0.
                      temprad = 0.
                      tempMaperr = 0.
                      tempraderr = 0.
                      tempnobj = -1
                  ENDELSE 

;print,tempMap,tempMaperr

                  binMap[idec,ira] = tempMap
                  binrad[idec,ira] = temprad
                  binMaperr[idec,ira] = tempMaperr
                  binraderr[idec,ira] = tempraderr
                  
                  nobj[idec,ira] = tempnobj

                  IF domap THEN BEGIN 
                      imMap[idec, ira] = tempMap
                      imrad[idec, ira] = temprad
                      imnum[idec, ira] = tnobj
                      imra[idec, ira]=(ra1+ra2)/2.
                      imdec[idec, ira]=(dec1+dec2)/2.
                  ENDIF 
              
              ENDFOR 
          ENDIF
      ENDFOR 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Output the images
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Indices to use
;      jj = index[0:idec-1, 0:ira-1] ;; -1 because of loop exit status

      
      jj = where(nobj NE -1, njj)
      IF njj EQ 0 THEN BEGIN
          print,'All bins are empty'
          return
      ENDIF ELSE BEGIN
          print
          print,'Used ',ntostr(njj),'/',ntostr(nb[ci]),'  bins'
          nb[ci] = njj
      ENDELSE 

      IF domap THEN BEGIN 
          nrastr = ntostr(nra[ci])
          ndecstr = ntostr(ndec[ci])
          fMap=imdir+'Mapim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          frad=imdir+'radim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fnum=imdir+'nim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fra=imdir+'raim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          fdec=imdir+'decim_'+nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          print
          print,'Creating images: ',imdir+'*_'+$
                nnm+'_'+colors[clr]+nrastr+'x'+ndecstr+'.fit'
          writefits, fMap,   imMap[jj]
          writefits, frad,   imrad[jj]
          writefits, fnum, imnum[jj]
          writefits, fra,   imra[jj]
          writefits, fdec, imdec[jj]
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Get averages from the individual bins
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
      vobj = variance(nobj[jj])
      tnobjavgerr = sqrt(vobj/nb[ci])
      tnobjavg    = median(nobj[jj])
      tnobjtot    = total(nobj[jj])

      wmom, binMap[jj], binMaperr[jj], avebinMap, wsig, avebinMaperr
      wmom, binrad[jj], binraderr[jj], avebinrad, wsig, avebinraderr

      out='   '+ntostr(nb[ci])+'     '+ntostr(binsize[ci])+'   '
      out=out+ntostr(ta[ci])+'   '
      out=out+ntostr(avebinMap)+'   '+ntostr(avebinMaperr)+'   '
      out=out+ntostr(avebinrad)+'   '+ntostr(avebinraderr)+'   '
      out=out+ntostr(tnobjavg) 
      printf, lun1, out
      flush, lun1

      nobjavg[ci] = tnobjavg
      nobjavgerr[ci] = tnobjavgerr
      nobjtot[ci] = tnobjtot

      Map[ci] = avebinMap           ; -average Map over all bins
      rad[ci] = avebinrad
      Maperr[ci] = avebinMaperr  ; -uncertainty in this average
      raderr[ci] = avebinraderr  ;

      taMaperr = spread/sqrt(tnobjavg)
      taraderr = spread/sqrt(tnobjavg)

      aMaperr[ci] = taMaperr
      araderr[ci] = taraderr

  ENDFOR 
      

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Sort them by binsize
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  s=sort(binsize)
  binsize = binsize(s)
  bina = bina(s)

  nobjavg = nobjavg(s)
  nobjavgerr = nobjavgerr(s)
  nobjtot = nobjtot(s)
  Map = Map[s]
  rad = rad[s]
  Maperr = Maperr[s]
  raderr = raderr[s]
  aMaperr = aMaperr(s)
  araderr = araderr(s)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Plot and print the results
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  openw, lun2, dataname, /get_lun
  !textunit = lun2

  mes=    '      Theta'
  mes1=mes+'          Map            Maperr           rad             raderr'
  mes2=mes+'        <nobj>      .32/sqrt(nobj)    nobjtot'
  
  ;; Plot verses radius = binsize/2.  change this later
  printf, lun2, mes1
  forprint, binsize, $
            Map, Maperr, rad, raderr, $
            TEXT=5, $  ;; prints to datafile ( I set !textunit = lun2)
            /silent
  printf, lun2, mes2
  forprint, binsize, $
            nobjavg, aMaperr, nobjtot, $
            TEXT=5, $
            /silent
  flush, lun2

  ;; Close these now in case plots mess up.
  close, lun1
  free_lun, lun1
  close, lun2
  free_lun, lun2

  title = 'Runs '+r1str+' and '+r2str+'  '+nnm+' in '+colors[clr]
  xt='Theta'

  ; create postscript file for output
  makeps, psname, /noland

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Log-Linear plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  yt=''
  ploterr, binsize, Map, Maperr, title=title+'  <M^2>', xtitle=xt, $
           ytitle=yt,psym=1, /xlog
  ploterr, binsize, rad, raderr, title='rad', xtitle=xt, $
           ytitle=yt,psym=1, /xlog
;  oplot, binsize, aMaperr
  legend, ['<Map^2>', '<rad^2>','.32/sqrt(N)'],psym=[1,3,0];,position=[3.,.035]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Log-Log plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  yt=''
  ploterr, binsize, Map, Maperr, title=title+'  <M^2>', xtitle=xt, $
            ytitle=yt,psym=1, /xlog, /ylog
  ploterr, binsize, rad, raderr,title='rad', xtitle=xt, ytitle=yt,psym=1, $
           /xlog, /ylog
;  oplot, binsize, aMaperr
;  legend, ['<Map^2>', '<rad^2>','.32/sqrt(N)'],psym=[1,3,0];,position=[.6,.06]
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Plots of <nobjavg> and nobjtot vs binsize
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  yt='<nobj>'
  ploterr, binsize, nobjavg, nobjavgerr, title=title, ytitle=yt, $
    xtitle=xt,psym=1

  yt = 'Total Objects Used'
  plot, binsize, nobjtot, title=title, ytitle=yt, xtitle=xt,psym=1

  ;close plot
  ep

  !p.multi=pold

  ptime, systime(1)-tt
  return 
END 































