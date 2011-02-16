PRO binfield2, clr, runs, get=get,stars=stars,addf=addf,log=log,$
             hist=hist

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


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: binfield2, clr, runs, get=get,stars=stars'
     print,''
     print,'Use doc_library,"bincol"  for more help.'  
     return
  ENDIF 


  pold=!p.multi
  !p.multi=[0,1,2]

  IF keyword_set(get) THEN get=1 ELSE get=0
  IF keyword_set(hist) THEN hist = 1 ELSE hist = 0

  colors=['u','g','r','i','z']
  ncol=6
  ncsize=12

  IF n_elements(addf) EQ 0 THEN addf=0
  IF keyword_set(stars) THEN stars=1 ELSE stars=0
  IF stars THEN nnm='star' ELSE nnm='gal'
  IF NOT keyword_set(log) THEN BEGIN
    ylog=0
    xlog=0
  ENDIF ELSE BEGIN
    ylog=1
    xlog=1
  ENDELSE 


  outfile='/sdss3/usrdevel/esheldon/tmp/run752_756_smooth_'+nnm+'.fit'
  dirone='/sdss3/usrdevel/philf/'
  rr=['run752/','run756/']
  

  IF get THEN BEGIN 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Read in each column and put all in one big struct
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    typ='rrun'
    FOR i=0,ncol-1 DO BEGIN
    
      istr=ntostr(i+1)
      name1=dirone+rr[0]+'adat'+istr+'c_smooth_'+nnm+'.fit'
      name2=dirone+rr[1]+'adat'+istr+'c_smooth_'+nnm+'.fit'
      IF i EQ 0 THEN print,'Reading from '+name1
      r52=mrdfits(name1,1,hdr,structyp=typ,/silent)
      r56=mrdfits(name2,1,hdr,structyp=typ,/silent)
      
      IF i EQ 0 THEN BEGIN
                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        minra=min(r52.ra)       ; Find overlap of 756, 752 not general now
        maxra=max(r56.ra)       ; 756 starts at smaller ra.  752 goes 
                                ; to a little bigger ra.
        print,'Minra: ',minra   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        print,'Maxra: ',maxra
        w56=where( r56.ra GT minra, n56)
        w52=where( r52.ra LT maxra, n52)

        IF n52 LT n56 THEN BEGIN
          nfield=n52
          w56=w56[0:nfield-1]
        ENDIF ELSE BEGIN 
          nfield=n56            
          w52=w52[0:nfield-1]
        ENDELSE 
        runs=replicate(r52[0], 2*ncol, nfield) 

      ENDIF 
      
      w56=where( r56.ra GT minra, n56)
      w52=where( r52.ra LT maxra, n52)
      
      r52=r52[w52[0:nfield-1]]
      r56=r56[w56[0:nfield-1]]
        
      runs[2*i,*] = r52
      runs[2*i+1,*] = r56
      
    ENDFOR 
    ;mwrfits,runs,outfile,/create  ;No point, won't be fits standard

  ENDIF ELSE nfield = n_elements(runs[0,*])

  ; square degrees in one field
  d1=.0344

  print
  print,'Number of overlapping fields: ',ntostr(nfield)
  print
  print,' binning #bins  sq.deg/bin  sq.deg.  e1err   e2err',$
    '  nobj/bin'
  print,'-------------------------------------------------------------'
    ;'-------'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bin by various sizes.  Will all be square for addf = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR csize=1, ncsize DO BEGIN 
    FOR fsize=csize, csize+addf DO BEGIN
      cstep=csize-1
      fstep=fsize-1

      is=ntostr(csize)
      js=ntostr(fsize)

      ;; sq. deg. per bin
      ba = csize*fsize*d1

      ;; total no. of bins.
      nc=2*ncol/csize
      nf=nfield/fsize
      nb = nc*nf

      ;; total area used
      ta = nb*ba

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Get all the individual bins
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR cc=0, nc-1 DO BEGIN
        FOR ff=0, nf-1 DO BEGIN 
                                
          IF csize EQ 1 AND fsize EQ 1 THEN BEGIN 
            ; Binsize is 1x1
            te1 = runs[cc,ff].rawe1[clr]
            te2 = runs[cc,ff].rawe2[clr]
            te1err = runs[cc,ff].rawe1err[clr]
            te2err = runs[cc,ff].rawe2err[clr]
            te1sig = runs[cc,ff].rawe1sig[clr]
            te2sig = runs[cc,ff].rawe2sig[clr]
            tnobj=runs[cc,ff].nobj[clr]

            te1sigerr = te1sig/sqrt(2*tnobj)
            te2sigerr = te2sig/sqrt(2*tnobj)
            te1errerr = te1err/sqrt(2*tnobj)
            te2errerr = te2err/sqrt(2*tnobj)

         ENDIF ELSE BEGIN 
            ; Binsize ne 1x1
            te1r = runs[cc*csize:cc*csize+cstep, $
                              ff*fsize:ff*fsize+fstep].rawe1[clr]
            te2r = runs[cc*csize:cc*csize+cstep, $
                              ff*fsize:ff*fsize+fstep].rawe2[clr]
            te1errr = runs[cc*csize:cc*csize+cstep, $
                              ff*fsize:ff*fsize+fstep].rawe1err[clr]
            te2errr = runs[cc*csize:cc*csize+cstep, $
                              ff*fsize:ff*fsize+fstep].rawe2err[clr]
            te1sigr = runs[cc*csize:cc*csize+cstep, $
                              ff*fsize:ff*fsize+fstep].rawe1sig[clr]
            te2sigr = runs[cc*csize:cc*csize+cstep, $
                              ff*fsize:ff*fsize+fstep].rawe2sig[clr]
            tnobj=total(runs[cc*csize:cc*csize+cstep, $
                              ff*fsize:ff*fsize+fstep].nobj[clr])

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ; estimate intrinsic spread for weighting
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

            v1 = (moment(te1r))[1]
            v2 = (moment(te2r))[1]

            tvar1 = te1errr^2
            tvar2 = te2errr^2

            v1int = (v1 - median(tvar1) ) > 0.
            v2int = (v2 - median(tvar2) ) > 0.

            wmom, te1r, sqrt( v1int + tvar1 ), te1, te1sig, te1err, $
              te1sigerr, te1errerr
            wmom, te2r, sqrt( v2int + tvar2 ), te2, te2sig, te2err, $
              te2sigerr, te2errerr


          ENDELSE 

          IF cc EQ 0 AND ff EQ 0 THEN BEGIN
            ; b is for bins
            be1 = te1
            be2 = te2
            be1err = te1err
            be2err = te2err
            be1sig = te1sig
            be2sig = te2sig
            be1sigerr = te1sigerr
            be2sigerr = te2sigerr
            be1errerr = te1errerr
            be2errerr = te2errerr

            nobj = tnobj
          ENDIF ELSE BEGIN 
            be1 = [be1, te1]
            be2 = [be2, te2]
            be1err = [be1err, te1err]
            be2err = [be2err, te2err]
            be1sig = [be1sig, te1sig]
            be2sig = [be2sig, te2sig]
            be1sigerr = [be1sigerr, te1sigerr]
            be2sigerr = [be2sigerr, te2sigerr]
            be1errerr = [be1errerr, te1errerr]
            be2errerr = [be2errerr, te2errerr]

            nobj = [nobj, tnobj]
          ENDELSE 
          
        ENDFOR 
      ENDFOR 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Get averages from the individual bins
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      title='binsize = '+is+'x'+js

      nobjavg=median(nobj)
      nobjtot=total(nobj)
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; find bin to bin spread
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ve1 = (moment(be1))[1]
      ve2 = (moment(be2))[1]
      s1=sqrt(ve1)
      s2=sqrt(ve2)

      var1 = be1err^2
      var2 = be2err^2
      
      v1int = (ve1 - median(var1) ) > 0.
      v2int = (ve2 - median(var2) ) > 0.
      IF v1int EQ 0 THEN print,'Dohh1'
      IF v2int EQ 0 THEN print,'Dohh2'
      
      wmom, be1, sqrt(v1int + var1), tbe1, tbe1sig, tt, tbe1sigerr
      wmom, be2, sqrt(v2int + var2), tbe2, tbe2sig, tt, tbe2sigerr
tbe1sigerr = tbe1sig/sqrt(nb)
tbe2sigerr = tbe2sig/sqrt(nb)
;print,tbe1, tbe2
      IF hist THEN BEGIN 
        plothist, be1, bin=.16*s1,xstyle=1, xrange=[-.2,.2],$
          xtitle='e1',title=title
        plothist, be2, bin=.16*s2,xstyle=1,  xrange=[-.2,.2],$
          xtitle='e2',title=title
        key=get_kbrd(20)
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Find mean uncertainty for each binning
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; try doing first order thing first.

      v1 = (moment(be1err))[1]
      v2 = (moment(be2err))[1]

      s1 = sqrt(v1)
      s2 = sqrt(v2) 

      ; Estimate intrinsic spread for weighting
      var1 = be1errerr^2
      var2 = be2errerr^2

      v1int = (v1 - median(var1) ) > 0.
      v2int = (v2 - median(var2) ) > 0.

      IF v1int EQ 0 THEN print,'Dohh2'
      IF v2int EQ 0 THEN print,'Dohh3'

      wmom, be1err, sqrt( v1int + var1 ), tbe1err, tbe1errsig, tbe1errerr
      wmom, be2err, sqrt( v2int + var2 ), tbe2err, tbe2errsig, tbe2errerr

;      tbe1errerr = sqrt( v1/nb )
;      tbe2errerr = sqrt( v2/nb )

      IF hist THEN BEGIN 
        plothist, be1err, bin=.16*s1,xrange=[0,.07],xstyle=1, $
          xtitle='e1err',title=title
        plothist, be2err, bin=.16*s2,xrange=[0,.07],xstyle=1, $
          xtitle='e2err',title=title
        key=get_kbrd(20)
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Approximate values (a for approximate)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tae1err = 0.32/sqrt(nobjavg)
      tae2err = 0.32/sqrt(nobjavg)

      out='   '+is+'x'+js
      out=out+'   '+ntostr(nb,6)+'     '+ntostr(ba,6)+'   '+ntostr(ta,6)
      out=out+'   '+ntostr(tbe1err,7)+'  '+ntostr(tbe2err,7)+'   '
      out=out+ntostr(nobjavg,6) ;+'   '+ntostr(nobjtot,6)
      print,out

      IF csize EQ 1 AND fsize EQ 1 THEN BEGIN 
        bina = ba

        e1 = tbe1               ;average e1 over all bins
        e2 = tbe2
        e1sig = tbe1sig         ; spread around average
        e2sig = tbe2sig
        e1sigerr = tbe1sigerr   ; Error in sigma.
        e2sigerr = tbe2sigerr

        e1err = tbe1err         ; Average uncertainty over all bins.
        e2err = tbe2err         ; Found uncertainty in each bin and averaged
        e1errerr = tbe1errerr
        e2errerr = tbe2errerr

        ae1err = tae1err
        ae2err = tae2err
      ENDIF ELSE BEGIN
        bina = [bina, ba]

        e1 = [e1, tbe1]
        e2 = [e2, tbe2]
        e1sig = [e1sig, tbe1sig]
        e2sig = [e2sig, tbe2sig]
        e1sigerr = [e1sigerr, tbe1sigerr]
        e2sigerr = [e2sigerr, tbe2sigerr]

        e1err = [e1err,tbe1err]
        e2err = [e2err, tbe2err]
        e1errerr = [e1errerr, tbe1errerr]
        e2errerr = [e2errerr, tbe2errerr]

        ae1err = [ae1err, tae1err]
        ae2err = [ae2err, tae2err]
      ENDELSE 

    ENDFOR 
  ENDFOR 
      

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Sort them by bin area
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  s=sort(bina)
  bina = bina(s)

  e1 = e1[s]
  e2 = e2[s]
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
  
  tt='Runs 752 and 756'+'  '+nnm+' in '+colors[clr]
  xt='square degrees'

  m = max([ 1.1*max(e1err), 1.1*max(e2err),1.1*max(e1sig), 1.1*max(e2sig) ])
;  yt = 'e1 Uncertainty'
  yt=''
  ploterr, bina, e1err, e1errerr, title=tt, xtitle=xt, ytitle=yt,psym=1, $
    ylog=ylog, xlog=xlog, yrange=[0,m], ystyle=1
  oploterr, bina, e1sig, e1sigerr, psym=3
  oplot, bina, ae1err
  legend, ['<e1err>/bin', 'StDev(e1)','.32/sqrt(N)'],psym=[1,3,0],$
    position=[3.,.035]

;  yt = 'e2 Uncertainty'
  yt=''
  ploterr, bina, e2err, e2errerr, title=tt, xtitle=xt, ytitle=yt,psym=1, $
    ylog=ylog, xlog=xlog, yrange=[0,m], ystyle=1
  oploterr, bina, e2sig, e2sigerr, psym=3
  oplot, bina, ae2err
  legend, ['<e2err>/bin', 'StDev(e2)','.32/sqrt(N)'],psym=[1,3,0],$
    position=[3.,.035]

  !p.multi=pold
  return 
END 































