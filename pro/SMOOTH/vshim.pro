PRO vshim, binra, bindec, clr, type, objtype, ngroup=ngroup

  IF n_params() NE 5 THEN BEGIN
      print,'-Syntax: vshim, clr, binra, bindec, type, objtype, ngroup=ngroup'
      print,' type=1 for e1  =2 for e2  =3 for e'
      print,' objtype=1 for gal =2 for star'
      return
  ENDIF 
  colors=['u','g','r','i','z']
  type = type > 1
  type = type < 3

  CASE type OF 
      1: add='e1im_'
      2: add='e2im_'
      3: add='eim_'
  ENDCASE 
  CASE objtype OF 
      1: mid='srcgal'
      2: mid='star'
  ENDCASE 

  IF n_elements(ngroup) EQ 0 THEN ngroup = 8
  dir = '/sdss4/data1/esheldon/CORRECTED/SHAPEIM/'+mid+'/'
  
  file=dir+add+mid+'_'+colors[clr]+ntostr(binra)+'x'+ntostr(bindec)+'.fit'

  t=mrdfits(file,0,hdr, /silent)
  t=t+10.0
  sigma_clip, t, mean, sig, nsig=3.5, niter=4

  gsize=long(binra)/ngroup
  left = binra MOD ngroup
  print,gsize

  up='[A'
  down='[B'
  pgup = '[5~'
  pgdown = '[6~'

  print,'use up and down arrow or page up and down Any other to quit'
  i=0
  WHILE 1 DO BEGIN 

      ra2=i+gsize-1
      IF ra2 LE binra-1 AND i NE -1 THEN BEGIN 
          tmp=t[*, i:ra2]
          rdis_setup, tmp, pls
          rdis, tmp, pls,/silent,/noframe,low=mean-2.5*sig,high=mean+3.5*sig
          axis,yaxis=0,yticks=1,ytickn=[ntostr(i),ntostr(ra2)]
      ENDIF ELSE i=iold
      go=1
      WHILE go DO BEGIN 
          getseq,key
          CASE key OF 
              up: BEGIN
                  iold=i && i=i+1 && go=0
              END 
              down: BEGIN 
                  iold=i && i=i-1 && go=0
              END
              pgup: BEGIN
                  iold=i && i=i+gsize-1 < (binra-gsize) && go=0
              END 
              pgdown: BEGIN
                  iold=i && i=i-(gsize-1) > 0 && go=0
              END 
              '': return
              ELSE: go=1
          ENDCASE 
      ENDWHILE 
  ENDWHILE 
  
  return
END 
