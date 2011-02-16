PRO makecand, rass, bcg, clusters=clusters, linear=linear, fchart=fchart

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    
;       
; PURPOSE:
;    
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
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
     print,'-Syntax: makecand, rass, bcg, clusters=clusters, linear=linear'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 
  
  zmax = 1.0
  zmin = 0.0

  w = where(rass.z LE zmax AND rass.z GE zmin AND $
            rass.use EQ 1 AND rass.l_x GT 0., ncand)

  IF ncand EQ 0 THEN BEGIN
      print,'No clusters found'
      return
  ENDIF 

  tags = tag_names(rass)
  wt=where(tags EQ 'L_O',nt)
  wgal=where(tags EQ 'NGALS',ngal)
  IF ngal EQ 0 THEN wgal=where(tags EQ 'NGAL',ngal)
  IF ngal EQ 0 THEN return

  cand = w[ reverse( sort(rass[w].l_x) ) ]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make the html document
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  

  home = '~/'
  HTML_DIR = home+'WWW/weaklens/clustercand/'
  HEADER_FILE=HTML_DIR+'header.html'
  FOOTER_FILE=HTML_DIR+'footer.html'
  htmlfile = HTML_DIR+'clusters.html'
  
  openw, outfile, htmlfile, /get_lun
  printf,outfile, '<!-- This file was generated by makecand.pro -->'

  openr, hdrfile, HEADER_FILE, /get_lun
  lin=''
  WHILE NOT EOF(hdrfile) DO BEGIN
      readf, hdrfile, lin
      printf, outfile, lin
  ENDWHILE 
  free_lun, hdrfile

  printf, outfile, "<P>This is a set of Jim Annis' clusters, matched to RASS.  They are all in the redshift range ["+ntostr(zmin,4)+', '+ntostr(zmax,4)+"].  They are ranked in order of X-ray luminosity. Here I have assumed H<SUB>0</SUB>=100 km/s/Mpc </P>"
  IF NOT keyword_set(fchart) THEN BEGIN 
      printf, outfile, '<P>Click on the rank to see an image of the field containing the cluster.  All clusters in the field are circled and labeled with z and ngals (the redshifts are uncorrected).  These images are from Steve Kents <a href="http://sdsslnx.fnal.gov:8015/">index of stripes</a><br><br>'
  ENDIF ELSE BEGIN
      printf, outfile, '<P>Click on the rank to see a color finding chart for the cluster. The finding chart is a 1 by 1 Mpc box at the redshift of the cluster (but it must be less than the size of a field, so boxes may be less than 1Mpc across for nearby clusters)<br><br></P>'
  ENDELSE 
;  printf, outfile, "<P>One way to tell bad candidates is looking at the image (click on rank).  Another is to compare X-ray luminosity and ngals.</P>"

  printf,outfile,''
  printf,outfile,'<TABLE border cellpadding=10>'
  printf,outfile,'<TR><TH> Rank            </TH>'+ $
                     '<TH> Z               </TH>'+ $
                     '<TH> L<SUB>X</SUB> (h<SUP>-2</SUP> ergs/s)     </TH>'+ $
                     '<TH> L<SUB>o</SUB> (h<SUP>-2</SUP> solar)      </TH>'+ $
                     '<TH> N<SUB>gals</SUB>           </TH>'+ $
                     '<TH> BCG RA (J2000) </TH>'+ $
                     '<TH> BCG DEC (J2000) </TH>'+ $
                 '</TR>'
  

  FOR i=0L, ncand-1 DO BEGIN
      index = cand[i]

      rank = ntostr(i+1)
;      z = ntostr( rnd(rass[index].z_orig, 2), 4)
      z = ntostr( rnd(rass[index].z, 2), 4)

;      IF keyword_set(linear) THEN BEGIN 
          xleft = ntostr( rass[index].l_x_bol ,4)
          xpower='44'
;      ENDIF ELSE BEGIN 
;          xpower = ntostr(long(rass[index].l_x_bol))
;          xleft = rass[index].l_x - long(rass[index].l_x)
;          xleft = ntostr( rnd(10.^xleft, 2), 4)
;      ENDELSE 
;      Opower = ntostr(long(rass[index].l_o))
;      Oleft = rass[index].l_o - long(rass[index].l_o)
;      Oleft = ntostr( rnd(10.^Oleft, 2), 4)

      IF nt NE 0 THEN BEGIN 
          Oleft = ntostr(rass[index].l_o/10., 4)
          Opower = '11'
      ENDIF ELSE BEGIN 
          Oleft = '1.0'
          Opower = '0.0'
      ENDELSE 

      ngals = ntostr( long(rass[index].(wgal[0])) )
      radecstr, rass[index].ra, rass[index].dec, rastr, decstr
      rastr = rastr+' ('+ntostr(rass[index].ra,7)+')'
      decstr = decstr+' ('+ntostr(rass[index].dec,6)+')'

      IF NOT keyword_set(fchart) THEN BEGIN 
          url = 'http://' + getkenturl(rass[index].run, $
                                       rass[index].camcol, $
                                       rass[index].field, clusters=clusters)
      ENDIF ELSE BEGIN
          ttt='RASS-'+run2string(rass[index].run)+'-'+ntostr(rass[index].camcol)+$
            '-'+ntostr(rass[index].rerun)+'-'$
            +field2string(rass[index].field)+'-'+ntostr(rass[index].id)+'.gif'
          url='./fchart/'+ttt
      ENDELSE 
      tmpstr = '<TR><TD nowrap>'+'<a href='+url+'>'+rank+'</A></TD>'+ $
                   '<TD nowrap>'+z+                  '</TD>' + $
                   '<TD nowrap>'+xleft+'e'+xpower+  '</TD>' + $
                   '<TD nowrap>'+Oleft+'e'+Opower+  '</TD>' + $
                   '<TD nowrap>'+ngals+             '</TD>' + $
                   '<TD nowrap>'+rastr+             '</TD>' + $
                   '<TD nowrap>'+decstr+            '</TD>' + $
               '</TR>'

      printf, outfile, tmpstr
  ENDFOR 

  printf,outfile,'</TABLE>'
  
  openr, ftrfile, FOOTER_FILE, /get_lun
  WHILE NOT EOF(ftrfile) DO BEGIN 
      readf,ftrfile,lin
      printf,outfile,lin
  ENDWHILE

  free_lun, ftrfile
  free_lun, outfile

  return 
END 
