
PRO mk_umidl_help_old,whichlib,list=list
;+
; NAME:
;   MK_UMIDL_HELP
; PURPOSE:
;   This program generates the large local help files (using RSI's
;   mk_html_help) for all installed software and also generates the master
;   index.html for the IDL page.
; CALLING SEQEUNCE:
;   mk_local_help 
; NOTES:
;   Modified from Eric Deutch's program. Erin Scott Sheldon ??/12/99
;-

  htmlfile='umich_idl.html'
  IDL_MAIN='/sdss3/products/idltools/um_idl/umsdss_idl/tools/'
  IDL_HTML_DIR='/sdss3/usrdevel/esheldon/WWW/weaklens/sdssidl/'
  HEADER_FILE='/sdss3/usrdevel/esheldon/WWW/weaklens/sdssidl/mk_local_help.header.html'
  FOOTER_FILE='/sdss3/usrdevel/esheldon/WWW/weaklens/sdssidl/mk_local_help.footer.html'
  TITLE='UMich IDL Homepage'

  libdirs $
    =[ $
       '', "General SDSS Library", 'General tools for SDSS data analysis',$
       'READPHOTO', "tsObj Reader Library", 'Tools for reading tsObj files', $
       'GETATLAS', "Atlas Image Reader Library", 'Tools for reading Atlas images', $
       'FLAGSELECT', "Flag Selection Library", 'Tools for selecting objects by their flags', $
       'UTIL',"Utilities", 'Utilities for general data manipulation', $
       'ROTSE', "ROTSE Routines",'Some helpful routines we borrowed from the ROTSE group at UM.  Contains image viewers, postscript makers and some good matching routines', $
       'ADMIN', "Administrative Routines", 'Some routines for setting up system variables', $
       'obsolete', "Obsolete Routines" , 'Routines no longer of general use'$
     ]


  IF (n_elements(list) NE 0) THEN BEGIN 
      FOR i=0,n_elements(libdirs)-2,2 DO BEGIN 
          print,i,'  ',libdirs(i),' ## ',libdirs(i+1)
      ENDFOR 
      return 
  ENDIF 

  stlib=0 & enlib=n_elements(libdirs)-2
  IF (n_elements(whichlib) NE 0) THEN BEGIN 
      stlib=whichlib & enlib=whichlib
  ENDIF 


  openw,outfile,IDL_HTML_DIR+htmlfile,/get_lun
  openw,outidx,'dirindex.dat',/get_lun
  printf,outfile,'<!-- This file was generated by mk_umidl_help.pro -->'
  
  openr,hdrfile,HEADER_FILE,/get_lun
  lin=''
  WHILE NOT EOF(hdrfile) DO BEGIN 
      readf,hdrfile,lin
      printf,outfile,lin 
  ENDWHILE 
  free_lun,hdrfile

  printf,outfile,'These are some procedures we and Dave Johnston of Chicago have written for retrieval and analysis of SDSS data.'
  printf,outfile,'You can get this code at Fermi lab through cvs.  The module name is sdssidl'

  printf,outfile,''
  printf,outfile,'<TABLE border cellpadding=10>'
 printf,outfile,'<TR> <TH> Library Name </TH><TH> Last Update </TH><TH> Comments </TH></TR>'
  
  FOR i=stlib, enlib,3 DO BEGIN 
      
      ctr=i/2+1
      libname='library'+strn(ctr,len=2,padchar='0')+'.html'

      printf,outidx,strn(ctr,len=2,padchar='0')+' '+IDL_MAIN+libdirs(i)
      print,'Building html file for library in ',+libdirs(i)

      mk_html_help,IDL_MAIN+libdirs(i),IDL_HTML_DIR+libname,title=libdirs(i+1)
      lastmod='-'

      IF (exist(IDL_MAIN+libdirs(i)+'/LAST_UPDATE')) THEN BEGIN 
          openr,5,IDL_MAIN+libdirs(i)+'/LAST_UPDATE'
          readf,5,lastmod
          close,5
      ENDIF 
;      libsource='-'
;      IF (exist(IDL_MAIN+libdirs(i)+'/SOURCE_PAGE')) THEN BEGIN 
;          openr,5,IDL_MAIN+libdirs(i)+'/SOURCE_PAGE'
;          readf,5,libsource
;          close,5
;      ENDIF 
      
      ;; link to the html help file.
      tmpstr='<TR> <TD nowrap> <A HREF="'+libname+'">'+libdirs(i+1)+'</A></TD> <TD nowrap> '+strn(lastmod)+'</TD><TD>'+libdirs(i+2)+' </TR>'
      
;      IF (strmid(libsource,0,4) EQ 'http') THEN $
;        tmpstr=tmpstr+' <TD> <A HREF="'+libsource+'">'+libsource+'</A> </TR>' $
;      ELSE IF (strmid(libsource,0,3) EQ 'ftp') THEN $
;        tmpstr=tmpstr+' <TD> <A HREF="'+libsource+'">'+libsource+'</A> </TR>' $
;      ELSE tmpstr=tmpstr+' <TD> '+libsource+' </TR>'

      printf,outfile,tmpstr

  ENDFOR 
  
  printf,outfile,'</TABLE>'
;  printf,outfile,'<HR size=3 noshade><BR>'
  dt=strmid(!stime,0,11)
  printf,outfile,'<H0><I>This page last updated on <strong>'$
    +dt+'</strong></I></H0><BR>'

  openr,ftrfile,FOOTER_FILE,/get_lun
  lin=''
  WHILE NOT EOF(ftrfile) DO BEGIN 
      readf,ftrfile,lin
      printf,outfile,lin
  ENDWHILE
  free_lun,ftrfile
  free_lun,outfile
  free_lun,outidx
end
