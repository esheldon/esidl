;+
; NAME:
;   MK_UMIDL_HELP
; PURPOSE:
;   This program generates help files for umidl procedure directories.
;   
; CALLING SEQEUNCE:
;   mk_local_help 
; NOTES:
;   Modified from Eric Deutch's program. Erin Scott Sheldon ??/12/99
;   Uses my program pro2html.pro to make a nice html version of each
;   procedure.
;-


PRO mk_umidl_help_header,lun

  printf,lun,'<HTML>'
  printf,lun,'<HEAD>'
  printf,lun,'<link rel="STYLESHEET" href="http://cheops1.uchicago.edu/weaklens/styles/default.css" type="text/css">'
  printf,lun,'<link rel="STYLESHEET" href="http://cheops1.uchicago.edu/weaklens/styles/home.css" type="text/css">'
  printf,lun,'<TITLE>UMich/UChicago IDL Libraries</TITLE>'
  printf,lun,'</HEAD>'
  printf,lun,'<BODY>'
  printf,lun
  printf,lun,'<table border="0">'
  printf,lun,'  <tr>'
  printf,lun,'    <!-- Navigation bar -->'
  printf,lun,'    <td rowspan="2" width="100" align="center" valign="top" nowrap>'
  printf,lun,'      <script language="JavaScript" src="http://cheops1.uchicago.edu/weaklens/javascript/navigation_bar.js">'
  printf,lun,'      </script>'
  printf,lun,'      <noscript>You have JavaScript turned off</noscript>'
  printf,lun,'      </td>'
  printf,lun,'      <!-- Gutter. Without some content, the width=20 is ignored -->'
  printf,lun,'      <td width="20" rowspan="2">'
  printf,lun,'        <img src="http://cheops1.uchicago.edu/weaklens/images/dot_clear.gif" width="20" height="1" alt="">'
  printf,lun,'      </td>'
  printf,lun,'        <!-- Main Body -->'
  printf,lun,'        <td valign="top">'
  printf,lun,'           <div align="left">'
  printf,lun,'               <H1>UMich/UChicago IDL Libraires'
  printf,lun,'                <HR size=3 noshade>'
  printf,lun,'                </H1>'
  printf,lun,'           </div>'
  
END 

PRO mk_umidl_help_footer, lun

  dt=systime()

  printf,lun,'      <hr size=3 noshade>'
  printf,lun,'      <!-- -------------- Signature --------------------- -->'
  
  printf,lun,'       <table border="0" width=100%>'
  printf,lun,'         <td>'
  printf,lun,'           Back to the <a href="../idlhelp.html">Umich/Uchicago IDL Help Archives</a>'
  printf,lun,'         </td>'
  printf,lun,'         <td>'
  printf,lun,'           <div class="modified"> Last modified: '+dt+'</div>'
  printf,lun,'         </td>'
  printf,lun,'       </table>'

  printf,lun,'      <script language="JavaScript" src="http://cheops1.uchicago.edu/weaklens/javascript/signature.js">'
  printf,lun,'      </script>'
  printf,lun,'      <noscript>You have JavaScript turned off</noscript>'
  printf,lun,'      <!-- Begin Nedstat Basic code -->'
  printf,lun,'      <!-- Title: Umich/Uchicago IDL libraries -->'
  printf,lun,'      <!-- URL: http://http://sdss4.physics.lsa.umich.edu:8080/~esheldon/idlhelp/sdssidl/umich_idl.html -->'
  printf,lun,'      <script language="JavaScript" src="http://m1.nedstatbasic.net/basic.js">'
  printf,lun,'      </script>'
  printf,lun,'      <script language="JavaScript">'
  printf,lun,'      <!--'
  printf,lun,'      nedstatbasic("ABy98AbyY55MdYOQ++/xS0auB+3g", 0);'
  printf,lun,'      // -->'
  printf,lun,'      </script>'
  printf,lun,'      <noscript>'
  printf,lun,'      <a target="_blank" href="http://v1.nedstatbasic.net/stats?ABy98AbyY55MdYOQ++/xS0auB+3g"><img src="http://m1.nedstatbasic.net/n?id=ABy98AbyY55MdYOQ++/xS0auB+3g" border="0" nosave width="18" height="18"></a>'
  printf,lun,'      </noscript>'
  printf,lun,'      <!-- End Nedstat Basic code -->'
  
  printf,lun,'    </td>'
  printf,lun,'  </tr>'
  printf,lun,'</table>'
  
  printf,lun,'</BODY>'
  printf,lun,'</HTML>'


END 

PRO mk_umidl_help,whichlib,list=list

;  IDL_MAIN='/sdss7/products/idltools/um_idl/umsdss_idl/tools/'
;  IDL_HTML_DIR='/sdss7/home/esheldon/WWW/idlhelp/sdssidl/'

  IDL_MAIN='~/umsdss_idl/pro/'
  IDL_HTML_DIR='/net/cheops1/home/www/html/idlhelp/sdssidl/'
  htmlfile='umich_idl.html'

  TITLE='UMich IDL Homepage'

  libdirs $
    =[ $
       'sdss',            "General SDSS Library",            'General tools for SDSS data analysis',$
       'read_tsobj',      "SDSS tsObj Reader Library",       'Tools for reading SDSS tsObj files', $
       'spec1d',          "SDSS SPECTRO routines",           'Tools for reading and displaying SDSS spectra',$
       'atlas',           "SDSS Atlas Image Reader Library", 'Tools for reading SDSS Atlas images', $
       'flag_select',     "SDSS Flag Selection Library",          'Tools for selecting objects by varyous types of flags', $
       'struct',          "IDL Structures",                  'Tools for manipulating IDL structures', $
       'idlstruct_files', "IDL structure i/o",               "Tools for reading/writing idl structures", $
       'fileio',          "File i/o",                        'File input/output', $
       'util',            "Utilities",                       'Utilities for general data manipulation', $
       'plotting',        "Plotting Routines",               'Generic Plotting Routines', $
       'probgal',         "Star-Galaxy Separation",          "Bayesian Star-Galaxy Separation Routines", $ 
       'admin',           "SDSS Administrative Routines",    'Some routines for SDSS config files and system variables', $
       'fits',            "Fits tools",                      'Modified versions of idlastron library tools' $
     ]

  latest_update_file = IDL_HTML_DIR+'umsdss_idl_links/LATEST_UPDATE'

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
  printf,outfile,'<!-- This file was generated by mk_umidl_help.pro -->'
  
  ;; read in latest update
  openr, updatelun, latest_update_file, /get_lun
  latest_update = ' '
  readf, updatelun, latest_update
  free_lun, updatelun
  print,'Latest Update: '+latest_update

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; print the header
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   mk_umidl_help_header, outfile

;  openr,hdrfile,HEADER_FILE,/get_lun
;  lin=''
;  WHILE NOT EOF(hdrfile) DO BEGIN 
;      readf,hdrfile,lin
;      printf,outfile,lin 
;  ENDWHILE 
;  free_lun,hdrfile

;  spawn,'date +%d-%b-%Y',dt
;  dt=dt[0]

printf,outfile,'This page contains documentaion for some procedures we have written for general data manipulation, as well as retrieval and analysis of SDSS data. You can download a recent CVS build:<br>&nbsp;&nbsp;<a href="umsdss_idl_links/umsdss_idl_'+latest_update+'.tar.gz">umsdss_idl_'+latest_update+'.tar.gz</a>.<br>  And the C source code for the get_atlas routines:<br>&nbsp;&nbsp;<a href="umsdss_idl_links/get_atlas_'+latest_update+'.tar.gz">get_atlas_'+latest_update+'.tar.gz</a><br>&nbsp;&nbsp;<a href="umsdss_idl_links/get_atlas_photo_v5.3_'+latest_update+'.tar.gz">get_atlas_photo_v5.3_'+latest_update+'.tar.gz</a><br><br>These routines can also be checked out from our CVS archive directly using<br><br>&nbsp;&nbsp;"cvs&nbsp;&nbsp;-d&nbsp;&nbsp;:pserver:cvsuser@sdss7.physics.lsa.umich.edu:/sdss7/products/cvs/cvsroot&nbsp;&nbsp;checkout&nbsp;&nbsp;umsdss_idl"<br><br>You can update with "update" instead of "checkout" from the umsdss_idl directory<br>For the C code, the process is analagous. <br><b>NOTE: </b> You make have to do a cvs login before this will work on your machine the first time: <br>"cvs&nbsp;&nbsp;-d&nbsp;&nbsp;:pserver:cvsuser@sdss7.physics.lsa.umich.edu:/sdss7/products/cvs/cvsroot&nbsp;&nbsp;login"<br><br>'

  printf,outfile,''
  printf,outfile,'<TABLE border cellpadding=10>'
  printf,outfile,'<TR> <TH> Library Name </TH><TH> Comments </TH></TR>'
  
  ntot = 0L
  FOR i=stlib, enlib,3 DO BEGIN 
      
      ctr=i/2+1
      libname='library'+strn(ctr,len=2,padchar='0')+'.html'

      print,'Building html file for library in ',+libdirs(i)

      ;; make a directory to put html files in
      HTML_LIB_DIR = IDL_HTML_DIR+libdirs[i]+'/'
      IF NOT exist(HTML_LIB_DIR) THEN spawn,'mkdir '+HTML_LIB_DIR,answer

      findstr = IDL_MAIN+libdirs[i]+'/*.pro'
      files = findfile(findstr)
      nfile = n_elements(files) & ntot=ntot+nfile
      delvarx, pronames, htmlnames
      FOR fi=0L, nfile-1 DO BEGIN 
          dirsep, files[fi], tmpdir, proname
          proname = ( str_sep(proname, '.pro') )[0]
          htmlname = proname + '.html'
          pro2html, files[fi], HTML_LIB_DIR+htmlname, /silent        
          add_arrval, strupcase(proname), pronames
          add_arrval, htmlname, htmlnames
      ENDFOR 

      npro = n_elements(pronames)

      libfile = IDL_HTML_DIR+libname
      openw, liblun, libfile, /get_lun
      printf, liblun, '<HTML>'
      printf, liblun, '<HEAD>'
      printf, liblun, '<TITLE>'+libdirs[i+1]+'</TITLE>'
      printf, liblun, '</HEAD>'
      printf, liblun, '<BODY bgcolor="white" link="#0066ff" vlink="#ff0000" text="#00000">'
      printf, liblun, '<H1>'+libdirs[i+1]+'</H1><P>'
      printf, liblun, '<HR>'
      printf, liblun, '<H1>List of Routines</H1><P>'
      printf, liblun, '<UL>'
      FOR np=0, npro-1 DO BEGIN 
          printf, liblun, '<LI><A href="./'+libdirs[i]+'/'+htmlnames[np]+'">'+pronames[np]+'</A>'
          get_pro_purpose, files[np], purp
          printf, liblun, purp
      ENDFOR 
      printf, liblun, '</UL><P>'
      printf, liblun, '<HR>'
      printf, liblun, 'This file was generated by mk_umidl_help.pro.<br>'
      printf, liblun, '<b>Email: esheldon at cfcp.uchicago.edu</b><br>'
      lastmod=systime()
      printf, liblun, '<b>Last modified</b> '+lastmod+'<br>'
      printf, liblun, '</BODY>'
      printf, liblun, '</HTML>'
      free_lun, liblun
            
      ;; link to the html help file.
      tmpstr='<TR> <TD nowrap> <A HREF="'+libname+'">'+libdirs(i+1)+'</A></TD><TD>'+libdirs(i+2)+' </TR>'
      
      printf,outfile,tmpstr

  ENDFOR 
  
  print
  print,'Processed ',ntostr(ntot),' .pro files'
  print

  printf,outfile,'</TABLE>'
;  printf,outfile,'<HR size=3 noshade><BR>'
;  dt=systime()
;  printf,outfile,'<H0><I>This page last updated on <strong>'$
;    +dt+'</strong></I></H0><BR>'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; print the footer
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mk_umidl_help_footer, outfile

;  openr,ftrfile,FOOTER_FILE,/get_lun
;  lin=''
;  WHILE NOT EOF(ftrfile) DO BEGIN 
;      readf,ftrfile,lin
;      printf,outfile,lin
;  ENDWHILE
;  free_lun,ftrfile

  free_lun,outfile

end

