;+
; NAME:
;   MK_MYIDL_HELP
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

PRO mk_myidl_help_header, lun

  printf,lun,'<HTML>'
  printf,lun,'<HEAD>'
  printf,lun,'<link rel="STYLESHEET" href="http://cheops1.uchicago.edu/weaklens/styles/default.css" type="text/css">'
  printf,lun,'<link rel="STYLESHEET" href="http://cheops1.uchicago.edu/weaklens/styles/home.css" type="text/css">'
  printf,lun,'<TITLE>My IDL Libraries</TITLE>'
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
  printf,lun,'               <H1>My IDL Libraires'
  printf,lun,'               <HR size=3 noshade>'
  printf,lun,'               </H1>'
  printf,lun,'           </div>'

END

PRO mk_myidl_help_footer, lun

  printf,lun,'      <p>Back to the <a href="../idlhelp.html">Umich IDL Help Archives</a></p>'
  printf,lun,'      <hr size=3 noshade>'
  printf,lun,'      <!-- -------------- Signature --------------------- -->'
  
  printf,lun,'       <script language="JavaScript" src="http://cheops1.uchicago.edu/weaklens/javascript/signature.js">'
  printf,lun,'      </script>'
  printf,lun,'      <noscript>You have JavaScript turned off</noscript>'
  printf,lun,'    </td>'
  printf,lun,'  </tr>'
  printf,lun,'</table>'
  
  printf,lun,'</BODY>'
  printf,lun,'</HTML>'


END


PRO mk_myidl_help,whichlib,list=list


  htmlfile='myidlhelp.html'
  IDL_MAIN='~/idl.lib/'
  IDL_HTML_DIR='/net/cheops1/home/www/html/idlhelp/myidl/'

  TITLE='UMich IDL Homepage'

  libdirs = $
    ['AD_MOM',                           "Adaptive Moment Code",   "IDL versions of the adaptive moment code", $
     'APERTURE',                         "Aperture Mass",          "IDL code for measuring the Aperture Mass Statistic", $
     'BOOT',                             "BOOT",                   "Bootstrap code", $
     'DATABASES',                        'DATABASES',              "Database Code", $
     'FIT',                              "FIT",                    "Some Fitting Routines", $
     'Galclass',                         'Galclass',               "Galaxy Classification Stuff", $
     'HTML',                             'HTML',                   "Code for building HTML documents", $
     'KAPPAMAP',                         'KAPPAMAP',               "Code for making mass maps", $
     'MASS_EXTRACT',                     "MASS_EXTRACT",           "Peak analysis code for Bhuv's mass map statistics", $
     'MATH',                             "MATH",                   "Math related routines", $
     'Myversions/sdss/',                 "Myversions/sdss",        "Developmental versions of SDSSIDL archive code", $
     'Myversions/ADMIN/',                'Myversions/ADMIN',       "",$
     'Myversions/FLAGSELECT/',           'Myversions/FLAGSELECT',  "",$
     'Myversions/GETATLAS/',             'Myversions/GETATLAS',    "",$
     'Myversions/NEW/',                  'Myversions/NEW',         "",$
     'Myversions/READPHOTO/',            'Myversions/READPHOTO',   "",$
     'Myversions/ROTSE/',                'Myversions/ROTSE',       "",$
     'Myversions/UTIL/',                 'Myversions/UTIL',        "",$
     'Myversions/plotting/',             'Myversions/plotting',    "",$
     'NED',                              "NED",                    "Ned related code", $
     'PHOTOZ',                           'PHOTOZ',                 "Polynomial Photo-z code (out of date DSE's versions are better)", $
     'PLOTTING',                         'PLOTTING',               "Plotting procedures", $
     'QSOSELECT',                        'QSOSELECT',              "Some crappy code for selecting QSO's", $
     'RADEC_RANGE',                      'RADEC_RANGE',            "Code for finding lambda-eta range of fields", $
     'REGRESSGAL',                       'REGRESSGAL',             "The regression GGL code", $
     'ROSAT',                            'ROSAT',                  "Some code for extracting ROSAT sources", $
     'SDSSEXTRACT',                      'SDSSEXTRACT',            "Sextractor pipeline for analysis of SDSS frames or SDSS-like simulations"]
  libdirs = $
    [libdirs, $
     'SELECTCLUST',                        'SELECTCLUST',                       "Some HTML making code for Jim's clusters", $
     'SHAPESIM',                           'SHAPESIM',                          "Code for making images with fake galaxies for testing the shape measurement code", $
     'SHEAR',                              'SHEAR',                             "Some out-of-date shear measurement stuff", $
     'SHEARSIM',                           'SHEARSIM',                          "Simulation code for making fake shear", $
     'SIGMA_CRIT',                         'SIGMA_CRIT',                        "Code for calculating SIGMA_CRIT based on r-band magnitude", $
     'SIGMA_CRIT/photoz',                  'SIGMA_CRIT/photoz',                 "Code for calculating SIGMA_CRIT based on photoz", $
     'SIMSDSS',                            'SIMSDSS',                           "Code for simulating SDSS fields", $
     'SMOOTH',                             'SMOOTH',                            "Code for doing shear variance (name misleading)", $
     'UTIL',                               'UTIL',                              "Utilities for general data manipulation", $
     'VORONOI',                            'VORONOI',                           "Code for making voronoi maps", $
     'WIDGET',                             'WIDGET',                            "Some widget based code",$
     'compea4',                            'compea4',                           "Code to call modified versions of Hirata's code", $
     'denscont_sim',                       'denscont_sim',                      "simple simulation to illustrate mass-sheet degeneracy", $
     'download',                           'download',                          "Find runs to get by comparing with Pitt and Fermi", $
     'etabounds',                          'etabounds',                         "Routines for determining positions and bounds of stripes", $
     'find_neighbors',                     'find_neighbors',                    "HTM and other routines for finding neighbors", $
     'htm',                                'HTM',                               "Tools for building and using HTM trees", $
     'lenstools/analysis',                 'lenstools/analysis',                "Lensing analysis code", $
     'lenstools/analysis/plots',           'lenstools/analysis/plots',          "For making lensing plots", $
     'lenstools/combine/',                 'lenstools/combine',                 "Tools for combining runs/stripes", $
     'lenstools/combine_stripes',          'lenstools/combine_stripes',         "Tools for combining multiple stripes", $
     'lenstools/combine/photo53',          'lenstools/combine/photo53',         "Tools for combining runs/stripes from photo 5.3 catalogs", $
     'lenstools/combine_stripes',          'lenstools/combine_stripes',         "Tools for combining multiple stripes", $
     'lenstools/correct/',                 'lenstools/correct',                 "Tools for correcting runs", $
     'lenstools/correct/admomatlas/',      'lenstools/correct/admomatlas',      "Old C pipeline", $
     'lenstools/correct/fitmom',           'lenstools/correct/fitmom',          "Phil's old psf fitting routine wrappers", $
     'lenstools/ellip_combine',            'lenstools/ellip_combine',           "Tools for combining ellipticities", $
     'lenstools/ellip_errors',             'lenstools/ellip_errors',            "Code for testing Dave's ellipticity errors", $
     'lenstools/get_cat',                  'lenstools/get_cat',                 "Tools for getting lensing catalogs", $
     'lenstools/likelihood',               'lenstools/likelihood',              "Likelihood code for lensing analysis (obsolete)", $
     'lenstools/make_lenscat/daves_spec',  'lenstools/make_lenscat/daves_spec', "Code for making lens catalogs", $
     'lenstools/masks',                    'lenstools/masks',                   "Tools for making and applying masks", $
     'lenstools/maxbcg_lensing',           'lenstools/maxbcg_lensing',          "lensing tools for MaxBCG", $
     'lenstools/meanscinv',                'lenstools/meanscinv',               "Calculating mean inverse critical density", $
     'lenstools/newweight',                'lenstools/newweight',               "Tools for doing new lens weighting", $
     'lenstools/photoz/photoz_cat',        'lenstools/photoz/photoz_cat',       "Dealing with the photoz outputs", $
     'lenstools/photoz/photoz_dist',       'lenstools/photoz/photoz_dist',      "Generate photometric redshift distributions", $
     'lenstools/photoz/photoz_errors',     'lenstools/photoz/photoz_errors',    "Test propogating photoz errors (obsolete)", $
     'lenstools/photoz/photoz_lrg',        'lenstools/photoz/photoz_lrg',       "Creating lrg samples, analysis thereof", $
     'lenstools/probcuts',                 'lenstools/probcuts',                "Find good probability cuts", $
     'lenstools/psf_rec',                  'lenstools/psf_rec',                 "Tools for reconstructing the PSF from psField files",$
     'lenstools/rsmear_cuts',              'lenstools/rsmear_cuts',             "Define cuts in Rsmear", $
     'lenstools/shear',                    'lenstools/shear',                   "Tools for calculating shear", $
     'lenstools/specgal_lensing',          'lenstools/specgal_lensing',         "galaxy-galaxy lensing tools", $
     'lenstools/structs',                  'lenstools/structs',                 "Programs for defining lens structures", $
     'lenstools/vagc',                     'lenstools/vagc',                    "Programs for creating lens samples from LSS/VAGC", $
     'lenstools/wtheta',                   'lenstools/wtheta',                  "Programs for calculating w(theta)",$
     'mpfit',                              'mpfit',                             "From Markwardt Library",$
     'ngc',                                'ngc',                               "Make finding charts for NGC objects", $
     'pixel_masks',                        'Pixel Masks',                       "Wrappers for some pixel mask C routines", $
     'setup',                              'Setup',                             "My setup routine", $
     'spherical_polygons',                 'Spherical Polygons',                "Wrappers for some spherical polygon C routines", $
     'sxread',                             'sxread',                            "Tools for reading sx outputs", $
     'tmp',                                'tmp',                               "Temporary Code"]

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
;  openw,outidx,'dirindex.dat',/get_lun
  printf,outfile,'<!-- This file was generated by mk_umidl_help.pro -->'
  
  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; print header
  ;;;;;;;;;;;;;;;;;;;;;;;

  mk_myidl_help_header, outfile

;  openr,hdrfile,HEADER_FILE,/get_lun
;  lin=''
;  WHILE NOT EOF(hdrfile) DO BEGIN 
;      readf,hdrfile,lin
;      printf,outfile,lin 
;  ENDWHILE 
;  free_lun,hdrfile

  printf,outfile,'This is documentaion for my IDL code.  Get a version of my idl.lib <a href="myumich_idldir.tar">here</a>'

  printf,outfile,''
  printf,outfile,'<TABLE border cellpadding=10>'
  printf,outfile,'<TR> <TH> Library Name </TH><TH> Comments </TH></TR>'
  
  ntot = 0L
  ctr = 0L
  FOR i=stlib, enlib,3 DO BEGIN 

      findstr = IDL_MAIN+libdirs[i]+'/*.pro'
      files = findfile(findstr, count=count)
      ntot = ntot+count

      IF count NE 0 THEN BEGIN 

          libname='library'+strn(ctr,len=3,padchar='0')+'.html'
          
;          printf,outidx,strn(ctr,len=3,padchar='0')+' '+IDL_MAIN+libdirs(i)
          print,'Building html file for library in ',+libdirs(i)
          
          ;; make a directory to put html files in
          HTML_LIB_DIR = IDL_HTML_DIR+libdirs[i]+'/'
          IF NOT exist(HTML_LIB_DIR) THEN spawn,'mkdir '+HTML_LIB_DIR,answer

          delvarx, pronames, htmlnames
          
          FOR fi=0L, count-1 DO BEGIN 
              dirsep, files[fi], tmpdir, proname
              proname = ( str_sep(proname, '.pro') )[0]
              htmlname = proname + '.html'
              pro2html, files[fi], HTML_LIB_DIR+htmlname, /silent        
              add_arrval, strupcase(proname), pronames
              add_arrval, htmlname, htmlnames
          ENDFOR 
          
          npro = n_elements(pronames)
          
          lastmod = systime()
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
          printf, liblun, 'This file was generated by mk_myidl_help.pro.<br>'
          printf, liblun, '<b>Email: esheldon at cfcp.uchicago.edu</b><br>'
          printf, liblun, '<b>Last modified</b> '+lastmod+'<br>'
          printf, liblun, '</BODY>'
          printf, liblun, '</HTML>'
          free_lun, liblun
            
          ;; link to the html help file.
          tmpstr='<TR> <TD nowrap> <A HREF="'+libname+'">'+libdirs(i+1)+'</A></TD><TD>'+libdirs(i+2)+' </TR>'
          
          printf,outfile,tmpstr

      ENDIF ELSE print,'Not found: '+findstr

      ctr = ctr + 1

  ENDFOR 

  print
  print,'Processed ',ntostr(ntot),' .pro files'
  print

  printf,outfile,'</TABLE>'
;  printf,outfile,'<HR size=3 noshade><BR>'
  dt=systime()
  printf,outfile,'<H0><I>This page last updated on <strong>'$
    +dt+'</strong></I></H0><BR>'

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; print footer
  ;;;;;;;;;;;;;;;;;;;;;;;

  mk_myidl_help_footer, outfile

;  openr,ftrfile,FOOTER_FILE,/get_lun
;  lin=''
;  WHILE NOT EOF(ftrfile) DO BEGIN 
;      readf,ftrfile,lin
;      printf,outfile,lin
;  ENDWHILE
;  free_lun,ftrfile

  free_lun,outfile
;  free_lun,outidx

end

