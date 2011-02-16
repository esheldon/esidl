;+
; NAME:
;   MK_SDSSIDL_HELP
; PURPOSE:
;   This program generates help files for sdssidl procedure directories.
;   
; CALLING SEQEUNCE:
;   mk_local_help 
; NOTES:
;   Modified from Eric Deutch's program. Erin Scott Sheldon ??/12/99
;   Uses my program pro2html.pro to make a nice html version of each
;   procedure.
;-

PRO mk_sdssidl_read_desc, directory, title, description

  file = filepath(root=directory,'DESCRIPTION') 

  IF NOT fexist(file) THEN BEGIN 
      title = 'NO TITLE'
      description = 'NO DESCRIPTION'
  ENDIF ELSE BEGIN

      nlines = numlines(file)

      openr, lun, file, /get_lun

      line = ' '
      FOR i=0L, nlines-1 DO BEGIN 
          readf, lun, line

          line=strtrim(line)

          tmp = strsplit(line, /extract)
          ntmp = n_elements(tmp)
          IF ntmp GT 1 THEN BEGIN
              str = strjoin(tmp[1:ntmp-1],' ')
              CASE strlowcase(tmp[0]) OF
                  'title': title = str
                  'description': description = str
                  ELSE: 
              ENDCASE 
          ENDIF 
      ENDFOR 

      free_lun, lun

  ENDELSE 
      
END 

PRO mk_sdssidl_help_header, lun, title

	common mk_sdssidl_block, sdssidl_url, local_url

    defstyle = filepath(root=local_url, 'styles/default.css')
    homestyle = filepath(root=local_url, 'styles/home.css')
    navurl = filepath(root=local_url, 'javascript/navigation_bar.js')
    doturl = filepath(root=local_url, 'images/dot_clear.gif')

    printf,lun,'<HTML>'
    printf,lun,'<HEAD>'
    printf,lun,'<link rel="STYLESHEET" href="'+defstyle+'" type="text/css">'
    printf,lun,'<link rel="STYLESHEET" href="'+homestyle+'" type="text/css">'
    printf,lun,'<TITLE>',title,'</TITLE>'
    printf,lun,'</HEAD>'
    printf,lun,'<BODY>'
    printf,lun
    printf,lun,'<table border="0">'
    printf,lun,'  <tr>'
    printf,lun,'    <!-- Navigation bar -->'
    printf,lun,'    <td rowspan="2" width="100" align="center" valign="top" nowrap>'


    ;printf,lun,'        <a href="http://www.sdss.org" ><img src="http://sdssidl.sourceforge.net/images/SDSS_logo.gif" alt="SDSS Public Site" border="0" width="100" height="125"></a>'
    ;printf,lun,'        <p><a style="color:white; font-weight:bold" href="http://sdssidl.sourceforge.net/">SDSSIDL</a>'

    printf,lun,'      <script language="JavaScript" src="'+navurl+'"></script>'
    printf,lun,'      <noscript>You have JavaScript turned off</noscript>'
    printf,lun,'      </td>'
    printf,lun,'      <!-- Gutter. Without some content, the width=20 is ignored -->'
    printf,lun,'      <td width="20" rowspan="2">'
    printf,lun,'        <img src="'+doturl+'" width="25" height="1" alt="">'
    printf,lun,'      </td>'
    printf,lun,'        <!-- Main Body -->'
    printf,lun,'        <td valign="top">'
    printf,lun,'           <div align="left">'
    printf,lun,'               <H1>',title
    printf,lun,'                <HR size=3 noshade>'
    printf,lun,'                </H1>'
    printf,lun,'           </div>'

END 

pro mk_sdssidl_help_dinfo, lun, version

  dt=systime()
  printf,lun,'       <table border="0" width=100%>'
  printf,lun,'         <td aligh=left>'
  printf,lun,'           Documentation for sdssidl '+version
  printf,lun,'         </td>'
  printf,lun,'         <td align=right>'
  printf,lun,'           <div class="modified"> Last modified: '+dt+'</div>'
  printf,lun,'         </td>'
  printf,lun,'       </table>'
end

PRO mk_sdssidl_help_footer, lun, version

	common mk_sdssidl_block, sdssidl_url, local_url
    sigurl = filepath(root=local_url,'javascript/signature.js') 

  printf,lun,'      <hr size=3 noshade>'
  printf,lun,'      <!-- -------------- Signature --------------------- -->'

  mk_sdssidl_help_dinfo, lun, version
  
  printf,lun,'      <script language="JavaScript" src="'+sigurl+'">'
  printf,lun,'      </script>'
  printf,lun,'      <noscript>You have JavaScript turned off</noscript>'  
  printf,lun,'    </td>'
  printf,lun,'  </tr>'
  printf,lun,'</table>'
  
  printf,lun,'</BODY>'
  printf,lun,'</HTML>'


END 


PRO mk_sdssidl_cinfo, filenames, funcnames, comments


  ;; Create the html versions
  filenames = ['fileio/ascii_writeIDL.cpp','fileio/ascii_readIDL.cpp',$
               'fileio/binary_readIDL.cpp',$
               'pgsql/pgsql_query.c',$
               'atlas/src/read_atlasIDL.c',$
               'htmIndexIDL/idl/htmIndexIDL.cpp','htmIndexIDL/idl/htmIntersectIDL.cpp',$
               'htmIndexIDL/idl/htmMatchIDL.cpp','sdsspixIDL/idl/applyPixelMaskIDL.c',$
               'sphpoly_masks/sphPolyCompIDL.c','total_int/total_int.c']
  funcnames = ['ascii_write', 'ascii_read','binary_read', 'pgsql_query','rdAtlas', $
               'htm_index','htm_intersect','htmMatchC','sdsspix_mask',$
               'sphpoly_comp','total_int']

  comments = ['Efficient C code to write an IDL structure to an ascii file.  Uses the IDL DLM mechanism.<br>IDL> ascii_write, struct, file, /append, delimiter=, status=',$
              $
              'Efficient C code to read from an ascii file into an IDL structure.  Can extract individual columns and rows, unlike built-in readf.  Uses the IDL DLM mechanism.<br>IDL> struct = binary_read(file/lun, structdef, nrows_in_file, rows=, columns=, status=)',$
              $
              'C code to read from a binary file into an IDL structure.  Can efficiently extract individual columns and rows, unlike the build in readu procedure.  Uses the IDL DLM mechanism.<br>IDL> struct = binary_read(file/lun, structdef, nrows_in_file, rows=, columns=, status=)',$
              $
              'Interface to the postgres database. Uses the IDL DLM mechanism.<br>IDL> struct = pgsql_query(query, connect_info=, file=, /append, /nointerrupt, /verbose, status=)',$
              $
              'Luptons atlas reader library linked to IDL with the DLM mechanism',$
              'Uses Peter Z. Kunszt HTM library for linking to IDL. Returns the HTM index of ra/dec pairs',$
              'Returns all HTM triangles that are within a circle centered at the given ra/dec',$
              'Match two sets of ra/dec points using HTM libraries.  Linked to IDL with DLM mechanism. The user should use the IDL wrapper htm_match.pro',$
              'Ryan Scrantons sdsspix library linked to IDL with DLM mechanism. Added code to check edges', $
              'Andreas Berlinds spherical polygon library, linked to IDL with DLM mechanism.  returns completeness for each input ra/dec',$
              'Like built in total(), but works on integers without converting to float. IDL versions greater than 6 also support total(x, /int)']



    return

  filenames = ['fileio/ascii_writeIDL.cpp','fileio/ascii_readIDL.cpp',$
               'fileio/binary_readIDL.cpp']

  funcnames=['ascii_write','ascii_read', 'binary_read']
  comments = ['Efficient C code to write an IDL structure to an ascii file.  Uses the IDL DLM mechanism.<br>IDL> ascii_write, struct, file, /append, delimiter=, status=',$
              $
              'Efficient C code to read from an ascii file into an IDL structure.  Can extract individual columns and rows, unlike built-in readf.  Uses the IDL DLM mechanism.<br>IDL> struct = binary_read(file/lun, structdef, nrows_in_file, rows=, columns=, status=)', $
              'C code to read from a binary file into an IDL structure.  Can efficiently extract individual columns and rows, unlike the build in readu procedure.  Uses the IDL DLM mechanism.<br>IDL> struct = binary_read(file/lun, structdef, nrows_in_file, rows=, columns=, status=)']

END 


PRO mk_sdssidl_help, version

	; THIS IS OBSOLETE!!!!!!!!!!!  WE NOW USE WIKIDOC TO GENERATE WIKI
	; PAGES!!!!
    if n_elements(version) eq 0 then begin
        version = 'trunk'
    endif


	common mk_sdssidl_block, sdssidl_url, local_url

	sdssidl_url = 'https://sdssidl.googlecode.com/svn/'
	local_irl = '.'

    tmpdir = expand_path('~/tmp')
    cd,tmpdir

    sname = 'sdssidl-'+version
    maindir = filepath(root=tmpdir, sname)
    dname = 'sdssidl-doc-'+version
    docdir=filepath(root=tmpdir,dname)

    tarname=filepath(root=tmpdir, sname+'.tar.gz')
    dtarname=filepath(root=tmpdir, dname+'.tar.gz')

    if not fexist(maindir) then begin
		url = path_join(sdssidl_url, version)
		command = 'svn export '+url+' '+sname

        spawn, command, exit_status=exit_status
        if exit_status ne 0 then message,'command failed: '+command

    endif
    if not fexist(docdir) then begin
        command = 'svn export http://howdy.physics.nyu.edu/svn/sdssidl-doc/trunk'+' '+dname
        spawn, command, exit_status=exit_status
        if exit_status ne 0 then message,'command failed: '+command
    endif

    print,'sdssidl exported to: ',maindir
    print,'creating docs in : ',docdir

    idl_main=filepath(root=maindir,'pro')
    idl_src=filepath(root=maindir,'src')
    htmlfile='index.html'


    ; copy in some files 
    readme = filepath(root=maindir, 'README.txt')
    cpreadme = filepath(root=docdir, 'README.txt')
    file_copy, readme, cpreadme, /overwrite

    authors = filepath(root=maindir, 'AUTHORS')
    cpauthors = filepath(root=docdir, 'AUTHORS')
    file_copy, authors, cpauthors, /overwrite

    changes = filepath(root=maindir, 'CHANGES.txt')
    cpchanges = filepath(root=docdir, 'CHANGES.txt')
    file_copy, changes, cpchanges, /overwrite



    TITLE='SDSSIDL Homepage'

    if not fexist(IDL_MAIN) then begin
        message,'IDL_MAIN: ',IDL_MAIN,' not found'
    endif

    cd,IDL_MAIN


    ;; sig for pro2html
    emailurl = filepath(root=sdssidl_url, '../images/email-14pt.png')
    signature = '<img src="'+emailurl+'">'


    libdirs = file_search('','*',/test_directory)
    wkeep = where(strmatch(libdirs,'*CVS*') EQ 0 and $
        strmatch(libdirs,'*obsolete*') eq 0, nlib)

    libdirs = libdirs[wkeep]

    openw,outfile,filepath(root=docdir,htmlfile),/get_lun
    printf,outfile,'<!-- This file was generated by mk_sdssidl_help.pro -->'

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; print the header
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    mk_sdssidl_help_header, outfile, 'Umich/UChicago/NYU SDSSIDL Libraries'

    printf,outfile,$
        'SDSSIDL is a set of IDL and C/C++ routines for general data manipulation. There are also routines specifically disigned for analysis of data from the <a href="http://www.sdss.org">Sloan Digital Sky Survey</a> (SDSS). You can download the latest release here: '
    printf,outfile,'<p>'
    printf, outfile, '<ul><li><a href="http://sourceforge.net/projects/sdssidl">Download page for sdssidl</a></ul>'
    printf,outfile,'The current version is '+version
    ;   printf,outfile,'<a href="umsdss_idl_links/umsdss_idl_'+latest_update+'.tar.gz">umsdss_idl_'+latest_update+'.tar.gz</a>'
    ;   printf,outfile,'<p>'
    printf,outfile,'Here are the <a href="./README.txt">README</a>, <a href="CHANGES.txt">CHANGELOG</a>, and <a href="./AUTHORS">AUTHORS</a> documents for this distribution.'
    printf,outfile,'Browse the documentation for version '+version+' below'
    printf,outfile,'<p>'
    printf,outfile,'<ul>'
    printf,outfile,'  <li><a href="#PRO">The IDL code</a>'
    printf,outfile,'  <li><a href="#CCODE">The C/C++ code</a> linked to IDL dynamically via the DLM mechanism'
    printf,outfile,'  <li>A tutorial for the postgres interface is <a href="http://howdy.physics.nyu.edu/index.php/Postgres#IDL_interface">here</a>'
    printf,outfile,'</ul>'

    printf,outfile,'<p>'
    printf,outfile,'<a name="PRO">'
    printf,outfile,'<TABLE border cellpadding=10>'
    printf,outfile,'<CAPTION><EM>SDSSIDL Documentation</EM></CAPTION>'
    printf,outfile,'<TR> <TH> Library Name </TH><TH> Comments </TH></TR>'

    ntot = 0L
    FOR i=0L, nlib-1 DO BEGIN 

        libDir = libDirs[i]

        libname='library'+strn(i,len=2,padchar='0')+'.html'

        print,'Building html file for library in ',libDir

        ;; make a directory to put html files in
        HTML_LIB_DIR = filepath(root=docdir,libDir)
        IF NOT fexist(HTML_LIB_DIR) THEN file_mkdir, HTML_LIB_DIR

        ;; Read the description file
        mk_sdssidl_read_desc, libDir, title, description

        findir = filepath(root=IDL_MAIN, libDir)
        findstr = filepath(root=findir,'*.pro')
        files = file_search(findstr)
        nfile = n_elements(files) & ntot=ntot+nfile

        delvarx, pronames, htmlnames
        FOR fi=0L, nfile-1 DO BEGIN 
            dirsep, files[fi], junkdir, proname
            proname = ( str_sep(proname, '.pro') )[0]
            htmlname = proname + '.html'
            pro2html, files[fi], filepath(root=HTML_LIB_DIR,htmlname), /silent, signature=signature
            add_arrval, strupcase(proname), pronames
            add_arrval, htmlname, htmlnames
        ENDFOR 

        npro = n_elements(pronames)

        libfile = filepath(root=docdir,libname)
        openw, liblun, libfile, /get_lun

        mk_sdssidl_help_header, liblun, title
        printf, liblun, '<H2>List of Routines</H2><P>'
        printf, liblun, '<UL>'
        FOR np=0, npro-1 DO BEGIN 
            tname = filepath(root=libDir, htmlnames[np])
            printf, liblun, '<LI><A href="'+tname+'">'+pronames[np]+'</A>'
            get_pro_purpose, files[np], purp
            printf, liblun, purp
        ENDFOR 
        printf, liblun, '</UL><P>'
        ;      printf, liblun, '<HR>'


        mk_sdssidl_help_footer, liblun, version
        free_lun, liblun

        ;; link to the html help file.
        title = libname + ': '+title
        title = libDir
        tmpstr='<TR> <TD nowrap> <A HREF="'+libname+'">'+title+'</A></TD><TD>'+description+' </TR>'

        printf,outfile,tmpstr

    ENDFOR 

    print
    print,'Processed ',string(ntot,f='(i0)'),' .pro files'
    print

    printf,outfile,'</TABLE>'


    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Documentation for the c functions
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



    printf,outfile,'<p>'
    printf,outfile,'<a name="CCODE">'
    printf,outfile,'<TABLE border cellpadding=10>'
    printf,$
        outfile,'<CAPTION><em>C/C++ source code, linked via DLM<br>NOTE: These links are to the MAIN program sections only!</em></CAPTION>'


    mk_sdssidl_cinfo, filenames, funcnames, comments

    c2html = 'c2html'
    srchtml = filepath(root=docdir, 'src')
    if not fexist(srchtml) then file_mkdir, srchtml

    filenames = filepath(root=IDL_SRC,filenames)

    nf = n_elements(filenames) 
    FOR i=0L, nf-1 DO BEGIN 

        cfile = filenames[i]

        dirsep, cfile, tdir, tfile

		; now do docs for fileio differently
		if i le 2 then begin
			tfile = 'doc/'+funcnames[i]+'.txt'
		endif else begin
			copy = filepath(root=srchtml,tfile)
			html = copy + '.html'
			spawn, ['cp','-f',cfile,copy], /noshell
			spawn, [c2html,copy], /noshell
			tfile = tfile + '.html'
		endelse


        printf,outfile,$
            '<TR><TD nowrap>'+$
            '<a href="./src/'+tfile+'>'+funcnames[i]+'</a>'+$
            '<td>'+comments[i]+'</td></tr>'

    ENDFOR 


    printf,outfile,'</TABLE>'

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; print the footer
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    mk_sdssidl_help_footer, outfile, version

    free_lun,outfile

    cd,tmpdir

    print,'Tarring to file: ',tarname
    if fexist(tarname) then file_delete, tarname
    command = 'tar cvfz '+tarname+' '+sname
    spawn, command, oput, exit_status=exit_status
    if exit_status ne 0 then message,'command failed: '+command



    print,'Tarring to file: ',dtarname
    if fexist(dtarname) then file_delete, dtarname
    command = 'tar cvfz '+dtarname+' '+dname
    spawn, command, oput, exit_status=exit_status
    if exit_status ne 0 then message,'command failed: '+command

    print,'tar files: '
    print,'  '+tarname
    print,'  '+dtarname

end

