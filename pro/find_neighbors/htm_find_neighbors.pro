;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    HTM_FIND_NEIGHBORS
;       
; PURPOSE:
;    Find the neighbors in radec list (ra2,dec2) around another set
;    of points (ra1,dec1) within radius searchrad.  For speed you should
;    make the first list the smaller of the two.  Uses HTM methods; all
;    objects in triangles that intersect the circle surrounding the position
;    in list 2 are kept, then objects outside the search radius are thrown out
;    (there will be some because intersecting triangles will be partly outside
;    the search radius).  This search method is N log(N) even for two
;    dimensional areas.
;
;
; CALLING SEQUENCE:
;    htm_find_neighbors,ra1,dec1,ra2,dec2,searchrad,ind1,ind2,
;                   outfile=outfile,leafids=leafids
;
; INPUTS: 
;         depth: level of HTM tree.  
;             depth   triangle area (square arcminutes)
;             -----------------------------------------
;               8     283
;               9      71
;              10      18
;              11       4.4
;              12       1.1
;              13       0.28 
;
;         ra1,dec1: list of ra's and dec's 
;         ra2,dec2: ""
;         searchrad: search radius for each object in list 1 in radians
;                    If all search radii will be 
;                    the same, the user may input just one number.
;
; OPTIONAL INPUTS: 
;    leafids: Input HTM leafids.  Saves time since leafids must be calculated.
;
; KEYWORD PARAMETERS:
;    outfile: set to named file for output rather than saving indices in
;             IDL variables.  This saves memory but is much slower.  
;             Of course, your going to write them to a file anyway....
;
;             Uses the writeu procedure to write unformatted binary data.
;             This is a very efficient way to write/read the data.     
;             One file is used for each set of indices.
;             The names are generated in the following way:
;
;                  frontfile = ( str_sep(outfile, '.fit') )[0]
;                  frontfile = ( str_sep(frontfile, '.fits') )[0]
;                  outfile1 = frontfile + '_ind1.bin'
;                  outfile2 = frontfile + '_ind2.bin'
;                  numfile = frontfile + '_num.dat'
;
;             If /convert is sent, then these files are read and output
;             as the first and second extensions of a .fit file called 
;             outfile.  
;
;   /convert: convert unformatted binary files to fits images, with ind1 and
;             ind2 in the 1st and 2nd extensions. This format is just as space
;             and speed efficient as those above, and is also platform
;             independent, so it is recommended that the user convert to fits.
;
;   /remove:   remove the unformatted binary files after converting to a fits
;              file. 
;   /output_dist: output the distance to each neighbor in radians. If
;                 outputting to a file, then it will have extension _dist.bin
;                 for the binary output, or will be in extension=2 in the fits
;                 file.
;
; OUTPUTS: (only output if outfile is not sent)
;    ind1: indices of matched list 1 
;    ind2: indices of matched list 2
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
;    htm_index() (C++ DLM)
;    htm_intersect() (C++ DLM)
;    (htm_bin2fits)
;
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    23-June-2002 Erin Scott Sheldon 
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO htm_find_neighbors, depth,ra1,dec1,ra2,dec2,srad,ind1,ind2,distance,$
                        outfile=outfile, convert=convert, remove=remove, $
                        step=step,leafids=leafids,output_dist=output_dist, $
                        numdefault=numdefault

  np = n_params()
  IF np LT 5 THEN BEGIN 
      print,'-Syntax: htm_find_neighbors, depth, ra1, dec1, ra2, dec2, srad, ind1, ind2, distance, outfile=outfile, convert=convert, remove=remove, step=step,leafids=leafids,/output_dist, numdefault=numdefault'
      return
  ENDIF 

  time=systime(1)
 
  nlist1 = n_elements(ra1)
  ;;pnlist1 = ntostr(nlist1+1)
  pnlist1 = ntostr(nlist1)
  imod = 1000
  dotmod = 100

  ;; should we return the distance?
  IF keyword_set(output_dist) THEN dodist = 1 ELSE dodist=0

  nrad = n_elements(srad)
  IF nrad EQ 1 THEN eachrad=0 ELSE BEGIN 
      IF nrad NE nlist1 THEN $
        message,'search radius must be scalar or same size as ra1'
      eachrad=1
  ENDELSE 

  IF n_elements(outfile) EQ 0 THEN BEGIN 
      ;; save the indices of the neighbors in 
      ;; a pointer (since we don't know how many there are ahead of time!)
      numlist = lonarr(nlist1)
      ptrlist1 = ptrarr(nlist1)
      ptrlist2 = ptrarr(nlist1)
      IF dodist THEN ptrlist_dist = ptrarr(nlist1)
      printout=0
  ENDIF ELSE BEGIN 
      ;; write to binary files
      ;; will convert to fits later
      printout=1

      frontfile = ( str_sep(outfile, '.fit') )[0]
      frontfile = ( str_sep(frontfile, '.fits') )[0]
      outfile1 = frontfile + '_ind1.bin'
      outfile2 = frontfile + '_ind2.bin'
      IF dodist THEN outfile3 = frontfile + '_dist.bin'
      numfile = frontfile + '_num.dat'

      print
      print,'Binary Files: '
      print,outfile1
      print,outfile2
      IF dodist THEN print,outfile3
      print,numfile
      IF keyword_set(convert) THEN BEGIN 
          print
          print,'Will convert to fits file: '+outfile
          IF keyword_set(remove) THEN print,'Will remove old files'
      ENDIF 
      print
      openw, lun1, outfile1, /get_lun
      openw, lun2, outfile2, /get_lun
      IF dodist THEN openw, lun3, outfile3, /get_lun

  ENDELSE 

  ntotal = 0L

  ;; find leaf id's for list 2
  IF n_elements(leafids) EQ 0 THEN BEGIN 
      print,'Looking up leaf ids for ra2/dec2 to depth = '+ntostr(depth)

      leafids = htm_index(ra2, dec2, depth)
      ;;htmlookupradec, ra2, dec2, depth, leafids

  ENDIF ELSE print,'Using Input leafids'
  print

  IF n_elements(leafids) EQ 1 THEN leafids = [leafids]

  ;; histogram leafids
  print,'Histograming leafids'
  minid = min(leafids, max=maxid)
  hist = histogram(leafids, min=minid, max=maxid, $
                   reverse_indices=rev_ind)
  nhist = n_elements(hist)
  print

  print,'Looping over lenses'
  FOR index=0L, nlist1-1 DO BEGIN 
              
      cenra = ra1[index]
      cendec = dec1[index]
      IF eachrad THEN searchrad = srad[index] ELSE searchrad=srad



      leaflist = htm_intersect(cenra, cendec, depth, searchrad)
      ;;htmintersectradec, cenra, cendec, searchrad, depth, leaflist, numdefault=numdefault

      nlist=n_elements(leaflist)

      ;; don't know how many objects match
      ;; need temporary pointers to deal with unknown
      tmp_ntotal=0L
      IF NOT printout THEN BEGIN 
          tmp_ptrlist = ptrarr(nlist)
          tmp_numlist = lonarr(nlist)
          IF dodist THEN tmp_ptrlist_dist = ptrarr(nlist)
      ENDIF 

      ;; see which leafids match
      FOR leaf=0L, nlist-1 DO BEGIN 
          
          binnum = leaflist[leaf] - minid
          IF (binnum LT nhist) AND (binnum GE 0) THEN BEGIN 

              ;; any matches?
              IF rev_ind[binnum] NE rev_ind[binnum+1] THEN BEGIN 
                  
                  keep = rev_ind[ rev_ind[binnum]:rev_ind[binnum+1]-1 ]
                  nkeep = n_elements(keep)

                  ;; Only keep those within radius (triangles that 
                  ;; overlap will most likely extend outside the circle)

                  mygcirc, cenra, cendec, ra2[keep], dec2[keep], dist, $
                           /radians_out
                  wtmp = where(dist LE searchrad, nkeep)
                  
                  ;; Are they within the circle?
                  IF nkeep NE 0 THEN BEGIN 

                      keep = keep[wtmp]
                      dist = dist[wtmp]
                      tmp_ntotal = tmp_ntotal+nkeep

                      IF NOT printout THEN BEGIN 
                          tmp_numlist[leaf] = nkeep
                          tmp_ptrlist[leaf] = ptr_new(keep,/no_copy)
                          IF dodist THEN tmp_ptrlist_dist[leaf] = ptr_new(dist, /no_copy)
                      ENDIF ELSE BEGIN 
                          rindex = replicate(index,nkeep)
                          writeu, lun1, rindex
                          writeu, lun2, keep
                          IF dodist THEN writeu, lun3, dist

                          keep=0
                          rindex=0
                          dist=0
                      ENDELSE 
                  ENDIF ;; in the circle?
              ENDIF ;; any in the leaf?
          ENDIF 
          
      ENDFOR 
      
      ntotal = ntotal + tmp_ntotal
      ;; collect the matches
      IF (NOT printout) AND (tmp_ntotal GT 0) THEN BEGIN 
          beg=0L
          tmp = lonarr(tmp_ntotal)
          IF dodist THEN tmp_dist = dblarr(tmp_ntotal)

          FOR leaf=0L, nlist-1 DO BEGIN 
              IF tmp_numlist[leaf] NE 0 THEN BEGIN 
                  tmp[beg:beg+tmp_numlist[leaf]-1] = *tmp_ptrlist[leaf]
                  IF dodist THEN tmp_dist[beg:beg+tmp_numlist[leaf]-1] = *tmp_ptrlist_dist[leaf]
              ENDIF 
              beg = beg+tmp_numlist[leaf]
          ENDFOR 
          ptrarr_free, tmp_ptrlist
          ptrarr_free, tmp_ptrlist_dist

          ;; now copy into bigger lists
          numlist[index] = tmp_ntotal
          ptrlist1[index] = ptr_new(replicate(index, tmp_ntotal), /no_copy)
          ptrlist2[index] = ptr_new(tmp, /no_copy)
          IF dodist THEN ptrlist_dist[index] = ptr_new(tmp_dist, /no_copy)

          tmp_numlist=0
      ENDIF 

      leaflist = 0
      IF ( index MOD imod ) EQ 0 THEN BEGIN
          print
          print,ntostr(index)+'/'+pnlist1
          print,'Total pairs: '+ntostr(ntotal)
          ;;help,/memory
      ENDIF 
      IF ( index MOD dotmod) EQ 0 THEN  print,'.',format='(a,$)'
  ENDFOR 
  print
  print,'Found '+ntostr(ntotal)+' Neighbors'
  print

  ;; free some memory
  rev_ind = 0

  ;; If we didn't print out the results, then convert pointer lists
  ;; into long arrays
  ind1 = -1L
  ind2 = -1L
  IF NOT printout THEN BEGIN 

      IF ntotal GT 0 THEN BEGIN 

          print,'Combining pointer arrays into match lists'
          print,'Match array 1'
          ind1 = lonarr(ntotal)
          ;;help,/memory
          beg = 0L
          FOR i=0L, nlist1-1 DO BEGIN 
              IF numlist[i] NE 0 THEN ind1[beg:beg+numlist[i]-1] = *ptrlist1[i]
              ptr_free, ptrlist1[i]
              beg = beg + numlist[i]
          ENDFOR 
          ptrarr_free, ptrlist1
          
          print,'Match array 2'
          ind2 = lonarr(ntotal)
          ;;help,/memory
          beg = 0L
          FOR i=0L, nlist1-1 DO BEGIN 
              IF numlist[i] NE 0 THEN ind2[beg:beg+numlist[i]-1] = *ptrlist2[i]
              ptr_free, ptrlist2[i]
              beg = beg + numlist[i]
          ENDFOR 
          ptrarr_free, ptrlist2
          
          IF dodist THEN BEGIN 
              print,'Distance array'
              distance = dblarr(ntotal)
              ;;help,/memory
              beg = 0L
              FOR i=0L, nlist1-1 DO BEGIN 
                  IF numlist[i] NE 0 THEN distance[beg:beg+numlist[i]-1] = *ptrlist_dist[i]
                  ptr_free, ptrlist_dist[i]
                  beg = beg + numlist[i]
              ENDFOR 
              ptrarr_free, ptrlist_dist
          ENDIF 
      ENDIF 

  ENDIF ELSE BEGIN

      free_lun, lun1, lun2
      IF dodist THEN free_lun, lun3
      openw, lun, numfile, /get_lun
      printf, lun, ntotal
      free_lun, lun

      IF keyword_set(convert) AND ntotal GT 0 THEN BEGIN 
          htm_bin2fits, outfile, remove=remove, output_dist=output_dist
      ENDIF 

  ENDELSE 
  print
  ptime,systime(1)-time
  ;;help,/memory

return
END 
