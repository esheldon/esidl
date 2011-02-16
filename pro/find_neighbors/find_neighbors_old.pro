;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    FIND_NEIGHBORS
;       
; PURPOSE:
;    Tool for searching through neighbor list to find hosts (i.e. source/lens).
;       It searches by hosts which must be the first list input.
;    You must use lambda/eta survey coordinates.  Convert from ra/dec to 
;       lambda/eta using eq2survey before running FIND_NEIGHBORS
;
; CALLING SEQUENCE:
;    find_neighbors,lambda1,eta1,lambda2,eta2,searchrad,ind1,ind2,
;                   issouth=issouth,outfile=outfile,step=step,dist=dist
;
; INPUTS: lambda/eta for each list where lambda1 and eta1 are for list to be
;         search by
;         searchrad: search radius for each host in degrees
;    
;
; OPTIONAL INPUTS: 
;    
;
; KEYWORD PARAMETERS:
;    issouth: use for southern survey stripes
;    outfile: set to named file for output rather than indices output. (saves
;             memory)
;    step: set step for lens group (default step=300)
;    dist: set if you want distance from each source to its lens in degrees
;
; OUTPUTS: 
;    ind1: indices of matched list 1 repeated for each lens
;    ind2: indices of matched list 2
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
;    ZOBJSHEAR_LAMBDA_GET_SOURCES
;    ARRSCL
;    ECHO
;    MYGCIRC_SURVEY
;    PTIME
;
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    17-JAN-2002 zobjshear_lambda_looplens by Erin Sheldon
;    17-MAY-2002 find_neigbors search algorithm extracted by Judith Racusin    
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO plot_find_neighbors,lambda1,eta1,lambda2,eta2,ind1,ind2

  plot,lambda2,eta2,psym=3,xrange=[min(lambda2),max(lambda2)],yrange=[min(eta2),max(eta2)]
  oplot,lambda1,eta1,psym=1,color=!red
  x=rem_dup(lambda1[ind1])
  oplot,lambda1[ind1[x]],eta1[ind1[x]],psym=5,color=!blue
;oplot,lambda2[ind2],eta2[ind2],psym=3,color=!blue

return
END 


PRO find_neighbors_old, lambda1,eta1,lambda2,eta2,srad,ind1,ind2,dist,issouth=issouth,outfile=outfile,step=step,output_dist=output_dist,doplot=doplot

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: find_neighbors,lambda1,eta1,lambda2,eta2,searchrad,ind1,ind2,dist,'
      print,'                        issouth=issouth,outfile=outfile,step=step,'
      print,'                        /output_dist,/doplot' 
      print,'         Use eq2survey to put ra/dec into lambda/eta'
      return
  ENDIF 

  time=systime(1)
 
  ;; save the indices of the neighbors in 
  ;; a pointer (since we don't know how many there are ahead of time!)

;;;;;;;;sorting of lambda's for both lists
  s1=sort(lambda1)
  s2=sort(lambda2)

  nlist1 = n_elements(lambda1)

  numlist = lonarr(nlist1)
  ptrlist1 = ptrarr(nlist1)
  ptrlist2 = ptrarr(nlist1)
  IF keyword_set(output_dist) THEN BEGIN 
      ptrlist3 = ptrarr(nlist1)
      get_dist=1
  ENDIF ELSE get_dist=0
  ntotal = 0L
  
  IF n_elements(step) EQ 0 THEN step= 300
  nstepOld = nlist1/step
  nstep=nstepOld
  left = nlist1 MOD step

  ;; To account for leftover stuff
  IF nstepOld EQ 0 THEN BEGIN
      nstepOld = -10
      step = left
      nstep = 1
  ENDIF ELSE BEGIN
      IF left NE 0 THEN nstep = nstepOld + 1
  ENDELSE 

  indices = lindgen(n_elements(lambda1))

  FOR group = 0L, nstep-1 DO BEGIN

      IF group EQ nstepOld THEN BEGIN
          ii = indices[ group*step : group*step+left-1  ]
          step = left
      ENDIF ELSE BEGIN
          ii = indices[ group*step : (group+1)*step -1 ]
      ENDELSE 

      ;; Choose sources around this lens group
      
      srad2 = max(srad[s1[ii]])    
      maxii = max(ii) & minii = min(ii)

      maxlambda1 = lambda1[s1[maxii]]+srad2
      minlambda1 = lambda1[s1[minii]]-srad2

      zobjshear_lambda_get_sources, lambda2[s2], minlambda1, maxlambda1, wsrc, nwsrc;,issouth=issouth

      IF nwsrc NE 0 THEN BEGIN 
           
          echo, ['Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep)],$
                color = ['cyan','none'], bold=[1,0],nonewline=[1,0]
          help,/memory

          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group

              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = s1[ii[gi]]
              ceneta = eta1[index]
              cenlambda  = lambda1[index]
              
              srad_i = srad[index]
               
              ;; choose sources around this lens
              maxlambda2 = cenlambda+srad_i
              minlambda2 = cenlambda-srad_i

              zobjshear_lambda_get_sources, lambda2[s2[wsrc]], minlambda2, $
                maxlambda2,twsrc, nwsrc2,issouth=issouth              
              
              ;; If there are any source left, save them
              nkeep=0L
              IF nwsrc2 GT 1 THEN BEGIN 
                  wsrc2 = s2[wsrc[twsrc]]

                  mygcirc_survey, cenlambda, ceneta, $
                    lambda2[wsrc2], eta2[wsrc2], $
                    R, theta
                  keep=where(R LE srad_i,nkeep)

                  IF nkeep GT 0 THEN BEGIN
                      R=R[keep]
                      keep=wsrc2[keep]
                      numlist[index] = nkeep
                      ptrlist1[index] = $
                        ptr_new(replicate(index,nkeep),/no_copy)
                      ptrlist2[index] = ptr_new(keep, /no_copy)
                      IF get_dist THEN ptrlist3[index] = ptr_new(R,/no_copy)
                      ntotal = ntotal + nkeep
                  ENDIF

;                  IF nkeep GT 0 THEN BEGIN
;                     keep=wsrc2[keep]
;                     add_arrval, replicate(index,nkeep), ind1
;                     add_arrval, keep, ind2
;                  ENDIF 

              ENDIF 
              IF (nwsrc2 EQ 0) OR (nkeep EQ 0) THEN BEGIN 
                  echo,'/',color='green',/bold,/nonewline
              ENDIF 
          ENDFOR 
      ENDIF ELSE BEGIN 
          print,'No objects found for group'
      ENDELSE 
      print
  ENDFOR 
;stop
  s1=0 & s2=0
;stop

  IF n_elements(outfile) NE 0 THEN BEGIN
      print,'Writing indices to file: ',outfile
      openw,lun,outfile,/get_lun
      FOR i=0L, nlist1-1 DO BEGIN 
          IF numlist[i] NE 0 THEN BEGIN
              IF get_dist THEN $
                FOR j=0L,numlist[i]-1 DO printf,lun,(*ptrlist1[i])[j],(*ptrlist2[i])[j],(*ptrlist3[i])[j],format='(3(I0,:,1X))' $
              ELSE $
                FOR j=0L,numlist[i]-1 DO printf,lun,(*ptrlist1[i])[j],(*ptrlist2[i])[j],format='(2(I0,:,1X))'
          ENDIF 
      ENDFOR
      ptrarr_free,ptrlist1
      ptrarr_free,ptrlist2
      IF get_dist THEN ptrarr_free,ptrlist3
      free_lun,lun
  ENDIF ELSE BEGIN 

      ind1 = lonarr(ntotal)
      beg = 0L
      FOR i=0L, nlist1-1 DO BEGIN 
          IF numlist[i] NE 0 THEN ind1[beg:beg+numlist[i]-1] = *ptrlist1[i]
          ptr_free, ptrlist1[i]
          beg = beg + numlist[i]
      ENDFOR 

      ind2 = lonarr(ntotal)
      beg = 0L
      FOR i=0L, nlist1-1 DO BEGIN 
          IF numlist[i] NE 0 THEN ind2[beg:beg+numlist[i]-1] = *ptrlist2[i]
          ptr_free, ptrlist2[i]
          beg = beg + numlist[i]
      ENDFOR 

      IF get_dist THEN BEGIN 
          dist = fltarr(ntotal)
          beg = 0L
          FOR i=0L, nlist1-1 DO BEGIN 
              IF numlist[i] NE 0 THEN dist[beg:beg+numlist[i]-1] = *ptrlist3[i]
              ptr_free, ptrlist3[i]
              beg = beg + numlist[i]
          ENDFOR 
      ENDIF 
  ENDELSE 


  print
  ptime,systime(1)-time

  IF keyword_set(doplot) THEN plot_find_neighbors,lambda1,eta1,lambda2,eta2,ind1,ind2

return
END 
