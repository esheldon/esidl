PRO maxbcg_reindex, bcg, neigh
  
;  on_error, 2
  nneigh = n_elements(neigh)
  nbcg = n_elements(bcg)
  IF nbcg EQ 0 OR nneigh EQ 0 THEN BEGIN 
      print,'-Syntax: maxbcg_reindex, bcg, neigh'
      print
      message,'Halting'
  ENDIF 

  keepneigh = lonarr(nneigh)

  ;; Ben indexed by stripe, then id within stripe
  ustripes = bcg[ rem_dup(bcg.stripe) ].stripe
  nstripe = n_elements(ustripes)
  
  FOR ist = 0L, nstripe-1 DO BEGIN 

      stripe = ustripes[ist]

      print,'//////////////////////////////////////////////////////////////'
      print,'Doing stripe: ',stripe


      ;; Must mach bcg stripe and bcg id
      wneighst = where( neigh.bcg_stripe EQ stripe, nneighst )
      wbcgst   = where( bcg.stripe EQ stripe, nbcgst)

      ;; histogram the neighbor bcgid
      print
      print,'Histogramming neighbor bcg_id'
      h = histogram(neigh[wneighst].bcg_id, $
                    min=0, max=max(bcg.bcg_id), $
                    rev=rev)

      print
      print,'Matching to bcgs'
      print

      FOR iibcg=0L, nbcgst-1 DO BEGIN 

          ibcg = wbcgst[iibcg]

          bcg_id = bcg[ibcg].bcg_id

          bcg[ibcg].bcg_id = ibcg
     
          ;; get neighbors
          IF rev[bcg_id] NE rev[bcg_id+1] THEN BEGIN 


              wmatch = rev[ rev[bcg_id]:rev[bcg_id+1]-1 ]
              wmatch = wneighst[wmatch]
              nmatch = n_elements(wmatch)
                  
              ;; re-index
              neigh[wmatch].bcg_id = ibcg
              keepneigh[wmatch] = 1
              
          ENDIF  

      ENDFOR ;; over bcgs in this stripe

  ENDFOR ;; over stripes

  keep = where(keepneigh NE 0, nkeep, comp=comp,  ncomp=ncomp)
  print
  print,'Neighbors that matched: '+ntostr(nkeep)+'/'+ntostr(nneigh)
  IF ncomp NE 0 THEN BEGIN 
      neigh[comp].bcg_id = -9999
  ENDIF 

END 
