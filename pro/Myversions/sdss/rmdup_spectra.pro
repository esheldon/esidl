PRO rmdup_spectra, spectra, good

  ;; remove duplicates, preferably the one with target selection
  ;; flags set

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: test3, spectra, good'
      return
  ENDIF 

  ns=n_elements(spectra)
  good = lonarr(ns)
  ss = lindgen(ns)

  ten=ulong64(10)
  super=ulong64(spectra.objid[4])+ulong64(spectra.objid[3])*ten^6 +ulong64(spectra.objid[2])*ten^11 $
    +ulong64(spectra.objid[1])*ten^13 + ulong64(spectra.objid[0])*ten^16

  s=sort(super)
  super = temporary(super[s])
  ss = temporary(ss[s])

  FOR i=0L, ns-2L DO BEGIN 
      nn=i+1
      
      WHILE super[nn] EQ super[i] DO BEGIN 
          
          CASE 1 OF
              (spectra[ss[i]].primtarget NE 0): good[ss[nn]] = -1
              ELSE: good[ss[i]] = -1
          ENDCASE 
          nn=nn+1
          IF (nn GE ns) THEN GOTO,jump
      ENDWHILE 
      jump:
  ENDFOR 

  good = where(good NE -1)

END 
