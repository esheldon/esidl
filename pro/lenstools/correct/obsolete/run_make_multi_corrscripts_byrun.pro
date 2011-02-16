PRO run_make_multi_corrscripts_byrun, cheops1=cheops1, cheops2=cheops2, $
                                      cheops3=cheops3, latest=latest

  ;; should be run on the correct machine since that's where the
  ;; crun stuff will be created

  input_index = where(!run_status.tsObj_photo_v GE 5.4)
  get_goodruns, runs, reruns, stripes, strips, indices, $
    /tsObj, /asTrans, input_index=input_index, /silent
  help,runs,reruns
  maxmag = 22.5

  nTotal = n_elements(runs)
  nSplit = 100

  IF keyword_set(cheops1) THEN BEGIN 

      outFront = "cheops1"
      disk = "data0"
      runs   = runs[0:nSplit-1]
      reruns = reruns[0:nSplit-1]

  ENDIF ELSE IF keyword_set(cheops2) THEN BEGIN 

      outFront = "cheops2"
      disk = "data0"
      runs   = runs[nSplit:2*nSplit-1]
      reruns = reruns[nSplit:2*nSplit-1]

  ENDIF ELSE IF keyword_set(cheops3) THEN BEGIN 

      outFront = "cheops3"
      disk = "data0"
      runs   = runs[2*nSplit:nTotal-1]
      reruns = reruns[2*nSplit:nTotal-1]

  ENDIF ELSE IF keyword_set(latest) THEN BEGIN 

      outFront = "latest"
      disk = "data5"
      w=where(!run_status[indices].adatc_photo_v EQ -1, nlatest)
      IF nlatest EQ 0 THEN BEGIN 
          print,'All runs are corrected'
          return
      ENDIF ELSE BEGIN 
          runs = runs[w]
          reruns = reruns[w]
          print
          print,'Correcting these runs: '
          forprint,runs,reruns

      ENDELSE 
  ENDIF ELSE BEGIN 
      message,'Choose either /cheops1 or /cheops2 or /cheops3'
  ENDELSE 

  help,runs,reruns

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; photoz not ready yet
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  noNewCrun = 1
  addBayes = 1
;  addPhotoz=1

  make_multi_corrscripts_byrun, runs, reruns, outfront, maxmag, $
    noNewCRun=noNewCRun, addBayes=addBayes, addPhotoz=addPhotoz, $
    disk=disk


END 
