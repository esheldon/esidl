PRO postgres_getruns, scriptnum, runs, reruns, new=new

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: postgres_getruns, scriptnum, runs, reruns, /new'
      return
  ENDIF 

  IF keyword_set(new) THEN BEGIN 

      readcol,'~/tmp/getruns_original.dat',runs,reruns,format='I,I'
      s = sort(runs)
      runs = runs[s]
      reruns = reruns[s]
      

  ENDIF ELSE BEGIN 

      ;; 5042 just has empty files, how stupid can you get
      w=where(!run_status.tsobj_photo_v GT 5.4 AND $
              !run_status.run NE 5042, nw)
      
      runs = !run_status[w].run
      reruns = !run_status[w].rerun

  ENDELSE 

  nrun = n_elements(runs)
  w = lindgen(nrun)

  eachdef = nrun/4
  remainder = nrun - (eachdef*4)
  CASE scriptnum OF
      1: w = w[0:eachdef-1]
      2: w = w[eachdef:2*eachdef-1]
      3: w = w[2*eachdef:3*eachdef-1]
      4: w = w[3*eachdef:nrun-1]
      ELSE: message,'unknown script number: '+ntostr(scriptnum)
  ENDCASE 

  runs = runs[w]
  reruns = reruns[w]

END 
