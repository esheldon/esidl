PRO testboot

  nchunk=10L
  nn=3000L
  
  ntot=nn*nchunk
  data = fltarr(ntot)

  FOR i=0L, nchunk-1 DO BEGIN 

      tdata=randomn(seed, nn)+(randomu(seed)-0.5)*10.0
      data[i*nn:(i+1)*nn-1] = tdata

  ENDFOR 

  plothist,data,bin=0.1
  oplot,[0,0],[0,100000]
  print,mean(data), sdev(data)/sqrt(ntot)


  runbootstrap, data, 1, 200, datamean, dataerr

  print,datamean, dataerr


END 
