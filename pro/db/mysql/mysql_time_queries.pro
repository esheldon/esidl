PRO mysql_time_queries, nrows_arr, time1_arr, time2_arr, time3_arr

  run = 3605
  rerun = 20
  camcol = 2

  fieldmin = 20
  fieldmax = [30,40,50,60,70,80,90,100,120,130,140,150,160,170,180]
  nf = n_elements(fieldmax)
  nrows_arr = lonarr(nf)
  time1_arr = dblarr(nf)
  time2_arr = dblarr(nf)
  time3_arr = dblarr(nf)

  FOR i=0L, nf-1 DO BEGIN 

      print,'----------------------------------------------------'
      minid = photoid(run,rerun,camcol,fieldmin,1)
      maxid = photoid(run,rerun,camcol,fieldmax[i],0)

      tt1=systime(1)
;      query1='select objID,run,rerun,camcol,field,id from srcgal53 WHERE objid BETWEEN 3605020200020000001 AND 3605020200200000000'
      
      query1='select objID,run,rerun,camcol,field,id from srcgal53 WHERE objid BETWEEN '+ntostr(minid)+' AND '+ntostr(maxid)

      st1=mysql_query("srcgal",query1, nrows=nrows)
      nrows_arr[i] = nrows
      time1_arr[i] = systime(1)-tt1
      print,time1_arr[i]
      
      print,'Nrows = ',nrows

      ;; Now construct a query to get them
      orStatement = 'objID = '+ntostr(st1.objID)
      orStatement = strjoin(orStatement,' OR ')
      ;;print,orStatement
      
      query2 = $
        'select objID,run,rerun,camcol,field,id FROM srcgal53 WHERE '+$
        '  ('+orStatement+')'
      
      tt2=systime(1)
      st2 = mysql_query("srcgal", query2)
      time2_arr[i]=systime(1)-tt2
      print,time2_arr[i]
      
      queries3 = $
        'select objID,run,rerun,camcol,field,id FROM srcgal53 WHERE '+$
        '  objID = '+ntostr(st1.objID)
      tt3=systime(1)
      FOR j=0L, nrows-1 DO BEGIN 
          st3 = mysql_query("srcgal",queries3[j], nrows=nrows)
      ENDFOR 
      time3_arr[i]=systime(1)-tt3
      print,time3_arr[i]
      

  ENDFOR 

  oneclr = (!d.name EQ 'X') ? !green : !blue
  betweenclr = !red

  aplot,!gratio,nrows_arr,time2_arr, $
    xtitle='# of Rows', ytitle='Time (sec)', title='Search by objID'
  oplot,nrows_arr,time3_arr,color=oneclr, line=2
  oplot,nrows_arr,time1_arr,color=betweenclr, line=3

  legend,['Big OR','One at a time','Between'],$
    line=[0,2,3],color=[!p.color,oneclr,betweenclr]

return

END 
