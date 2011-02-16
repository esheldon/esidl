pro collate_plate,pls,platef=platef,platen=platen,ln=ln, plot=plot

  IF NOT keyword_set(platef) then begin
      platef='/data1/spectra/1d_15/tsObj/tsObj-0'+strtrim(ntostr(platen),2)+$
        '.fit'
  ENDIF 

    pl=mrdfits(platef,1,/silent)

    newtags = ['z','zerr',$
               'e1','e2','momerr','r',$
               'ixx','iyy','ixy','rho4',$
               'seeing','rotation','matchrerun']

    edef = '[-1.e10,-1.e10,-1.e10,-1.e10,-1.e10]'
    tagvals = ['0.0','0.0',$
               edef,edef,edef,edef,$
               edef,edef,edef,edef,$
               'fltarr(5)','fltarr(5)','-1L']
    add_tags, pl, newtags, tagvals, pls

    IF keyword_set(plot) THEN plot,pl.ra,pl.dec,psym=1,/ynozero,title=platef

    runs=pl(rem_dup(pl.run)).run
    nruns=n_elements(runs)
    print
    print,'Runs in plate: ',runs

    make_aligntags,tl
 
    make_runstatus_struct,rs
    rs.adatc_exist = 'Y'
    runstatus_select, rs, good

    ;; match reruns to !run_status

    for i=0,nruns-1 do BEGIN
        print
        print,'Processing run: ',ntostr( runs[i] )
        inr=where(pl.run eq runs(i))

        ;; Look for run
        wr=where(!run_status[good].run eq runs(i), nw)
        if (nw ne 0) then begin
            wr=good[wr]

            ;; now see what reruns are in plate (should only be one)
            reruns = pl[rem_dup(pl[inr].rerun)].rerun
            nrerun = n_elements(reruns)
            FOR ir=0L, nrerun-1 DO BEGIN 
                ;; do we have this rerun in !run_status?
                rerun = reruns[ir]
                wrer=where(!run_status[wr].rerun EQ rerun, nrer)
                IF nrer NE 0 THEN BEGIN 
                    matchrerun = rerun
                ENDIF ELSE BEGIN 
                    ;; no match, so just close_match with the
                    ;; highest available rerun
                    matchrerun = fix( max(!run_status[wr].rerun) )
                ENDELSE 

                print,'Rerun: ',ntostr(rerun),' Match rerun: ',ntostr(matchrerun)
                ;; which columns do we have?
                cols=pl[inr[rem_dup(pl[inr].camcol)]].camcol
                ncols=n_elements(cols)
                            
                for j=0,ncols-1 do begin
                    print,'Processing column: ',ntostr(cols[j])
                    inc=where(pl(inr).camcol eq cols(j))
                    inc=inr(inc)
                    fields=pl(inc(rem_dup(pl(inc).field))).field
                    nfields=n_elements(fields)
                
                    fetch_dir,runs(i),cols(j),matchrerun,$
                              dir,atldir,corrdir=corrdir,/check
                    if (dir ne '' and corrdir ne '') then begin
                        for k=0,nfields-1 do begin
                            
                            read_tsobj,corrdir,ln,start=fields(k),nframes=1,$
                                       front='adatc',tsobjstr=tsobjstr,$
                                       taglist=tl,verbose=0,status=status
                            if (status eq 0) then begin
                                if keyword_set(plot) then $
                                  oplot,ln.ra,ln.dec,psym=3
                                infield=where(pl.field eq fields(k) and $
                                              pl.camcol eq cols(j) $
                                              and pl.run eq runs(i))
                                if (matchrerun EQ rerun) then begin
                                    sphoto_match,pl(infield),ln,m1,m2
                                    IF ((m2(0) NE -1) AND $
                                        keyword_set(plot) ) THEN BEGIN 
                                        oplot,[ln(m2).ra],[ln(m2).dec],psym=4,$
                                              color=!green
                                    ENDIF 
                                endif else begin
                                    close_match_radec,pl(infield).ra,$
                                                      pl(infield).dec,$
                                                      ln.ra,ln.dec,$
                                                      m1,m2,0.0004,1,/silent
                                    IF ((m2(0) NE -1) AND $
                                        keyword_set(plot) ) THEN BEGIN
                                        oplot,[ln(m2).ra],[ln(m2).dec],psym=7,$
                                              color=!red
                                    ENDIF 
                                endelse
                                if (m1(0) ne -1) then begin
                                    m1=infield(m1)
                                    pls(m1).e1=ln(m2).e1
                                    pls(m1).e2=ln(m2).e2
                                    for ind=0,4 do begin
                                        bad=where(abs(pls.e1(ind)) gt 1.0)
                                        if (bad(0) ne -1) then $
                                          pls(bad).e1(ind) = 1.e10
                                    endfor
                                    for ind=0,4 do begin
                                        bad=where(abs(pls.e2(ind)) gt 1.0)
                                        if (bad(0) ne -1) then $
                                          pls(bad).e2(ind) = 1.e10
                                    endfor
                                    pls(m1).momerr=ln(m2).momerr
                                    pls(m1).rho4=ln(m2).rho4
                                    pls(m1).r=ln(m2).r
                                    pls(m1).ixx=ln(m2).ixx
                                    pls(m1).iyy=ln(m2).iyy
                                    pls(m1).ixy=ln(m2).ixy
                                    pls(m1).rotation=ln(m2).rotation
                                    ;; E.S.S.
                                    pls[m1].seeing = ln[m2].seeing
                                    pls[m1].matchrerun = matchrerun
                                ENDIF 
                            ENDIF 
                        ENDFOR 
                    ENDIF 
                ENDFOR 
            ENDFOR 
        endif else begin
            print,'No data available for run:'+ntostr(runs(i))
        endelse
    endfor
    
    pls.z    = pl.zfinal
    ;; E.S.S.
    pls.zerr = pl.zerrfinal
    
    return
    end
