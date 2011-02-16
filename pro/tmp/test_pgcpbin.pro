    
pro test_pgcpbin, struct, $
        maketables=maketables, check=check, clean=clean, $
        copy=copy

    pg=obj_new('postgres')
   
    if keyword_set(copy) then begin
        ; compare binary and ascii speeds on copy
        n=100000
        stdef = {id:0L, i:0, l:0L, sarr:['',''],larr:lonarr(2),f:0.0, d:0d}
        ;stdef = {id:0L, i:0, l:0L, larr:lonarr(2),f:0.0, d:0d}
        ;stdef = {sarr:['','']}
        st = replicate(stdef, n)
        st.id = lindgen(n)
        st.i = fix( 10000*randomu(seed,n) )
        st.l = long( 10000*randomu(seed,n) )

        st.sarr[0] = 'first-'+ntostr( long(10000*randomu(seed,n)) )
        st.sarr[1] = 'sec-'+ntostr( long(10000*randomu(seed,n)) )

        st.f = randomu(seed, n)
        st.d = randomu(seed, n, /double)

        st.larr[0] = long( 10000*randomu(seed,n) )
        st.larr[1] = long( 10000*randomu(seed,n) )

        pg->query, 'Drop table tcopy_binary'
        pg->query, 'Drop table tcopy_ascii'
        pg->struct2table, st, 'tcopy_binary', /binary
        pg->struct2table, st, 'tcopy_ascii'

    endif else if keyword_set(check) then begin
        ;pg->read_binary, '/tmp/testscalar.bin', n=1
        print,'------------------------------------------------'
        pg->read_binary, '/tmp/testlarr.bin', n=9
        print,'------------------------------------------------'
        pg->read_binary, '/tmp/testiarr.bin', n=3
        print,'------------------------------------------------'
        pg->read_binary, '/tmp/testbiarr.bin', n=9
        print,'------------------------------------------------'
        pg->read_binary, '/tmp/testfarr.bin', n=9
        print,'------------------------------------------------'
        pg->read_binary, '/tmp/testdarr.bin', n=9
        print,'------------------------------------------------'
        pg->read_binary_char, '/tmp/testchar.bin',1
        print,'------------------------------------------------'
        pg->read_binary_char, '/tmp/testvarchar.bin',1
        print,'------------------------------------------------'
        pg->read_binary_char, '/tmp/testtext.bin',1
        print,'------------------------------------------------'
        pg->read_binary_char, '/tmp/testtarr.bin',2,/array
        ;print,'------------------------------------------------'
        ;pg->read_binary, '/tmp/testarr2.bin', n=9
        ;print,'------------------------------------------------'
        ;pg->read_binary, '/tmp/testmatrix.bin', n=19
        ;print,'------------------------------------------------'
        ;pg->read_binary, '/tmp/pgtest.bin', n=19
    endif else if keyword_set(maketables) then begin


        test_pgcpbin, /clean
        ;    n=1

        ;    st1 = {l:5L}
        ;    pg->struct2table, st1, 'testscalar'

        table='testiarr'
        pg->query,"CREATE TABLE "+table+" (i smallint[2])"
        pg->query,"INSERT INTO "+table+" values('{5,5}')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

        table='testlarr'
        pg->query,"CREATE TABLE "+table+" (l int[2])"
        pg->query,"INSERT INTO "+table+" values('{5,5}')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

        table='testbiarr'
        pg->query,"CREATE TABLE "+table+" (l bigint[2])"
        pg->query,"INSERT INTO "+table+" values('{5,5}')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'


        table='testfarr'
        pg->query,"CREATE TABLE "+table+" (f real[2])"
        pg->query,"INSERT INTO "+table+" values('{5.0,5.0}')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

        table='testdarr'
        pg->query,"CREATE TABLE "+table+" (d double precision[2])"
        pg->query,"INSERT INTO "+table+" values('{5.0,5.0}')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

        table='testchar'
        pg->query,"CREATE TABLE "+table+" (c character(6))"
        pg->query,"INSERT INTO "+table+" values('hi')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

        table='testvarchar'
        pg->query,"CREATE TABLE "+table+" (c varchar)"
        pg->query,"INSERT INTO "+table+" values('hi')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

        table='testtext'
        pg->query,"CREATE TABLE "+table+" (c text)"
        pg->query,"INSERT INTO "+table+" values('hi')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

        table='testtarr'
        pg->query,"CREATE TABLE "+table+" (c text[2])"
        values = '{"hi","there"}'
        pg->query,"INSERT INTO "+table+" values('"+values+"')"
        values = '{"stuff","876543210"}'
        pg->query,"INSERT INTO "+table+" values('"+values+"')"
        pg->query, "COPY BINARY "+table+" to '/tmp/"+table+".bin'", $
            conn='user=postgres'

    endif else if keyword_set(clean) then begin
        ;pg->query, 'Drop table testscalar'
        pg->query, 'Drop table testiarr'
        pg->query, 'Drop table testlarr'
        pg->query, 'Drop table testbiarr'
        pg->query, 'Drop table testfarr'
        pg->query, 'Drop table testdarr'
        pg->query, 'Drop table testchar'
        pg->query, 'Drop table testvarchar'
        pg->query, 'Drop table testtext'
        pg->query, 'Drop table testtarr'
        pg->query, 'Drop table tcopy_binary'
        pg->query, 'Drop table tcopy_ascii'

        ;pg->query, 'Drop table testmatrix'
    endif

    obj_destroy, pg
end
