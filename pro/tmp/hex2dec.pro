pro hex2dec,inp,out,quiet=quiet

;
;  trap invalid input
;
if datatype(inp) ne 'STR' then begin
   print,'Error: input must be string.'
   return
endif

;  
;  initialise output etc
;
out = 0L
n = strlen(inp)

;
;  convert each character in turn
;
for i=n-1,0,-1 do begin
  c = strupcase(strmid(inp,i,1))
  case c of
   'A': c = 10
   'B': c = 11
   'C': c = 12
   'D': c = 13
   'E': c = 14
   'F': c = 15
  else: begin
         if not valid_num(c,/integer) then begin
           print,'Invalid character **',c,'**'
           out = 0
           return
         endif
        end
  endcase
  out = out + long(c)*16L^long(n-1-i)
endfor

;
;  if not silenced, print result
;
if not keyword_set(quiet) then print,out

end
