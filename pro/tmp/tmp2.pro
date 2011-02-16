;+
;
; NAME:
;   match_multi 
;       
; PURPOSE:
;   Find the elements of an array of integers (arr2, which may contain
;   duplicates) which have a match in another array (arr1, which may not
;   contain duplicates).
;
; CALLING SEQUENCE:
;    match_multi, arr1, arr2, match2, numlist=, hist=, $
;       match1=, dupmatch1=, reverse_indices=
;
; INPUTS: 
;    arr1: integer array containing unique integers
;    arr2: integer array which may contain duplicates.
;
; OUTPUTS: 
;    match2: match indices. The indices of the second list that have matches
;          in the first. 
;
; OPTIONAL OUTPUTS:
;   numlist: the number that matched each element in arr1
;   hist: the histogram of arr2, for numbers 0 to max(arr1)
;   reverse_indices: the reverse indices from the histogram function on arr2
;   match1: The array of match indices from array 1.  This may not be the same
;       length as the match array for array 2.  See dupmatch1 for that.
;   dupmatch1: An array of matches from array 1 that correspond direclty
;       to all matches from array2.  E.g.
;           arr1 = [0,2,3,5,6]
;           arr2 = [2,3,3,4,5,6,6,7]
;           match_multi, arr1, arr2, m2, match1=m1, dupmatch1=dm1
;             1           0
;             2           1
;             2           2
;             3           4
;             4           5
;             4           6
;
; CALLED ROUTINES:
;    histogram
; 
; PROCEDURE: 
;    see code
;	
;
; REVISION HISTORY:
;    Jun-18-2002: Erin Scott Sheldon UofMich
;       Adapted from the old match2rand program
;                                      
;-                                       
;  Copyright (C) 2005  Erin Sheldon, NYU.  erin dot sheldon at gmail dot com
;
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, version 2.
;

pro match_multi, arr1, arr2, match2, match1=match1, dupmatch1=dupmatch1, numlist=numlist, hist=hist, reverse_indices=reverse_indices

    if n_params() lt 2 then begin 
        print,'-Syntax: match_multi, arr1, arr2, match2, numlist=numlist, hist=hist, reverse_indices=reverse_indices'
        print,'This works even for subsamples of subsamples'
        return
    endif 

    n1 = n_elements(arr1)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; using histogram this way, with min=0, forces a bin for every
    ;; integer from 0 to the max. that way, we can use
    ;; arr1[i] as subscript for reverse_indices! then it
    ;; works for sub-samples
    ;;
    ;; can use lots of memory, though, if the numbers in arr1 span
    ;; a large range, since memory usage is proportional to 
    ;; sizeof(long)*span
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    hist=histogram(arr2, reverse_indices=reverse_indices, $
        min=0, max=max(arr1) )

    ptrlist = ptrarr(n1)
    numlist = lonarr(n1)
    nmatch = 0l

    if arg_present(dupmatch1) then begin
        ptrlist1 = ptrarr(n1)
    endif
    if arg_present(match1) then begin
        match1 = lindgen(n1) 
    endif

    for i=0l, n1-1 do begin 

        index1 = arr1[i]

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Use of histogram with min=0 is why this works for 
        ;; subsamples of subsamples: when you do histogram
        ;; on integers with bin=1, min=0, bins are created
        ;; for each integer from 0 to max(zindex), and 
        ;; we can just subscript reverse_indices with zi. 
        ;; this could be bad for memory if all the integers
        ;; in arr1 are large
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        if reverse_indices[index1] ne reverse_indices[index1+1] then begin 

            s2 = reverse_indices[ $
                reverse_indices[index1]:reverse_indices[index1+1] -1 $
                ]
            ns2 = n_elements(s2)
            ptrlist[i] = ptr_new(s2)
            numlist[i] = ns2
            nmatch = nmatch + ns2

            if arg_present(dupmatch1) then begin
                ptrlist1[i] = ptr_new(replicate(i, ns2))
            endif
        endif else begin
            if arg_present(match1) then match1[i] = -1
        endelse

    endfor 

    match2 = lonarr(nmatch)
    if arg_present(dupmatch1) then begin
        dupmatch1 = lonarr(nmatch)
    endif
    beg=0l
    for i=0l, n1-1 do begin 
        if numlist[i] ne 0 then begin 
            match2[beg:beg+numlist[i]-1] = *ptrlist[i]
            if arg_present(dupmatch1) then begin
                dupmatch1[beg:beg+numlist[i]-1] = *ptrlist1[i]
            endif
        endif 
        ptr_free, ptrlist[i]
        if arg_present(dupmatch1) then ptr_free, ptrlist1[i]
        beg = beg + numlist[i]
    endfor 

    if arg_present(match1) then begin
        w=where(match1 ne -1, nw)
        if nw ne 0 then match1 = match1[w] else match1=-1
    endif

    return
end 



pro tmp2

arr1 = [0,2,3,5,6]
arr2 = [2,3,3,4,5,6,6,7]

match_multi_test, $
    arr1, arr2, m2, match1=m1, dupmatch1=dm1

colprint, dm1, m2
colprint, m1

;Example output:
;m1 = [1,2,2,3,4,4]
;m2 = [0,1,2,4,5,6]

;Current output:
;m1 = [1,2,3,4]
;m2 = [0,1,2,4,5,6]

end
