pro  red_allgood, oldstruct, newstruct


        IF N_PARAMS() eq 0 THEN BEGIN

                PRINT,"Syntax - red_allgood, oldstruct, newstruct"
                return
        END

only = where(oldstruct.m[1,*] ne -1 and oldstruct.m[2,*] ne -1 and oldstruct.m[3,*] ne -1 and $
		oldstruct.m[4,*] ne -1 and oldstruct.m[5,*] ne -1 and oldstruct.cz ne -1)
tot = n_elements(only)

newstruct = create_struct("M", dblarr(6,tot), "CZ", dblarr(tot),"CZERR",dblarr(tot) )

newstruct.m[0,*]=1.0
FOR i=1,5 DO BEGIN

	newstruct.m[i,*] = oldstruct.m[i,only]

ENDFOR

newstruct.cz[*] = oldstruct.cz[only]
newstruct.czerr[*] = oldstruct.czerr[only]

return
END
