PRO red_data_two,list,cZ,X

   IF N_PARAMS() eq 0 THEN BEGIN

         PRINT,"Syntax - red_data_two, list, cZ, X"
           return
           END


   only = where(list.r ne -1 and list.g ne -1 and list.u ne -1 and $
                list.i ne -1 and list.z ne -1 and list.cz ne -1)
   b = size(list.r[only])
   X = dblarr(b[1],6)
   cZ = dblarr(b[1])   

   X[*,0]=1.0
   X[*,1]=list.r[only]
   X[*,2]=list.g[only]
   X[*,3]=list.u[only]
   X[*,4]=list.i[only]
   X[*,5]=list.z[only]

   cZ=list.cz[only]/3.0e8

return
END
