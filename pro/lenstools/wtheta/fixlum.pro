PRO fixlum, struct, clr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: fixlum, struct, clr'
      return
  ENDIF 

  lumfix = [1.314430862, 1.235023618, 1.194514446, 1.227677096, 1.220581657]
  
  lowdiff = abs(struct.norm[clr] - struct.normlow[clr])
  highdiff = abs(struct.normhigh[clr] - struct.norm[clr])

  struct.norm[clr] = struct.norm[clr]*lumfix[clr]
  struct.normlow[clr] = struct.norm[clr] - lowdiff*lumfix[clr]
  struct.normhigh[clr] = struct.norm[clr] + highdiff*lumfix[clr]
  CASE clr OF 
      1: BEGIN 
          struct.glumdense = struct.glumdense*lumfix[clr]
          struct.glumdenserr = struct.glumdenserr*lumfix[clr]
          struct.gtlumdense = struct.gtlumdense*lumfix[clr]
          struct.gtlumdenserr = struct.gtlumdenserr*lumfix[clr]
      END 
      2: BEGIN 
          struct.rlumdense = struct.rlumdense*lumfix[clr]
          struct.rlumdenserr = struct.rlumdenserr*lumfix[clr]
          struct.rtlumdense = struct.rtlumdense*lumfix[clr]
          struct.rtlumdenserr = struct.rtlumdenserr*lumfix[clr]
      END 
      3: BEGIN 
          struct.ilumdense = struct.ilumdense*lumfix[clr]
          struct.ilumdenserr = struct.ilumdenserr*lumfix[clr]
          struct.itlumdense = struct.itlumdense*lumfix[clr]
          struct.itlumdenserr = struct.itlumdenserr*lumfix[clr]
      END 
      4: BEGIN 
          struct.zlumdense = struct.zlumdense*lumfix[clr]
          struct.zlumdenserr = struct.zlumdenserr*lumfix[clr]
          struct.ztlumdense = struct.ztlumdense*lumfix[clr]
          struct.ztlumdenserr = struct.ztlumdenserr*lumfix[clr]
      END 
  ENDCASE 
END 
