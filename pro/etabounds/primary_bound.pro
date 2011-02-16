PRO primary_bound, stripe, bound

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: primary_bound, stripe, bound'
      return
  ENDIF 

;; Convert to LAMBDA/ETA coordinates from RA/DEC if necessary


bound = create_struct(name=type,'lammin',double(0.0),'lammax',double(0.0),$
			'etamin',double(0.0),'etamax',double(0.0))

stripe_inclination,stripe,inc

etamin = inc - 32.5 - 1.25
etamax = inc - 32.5 + 1.25

SWITCH stripe OF
    9: BEGIN
        lammin = -58.8
        lammax = 53.5
        BREAK
    END
    10: BEGIN
        etamin = -33.75
        etamax = -31.25
        lammin = -63.0
        lammax = 64.95
        BREAK
    END
    11: BEGIN
        lammin = -60.4
        lammax = 55.2
        BREAK
    END
    12: BEGIN
        lammin = -64.1
        lammax = 56.6
        BREAK
    END
    13: BEGIN
        lammin = -62.15
        lammax = 57.8
        BREAK
    END
    14: BEGIN
        lammin = -62.4
        lammax = 58.9
        BREAK
    END
    15: BEGIN
        lammin = -64.95
        lammax = 59.8
        BREAK
    END
    16: BEGIN
        lammin = -63.1
        lammax = 60.6
        BREAK
    END
    17: BEGIN
        lammin = -63.4
        lammax = 61.2
        BREAK
    END
    18: BEGIN
        lammin = -63.6
        lammax = 61.8
        BREAK
    END
    19: BEGIN
        lammin = -63.7
        lammax = 62.3
        BREAK
    END
    20: BEGIN
        lammin = -63.8
        lammax = 62.7
        BREAK
    END
    21: BEGIN
        lammin = -63.7
        lammax = 63.1
        BREAK
    END
    22: BEGIN
        lammin = -63.7
        lammax = 63.3
        BREAK
    END
    23: BEGIN
        lammin = -63.5
        lammax = 63.5
        BREAK
    END
    24: BEGIN
        lammin = -63.3
        lammax = 63.7
        BREAK
    END
    25: BEGIN
        lammin = -63.1
        lammax = 63.7
        BREAK
    END
    26: BEGIN
        lammin = -62.7
        lammax = 63.8
        BREAK
    END
    27: BEGIN
        lammin = -64.75
        lammax = 63.7
        BREAK
    END
    28: BEGIN
        lammin = -65.55
        lammax = 63.6
        BREAK
    END
    29: BEGIN
        lammin = -62.8
        lammax = 63.4
        BREAK
    END
    30: BEGIN
        lammin = -63.0
        lammax = 63.1
        BREAK
    END
    31: BEGIN
        lammin = -60.4
        lammax = 62.8
        BREAK
    END
    32: BEGIN
        lammin = -60.0
        lammax = 63.35
        BREAK
    END
    33: BEGIN
        lammin = -60.0
        lammax = 61.9
        BREAK
    END
    34: BEGIN
        lammin = -59.25
        lammax = 61.8
        BREAK
    END
    35: BEGIN
        lammin = -55.2
        lammax = 60.95
        BREAK
    END
    36: BEGIN
        lammin = -53.6
        lammax = 61.65
        BREAK
    END
    37: BEGIN
        lammin = -52.3
        lammax = 58.8
        BREAK
    END
    38: BEGIN
        lammin = -49.4
        lammax = 57.6
        BREAK
    END
    39: BEGIN
        lammin = -46.8
        lammax = 56.1
        BREAK
    END
    40: BEGIN
        lammin = -43.6
        lammax = 54.6
        BREAK
    END
    41: BEGIN
        lammin = -39.6
        lammax = 52.8
        BREAK
    END
    42: BEGIN
        lammin = -34.7
        lammax = 50.4
        BREAK
    END
    43: BEGIN
        lammin = -28.3
        lammax = 47.2
        BREAK
    END
    44: BEGIN
        lammin = -19.8
        lammax = 42.8
        BREAK
    END
    45: BEGIN
        lammin = -7.1
        lammax = 35.5
        BREAK
    END
    76: BEGIN
        etamin = 131.25
        etamax = 133.75
        lammin = -27.95
        lammax = 48.5
        BREAK
    END
    82: BEGIN
        etamin = 146.25
        etamax = 148.75
        lammin = -60.0
        lammax = 60.0
        BREAK
    END
    86: BEGIN
        etamin = 156.25
        etamax = 158.75
        lammin = -61.8
        lammax = 55.7
        BREAK
    END

ENDSWITCH

bound.lammin = lammin
bound.lammax = lammax
bound.etamin = etamin
bound.etamax = etamax

return

;; old is below

SWITCH stripe OF
    9: BEGIN
        lammin = -58.8
        lammax = 51.7
        BREAK
    END
    10: BEGIN
        etamin = -33.75
        etamax = -31.25
        lammin = -59.6
        lammax = 53.6
        BREAK
    END
    11: BEGIN
        lammin = -60.4
        lammax = 55.2
        BREAK
    END
    12: BEGIN
        lammin = -61.2
        lammax = 56.6
        BREAK
    END
    13: BEGIN
        lammin = -61.9
        lammax = 57.8
        BREAK
    END
    14: BEGIN
        lammin = -62.4
        lammax = 58.9
        BREAK
    END
    15: BEGIN
        lammin = -62.8
        lammax = 59.8
        BREAK
    END
    16: BEGIN
        lammin = -63.1
        lammax = 60.6
        BREAK
    END
    17: BEGIN
        lammin = -63.4
        lammax = 61.2
        BREAK
    END
    18: BEGIN
        lammin = -63.6
        lammax = 61.8
        BREAK
    END
    19: BEGIN
        lammin = -63.7
        lammax = 62.3
        BREAK
    END
    20: BEGIN
        lammin = -63.8
        lammax = 62.7
        BREAK
    END
    21: BEGIN
        lammin = -63.7
        lammax = 63.1
        BREAK
    END
    22: BEGIN
        lammin = -63.7
        lammax = 63.3
        BREAK
    END
    23: BEGIN
        lammin = -63.5
        lammax = 63.5
        BREAK
    END
    24: BEGIN
        lammin = -63.3
        lammax = 63.7
        BREAK
    END
    25: BEGIN
        lammin = -63.1
        lammax = 63.7
        BREAK
    END
    26: BEGIN
        lammin = -62.7
        lammax = 63.8
        BREAK
    END
    27: BEGIN
        lammin = -62.3
        lammax = 63.7
        BREAK
    END
    28: BEGIN
        lammin = -61.8
        lammax = 63.6
        BREAK
    END
    29: BEGIN
        lammin = -61.2
        lammax = 63.4
        BREAK
    END
    30: BEGIN
        lammin = -60.6
        lammax = 63.1
        BREAK
    END
    31: BEGIN
        lammin = -59.8
        lammax = 62.8
        BREAK
    END
    32: BEGIN
        lammin = -58.9
        lammax = 62.4
        BREAK
    END
    33: BEGIN
        lammin = -57.8
        lammax = 61.9
        BREAK
    END
    34: BEGIN
        lammin = -56.6
        lammax = 61.2
        BREAK
    END
    35: BEGIN
        lammin = -55.2
        lammax = 60.4
        BREAK
    END
    36: BEGIN
        lammin = -53.6
        lammax = 59.6
        BREAK
    END
    37: BEGIN
        lammin = -51.7
        lammax = 58.8
        BREAK
    END
    76: BEGIN
        etamin = 131.25
        etamax = 133.75
        lammin = -60.0
        lammax = 60.0
        BREAK
    END
    82: BEGIN
        etamin = 146.25
        etamax = 148.75
        lammin = -60.0
        lammax = 60.0
        BREAK
    END
    86: BEGIN
        etamin = 156.25
        etamax = 158.75
        lammin = -60.0
        lammax = 60.0
        BREAK
    END

ENDSWITCH

bound.lammin = lammin
bound.lammax = lammax
bound.etamin = etamin
bound.etamax = etamax

RETURN

END
