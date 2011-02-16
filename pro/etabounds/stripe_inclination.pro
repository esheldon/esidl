PRO stripe_inclination, stripe_in, inc

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: stripe_inclination, stripe, inc'
      return
  ENDIF 

  y_offset = -25d
  slope = 2.5d

  stripe = long(rnd(stripe_in))
  CASE stripe OF
      10: inc = 0.005888d
      76: inc = -15.0d
      82: inc = 0.00836007d
      86: inc = 10.0d
      ELSE: BEGIN
          IF (stripe GT 45) OR (stripe LT 1) THEN BEGIN
              message,'No such stripe: '+ntostr(stripe)
          ENDIF 
          inc = y_offset + slope*stripe
      END 
  ENDCASE
return

  CASE long(rnd(stripe)) OF
       9: inc = -2.5d
      10: inc = 0.005888d
      11: inc = 2.5d
      12: inc = 5.0d
      13: inc = 7.5d
      14: inc = 10.0d
      15: inc = 12.5d
      16: inc = 15.0d
      17: inc = 17.5d
      18: inc = 20.0d
      19: inc = 22.5d
      20: inc = 25.0d
      21: inc = 27.5d
      22: inc = 30.0d
      23: inc = 32.5d
      24: inc = 35.0d
      25: inc = 37.5d
      26: inc = 40.0d
      27: inc = 42.5d
      28: inc = 45.0d
      29: inc = 47.5d
      30: inc = 50.0d
      31: inc = 52.5d
      32: inc = 55.0d
      33: inc = 57.5d
      34: inc = 60.0d
      35: inc = 62.5d
      36: inc = 65.0d
      37: inc = 67.5d
      38: inc = 70.0d
      39: inc = 72.5d
      40: inc = 75.0d
      41: inc = 77.5d
      42: inc = 80.0d
      43: inc = 82.5d
      44: inc = 85.0d
      45: inc = 87.5d
      76: inc = -15.0d
      82: inc = 0.00836007d
      86: inc = 10.0d
      ELSE: message,'No info on stripe: '+ntostr(long(stripe))
  ENDCASE 

END 
