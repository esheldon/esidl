FUNCTION quadrant_integral, d, R
  
  return, 0.5*d*sqrt(R^2-d^2) + 0.5*R^2*atan(d, sqrt(R^2-d^2) )

END 

FUNCTION quadrant_area, d1, d2, radius,hit=hit

  hit=1
  dfromcorner = sqrt(d1^2 + d2^2)
  CASE 1 OF
      dfromcorner GE radius: BEGIN ;; circle not outside corner
          CASE 1 OF
              (d1 LT radius) AND (d2 LT radius): BEGIN 
                  ;; just corner outside circle, both sides crossed
                  area = quadrant_integral(d1, radius) - $
                    (0.25*!pi*radius^2 - quadrant_integral(d2, radius))
              END 
              (d1 LT radius): BEGIN 
                  ;; only one side crossed
                  area=quadrant_integral(d1, radius)
              END 
              (d2 LT radius): BEGIN 
                  ;; only one side crossed
                  area=quadrant_integral(d2, radius)
              END 
              ELSE: BEGIN 
                  ;; neither side crossed
                  hit=0
                  area=0.25*!pi*radius^2
              END 
          ENDCASE 
      END 
      dfromcorner LT radius: BEGIN ;; circle passes ouside corner
          area = d1*d2
      END 
  ENDCASE 

  return,area

END 

PRO circ_rect_intersect, dfromtop, dfrombottom, dfromright, dfromleft, $
                              radius, area, hit=hit

  IF n_params() LT 5 THEN BEGIN 
      print,'Returns the area of intersection between circle and rectangle'
      print,'-Syntax: circ_rect_intersect,dfromtop, dfrombottom, dfromright, dfromleft, radius, area'
      return
  ENDIF 

  ;; calculate the area of intersection between a circle and
  ;; a rectangle given the distances from center to edges
  ;; center of circle must be within the rectangle

  ;; quadrant 1 (top right)
  ;; since center is within rectangle only need to consider distance
  ;; to right side and top side. similar for the others.
  
  area1 = quadrant_area(dfromtop, dfromright, radius,hit=hit1)
  area2 = quadrant_area(dfromtop, dfromleft, radius,hit=hit2)
  area3 = quadrant_area(dfrombottom, dfromleft, radius,hit=hit3)
  area4 = quadrant_area(dfrombottom, dfromright, radius,hit=hit4)

  hit = hit1+hit2+hit3+hit4 < 1
  area= area1+area2+area3+area4

  return

END 
