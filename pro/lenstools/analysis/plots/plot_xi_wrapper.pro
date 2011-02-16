PRO plot_xi_wrapper

  FOR plotfits = 0, 1 DO BEGIN 
      FOR color=0, 1 DO BEGIN 
          FOR encap = 0, 1 DO BEGIN 
;              plot_typexi, /dops, encap=encap, color=color, plotfits=plotfits
              plot_allband_lumxi, 'all', /dops, encap=encap, color=color, plotfits=plotfits
              plot_allband_lumxi, 'redthree', /dops, encap=encap, color=color, plotfits=plotfits
          ENDFOR 
      ENDFOR 
  ENDFOR 
      



END 
