PRO ps_alpha, n, beta, alpha_implied

; without differential thing
alpha_implied = float(beta)*(float(n) - 9.)/6.

; with differential thing
;alpha_implied = float(beta)*(float(n) - 3.)/6. - 1.

return 
END 

PRO compare_alpha, beta, alpha_n1, alpha_n2, alpha_lit, higherr, lowerr

;; Look at power law part of schecter function

;; fit alpha from literature

alpha_lit = -1.52
higherr=0.11
lowerr=0.11
alpha_lit_high = alpha_lit + higherr
alpha_lit_low = alpha_lit - lowerr

;; implied alpha from Press Schecter formalism

print
print,'--------------------------------------------'
print,'alpha from literature: ',ntostr(alpha_lit),' (',ntostr(alpha_lit_high),', ',$
  ntostr(alpha_lit_low),')  68%'
n=-1.5
ps_alpha, n, beta, alpha_n1
print,'alpha (beta='+ntostr(beta,5)+', n=-1.5):  ',ntostr(alpha_n1,5)
n=-2.
ps_alpha, n, beta, alpha_n2
print,'alpha (beta='+ntostr(beta,5)+', n=-2.0):  ',ntostr(alpha_n2,5)


print,'--------------------------------------------'
print

return
END


PRO press_schecter_test

;; our power law M vs L^beta

beta = 0.9
beta_low = 0.4
beta_high = 1.5

print,'######## our best fit beta
compare_alpha, beta, m_alpha_n1, m_alpha_n2, alpha_lit, higherr, lowerr
print,'######## our low beta
compare_alpha, beta_low, m_alpha_n1_high, m_alpha_n2_high ;low-beta -> high-alpha
print,'######## our high beta'
compare_alpha, beta_high, m_alpha_n1_low, m_alpha_n2_low

;; expected beta

beta_exp = 0.75
print,'######## Expected beta from scaling laws'
compare_alpha, beta_exp, exp_alpha_n1, exp_alpha_n2

simpctable
!p.color=!black
!p.background = !white
plot,[0],/nodata, xrange=[0,5], yrange=[-3.0, -0.5]

;; plot literature result
xmin = 1.0
xmax = 3.0
ymax = alpha_lit + higherr
ymin = alpha_lit - lowerr
simpctable
;plot_box, xmin, xmax, ymin, ymax
polyfill, [xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax], color=!magenta, orientation=90.

;; plot allowed by our results
;xmin = 1.5
;xmax = 2.5
ymax = m_alpha_n2_high
ymin = m_alpha_n1_low
polyfill, [xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax],color=!red, orientation=0.

;; plot expected from scaling relations
;xmin = 1.0
;xmax = 3.0
ymin = exp_alpha_n1
ymax = exp_alpha_n2
polyfill, [xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax], orientation=45.

mess = ['lit','ours','exp']
legend,mess,line=[0,0,0],colors=[!magenta,!red,!black],/right,charsize=1.5

return
END 
