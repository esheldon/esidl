pro maxbcg_m2l_paperplots

    c2 = obj_new('maxbcg_correlate', 2)
    c4 = obj_new('maxbcg_correlate', 4)
    m2l = obj_new('maxbcg_m2l', [21,22], 4)

    c4->plot_usedarea, 'matchdr', sub='ngals200_12', bin=7, $
        yrange=[0.8,1.02],ystyle=3,xrange=[0.01,40],xstyle=3, $
        /dops

    c2->plot_lumcolor_vs_rad, 'corrected', 'ngals200_12', 7, 0.25, /dops

    c2->plot_lumfunc, 'corrected', 'ngals200_12', 0,0.25, /linear,/dops
    c2->plot_lumfunc, 'corrected', 'ngals200_12', 0,0.25, /linear,/dops,/color

    yt = textoidl('#/Area [h^2 Mpc^{-2}]')
    c4->plot_profile_over, 'jackknife', 'ngals200_12', 'radnumdens',/reverse,$
        yt=yt, charsize=2.5, lcharsize=1.5, /dops

    c4->plot_profile_over, 'jackknife', 'ngals200_12', 'radnumdens',/reverse,$
        yt=yt, charsize=2.5, lcharsize=1.5, /dops, /color

    yt=textoidl('L_{0.25i}/Area [L_{\odot} Mpc^{-2}]')
    c4->plot_profile_over, 'jackknife', 'ngals200_12', 'radilumdens',/reverse,$
        yt=yt, charsize=2.5, lcharsize=1.5, /dops
    c4->plot_profile_over, 'jackknife', 'ngals200_12', 'radilumdens',/reverse,$
        yt=yt, charsize=2.5, lcharsize=1.5, /dops,/color

    m2l->plot_luminosities_over, 'ngals200_12', /both
    m2l->plot_luminosities_over, 'ngals200_12', /both,/color

    m2l->plot_masses, 'ngals200_12'
    m2l->plot_masses, 'ngals200_12',/color

    m2l->plot_m2lfits, 'ngals200_12'
    m2l->plot_m2lfits, 'ngals200_12',/color

    m2l->plot_m2l200_vs_ngals200, 'ngals200_12'
    m2l->plot_m2llast_vs_ngals200, 'ngals200_12'


    obj_destroy, c2, c4, m2l



end
