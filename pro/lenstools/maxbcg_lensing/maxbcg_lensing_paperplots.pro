pro maxbcg_lensing_paperplots

    ; make all but the scat weight plot

    common maxbcg_lensing_paperplots_common, lsamp, bcg

    m21=obj_new('maxbcg_lensing', 21)
    m22=obj_new('maxbcg_lensing', 22)
    zw2=obj_new('zweight', 2)
    p=obj_new('percolate',2)

    m21->plot_aitoff, l=lsamp, sample=21, /dops 
    m21->plot_aitoff, l=lsamp, sample=21, /dops, /color

    zw2->plothist_zbins, [2,6,14,20,24]
    zw2->plothist_zbins, [2,6,14,20,24],/color

    m21->plothist_ngals200, bcg=bcg, /dops

    m21->plot_zdensity, 16, bcg=bcg
    m21->plot_zdensity, 16, bcg=bcg, /color

    m22->compare_real_rand_two, sub='ngals200_12', last=7, /dops

    m21->plot_clustcorr_over, 'ngals200_12',/dops
    m21->plot_clustcorr_over, 'ngals200_12',/dops, /color
    
    m21->plot_profile, 'jackknife', sub='ngals200_12', sample=[21,22], /dops

    m21->powerlaw_fits, sample=[21,22]
    m21->powerlaw_fits, sample=[21,22],/color

    m21->plot_ngals_lumsplit, sample=[21,22], /dops
    m21->plot_ngals_lumsplit, sample=[21,22],/color, /dops

    m21->plot_profile, 'jackknife', sub='ilum200_16', sample=[21,22], /dops

    p->plot_intrinsic,/dops

    m21->compare_intrinsic, 'ngals200_12', /cumulative, /dops
    m21->compare_intrinsic, 'ngals200_12', /cumulative, /dops, /color

    obj_destroy, m21, m22, zw2, p
end
