function estitle, tname, tex=tex

	case strlowcase(tname) of
		'kpcxtitle': title = 'R [h^{-1} kpc]'
		'kpcxtitle2': return, estitle('kpcxtitle')
		'mpcxtitle': title = 'R [h^{-1} Mpc]'
		'mpcxtitle2': return, estitle('mpcxtitle')
		'mpcxtitle32': title = 'r [h^{-1} Mpc]'
		'shytitle': title = '\gamma_T'
		'orthoytitle': title = 'gamma_{\times}'
		'deltaytitle': title = '\Delta\Sigma [h M_{\sun} pc^{-2}]'
		'rdeltaytitle': title = '\Delta\Sigma_{\times} [h M_{\sun} pc^{-2}]'
		'randdeltaytitle': title = '\Delta\Sigma_{rand} [h M_{\sun} pc^{-2}]'

		'lumytitle': title = 'L/Area [10^{10} L_{\sun} Mpc^{-2}]'

		'wgmytitle': title = 'w_{gm} [h^{-1} Mpc]'

		'xigmxtitle': title = 'r [h^{-1} Mpc]'
		'xigmytitle': title = '\xi_{gm}(r) \times \Omega_{m}/0.27'

		'rhoytitle': title = '\rho(r) - !S\rho!R!U-!N} [h^{2} 10^{12} M_{\sun} Mpc^{-3}]'
		'rhoytitle2': title = '\rho(r) - !S\rho!R!U-!N} [h^{2} 10^{12} M_{\sun} pc^{-3}]'

		'mytitle': title = 'M(<r) [10^{12} h^{-1} M_{\sun}]'
		'mytitle2': title = 'M(<r) [10^{14} h^{-1} M_{\sun}]'

		'm200ytitle': title = 'M_{200} [10^{14} h^{-1} M_{\sun}]'
		'r200ytitle': title = 'r_{200} [h^{-1} Mpc]'
		'mVirYtitle': title = 'M_{Vir} [10^{14} h^{-1} M_{\sun}]'
		'rVirYtitle': title = 'r_{Vir} [h^{-1} Mpc]'

		else:begin
			on_error, 2
			message,'Unknown title name: '+string(tname)
		end
	endcase

	if keyword_set(tex) then begin
		return, title
	endif else begin
		return, textoidl(title)
	endelse
end
