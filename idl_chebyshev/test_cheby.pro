@cheby_routines.pro
PRO test_cheby, nx = nx
	if (not keyword_set(nx)) then nx = 32
	gen_coloc,nx,x

	print, x
	f = x
	to_spectral, f, c
	print, c

	f = 2.0d0*x^2-1.0d0
	to_spectral, f, c
	print, c

	f = 4.0d0*x^3-3.0d0*x
	to_spectral, f, c
	print, c

	!p.multi = [0,1,3]
	f = exp(-4.0d*x^2)
	f = f+ exp(-4.0d*(x-0.25)^2)
	f = f- exp(-4.0d*(x+0.25)^2)
	
	to_spectral,f,c
	print, 'f    : ', f
	print, 'Coefs: ', c
	from_spectral,c, f2
	plot, x, f, psym = 2, title = 'The Function'
	oplot, x, f2, color = 250

	;cheby_deriv,c,cprime
	gen_deriv_weights,nx,alpha
	cprime = alpha#c
	from_spectral, cprime, fprime
	ifprime = deriv(x,f)
	plot, x, ifprime, title = 'The Derivative', psym = 2
	oplot, x, fprime, color = 250



END
