PRO GEN_COLOC, nx, x
	x = dblarr(nx)
	FOR k = 0,  nx -1 DO BEGIN
		arg = !pi*(k+0.5d0)/double(nx)
		x[k] = cos(arg)
	ENDFOR
END

PRO TO_SPECTRAL, f, c
	n = N_ELEMENTS(f)
	ndouble = DOUBLE(n)
	c = DBLARR(n)
	FOR j = 0, n-1 DO BEGIN
		c[j] = 0.0d0
		FOR k = 0, n -1 DO BEGIN
			arg = !pi*j*(k+0.5d0)/ndouble
			c[j] = c[j]+f[k]*cos(arg)
		ENDFOR
	ENDFOR
	c = c*(2.0d0)/ndouble
END

PRO FROM_SPECTRAL, c, f
	n = N_ELEMENTS(c)
	f = DBLARR(n)
	ndouble = DOUBLE(n)
	FOR j = 0, n-1 DO BEGIN
		f[j] = -0.5d0*c[0]
		FOR k = 0, n -1 DO BEGIN
			arg = !pi*k*(j+0.5d0)/ndouble
			f[j] = f[j]+c[k]*cos(arg)
		ENDFOR
	ENDFOR

END

PRO CHEBY_DERIV, c, cprime, cprime2 = cprime2, cprime3 = cprime3
	n = N_ELEMENTS(c)
	cprime = DBLARR(n)
	cprime[n-1] = 0.0d0
	cprime[n-2] = 2.0d0*(n-1)*c[n-1]
	FOR i = n-3, 0,-1 DO BEGIN
		cprime[i] = cprime[i+2]+2.0d0*(i+1)*c[i+1]
	ENDFOR



END 

PRO GEN_DERIV_WEIGHTS, nmax, alpha
	alpha = LONARR(nmax,nmax)
	alpha[*,*] = 0L
	alpha[nmax-1,*] = 0L
	alpha[nmax-2,nmax-1] = 2L*(nmax-1L)
	FOR n = nmax-3,0L,-1 DO BEGIN
		alpha[n,n+1] = 2L*(n+1L)
		alpha[n,*] = alpha[n,*]+alpha[n+2,*] 
	ENDFOR
END


PRO cheby_routines
	print, 'Void function'
END
