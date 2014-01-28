PRO READ_SHELL, file, res
	dims = lonarr(4)
	CLOSE, 13
	OPENR, 13, file, /f77_unformatted
	READU,13, dims
	nphi = dims[0]
	ntheta = dims[1]
	nr = dims[2]
	nq = dims[3]
	radius = dblarr(nr)
	qvals = lonarr(nq)
	sintheta = dblarr(ntheta)
	costheta = dblarr(ntheta)
	vals = dblarr(nphi,ntheta,nr,nq)

	readu,13,qvals
	readu,13,radius
	readu,13,sintheta
	readu,13, vals
	readu,13,costheta
	res = {vals:vals, radius:radius, sintheta:sintheta,q:qvals, costheta:costheta}

	close, 13
END
