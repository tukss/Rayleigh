PRO READ_AZAVG, file, res
	dims = lonarr(3)
	CLOSE, 13
	OPENR, 13, file, /f77_unformatted
	READU,13, dims
	nr = dims[0]
	nt = dims[1]
	nq = dims[2]
	radius = dblarr(nr)
	sintheta = dblarr(nt)
	qvals = lonarr(nq)
	vals = dblarr(nr,nt,nq)

	readu,13,qvals
	readu,13,radius
	readu,13, sintheta
	readu,13, vals
	res = {vals:vals, radius:radius, sintheta:sintheta, q:qvals}
	close, 13
END
