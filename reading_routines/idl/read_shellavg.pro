PRO READ_SHELLAVG, file, res
	time = 0.0d0
	dims = lonarr(2)
	CLOSE, 13
	OPENR, 13, file, /f77_unformatted
	READU,13, dims
	nr = dims[0]
	nq = dims[1]
	radius = dblarr(nr)
	qvals = lonarr(nq)
	vals = dblarr(nr,nq)

	readu,13,qvals
	readu,13,radius
	readu,13, vals
	readu,13, time
	res = {vals:vals, radius:radius, q:qvals, time:time}
	close, 13
END
