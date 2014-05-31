PRO READ_SHELLAVG, file, res, SWAP = swap
	time = 0.0d0
	dims = lonarr(2)
	CLOSE, 13
	IF (KEYWORD_SET(SWAP)) THEN BEGIN
		OPENR, 13, file, /f77_unformatted, /swap_if_little_endian
	ENDIF ELSE BEGIN
		OPENR, 13, file, /f77_unformatted
	ENDELSE

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
