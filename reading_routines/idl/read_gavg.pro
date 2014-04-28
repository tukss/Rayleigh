pro read_gavg, file, res
	time = 0.0d0
	openr, 13, file, /f77_unformatted
	nq = 0L
	readu,13,nq
	qvals = lonarr(nq)
	readu,13,qvals
	vals = dblarr(nq)
	readu,13, vals
	readu,13, time
	close, 13
	res = {nq:nq, vals:vals, qvals:qvals, time:time}
end
