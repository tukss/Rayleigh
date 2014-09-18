PRO READ_SPECTRA, file, res
	endian_tag = 0L
	version = 0L
	nrec = 0L
	lmax = 0L
	nr = 0L
	nq = 0L
	CLOSE, 13
	OPENR, 13, file
	READU,13, endian_tag
	If (endian_tag ne 314) THEN BEGIN
		endian1 = endian_tag
		CLOSE,13
		OPENR,13,file, /swap_if_little_endian
		READU,13,endian_tag
		IF (endian_tag ne 314) THEN BEGIN
			print, 'Unable to discern endianess of file!'
			print, "Expected integer value 314 in first 4 bytes"
			print, "Found : ", endian1
			print, "And   : ", endian_tag
			STOP
		ENDIF

	Endif
	; Read the data dimensions
	READU,13, version
	READU,13, nrec
	READU,13, lmax
	READU,13, nr
	READU,13, nq


	; Read the header arrays
	qvals      = LONARR(nq)
	radius     = DBLARR(nr)
	shell_inds = LONARR(nr)

	READU,13,qvals
	READU,13,radius
	READU,13,shell_inds


	; Read the individual records
	vals = DCOMPLEXARR(lmax+1,lmax+1,nr,nq,nrec)
	tmp  = DBLARR(lmax+1,lmax+1)
	time = DBLARR(nrec)
	iter = LONARR(nrec)
	it = 0L
	tm = 0.0d0
	FOR i = 0, nrec-1 DO BEGIN
		FOR im = 0, 1 DO BEGIN
		FOR k = 0, nq -1 DO BEGIN
			FOR j = 0, nr -1 DO BEGIN
				READU,13,tmp
				print,min(tmp),max(tmp)
				if (im eq 0) THEN BEGIN
					vals[*,*,j,k,i] = dcomplex(tmp,tmp*0)					
				endif else begin
					vals[*,*,j,k,i] =vals[*,*,j,k,i]+ dcomplex(0*tmp,tmp)
				endelse
			ENDFOR
		ENDFOR
		ENDFOR
		READU,13,tm
		READU,13,it
		time[i] = tm
		iter[i] = it
	ENDFOR
	CLOSE, 13

	; Build a lookup table for the quantity codes
        qmax = 400L
        lut = LONARR(qmax+1)
        lut[*] = qmax*2
        FOR i =	0, nq -1 DO BEGIN
                lut[qvals[i]] =	 i
        ENDFOR


	tmp = abs(vals)^2
	pow = total(tmp,2)	
	res = {vals:vals, radius:radius, qvals:qvals, shell_inds:shell_inds, time:time, $
			iter:iter, version:version, lut:lut,pow:pow}

END
