PRO READ_ARRAY, file, arr, COMPLEX = complex
	OPENU, 10, file, /f77_unformatted
	ndim = 0L
	READU, 10, ndim
	dim = lonarr(ndim)
	FOR i = 0, ndim -1 DO BEGIN
		tmp = 0L
		READU, 10, tmp
		dim[i] = tmp
	ENDFOR
	IF (KEYWORD_SET(COMPLEX)) THEN BEGIN
		arr = DCOMPLEXARR(dim)

	ENDIF ELSE BEGIN
		arr = DBLARR(dim)

	ENDELSE
	READU,10, arr
	CLOSE, 10
END
