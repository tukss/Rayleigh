pro read_reference, file, res
	openr, 13, file, /f77_unformatted
	nr = 0L
	readu,13,nr

	radius = dblarr(nr)
	density = dblarr(nr)
	dlnrho = dblarr(nr)
	d2lnrho = dblarr(nr)
	pressure = dblarr(nr)
	temperature = dblarr(nr)
	dlnt = dblarr(nr)
	dsdr = dblarr(nr)
	entropy = dblarr(nr)
	gravity = dblarr(nr)
	readu,13, radius
	readu,13, density
	readu,13, dlnrho
	readu,13, d2lnrho
	readu,13, pressure
	readu,13, temperature
	readu,13, dlnt
	readu,13, dsdr
	readu,13, entropy
	readu,13, gravity

	close, 13
	res = {radius:radius, density:density,pressure:pressure,temperature:temperature,dlnro:dlnrho,  $
		gravity:gravity, dlnT:dlnT, d2lnrho:d2lnrho, dsdr:dsdr,entropy:entropy}
end
