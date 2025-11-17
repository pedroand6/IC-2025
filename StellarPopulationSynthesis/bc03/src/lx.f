	REAL*4 FUNCTION LX(X,Y,N)

	! Compute X ray luminosity in ergs/sec
	! Assumes X is in Angstroms and Y in Lo/Angstroms
	! Lx: between 0.5 and 8 keV (i.e. between 1.5498 and 24.7969 Å)

	! Variables
	real*4 x(n),y(n)
	real*8 lsun
	data x1/1.5498/, x2/24.7969/, lsun/3.846D33/

	! Check if SED extends to X ray region SED extends to X ray region (0.1 A)
	if (x(1) > 0.1) then
		lx = 0.
	else
		lx = trapz2(x,y,n,x1,x2,ierr)
	endif
	return
	end
