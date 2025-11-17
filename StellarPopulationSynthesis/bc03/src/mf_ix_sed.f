	REAL FUNCTION MF_IX_SED(J,X,F,N)

	! Computes value of Fanelli et al. UV index J directly from sed
	! Follows procedure outlined by Trager at el. (1998, ApJS 116, 1)
	! See their eq. 1-3, pag. 5.

	! J = Index identification
	! X = wavelength scale
	! F = flux scale (sed)
	! N = dimension of X, F

	! Array declaration
	character*40 idxlbl
	real x(n),f(n)
	real w(6),u(10000),z(10000)

	! Equation of straight line joining (wb,fb) and (wr,fr)
	fc(v) = ( (wr-v)*fb + (v-wb)*fr ) / (wr-wb)

	! Find wavelength points everytime to allow calls to this routine in
	! the same program using different wavelength arrays.
	! Get Parameters
	call mf_dat(j,w,idxlbl)

	! Compute mean height inside blue pseudo continuum band
       	fb=trapz3(x,f,n,w(1),w(2),ierr)/(w(2)-w(1))
	wb=(w(1)+w(2))/2.
	if (ierr.ne.0) then
		mf_ix_sed = -99.
		return
	endif

	! Compute mean height inside red pseudo continuum band
       	fr=trapz3(x,f,n,w(5),w(6),ierr)/(w(6)-w(5))
	wr=(w(5)+w(6))/2.
	if (ierr.ne.0) then
		mf_ix_sed = -99.
		return
	endif

	! Compute equivalent width according to eq. (2) of
	! Trager at el. (1998, ApJS 116, 1)
	! Locate central band in array x: w(3) and w(4)
       	call locate (x,n,w(3),i1)
       	call locate (x,n,w(4),i2)
	! Fill in auxiliary array from i1-5 to i2+5 with ratio f/fc
	nu=0
	do i=i1-5,i2+5
	nu=nu+1
	u(nu)=x(i)
	if (fc(x(i)) <= 0.) then
		mf_ix_sed = -98.
		return
	endif
	z(nu)=f(i)/fc(x(i))
	enddo
	! Compute integral of ratio f(i)/fc(x(i)) from w(3) to w(4)
       	fm=trapz3(u,z,nu,w(3),w(4),ierr)
	if (ierr.ne.0) then
		mf_ix_sed = -99.
		return
	endif

	! Compute index
	if (fm > 1000.) then
		mf_ix_sed = -97.
	elseif (fm < 0.) then
		mf_ix_sed = -96.
	elseif (fm == 0.) then
		mf_ix_sed = -95.
	else
		mf_ix_sed = w(4) - w(3) - fm
	endif
	return
	end
