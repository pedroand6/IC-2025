	PROGRAM SEEFILTER

c	Writes an ASCII file with the response function of the
c	selected filters.

	include 'filter.dec'
	real x(5000),y(5000),z(5000)

	call FILTER0

	l1=0
	l2=0
1	write (6,100) 'Extract filter number = '
	read (5,101,err=1,end=10) n
	if (n <= 0) stop
100	format (1x,a,$)
101	format (i10)
	nt=np(n)
	l1=l1+3
	l2=l2+nt+2
	write (28,'(''# '',a,i4,2a)') 'Filter',n,': ',fid(n)(1:largo(fid(n)))
	write (28,'(''# '',i4,a,i4,a,i4)') nt,' data points. Lines',l1,' to',l2
	nt=0
	do i=ni(n),nl(n)
	nt=nt+1
	x(nt)=r(i,1)
	y(nt)=r(i,2)
	z(nt)=y(nt)*x(nt)
	write (28,102) x(nt),y(nt)
102	format (1p2e14.5)
	enddo
	area=trapz1(x,y,nt)
	weff=trapz1(x,z,nt)/area

	! Compute zeropoint
	Rvega = 1.6432E6 	! Radiis of Vega in km
	Dvega = 7.76 * 3.086E13	! Distance to Vega in km = 7.76 pc
	Dfact = (Rvega/Dvega)**2
	pi    = 4.*atan(1.)
	f1    = 3.631E-9/6.4829E-10
	f2    = 3.631E-9/6.4061E-10
	write (6,*) dfact,pi
	if (jread == 0) call readA0V
	do i=1,jread
	if (xa0v(i) > weff) then
		write (6,*) i,xa0v(i),ya0v(i),ya0v(i)*Dfact*pi*f1,ya0v(i)*Dfact*pi*f2,f1,f2
		goto 2
	endif
	enddo
2	f0 = filter(n,xa0v,ya0v,jread,0.,kerr) * dfact * pi / area * f1
	g0 = f_meanl(n,xa0v,ya0v,jread,0.,kerr) * dfact * pi * f2
	write (6,'(x,a,i4,2a)')      'Filter',n,': ',fid(n)
	write (6,'(x,a,1p2e12.4)')   'Range covered        = ',x(1),x(nt)
	write (6,'(x,a,1pe12.4)')    'Area below filter    = ',area
	write (6,'(x,a,1pe12.4)')    'Effective wavelength = ',weff
	write (6,'(x,a,1pe12.4,a )') 'Zero point-old       = ',f0,' ergs/s/cm**2/A'
	write (6,'(x,a,1pe12.4,a/)') 'Zero point-new       = ',g0,' ergs/s/cm**2/A'
	write (6,*) '------------------------------------------------------'
	write (9,'(x,a,i4,2a)')      'Filter',n,': ',fid(n)
	write (9,'(x,a,1p2e12.4)')   'Range covered        = ',x(1),x(nt)
	write (9,'(x,a,1pe12.4)')    'Area below filter    = ',area
	write (9,'(x,a,1pe12.4)')    'Effective wavelength = ',weff
	write (9,'(x,a,1pe12.4,a )') 'Zero point-old       = ',f0,' ergs/s/cm**2/A'
	write (9,'(x,a,1pe12.4,a/)') 'Zero point-new       = ',g0,' ergs/s/cm**2/A'
	write (9,*) '------------------------------------------------------'
	l1=l2
	goto 1
10	write (6,*) 'Output in files fort.28 and fort.9'
	end

