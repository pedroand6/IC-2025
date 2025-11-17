	SUBROUTINE SDSS_COLOR(T,X,Y,INW,LUN,BOLFLUX,SISSA20)

	! computes magnitude and colors in AB system

	! Array declarations
	parameter (nc=70,lb=2*nc,mxft=300)
	character genfile*96,envfile*96,rfcolorfile*96,fid(0:mxft)*64
	integer no(lb),n1(nc),n2(nc),kerr(nc)
	real x(inw),y(inw),col(nc),fx(lb),ab(0:nc)
	real gmag
	common /kpercent/ ipcall,iplast,iphead,jc
	common /r_filter/ iread,jread,ireset
	data icall/0/,jcall/0/,last/0/,nb/0/
	data nb,no,n1,n2,ng,nv,mc,kb,v0/0,lb*0,nc*0,nc*0,4*0,0./
	data kerr/nc*0/
	! Add V and g filters to list of desired filters
	logical sissa20
	integer glog
	data glog/121/			! SDSS g filter
	common /f1001/ mc,ab

	! Check for zero flux sed
	if (t.eq.0..or.bolflux.le.0.) return

	! Check if number of points has changed
	if (inw.ne.last) then
		! icall=0
		! write (6,*) 'Resetting sdss:',inw,last
		! ireset=1
		last=inw
	endif

	if (icall.eq.0) then
		icall=1

		! Format of file modified 07/01/2003:
		! To avoid recompiling this routine the arrays n1,n2 of up to nc elements each are now read from a file.
		! Get file name from environment variable RF_COLORS_ARRAYS
		envfile='RF_COLORS_ARRAYS'
		call getenv(envfile,rfcolorfile)
		close (1)
		open (1,file=rfcolorfile,status='old',form='formatted',err=2)
		! write (6,*)
		! write (6,*) 'List of filters in file: ',rfcolorfile(1:largo(rfcolorfile))
		! Skip first lines
		do i=1,100
		read (1,'(a)') genfile
		if (index(genfile,'AB system').gt.0) goto 5
		enddo
5		kb = 0
		jb = 0
		! Read filter pairs
		do i=1,100
		read (1,'(a)',end=10) genfile
		if (genfile(1:1) /= '#') then
			jb = jb+1
			read (genfile,*) n1(jb),n2(jb)
			kb=kb+1
			fx(kb)=n1(jb)
			kb=kb+1
			fx(kb)=n2(jb)
		endif
		enddo
10		mc=jb
		close (1)
		! Make sure that SDSS g filter = filter glog in filter file, are included
		kb=kb+1
		fx(kb)=glog
		! Sort filters in numerical order
		call sort(kb,fx)
		! Find independent filters in arrays n1 and n2 and store in array no
		nb=1
		no(1)=fx(1)
		do i=2,kb
		if (fx(i).gt.fx(i-1)) then
			nb=nb+1
			no(nb)=nint(fx(i))
		endif
		enddo
		! Write filter ID
		  call filterid(fid)
		! nun = 1002
		! nun = 6
		! write (nun,'(/x,a)') 'Selected filters (AB mag):'
		! do i=1,nb
		! write (nun,'(i4,'': '',i4,3x,a)') i,no(i),fid(no(i))
		! enddo
		! write (nun,'(i3,a)') mc,' colors selected (AB mag):'
		! do i=1,mc
		! write (nun,'(3i4,3x,3a)') i,n1(i),n2(i),trim(fid(n1(i))),' - ',trim(fid(n2(i)))
		! enddo

		! Fill arrays with filter numbers
		do i=1,mc
		do j=1,nb
		if (n1(i).eq.no(j)) n1(i) = j
		if (n2(i).eq.no(j)) n2(i) = j
		if (no(j).eq.glog)  ng    = j
		enddo
		enddo
	endif

	! Compute flux through each of nb filters
	do i=1,nb
	if (kerr(i).eq.0) then
		fx(i)=f_mean(no(i),x,y,inw,0.,kerr(i))
	else
		fx(i)=0.
	endif
	enddo

	! Compute colors in file RF_COLORS in the AB mag system
	do i=1,mc
	col(i)=-2.5*alog10(fx(n1(i))/fx(n2(i)))
	enddo

	! Compute absolute g AB magnitude for a 1 Mo galaxy
	! It is -27.5 magnitudes brighter for a 1E11 Mo galaxy
	! 10 pc in units of Mpc
	dl=1.e-5
	gmag=-2.5*alog10(fx(ng))-48.6
	gmag=gmag+(5. * alog10(1.7684e+08 * dl))	! this factor is SQRT(4*pi*(3.0856E24)^2/Lsun)

	! Compute flux inside GALEX filters
	nfuv = 139	! FUV Galex filter
	ffuv = f_mean(nfuv,x,y,inw,0.,kerruv)
	nnuv = 140	! NUV Galex filter
	fnuv = f_mean(nnuv,x,y,inw,0.,kerruv)

	! Compute flux inside a 100 A square filtered centered at 1500 A
	n1500 = 235     ! filter number
	f1500 = f_mean(n1500,x,y,inw,0.,kerruv)

	! Write magnitudes to file
	! Bands 01:18 go to file lun+13 = *.1ABmag
	! Bands 39:50, 30:38, 19:29 go to file lun+88 = *.2ABmag

	! Write results
	if (jcall == 0 .and. .not. sissa20) then
		! Write a record corresponding to t = 0 (use tl = 0 => 1 yr) for compatibility with number of time steps on SED
		! This record is a duplicate of next record
		jcall = 1
		write (lun+13,101) 0.,gmag-col(1),gmag,(gmag-col(i),i=2,18),ffuv,fnuv,f1500
		write (lun+88,102) 0.,(gmag-col(i),i=39,50),(gmag-col(i),i=30,38),(gmag-col(i),i=19,29)
	endif
	tl=alog10(t)
	write (lun+13,101) tl,gmag-col(1),gmag,(gmag-col(i),i=2,18),ffuv,fnuv,f1500
	write (lun+88,102) tl,(gmag-col(i),i=39,50),(gmag-col(i),i=30,38),(gmag-col(i),i=19,29)
101	format (f10.6,5f10.4,3x,10f10.4,3x,2f10.4,3x,2f10.4,1p3e13.4)
102	format (f10.6,40f10.4)

	! Store AB mags
	ab(0) = gmag-col(1)
	ab(1) = gmag
	do i=2,mc
	ab(i) = gmag-col(i)
	enddo

	return
2	write (6,*) 'Error opening file: ',rfcolorfile(1:largo(rfcolorfile))
	stop
	end

	! Prepare filters for next call to filter_n
	! ireset=1

	! Bands 19:29 go to file lun+86 = *.acs_wfc_ABmag
	! Bands 30:38 go to file lun+84 = *.wfc3_ABmag
	! Bands 39:50 go to file lun+88 = *.wfc3_uvis1_legus_ABmag
	! Bands 51:64 go to file lun+85 = *.wfc3_uvis1_ABmag

		! write (lun+86,102) 0.,(gmag-col(i),i=19,29)
		! write (lun+84,102) 0.,(gmag-col(i),i=30,38)
		! write (lun+88,102) 0.,(gmag-col(i),i=39,50)
		! write (lun+85,102) 0.,(gmag-col(i),i=51,64)

	! write (lun+86,102) tl,(gmag-col(i),i=19,29)
	! write (lun+84,102) tl,(gmag-col(i),i=30,38)
	! write (lun+88,102) tl,(gmag-col(i),i=39,50)
	! write (lun+85,102) tl,(gmag-col(i),i=51,64)
