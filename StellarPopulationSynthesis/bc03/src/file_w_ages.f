	SUBROUTINE FILE_W_AGES(IO,NAME,Z,T,X,Y,N,UFWA,GFWA,RFWA,IFWA,ZFWA,KFWA,MWA,UFWLA,GFWLA,RFWLA,IFWLA,ZFWLA,KFWLA,MWLA)

c	Open file to write flux weighted age

c	Variables
	character*(*) name		!	,aux*10
	real x(n),y(n),ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,mwa,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,mwla
	data ifile/0/

c	Check kind of CSP
	if (io <= 0) then
		ifile=0
		close (2002)
		return
	endif

c	Check first time step
	if (t <= 0) then
		! ifile=0
		return
	endif

c	Build file name
	if (ifile == 0) then
c		Modify filename if z > 0.
		if (z > 0) then
			name = trim(name) // '.w_age_of'
		else
			name = trim(name) // '.w_age_rf'
		endif
		open (2002,file=name,form='formatted')
		call file_header(2002,name,0)
		write (2002,'(a)') '#            <---- z = 0 ---->  <---------------------------------------------------------------- computed after red shifting sed''s --------------------------------------------------------------------->'
		write (2002,'(a)') '# log t(yr)  B4000_N  Hdelta_A   MWA (yr)    uFWA(yr)    gFWA(yr)    rFWA(yr)    iFWA(yr)    zFWA(yr)    KFWA(yr)    MWLA(yr)  uFWLA(yr)  gFWLA(yr)  rFWLA(yr)  iFWLA(yr)  zFWLA(yr)  KFWLA(yr)      z'
		write (2002,'(a)') '#   (1)        (2)       (3)       (4)         (5)         (6)         (7)         (8)         (9)        (10)         (11)       (12)       (13)       (14)       (15)       (16)       (17)       (18)'
		ifile=1
		return
	endif

c	Compute B_4000 and Hdelta_A
        b4=b4000vn(x,y,n)
        hd=gw_ix_sed(22,x,y,n,0)

c	Report flux and mass weighted age for the composite population
	tl=alog10(t)
        write (2002,100) tl,b4,hd,mwa,ufwa,gfwa,rfwa,ifwa,zfwa,kfwa,mwla,ufwla,gfwla,rfwla,ifwla,zfwla,kfwla,z
100    	format (f10.6,2f10.5,1p7e12.4,0p8f11.5)
	return
	end
