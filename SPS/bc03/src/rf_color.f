	SUBROUTINE RF_COLOR(IO,T,X,Y,INW,LUN,BOLFLUX,STRMASS,SF,EVFLUX,SNBR,PNBR,BH,SN,WD,RM,XM,BOLMS,GASMASS,GALMASS,TMLR,TDPR,SISSA20)

	! Array declarations
	parameter (nc=100,lb=2*nc,jid=35,its=400)
	logical sissa20
	character genfile*96,envfile*96,rfcolorfile*96,h(0:200)*15	! ,h(0:200)*24
	! Negative indices allow for mandatory colors
	integer p,n1(-4:nc),n2(-4:nc),no(lb),ly(6),kerr(nc)
	real zp(-4:nc),col(-4:nc),phot(0:200),f(200),w(200),s(200),snbr(0:4)
	real umag,bmag,rmag,jmag,kmag,bsun,vsun,ksun,blr1,vlr1,klr1,blr,vlr,klr,lx
	real x(inw),y(inw),fx(lb),balm(3),gwx(jid),gwl(jid),ab(0:70)
	real m_pn,m_hm,m_bh,m_ns,m_wd
	save no,ly,zp,sl,mb,nv,mc,kb,nb,v0,p,w,s
	include 'stelib.dec'
	include 'filter.dec'
	common /kpercent/ ipcall,iplast,iphead,jc
	! Common to store fluxes used to compute indices
	integer iim(jid)
	real ffb(jid),ffr(jid),ffc(jid),ffl(jid),mf_ix(21),mf_ix_sed
	common /fluxes/ iim,ffb,ffr,ffc,ffl
	common /massrem/ m_pn,m_hm,m_bh,m_ns,m_wd
	real*8 mabo(its),mbel(its),mcut,rmup(its),rmlo(its)
	common /topimf/ mabo,mbel,rmup,rmlo,mcut,ic
	data icall/0/,jcall/0/,ly/6*0/,last/0/
	! Define mandatory colors
	data n1/ 15,15,15,12,14, nc*0/
	data n2/125,57,84,13,15, nc*0/
        data j1,j2/1,21/
	data kerr/nc*0/,jfits/0/
	common /f1001/ mab,ab

	! Return if t = 0 (log t undefined)
	if (t <= 0.) return
	tl = alog10(t)

	! Check for zero flux sed
	if (bolflux <= 0.) then
		bolflux=trapz1(x,y,inw)
		if (bolflux <= 0.) return
	endif
	izero=0
	do i=1,inw
	if (y(i) > 0) then
		izero=izero+1
	endif
	enddo
	if (izero == 0) then
		! write (6,*) 'Exiting rf_color at t =',t,' because of zero flux'
		! return
		! To allow using .?color files with a single number of lines
		do i=1,inw
		y(i) = 1.E-33
		enddo
	endif

	! Compute galaxy mass and mass in gas (except for add_bursts io = 5)
	! if (io.ne.10) then
	! 	galmass=gal_mass(io,t,sf)
	! endif

	! Report what we are doing
	if (iphead == 0) then
		write (6,*)
		write (6,'(x,a/ )') 'Computing magnitudes, colours and indices'
		iphead=1
	endif

	! Compute Worthey and other indices directly from sed and in the Lick system
	ism=0
	do j=1,jid
	gwx(j)=gw_ix_sed(j,x,y,inw,0)
	gwl(j)=gw_ix_sed_lick_system(j,x,y,inw,0,ism)
	enddo

	! Compute Fanelli et al. indices
	do j=j1,j2
	! Indices 1 to 11 are below wavelength coverage of HNGSL
	mf_ix(j)=mf_ix_sed(j,x,y,inw)
	enddo

	! Check if number of points has changed
	if (inw <= last) then
	      ! write (6,*) 'Resetting:',inw,last
	      !	ireset=1
		last=inw
	endif

	if (icall.eq.0) then
		icall=1

		! Format of file modified 07/01/2003:
		! To avoid recompiling this routine the arrays n1,n2 of
		! up to nc elements each are now read from a file.
		! Get file name from environment variable RF_COLORS_ARRAYS
		envfile='RF_COLORS_ARRAYS'
		call getenv(envfile,rfcolorfile)
		close (1)
		open (1,file=rfcolorfile,status='old',form='formatted',err=2)
		! write (6,*) 'List of filters in file: ',rfcolorfile(1:largo(rfcolorfile))
		! Read one line of header
		read (1,'(a)') genfile
		! write (6,*) genfile(1:largo(genfile))
		kb=0
		! Read filter pairs
		do i=1,nc
		read (1,'(2i4)',err=10) n1(i),n2(i)
		kb=kb+1
		fx(kb)=n1(i)
		kb=kb+1
		fx(kb)=n2(i)
		enddo
10		mc=i-1
		close (1)
		! Add mandatory filters
		do i=-4,0
		kb=kb+1
		fx(kb)=n1(i)
		kb=kb+1
		fx(kb)=n2(i)
		enddo
		! Sort filters in numerical order
		call sort(kb,fx)
		! Find independent filters in arrays n1 and n2, and store in no
		nb=1
		no(1)=fx(1)
		do i=2,kb
		if (fx(i).gt.fx(i-1)) then
			nb=nb+1
			no(nb)=nint(fx(i))
		endif
		enddo
		! Get filter ID
		call filterid(fid)
		! nun=1002
		! nun=6
		! write (nun,'(x,a)') 'Selected filters (Vega mag):'
		! do i=1,nb
		! write (nun,'(i4,'': '',i4,3x,a)') i,no(i),fid(no(i))
		! enddo
		! write (nun,'(i3,a)') mc,' colors selected (Vega mag):'
		! do i=1,mc
		! write (nun,'(3i4,3x,3a)') i,n1(i),n2(i),trim(fid(n1(i))),' - ',trim(fid(n2(i)))
		! enddo

		! Log of solar luminosity
		sl=33.+alog10(3.826)

		! Read filter file and compute zero points
		! write (6,*)
		! write (6,*) 'Computing Zero Points:'
		! Compute V magnitude Vega zero point
		v0=vega_0p_n(15)
		do i=-4,mc
		zp(i)=zerop_n(n1(i),n2(i))

		! Fill arrays with filter numbers
		do j=1,nb
		if (n1(i).eq.no(j)) n1(i)=j
		if (n2(i).eq.no(j)) n2(i)=j
		if (no(j).eq.14) mb=j
		if (no(j).eq.15) nv=j
		enddo
		enddo

		! Find in array x the position of points used to define the continuum at Lyman alpha
		do i=1,inw
		if (x(i).le.1120.) ly(1)=i
		if (x(i).le.1140.) ly(2)=i
		if (x(i).le.1160.) ly(3)=i
		if (x(i).le.1280.) ly(4)=i
		if (x(i).le.1300.) ly(5)=i
		if (x(i).ge.1320.) then
			ly(6)=i
			goto 1
		endif

		enddo
	endif
1	continue

	! Compute AB magnitudes for SDSS and other bands
	! ireset=1
	call sdss_color(t,x,y,inw,lun,bolflux,sissa20)

	! Compute flux through each of nb filters
	do i=1,nb
	kerr(i)=0	! This is useful just in case the first sed in a model has zero fluxes for a given filter
	! 		! otherwise the flux through this filter is not computed at later ages.
	if (kerr(i).eq.0) then
		fx(i)=filter_n(no(i),x,y,inw,0.,kerr(i))
	else
		fx(i)=0.
	endif
	enddo

	! Compute colors in Vega system
	do i=-4,mc
	if (fx(n1(i)).gt.0..and.fx(n2(i)).gt.0.) then
		col(i)=zp(i)-2.5*alog10(fx(n1(i))/fx(n2(i)))
	! 	write (6,*) i,n1(i),kerr(n1(i)),n2(i),kerr(n2(i))
	else
		col(i)=-99.99
	! 	write (6,*) i,n1(i),kerr(n1(i)),fx(n1(i)),n2(i),kerr(n2(i)),fx(n2(i))
	endif
	enddo

	! Compute bolometric magnitude
	bolmag=4.75-2.5*alog10(bolflux)
	bolpms=bolflux-bolms
	if (bolms > 0.) then
		bolrat=bolpms/bolms
	else
		bolrat=0.
	endif

	! Compute V magnitude for a 1 Mo galaxy
	! It is -27.5 magnitudes brighter for a 1E11 Mo galaxy
	vmag=v0-2.5*alog10(fx(nv))

	! Compute U, B, R, K, and J2MASS magnitudes
	bmag=vmag+col( 0)
	umag=bmag+col(-1)
	rmag=vmag-col(-2)
	kmag=vmag-col(-3)
	jmag=vmag-col(-4)

	! Compute mass-to-visual-light ratio in solar units SUPERSEDED (Feb. 27, 2004).
	! Using a G2 V sed and the filters number 14 and 15
	! in the filter file, one derives:
	! 	fblue(sun) = 0.138Lo
	! 	fvis (sun) = 0.113Lo
	! This numbers apply only to the filters in this filter library.
	! Express blue and visual flux of model galaxy (also measured in
	! Lo) in units of the blue and visual flux of the sun:
	! fblu=fx(mb)/0.138
	! fvis=fx(nv)/0.113
	! fblu=filter(14,x,y,inw,0.)/0.138
	! fvis=filter(15,x,y,inw,0.)/0.113
	! Total mass in galaxy  = 1 Mo
	! Compute mass-to-visual-light ratio
	! blr=1./fblu
	! vlr=1./fvis
	! Compute stellar-mass-to--light ratios
	! blr=strmass/fblu
	! vlr=strmass/fvis

	! The solar absolute magnitudes for U,B,V,R,I,J,H,K were calibrated against the values
	! of Binney and Merrifield 1998, Galactic Astronomy, Table 2.1 (page 53), assuming
	! Bessell filters, and the offsets used to calibrate the entire set of filters.
	! Some values need checking - particularly those using UV filters (FOCA and Galex).
	! Taken from: http://mips.as.arizona.edu/~cnaw/sun.html

	! 	Filter	 B&M	 here	 difference
	! 	U	 5.61	 5.55	 0.06
	! 	B	 5.48	 5.45	 0.03
	! 	V	 4.83	 4.80	 0.03
	! 	R	 4.42	 4.46	 -0.04
	! 	I	 4.08	 4.11	 -0.03
	! 	J	 3.64	 3.67	 -0.02
	! 	H	 3.32	 3.33	 0.01
	! 	K	 3.28	 3.29	 0.01


	! Compute mass-to-visual-light ratio in solar units. Improved definition (Feb. 27, 2004).
	! Use solar absolute V and B magnitudes and (B-V)sun = 0.65
	vsun=4.80
	bsun=5.45
	ksun=3.29
	! M/L using the total mass in stars (M*), i.e. no remmnants
	blr1=strmass*10.**(0.4*(bmag-bsun))
	vlr1=strmass*10.**(0.4*(vmag-vsun))
	klr1=strmass*10.**(0.4*(kmag-ksun))
	! Modified Jan. 2011 to add mass of remmnants (rm = mBH + mNS + mWD)
	blr=(strmass+rm)*10.**(0.4*(bmag-bsun))
	vlr=(strmass+rm)*10.**(0.4*(vmag-vsun))
	klr=(strmass+rm)*10.**(0.4*(kmag-ksun))

	! Number of Lyman Continuum photons = Cly (log)
	! Flux in Lyman alpha from recombination theory
	! E(Lalpha) = 4.78E-13 * 33.1 * Nuv
	! log E = log(Nuv) -10.8 = cly -10.8
	! Number of Lyman continuum photons
	phly=clyman(x,y,inw)
	if (phly > 0.) then
		cly=sl+alog10(phly)
		fa=cly-10.8
	else
		cly=0.
		fa=0.
	endif

	! Number of Helium ionizing photons
	phe=chelium(x,y,inw,phe2)
	if (phe > 0.) then
		che=sl+alog10(phe)
	else
		che=0.
	endif
	if (phe2 > 0.) then
		che2=sl+alog10(phe2)
	else
		che2=0.
	endif
	! Ratio de He II to H I ionizing photons
	cher=che2-cly

	! Stellar continuum at Lyman alpha
	if (x(1).ge.1320.) then
		scly=0.
	else
		scly=(y(ly(1))+y(ly(2))+y(ly(3))+y(ly(4))+y(ly(5))+y(ly(6)))/6.
	endif
	if (scly.gt.0.) then
		fc=sl+alog10(scly)
	else
		fc=0.
	endif

	! Ly alpha equivalent width assuming that the continuum is the stellar continuum
	if (fa.gt.0..and.fc.gt.0.) then
		ew=10.**(fa-fc)
		ew2=phly/scly/10.**(10.8)
	else
		ew=0.
		ew2=0.
	endif

	! X ray luminosity (> 0 in models with X ray binaries)
	xlx = lx(x,y,inw)

	! Compute Mg2 index
	ymg2=ymag2(x,y,inw)

	! Compute 912 A break
	b9=b912(x,y,inw)

	! Compute 4000 A break
	b4=b4000(x,y,inw)

	! Compute narrow version of D4000
	b4_n=b4000vn(x,y,inw)

	! Compute SDSS version of D4000
	b4_s=b4000_sdss(x,y,inw)

	! Compute equivalent width of Balmer lines (Hgamma, Hdelta, Hbeta)
	ewbl=ew_balmer(x,y,inw,balm)

	! Compute bolometric magnitude
	bolmag=4.75-2.5*alog10(bolflux)

	! Compute specific flux, snbr, pnbr
	evf=evflux/bolflux
	! write (6,*) t,evflux,evf,xm

	! Compute SNR and PNR per unit Luminosity (superseded, Sept. 2021)
	! snr=snbr/bolflux
	! pnr=pnbr/bolflux

	! SFR/year
	! sf=sfr(t)

	! Compute quantities requested by C. Popescu
	! call popescu(tl,sf,x,y,inw)

	! Compute flux from 1500 to 2800A
	! it=0
	! do i=1,inw
	! if (x(i).ge.1500.0.and.x(i).le.2800.0) then
	! 	it=it+1
	! 	xx(it)=x(i)
	! 	yy(it)=y(i)
	! endif
	! enddo
	! fuv=trapz1(x,y,it)
	! write (9,*) tl, sf,fuv

	! Write results
	if (jcall == 0 .and. .not. sissa20) then
		! Write a record corresponding to t = 0 (use tl = 0 => 1 yr) for compatibility with number of time steps on SED
		! This record is a duplicate of next record
		jcall = 1
		write (lun+1 ,108) 0.,bolmag,(vmag-col(i),i=1,3),vmag,(vmag-col(i),i=4,14)
		write (lun+2 ,108) 0.,bolmag,vmag,(vmag-col(i),i=15,28)
		!rite (lun+3 ,102) 0.,b4,b4_n,b4_s,b9,cly,che,che2,cher,bolmag,bolflux,snr,snt,bh,sn,pnr,wd,rm,xlx
		!rite (lun+3 ,102) 0.,b4,b4_n,b4_s,b9,cly,che,che2,cher,bolmag,bolflux,snbr,bh,sn,pnbr,wd,rm,xlx
		!rite (lun+3 ,102) 0.,b4,b4_n,b4_s,b9,cly,che,che2,cher,snbr,bh,sn,pnbr,wd,rm,xlx
		write (lun+3 ,102) 0.,b4,b4_n,b4_s,b9,cly,che,che2,cher,snbr,pnbr,bh,sn,wd,xlx
		write (lun+4 ,103) 0.,bolmag,bmag,vmag,kmag,strmass,rm,gasmass,galmass,sf,strmass+rm,blr,vlr,klr,blr1,vlr1,klr1
		write (lun+5 ,104) 0.,bolmag,evflux,evf,xm,bolrat,tmlr,tdpr
		write (lun+8 ,105) 0.,(gwx(j),j= 1,21)
		write (lun+9 ,106) 0.,(gwx(j),j=22,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
		write (lun+10,110)(0.,j,iim(j),ffb(j),ffr(j),ffc(j),ffl(j),gwx(j),j=1,jid)			! Store fluxes used to compute spectral indices
		write (lun+11,105) 0.,(gwl(j),j= 1,21)
		write (lun+12,106) 0.,(gwl(j),j=22,25),(gwl(j),j=30,31),(gwl(j),j=26,29),(gwl(j),j=32,jid)
		write (lun+87,107) 0.,(mf_ix(i),i=j1,j2)
	endif
	write (lun+1 ,108) tl,bolmag,(vmag-col(i),i=1,3),vmag,(vmag-col(i),i=4,14)
	write (lun+2 ,108) tl,bolmag,vmag,(vmag-col(i),i=15,28)
	!rite (lun+3 ,102) tl,b4,b4_n,b4_s,b9,cly,che,che2,cher,bolmag,bolflux,snr,snt,bh,sn,pnr,wd,rm,xlx
	!rite (lun+3 ,102) tl,b4,b4_n,b4_s,b9,cly,che,che2,cher,bolmag,bolflux,snbr,bh,sn,pnbr,wd,rm,xlx
	!rite (lun+3 ,102) tl,b4,b4_n,b4_s,b9,cly,che,che2,cher,snbr,bh,sn,pnbr,wd,rm,xlx
	write (lun+3 ,102) tl,b4,b4_n,b4_s,b9,cly,che,che2,cher,snbr,pnbr,bh,sn,wd,xlx
	write (lun+4 ,103) tl,bolmag,bmag,vmag,kmag,strmass,rm,gasmass,galmass,sf,strmass+rm,blr,vlr,klr,blr1,vlr1,klr1
	write (lun+5 ,104) tl,bolmag,evflux,evf,xm,bolrat,tmlr,tdpr
	write (lun+8 ,105) tl,(gwx(j),j= 1,21)
	write (lun+9 ,106) tl,(gwx(j),j=22,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
	write (lun+10,110)(tl,j,iim(j),ffb(j),ffr(j),ffc(j),ffl(j),gwx(j),j=1,jid)				! Store fluxes used to compute spectral indices
	write (lun+11,105) tl,(gwl(j),j= 1,21)
	write (lun+12,106) tl,(gwl(j),j=22,25),(gwl(j),j=30,31),(gwl(j),j=26,29),(gwl(j),j=32,jid)
	write (lun+87,107) tl,(mf_ix(i),i=j1,j2)
!02	format (f10.6,3f10.4,1pe12.3,1x,0p4f10.4,2x,f10.4,1p12e12.4)
102	format (f10.6,3f9.4,1pe11.3,0p4f9.4,1p11e12.4)
103	format (f10.6,4f10.4,1p6e13.4,1x,6e12.4)
104	format (f10.6,f10.4,1p6e13.4)
105     format (f10.6,21f9.4)
106     format (f10.6,1x,6f9.4,3f11.4,4f10.4,f10.4)
107	format (f10.6,21f11.4)
108	format (f10.6,18f10.4)
110	format (f10.6,i4,i4,1p4e12.3,0pf12.4)

	! return					 ! comment this statement to get file fort.1000 written

	! Write temporary file fort.1000 with photometric magnitudes
	! Organize photometric magnitudes in order of wavelength to be included in fits file

	! Bands listed in file *.1VEGAmag: Mbol, U, B2, B3, V, Rc, Ic, R, I, J, K, L, PalJ, PalH, PalK, K'
	! Johnson UBV filters
	p =   0 ; phot(p) = bolmag         ;                h(p) = '           Mbol'  !  Mbol
	p = p+1 ; phot(p) = vmag - col( 1) ; f(p) =  12   ; h(p) = '      U_Johnson'  !  U
	p = p+1 ; phot(p) = vmag - col( 2) ; f(p) =  13   ; h(p) = '     B2_Johnson'  !  B2
	p = p+1 ; phot(p) = vmag - col( 3) ; f(p) =  14   ; h(p) = '     B3_Johnson'  !  B3
	p = p+1 ; phot(p) = vmag           ; f(p) =  15   ; h(p) = '      V_Johnson'  !  V
	! Cousins RI filters
	p = p+1 ; phot(p) = vmag - col( 4) ; f(p) =  84   ; h(p) = '      R_Cousins'  !  Rc
	p = p+1 ; phot(p) = vmag - col( 5) ; f(p) =  85   ; h(p) = '      I_Cousins'  !  Ic
	! Johnson RIJKL filters
	p = p+1 ; phot(p) = vmag - col( 6) ; f(p) =  32   ; h(p) = '      R_Johnson'  !  R
	p = p+1 ; phot(p) = vmag - col( 7) ; f(p) =  33   ; h(p) = '      I_Johnson'  !  I
	p = p+1 ; phot(p) = vmag - col( 8) ; f(p) =  34   ; h(p) = '      J_Johnson'  !  J
	p = p+1 ; phot(p) = vmag - col( 9) ; f(p) =  35   ; h(p) = '      K_Johnson'  !  K
	p = p+1 ; phot(p) = vmag - col(10) ; f(p) =  36   ; h(p) = '      L_Johnson'  !  L
	! Palomar JHK filters
	p = p+1 ; phot(p) = vmag - col(11) ; f(p) =  55   ; h(p) = '      J_Palomar'  !  PalJ
	p = p+1 ; phot(p) = vmag - col(12) ; f(p) =  56   ; h(p) = '      H_Palomar'  !  PalH
	p = p+1 ; phot(p) = vmag - col(13) ; f(p) =  57   ; h(p) = '      K_Palomar'  !  PalK
	! Cowie K' filter
	p = p+1 ; phot(p) = vmag - col(14) ; f(p) =  88   ; h(p) = '   Kprime_Cowie'  !  K'

	! Bands listed in file *.2VEGAmag: Mbol, Vmag, 2MASSJ, 2MASSH, 2MASSKs, IRAC3.5, IRAC4.5, IRAC5.7, IRAC7.9, IRAS12, IRAS25, IRAS60, IRAS100, MIPS24, MIPS70, MIPS160
	! 2Mass JHKs filters
	p = p+1 ; phot(p) = vmag - col(15) ; f(p) = 125   ; h(p) = '        J_2Mass'  !  J2MASS
	p = p+1 ; phot(p) = vmag - col(16) ; f(p) = 126   ; h(p) = '        H_2Mass'  !  H2MASS
	p = p+1 ; phot(p) = vmag - col(17) ; f(p) = 127   ; h(p) = '       Ks_2Mass'  ! Ks2MASS
	! IRAC filters
	p = p+1 ; phot(p) = vmag - col(18) ; f(p) = 128   ; h(p) = '      I3p6_IRAC'  ! IRAC3.5
	p = p+1 ; phot(p) = vmag - col(19) ; f(p) = 129   ; h(p) = '      I4p5_IRAC'  ! IRAC4.5
	p = p+1 ; phot(p) = vmag - col(20) ; f(p) = 130   ; h(p) = '      I5p7_IRAC'  ! IRAC5.7
	p = p+1 ; phot(p) = vmag - col(21) ; f(p) = 131   ; h(p) = '      I7p9_IRAC'  ! IRAC7.9
	! IRAS filters
	p = p+1 ; phot(p) = vmag - col(22) ; f(p) =  71   ; h(p) = '       I12_IRAS'  ! IRAS12 
	p = p+1 ; phot(p) = vmag - col(23) ; f(p) =  72   ; h(p) = '       I25_IRAS'  ! IRAS25 
	p = p+1 ; phot(p) = vmag - col(24) ; f(p) =  73   ; h(p) = '       I60_IRAS'  ! IRAS60 
	p = p+1 ; phot(p) = vmag - col(25) ; f(p) =  74   ; h(p) = '      I100_IRAS'  ! IRAS100
	! MIPS filters
	p = p+1 ; phot(p) = vmag + col(26) ; f(p) = 132   ; h(p) = '       M24_MIPS'  ! MIPS24 
	p = p+1 ; phot(p) = vmag + col(27) ; f(p) = 133   ; h(p) = '       M70_MIPS'  ! MIPS70 
	p = p+1 ; phot(p) = vmag - col(28) ; f(p) = 134   ; h(p) = '      M160_MIPS'  ! MIPS160

	! Bands listed in file *.1ABmag: u, g, r, i, z, u1, u3, g1, g3, r1, r3, i2, i3, z1, z3, H, Ks, FUV, NUV
	! SDSS filters (AB mag)
	p = p+1 ; phot(p) = ab( 0)         ; f(p) = 120   ; h(p) = '      u_SDSS_AB'  ! SDSS u AB mag
	p = p+1 ; phot(p) = ab( 1)         ; f(p) = 121   ; h(p) = '      g_SDSS_AB'  ! SDSS g AB mag
	p = p+1 ; phot(p) = ab( 2)         ; f(p) = 122   ; h(p) = '      r_SDSS_AB'  ! SDSS r AB mag
	p = p+1 ; phot(p) = ab( 3)         ; f(p) = 123   ; h(p) = '      i_SDSS_AB'  ! SDSS i AB mag
	p = p+1 ; phot(p) = ab( 4)         ; f(p) = 124   ; h(p) = '      z_SDSS_AB'  ! SDSS z AB mag
	! CFHT MegaCam filters (AB mag )
	p = p+1 ; phot(p) = ab( 5)         ; f(p) = 237   ; h(p) = '  u1_CFHT_MC_AB'  ! CFHT MC u1 AB mag
	p = p+1 ; phot(p) = ab( 6)         ; f(p) = 266   ; h(p) = '  u3_CFHT_MC_AB'  ! CFHT MC u3 AB mag
	p = p+1 ; phot(p) = ab( 7)         ; f(p) = 238   ; h(p) = '  g1_CFHT_MC_AB'  ! CFHT MC g1 AB mag
	p = p+1 ; phot(p) = ab( 8)         ; f(p) = 267   ; h(p) = '  g3_CFHT_MC_AB'  ! CFHT MC g3 AB mag
	p = p+1 ; phot(p) = ab( 9)         ; f(p) = 239   ; h(p) = '  r1_CFHT_MC_AB'  ! CFHT MC r1 AB mag
	p = p+1 ; phot(p) = ab(10)         ; f(p) = 268   ; h(p) = '  r3_CFHT_MC_AB'  ! CFHT MC r3 AB mag
	p = p+1 ; phot(p) = ab(11)         ; f(p) = 255   ; h(p) = '  i2_CFHT_MC_AB'  ! CFHT MC i2 AB mag
	p = p+1 ; phot(p) = ab(12)         ; f(p) = 269   ; h(p) = '  i3_CFHT_MC_AB'  ! CFHT MC i3 AB mag
	p = p+1 ; phot(p) = ab(13)         ; f(p) = 241   ; h(p) = '  z1_CFHT_MC_AB'  ! CFHT MC z1 AB mag
	p = p+1 ; phot(p) = ab(14)         ; f(p) = 270   ; h(p) = '  z3_CFHT_MC_AB'  ! CFHT MC z3 AB mag
	p = p+1 ; phot(p) = ab(15)         ; f(p) = 126   ; h(p) = '     H_2MASS_AB'  ! H_2MASS AB mag
	p = p+1 ; phot(p) = ab(16)         ; f(p) = 256   ; h(p) = '  Ks_CFHT_WC_AB'  ! CFHT WIRCam Ks AB mag
	! GALEX filters (AB mag)
	p = p+1 ; phot(p) = ab(17)         ; f(p) = 139   ; h(p) = '   FUV_GALEX_AB'  ! GALEX 1500 AB mag
	p = p+1 ; phot(p) = ab(18)         ; f(p) = 140   ; h(p) = '   NUV_GALEX_AB'  ! GALEX 2300 AB mag

	! Bands listed in file *.2ABmag: F225w, F275w, F336w, F438w, F547m, F555w, F606w, F625w, F656n, F657n, F658n, F814w,   
	!                                F110W, F125W, F160W, F225W, F336W, FR388N, F438W, F555W, F814W,
	!                                F220w, F250w, F330w, F410w, F435w, F475w, F555w, F606w, F625w, F775w, F814w
	! HST WFC3 UVIS1 LEGUS filters (AB mag)
	p = p+1 ; phot(p) = ab(39)	   ; f(p) = 242   ; h(p) = ' UVIS1_f225w_AB'  ! UVIS1 LEGUS f225w AB mag
	p = p+1 ; phot(p) = ab(40)	   ; f(p) = 243   ; h(p) = ' UVIS1_f275w_AB'  ! UVIS1 LEGUS f275w AB mag
	p = p+1 ; phot(p) = ab(41)	   ; f(p) = 244   ; h(p) = ' UVIS1_f336w_AB'  ! UVIS1 LEGUS f336w AB mag
	p = p+1 ; phot(p) = ab(42)	   ; f(p) = 245   ; h(p) = ' UVIS1_f438w_AB'  ! UVIS1 LEGUS f438w AB mag
	p = p+1 ; phot(p) = ab(43)	   ; f(p) = 246   ; h(p) = ' UVIS1_f547m_AB'  ! UVIS1 LEGUS f547m AB mag
	p = p+1 ; phot(p) = ab(44)	   ; f(p) = 247   ; h(p) = ' UVIS1_f555w_AB'  ! UVIS1 LEGUS f555w AB mag
	p = p+1 ; phot(p) = ab(45)	   ; f(p) = 248   ; h(p) = ' UVIS1_f606w_AB'  ! UVIS1 LEGUS f606w AB mag
	p = p+1 ; phot(p) = ab(46)	   ; f(p) = 249   ; h(p) = ' UVIS1_f625w_AB'  ! UVIS1 LEGUS f625w AB mag
	p = p+1 ; phot(p) = ab(47)	   ; f(p) = 250   ; h(p) = ' UVIS1_f656n_AB'  ! UVIS1 LEGUS f656n AB mag
	p = p+1 ; phot(p) = ab(48)	   ; f(p) = 251   ; h(p) = ' UVIS1_f657n_AB'  ! UVIS1 LEGUS f657n AB mag
	p = p+1 ; phot(p) = ab(49)	   ; f(p) = 252   ; h(p) = ' UVIS1_f658n_AB'  ! UVIS1 LEGUS f658n AB mag
	p = p+1 ; phot(p) = ab(50)	   ; f(p) = 253   ; h(p) = ' UVIS1_f814w_AB'  ! UVIS1 LEGUS f814w AB mag

	! HST WFC3 filters (AB mag)
	p = p+1 ; phot(p) = ab(30)	   ; f(p) = 197   ; h(p) = '  WFC3_F110W_AB'  ! WFC3 F110W AB mag
	p = p+1 ; phot(p) = ab(31)	   ; f(p) = 198   ; h(p) = '  WFC3_F125W_AB'  ! WFC3 F125W AB mag
	p = p+1 ; phot(p) = ab(32)	   ; f(p) = 199   ; h(p) = '  WFC3_F160W_AB'  ! WFC3 F160W AB mag
	p = p+1 ; phot(p) = ab(33)	   ; f(p) = 200   ; h(p) = '  WFC3_F225W_AB'  ! WFC3 F225W AB mag
	p = p+1 ; phot(p) = ab(34)	   ; f(p) = 201   ; h(p) = '  WFC3_F336W_AB'  ! WFC3 F336W AB mag
	p = p+1 ; phot(p) = ab(35)	   ; f(p) = 202   ; h(p) = ' WFC3_FR388N_AB'  ! WFC3 FR388N AB mag
	p = p+1 ; phot(p) = ab(36)	   ; f(p) = 203   ; h(p) = '  WFC3_F438W_AB'  ! WFC3 F438W AB mag
	p = p+1 ; phot(p) = ab(37)	   ; f(p) = 204   ; h(p) = '  WFC3_F555W_AB'  ! WFC3 F555W AB mag
	p = p+1 ; phot(p) = ab(38)	   ; f(p) = 205   ; h(p) = '  WFC3_F814W_AB'  ! WFC3 F814W AB mag

	! HST ACS wide filters (AB mag)
	p = p+1 ; phot(p) = ab(19)	   ; f(p) = 219   ; h(p) = 'ACSWFC_F220w_AB'  ! ACS WFC F220w AB mag
	p = p+1 ; phot(p) = ab(20)	   ; f(p) = 220   ; h(p) = 'ACSWFC_F250w_AB'  ! ACS WFC F250w AB mag
	p = p+1 ; phot(p) = ab(21)	   ; f(p) = 221   ; h(p) = 'ACSWFC_F330w_AB'  ! ACS WFC F330w AB mag
	p = p+1 ; phot(p) = ab(22)	   ; f(p) = 222   ; h(p) = 'ACSWFC_F410w_AB'  ! ACS WFC F410w AB mag
	p = p+1 ; phot(p) = ab(23)	   ; f(p) = 223   ; h(p) = 'ACSWFC_F435w_AB'  ! ACS WFC F435w AB mag
	p = p+1 ; phot(p) = ab(24)	   ; f(p) = 224   ; h(p) = 'ACSWFC_F475w_AB'  ! ACS WFC F475w AB mag
	p = p+1 ; phot(p) = ab(25)	   ; f(p) = 225   ; h(p) = 'ACSWFC_F555w_AB'  ! ACS WFC F555w AB mag
	p = p+1 ; phot(p) = ab(26)	   ; f(p) = 226   ; h(p) = 'ACSWFC_F606w_AB'  ! ACS WFC F606w AB mag
	p = p+1 ; phot(p) = ab(27)	   ; f(p) = 227   ; h(p) = 'ACSWFC_F625w_AB'  ! ACS WFC F625w AB mag
	p = p+1 ; phot(p) = ab(28)	   ; f(p) = 228   ; h(p) = 'ACSWFC_F775w_AB'  ! ACS WFC F775w AB mag
	p = p+1 ; phot(p) = ab(29)	   ; f(p) = 229   ; h(p) = 'ACSWFC_F814w_AB'  ! ACS WFC F814w AB mag

	! Write temporary file fort.1000 with photometric magnitudes to be inserted in fits file
	if (jfits == 0) then
		jfits = 1
		! Compute effective wavelength for each of the p filters
		do i=1,p
		s(i) = i
		area = areafilt(int(f(i)),w(i))
		enddo
		! Sort in increasing order of effective wavelength
		call sort2(p,w,s)
		write (1000, '(a,200(a16))') '# logage_yr',h(0),(h(int(s(i))),i=1,p)
		! Write a record corresponding to t = 0 (use tl = 0 => 1 yr) for compatibility with number of time steps on SED
		! This record is a duplicate of next record
		if (.not. sissa20) then
			write (1000,109) 0.,phot(0),(phot(int(s(i))),i=1,p)
		else
			if (t == 1.5E4) then
				write (1000,109) 0.,phot(0),(phot(int(s(i))),i=1,p)
			endif
		endif
		! do i=1,p
		! write (1001,112) 0.,phot(int(s(i))),w(int(s(i)))
		! enddo
	endif
	write (1000,109) tl,phot(0),(phot(int(s(i))),i=1,p)
	! do i=1,p
	! write (1001,112) tl,phot(int(s(i))),w(int(s(i)))
	! enddo
	return
2	write (6,*) 'Error opening file: ',rfcolorfile(1:largo(rfcolorfile))
109	format (f11.6,90f16.4)
!112	format (f10.6,f12.4,f12.1)
	stop
	end
