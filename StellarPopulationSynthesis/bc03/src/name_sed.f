	SUBROUTINE NAME_SED(LUN,JUN,KDIST,IHRD,IOPTION,JDEF,LDEF,ML,MU,EVOFILE,ATLAS)

!	Selects options and builds output file names

!	Array declaration
	parameter (nb=20,nph=20)
	logical done
	real ml,mu
	include 'stelib.dec'
        common /vazdekis/ jvaz,xmu
	character atlas,evofile*(*),phase(nph)*7,fract(nb)*4,stage(nb)*5,savefile*512,imf(-1:24)*7,tempfile*512,subdir*512,auxname*16

!	data phase/'ms','sgb','rgb','hb','agb','dcspn','cspn', 'wd',     'ccb', 'tpagb','coolg','coold'/
!	data phase/'ms','sgb','rgb','hb','agb','Orich','Crich','SWphase','cspn','wd',   'ccb',  'other'/
!	Modified to include WR stars (Dec 2013) Corrected Apr. 2017.
	data phase/'ms','sgb','rgb','hb','agb','Orich','Crich','Orich2','Crich2','SWphase','cspn','wd','ccb','other','WNL','WNE','WC','WO','WR','hotsed'/

!	data fract/'H'  , 'Ks','22','27','U','B','V','R','I','J','K','L','BOL'/
!	data fract/'I25','I60','22','27','U','B','V','R','I','J','K','L','BOL'/
!	data fract/'14' , '17','22','27','U','B','V','R','I','J','K','L','BOL'/
!	Modified to include IRAC and MIPS filters (July 2013)
	data fract/'14' , '17','22','27','U','B','V','R','I','J','K','L','I3p5','I4p5','I5p7','I7p9','M24','M70','M160','BOL'/

!	data stage/'nH'  ,'nKs' ,'n22','n27','nU','nB','nV','nR','nI','nJ','nK','nL','nBOL'/
!	data stage/'nI25','nI60','n22','n27','nU','nB','nV','nR','nI','nJ','nK','nL','nBOL'/
!	data stage/'n14' ,'n17' ,'n22','n27','nU','nB','nV','nR','nI','nJ','nK','nL','nBOL'/
!	Modified to include IRAC and MIPS filters (July 2013)
	data stage/'n14' ,'n17' ,'n22','n27','nU','nB','nV','nR','nI','nJ','nK','nL','nI3p5','nI4p5','nI5p7','nI7p9','nM24','nM70','nM160','nBOL'/

!	Various IMF's options
	data imf /'chab','chab','kroup','KROUP','salp','scalo','SCALO','mi_sc','ktg93','kenn', 'leque','ferra','piott','carig','thIMF','myIMF','vxpxx','Schab','Skroup','SKROUP','Ssalp',5*'NumIMF'/

!	Use to file name of models with IMF populated stochastically
	include 'stochnam.dec'

!	Check for Vazdekis IMF
	if (jvaz > 0) then
		write (imf(15)(2:),'(f4.2)') xmu
		imf(15)(3:3) = 'p'
	endif

!	Check for subdirectory at begining of name
	if (lun.eq.-99)	then
		do i=largo(evofile),1,-1
		if (evofile(i:i).eq.'/') then
			subdir=evofile(1:i)
			goto 2
		endif
		enddo
	endif

!	Ask for file name (if not entered in command line)
2	if (lun.lt.0) then
		evofile=' '
		call getenv('PREFIX',tempfile)
		if (largo(tempfile).eq.0) tempfile='op'
		if (stelib) then
			if (ngsl > 0) then
				tempfile=tempfile(1:largo(tempfile)) // '_lrs_hngsl'
			elseif (hires == 3 .and. qran) then
				tempfile=tempfile(1:largo(tempfile)) // '_hr_rmiles'	! Miles only (or pure Miles) including random errors in stellar parameters
			elseif (hires == 3) then
				tempfile=tempfile(1:largo(tempfile)) // '_hr_pmiles'	! Miles only (or pure Miles)
			elseif (hires == 4) then
				tempfile=tempfile(1:largo(tempfile)) // '_hr_smiles'	! Synthetic Miles
			elseif (hires == 5) then
				tempfile=tempfile(1:largo(tempfile)) // '_hrs_indous'	! IndoUS only (or pure IndoUS)
			elseif (hires == 6) then
				tempfile=tempfile(1:largo(tempfile)) // '_hrx_mastar'	! Miles + Tlusty UV
			! elseif (hires == 6) then
			! 	tempfile=tempfile(1:largo(tempfile)) // '_hrt_miles'	! Miles + Tlusty UV
			! elseif (hires == 7) then
			! 	tempfile=tempfile(1:largo(tempfile)) // '_hrl_miles'	! Miles + Tlusty + Leitherer et al. UV
			! elseif (hires == 8) then
			! 	tempfile=tempfile(1:largo(tempfile)) // '_hrw_miles'	! Miles + Tlusty + Smith et al. WR models
			elseif (hires == 9) then
				tempfile=tempfile(1:largo(tempfile)) // '_hrx_miles'	! Extended Miles=Miles+Tlusty+Leitherer+WR
			else
				tempfile=tempfile(1:largo(tempfile)) // '_hrs_stelib'
			endif
		else
			if (abs(ioption).ge.15.and.abs(ioption).le.19) then
				if (ioption.gt.0) then
					tempfile=tempfile(1:largo(tempfile)) // '_lr_bpgs'
				else
					tempfile=tempfile(1:largo(tempfile)) // '_lr_jacoby'
				endif
				if (abs(ioption).eq.15) then
					ioption = -62
				elseif (abs(ioption).eq.16) then
					ioption = -162
				elseif (abs(ioption).eq.17) then
					ioption = -96
					tempfile(1:6)='cb2007'
				elseif (abs(ioption).eq.18) then
					ioption = -107
					tempfile(1:6)='cb2009'
				elseif (abs(ioption).eq.19) then
					ioption = -117
					tempfile(1:6)='cb2009'
				endif
			elseif (abs(ioption).ge.5.and.abs(ioption).le.10) then
				if (ioption.gt.0) then
					tempfile=tempfile(1:largo(tempfile)) // '_lr_pickles_fir'
				else
					if (atlas.eq.'f') then
						tempfile=tempfile(1:largo(tempfile)) // '_lr_fanelli_pickles_fir'
					elseif (atlas.eq.'F') then
						tempfile=tempfile(1:largo(tempfile)) // '_lr_fanelli_pickles_irtf_fir'
					else
						write (6,*) 'Unknown option for Pickles extended atlas: ',atlas
						stop
					endif
				endif
				if (abs(ioption).eq.5) then
					ioption = -62
				elseif (abs(ioption).eq.6) then
					ioption = -162
				elseif (abs(ioption).eq.7) then
					ioption = -96
					tempfile(1:6)='cb2007'
!				elseif (abs(ioption).eq.8) then
!					ioption = -107
!					tempfile(1:6)='cb2009'
!				elseif (abs(ioption).eq.9) then
!					ioption = -117
!					tempfile(1:6)='cb2009'
!				elseif (abs(ioption).eq.10) then
!					ioption = -217
!					tempfile(1:6)='cb2010'
!				elseif (abs(ioption).eq.11) then
!					ioption = -317
!					tempfile(1:6)='cb2013'
				elseif (abs(ioption).eq.8) then
					ioption = -410
					tempfile(1:6)='cb2018'
				elseif (abs(ioption).eq.9) then
					ioption = -411
					tempfile(1:6)='cb2018'
				elseif (abs(ioption).eq.10) then
					ioption = -412
					tempfile(1:6)='cb2018'
				endif
			else
				if (galaxev) tempfile=tempfile(1:largo(tempfile)) // '_lr_BaSeL'
			endif
		endif

     	  if (ioption.gt.99) then
	  	tempfile=tempfile(1:largo(tempfile)) // '_p000_'
     	  	ii=index(tempfile,'_p000_')+2
     	  	kk=index(tempfile,    '0_')
		write (tempfile(ii:ii+2),'(i3)') ioption
     	  elseif (ioption.ge.91.and.ioption.le.96) then
!		Solution to handle Garching and BaSTI tracks
		if (.not.basti) then
	  		tempfile=tempfile(1:largo(tempfile)) // '_M0_'
     	  		ii=index(tempfile,'_M0_')+2
		else
	  		tempfile=tempfile(1:largo(tempfile)) // '_B0_'
     	  		ii=index(tempfile,'_B0_')+2
		endif
		kk=index(tempfile,    '0_')
		write (tempfile(ii:ii),'(i1)') abs(ioption)-90
		tempfile(1:6) = 'cb2007'
     	  elseif (ioption.gt.10) then
	  	tempfile=tempfile(1:largo(tempfile)) // '_p00_'
     	  	ii=index(tempfile,'_p00_')+2
     	  	kk=index(tempfile,   '0_')
		write (tempfile(ii:ii+1),'(i2)') ioption
     	  elseif (ioption.gt.0) then
	  	tempfile=tempfile(1:largo(tempfile)) // '_p0_'
     	  	ii=index(tempfile,'_p0_')+2
     	  	kk=index(tempfile,  '0_')
		write (tempfile(ii:ii),'(i1)') ioption
     	  elseif (ioption.gt.-10) then
	  	tempfile=tempfile(1:largo(tempfile)) // '_m0_'
     	  	ii=index(tempfile,'_m0_')+2
     	  	kk=index(tempfile,  '0_')
		write (tempfile(ii:ii),'(i1)') abs(ioption)
     	  elseif (ioption.gt.-99) then
	  	tempfile=tempfile(1:largo(tempfile)) // '_m00_'
     	  	ii=index(tempfile,'_m00_')+2
     	  	kk=index(tempfile,   '0_')
		write (tempfile(ii:ii+1),'(i2)') abs(ioption)
     	  elseif (ioption.gt.-800) then
	  	tempfile=tempfile(1:largo(tempfile)) // '_m000_'
     	  	ii=index(tempfile,'_m000_')+2
     	  	kk=index(tempfile,    '0_')
		write (tempfile(ii:ii+2),'(i3)') abs(ioption)
     	  endif

!	  Check for mass loss option
	  if     (massloss.eq.-2) then
	  		diewssvl=' '
			massloss=0				! Distorted Schultheis models
	  elseif (massloss.eq.-1) then
			tempfile=tempfile(1:kk) // 't_'		! Distorted IRTF + Aringer et al. C-star models
			diewssvl=tempfile(ii-1:kk+1)
			massloss=0
	  elseif (massloss.eq.10) then
			tempfile=tempfile(1:kk) // 'f_'		! IRTF + Aringer et al. C-star model, dustfree sed''s used as input for Dusty code (added in CB10 and CB13 models)
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.1) then
			tempfile=tempfile(1:kk) // 'n_'		!    1 x instantaneous normal mass loss rate from tracks
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.2) then
			tempfile=tempfile(1:kk) // 'p_'		!   10 x instantaneous normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.3) then
			tempfile=tempfile(1:kk) // 'q_'		!    5 x instantaneous normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.4) then
			tempfile=tempfile(1:kk) // 'r_'		!    2 x instantaneous normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.5) then
			tempfile=tempfile(1:kk) // 'm_'		! 1/10 x instantaneous normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.6) then
			tempfile=tempfile(1:kk) // 'i_'		! 2/10 x instantaneous normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.7) then
			tempfile=tempfile(1:kk) // 'j_'		! 5/10 x instantaneous normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+1)
	  elseif (massloss.eq.11) then
			tempfile=tempfile(1:kk) // 'na_'	!    1 x average normal mass loss rate from tracks
			diewssvl=tempfile(ii-1:kk+2)
	  elseif (massloss.eq.12) then
			tempfile=tempfile(1:kk) // 'pa_'	!   10 x average normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+2)
	  elseif (massloss.eq.13) then
			tempfile=tempfile(1:kk) // 'qa_'	!    5 x average normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+2)
	  elseif (massloss.eq.14) then
			tempfile=tempfile(1:kk) // 'ra_'	!    2 x average normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+2)
	  elseif (massloss.eq.15) then
			tempfile=tempfile(1:kk) // 'ma_'	! 1/10 x average normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+2)
	  elseif (massloss.eq.16) then
			tempfile=tempfile(1:kk) // 'ia_'	! 2/10 x average normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+2)
	  elseif (massloss.eq.17) then
			tempfile=tempfile(1:kk) // 'ja_'	! 5/10 x average normal low mass loss rate extrapolated from refereed paper figure
			diewssvl=tempfile(ii-1:kk+2)
	  elseif (massloss.eq.21) then
			tempfile=tempfile(1:kk) // 'nh_'	!    1 x hybrid normal mass loss rate from tracks
			diewssvl=tempfile(ii-1:kk+2)
!	Superseded
!	  elseif (massloss.eq.-2) then
!			tempfile=tempfile(1:kk) // 'p2_'  !   10 x normal low mass loss rate extrapolated from non-refereed paper Fig 2
!			diewssvl=tempfile(ii-1:kk+1)
!	  elseif (massloss.eq.-3) then
!			tempfile=tempfile(1:kk) // 'q2_'  !    5 x normal low mass loss rate extrapolated from non-refereed paper Fig 2
!			diewssvl=tempfile(ii-1:kk+1)
!	  elseif (massloss.eq.-4) then
!			tempfile=tempfile(1:kk) // 'r2_'  !    2 x normal low mass loss rate extrapolated from non-refereed paper Fig 2
!			diewssvl=tempfile(ii-1:kk+1)
!	  elseif (massloss.eq.-5) then
!			tempfile=tempfile(1:kk) // 'm2_'  ! 1/10 x normal low mass loss rate extrapolated from non-refereed paper Fig 2
!			diewssvl=tempfile(ii-1:kk+1)
!	  elseif (massloss.eq.-6) then
!			tempfile=tempfile(1:kk) // 'i2_'  ! 2/10 x normal low mass loss rate extrapolated from non-refereed paper Fig 2
!			diewssvl=tempfile(ii-1:kk+1)
!	  elseif (massloss.eq.-7) then
!			tempfile=tempfile(1:kk) // 'j2_'  ! 5/10 x normal low mass loss rate extrapolated from non-refereed paper Fig 2
!			diewssvl=tempfile(ii-1:kk+1)
!
!	  elseif (massloss.eq.-8) then
!			tempfile=tempfile(1:kk) // 's_'  ! No mass loss, Schultheis sed''s at all phases positions 261:275
!			diewssvl=tempfile(ii-1:kk+1)
!	  elseif (massloss.eq.8) then
!			tempfile=tempfile(1:kk) // 'b_'  ! black body seds (only March 07 option)
!			diewssvl=tempfile(ii-1:kk+1)
!	  elseif (massloss.eq.-1) then
!			tempfile=tempfile(1:kk) // 'k_'  ! use Basel (Kurucz) for TP-AGB stars (only March 07 option)
!	  		diewssvl=' '
!			massloss=0
	  endif

	  if (ioption.eq.99) then
	     ioption=0
	  else
	     tempfile=tempfile(1:largo(tempfile)) // imf(jdef)(1:largo(imf(jdef))) // popprm(1:largo(popprm)) //'_ssp'
!	     Add chi2 subfix to stelib model if kfit=2
	     if (kfit.eq.2) then
		tempfile=tempfile(1:largo(tempfile)) // '_chi2'
	     endif
	  endif

!	  Check if Padova 2007+ tracks are selected ( = Padova 1994/2000/2007/2008 + Marigo 2007/2008 TP-AGB)
	  call check2007(tempfile,cb07)

!	  Write file name in 2016 style
	  call name2016(tempfile)

!	  Store tempfile in variable popprm in case option to write "cmd" file (Theo B.) has been requested
!	  popprm = tempfile
	  if (largo(popprm) > 0) popprm = tempfile

!	  Check if IMF mass limits have been changed from default
	  if (ldef.ne.0) then
		auxname='_ML0000_MU00000'
     	  	ii=index(auxname,'_ML0000_')+3
!		Write value of ML
		if (ml.lt.10) then
			write (auxname(ii:ii+3),'(f4.2)') ml
		elseif (ml.lt.100) then
			write (auxname(ii:ii+3),'(f4.1)') ml
		else
			write (auxname(ii:ii+3),'(f4.0)') ml
		endif
!		Write value of MU
		if (mu.ge.100) then
			write (auxname(ii+7:ii+11),'(f5.1)') mu
		elseif (mu.ge.10) then
			write (auxname(ii+8:ii+11),'(f4.1)') mu
		elseif (mu.ge.1) then
			write (auxname(ii+9:ii+11),'(f3.1)') mu
		endif
!		write (auxname(ii+7:ii+9),'(i3)') int(mu)
11     	  	ii=index(auxname(1:largo(auxname)),' ')
		if (ii.gt.0) then
			auxname(ii:ii)='0'
			goto 11
		endif
12		ii=index(auxname(1:largo(auxname)),'.')
		if (ii.gt.0) then
			auxname(ii:ii)='p'
			goto 12
		endif
		if (auxname /= '_ML0p10_MU100p0') then
			if (auxname(:7) == '_ML0p10') then
				auxname = auxname(8:)
				ii = index(auxname,'p0')
				if (ii > 0) auxname = auxname(:ii-1)
			endif
!			tempfile=tempfile(1:largo(tempfile)) // auxname    ! old style = cb2018_z004_chab_hr_xmilesi_ssp_ML0p10_MU600p0.ised
			ii = index(tempfile,'_')
			ii = index(tempfile(ii+1:),'_') + ii
			ii = index(tempfile(ii+1:),'_') + ii
			tempfile = tempfile(:ii-1) // trim(auxname) // tempfile(ii:)
		endif
		if (jvaz == 0) then
			ii = index(tempfile,'p00')
			if (ii > 0) tempfile = tempfile(:ii-1) // tempfile(ii+3:)
			ii = index(tempfile,'p0')
			if (ii > 0) tempfile = tempfile(:ii-1) // tempfile(ii+2:)
		endif
	  endif
	  if (lun.eq.-2) then
!		Option number entered as first argument of command line
		write (6,'(/1x,a,i4)') 'Using Option: ',ioption
		evofile=tempfile
	  elseif (lun.ne.-99) then
	      if (ioption.ne.0) then
		  if (index(tempfile,'NumIMF') > 0) then
!			Correct file name for Stochastic IMF read from file:
			ii = index(tempfile,'NumIMF') + 5
			jj = index(tempfile,'yr') + 1
			ll = index(tempfile,'single') + 5
			if (ll > 5) then
				tempfile = tempfile(:ii) // tempfile(jj+1:ll) // tempfile(ii+1:jj) // tempfile(ll+1:)
			endif
			! if (jdef == 24) then
			! 	tempfile = tempfile(:ii-6) // trim(bpimf) // tempfile(ii+1:)
			! 	call no_bl(tempfile)
			! elseif (jdef >= 20) then
			if (jdef >= 20) then
				tempfile(ii-largo(bpimf)+1:ii) = bpimf
				!ii=index(tempfile,'BF') + 10
				ii=index(tempfile,'_lr_') - 1
				if (ii <= 0) ii=index(tempfile,'_hr_') - 1
				if (ii <= 0) ii=index(tempfile,'_hrx_') - 1
				!all system ('/bin/mv -f fort.668      ' // tempfile(:ii) // '.binned_imf')
				call system ('/bin/rm -f fort.668')
				call system ('/bin/mv -f fort.776      ' // tempfile(:ii) // '.stat')
				inquire (file='fort.777',exist=done)
				if (done) then
					if (tempfile(ii-3:ii) == 'BF_0') then
						call system ('/bin/rm -f fort.777')
					else
						call system ('/bin/mv -f fort.777      ' // tempfile(:ii) // '.bp')
					endif
					call system ('/bin/mv -f fort.778      ' // tempfile(:ii) // '.allstrs')
					call system ('/bin/mv -f fort.779      ' // tempfile(:ii) // '.allsngl')
					call system ('/bin/cp -f random.seed12 ' // tempfile(:ii) // '.random.seed12')
				endif
			endif
		  endif
		  write (6,'(/1x,a,i4,3a,$)') 'Using Option: ',ioption,'. Output file name (',tempfile(1:margo(tempfile)),') = '
		  if (largo(popprm) > 0 .and. jdef < 20) then
!			Build file name to store stochastic imf
			ii = index(tempfile,'_z') + 1
			jj = index(tempfile,'yr') + 1
			if (index(tempfile,'IMF') > 0) then
				if (index(tempfile,'NumIMF') > 0) then
					kk = index(tempfile,'IMF') - 3
				elseif (index(tempfile,'S2plIMF') > 0) then
					kk = index(tempfile,'IMF') - 4
				endif
				subdir = tempfile(kk:jj) // '.binned_imf'
			else
				kk = index(tempfile,'_S') + 1
				if (kk == 1 ) kk = index(tempfile,'_Num') + 1
				auxname = tempfile(ii:kk-2)
				subdir  = tempfile(kk:jj)
				ii = index(subdir,'_')
				subdir = subdir(:ii) // trim(auxname) // subdir(ii:largo(subdir)) // '.binned_imf'
			endif
			open (668,file=subdir,status='unknown',form='unformatted')
			if (wrimf) then
				ii = index(subdir,'.binned_imf')
				subdir = subdir(:ii) // 'unbinned_imf'
				open (669,file=subdir,status='unknown',form='formatted')
				! write file name in temporary file
				write (679,'(a)') trim(subdir)
			endif
		  endif
	      else
		  write (6,'(/1x,a,$)') 'Output file name = '
	      endif
	      read (5,'(a)',end=1) evofile
              if (largo(evofile) == 0) then
			evofile=tempfile
	      elseif (largo(evofile) == 1) then
			evofile='junk'
	      endif
	      call mixfilenames(tempfile,evofile,ioption)
	  else
	      evofile=subdir(1:largo(subdir)) // tempfile
	  endif
	  write (6,'(/1x,2a )') 'Generic file name for model = ',evofile(1:margo(evofile))
        else
          	write (6,'(/1x,2a )') 'Generic file name for model = ',evofile(1:margo(evofile))
        endif
	savefile=evofile
	if (jun.eq.-987) return
	lun=7
	jun=lun+14

!	File to write sed
	nm = -1 + index(evofile,'.ised')
	if (nm < 0) then
		evofile = evofile(1:largo(evofile)) // '.ised'
	else
		call chaext(evofile,'ised',nm)
	endif
	close (lun)
	open (lun,file=evofile,form='unformatted',status='unknown')

!	Open files to write model colors (14 files: lun+1 - lun+14)
	call open_color_files(lun,evofile)
	if (wrimf) then
		write (679,'(a)') trim(evofile)
		close (679)
	endif

!	Check for more files to open
	if (kdist.eq.0) then
		evofile=savefile
		return
	endif

!	Open files to write sed''s of each stellar group (20 files: lun+15 - lun+34)
	if (kdist == 1 .or. kdist == 3) then
		do i=1,nph
		evofile=savefile(1:largo(savefile)) // '_' // phase(i)
		evofile= evofile(1:largo(evofile))  // '.ised'
!		call chaext(evofile,'ised',nm)
		close (lun+14+i)
!		write (6,*) 'i,lun,lun+14+i =',i,lun,lun+14+i,' ',evofile(1:largo(evofile))
		open (lun+14+i,file=evofile,form='unformatted',status='unknown')
		enddo
		evofile=savefile
	endif

!	Open files to write % contribution of each stellar phase (26 files: lun+31 - lun+56)
	if (kdist.eq.2.or.kdist.eq.3) then
		do i=1,nb
		call chaext(evofile,fract(i),nm)
		close (lun+30+i)
!		write (6,*) 'i,lun,lun+30+i =',i,lun,lun+30+i,' ',evofile(1:largo(evofile))
		open (lun+30+i,file=evofile,status='unknown')
       		call file_header(lun+30+i,evofile,0)
		write (lun+30+i,105) fract(i)(1:4)
105		format ('#log-age-yr  %MS        %SGB       %RGB       %CHeB      %AGB       %TP-AGB    %CSPN      %WD        %CCB       %NO-NAME   TOTAL-',a)   
		enddo
		do i=1,nb
		call chaext(evofile,stage(i),nm)
		close (lun+30+nb+i)
!		write (6,*) 'i,lun,lun+30+nb+i =',i,lun,lun+30+nb+i,' ',evofile(1:largo(evofile))
		open (lun+30+nb+i,file=evofile,status='unknown')
       		call file_header(lun+30+nb+i,evofile,0)
		write (lun+30+nb+i,106) stage(i)(2:4)
106		format ('#log-age-yr    nMS   nSGB   nRGB   nCHeB   nAGB  nCSPN    nWD   nCCB  nTRPH   NMS   NSGB   NRGB   NCHeB   NAGB  NCSPN    NWD   NCCB  NTRPH ',1x,a)
		enddo
	endif
	l=index(evofile,'.legus') - 1
	if (kdist.eq.4) then
!		l=index(evofile,'.') - 1
		evofile=evofile(1:l) // '_msOB.ised'
!		write (6,*) 'jun+1 =',jun+1,' ',evofile(1:largo(evofile))
		open (jun+1,file=evofile,form='unformatted',status='unknown')
		evofile=evofile(1:l) // '_msOBA.ised'
!		write (6,*) 'jun+2 =',jun+2,' ',evofile(1:largo(evofile))
		open (jun+2,file=evofile,form='unformatted',status='unknown')
		evofile=evofile(1:l) // '_rest.ised'
!		write (6,*) 'jun+3 =',jun+3,' ',evofile(1:largo(evofile))
		open (jun+3,file=evofile,form='unformatted',status='unknown')
	elseif (kdist.eq.5) then
!		l=index(evofile,'.') - 1
!		evofile=evofile(1:l) // '_ms.ised'
		evofile=evofile(1:l) // '_hm.ised'
!		write (6,*) 'jun+1 =',jun+1,' ',evofile(1:largo(evofile))
		open (jun+1,file=evofile,form='unformatted',status='unknown')
!		evofile=evofile(1:l) // '_mms.ised'
		evofile=evofile(1:l) // '_lm.ised'
!		write (6,*) 'jun+2 =',jun+2,' ',evofile(1:largo(evofile))
		open (jun+2,file=evofile,form='unformatted',status='unknown')
		evofile=evofile(1:l) // '.mass_contr'
!		write (6,*) 'lun+17 =',lun+17,' ',evofile(1:largo(evofile))
		open (lun+17,file=evofile,form='formatted',status='unknown')
		call file_header(lun+17,evofile,0)
		write (lun+17,'(a)') '#    (1)       (2)         (3)         (4)         (5)         (6)         (7)'
		write (lun+17,'(a)') '#'
		write (lun+17,'(a)') '#log-age-yr  Mass_high   Mass_low    Mass_tot      Mcut        M*_up       M*_low'
	elseif (kdist.eq.6) then
!		l=index(evofile,'.') - 1
		evofile=evofile(1:l) // '_coll_merged_3_5.ised'
!		write (6,*) 'jun+1 =',jun+1,' ',evofile(1:largo(evofile))
		open (jun+1,file=evofile,form='unformatted',status='unknown')
		evofile=evofile(1:l) // '_coll_merged_6_7.ised'
!		write (6,*) 'jun+2 =',jun+2,' ',evofile(1:largo(evofile))
		open (jun+2,file=evofile,form='unformatted',status='unknown')
		evofile=evofile(1:l) // '_coll_merged_3_5.indx'
!		write (6,*) 'jun+3 =',jun+3,' ',evofile(1:largo(evofile))
		open (jun+3,file=evofile,form='unformatted',status='unknown')
		evofile=evofile(1:l) // '_coll_merged_6_7.indx'
!		write (6,*) 'jun+4 =',jun+4,' ',evofile(1:largo(evofile))
		open (jun+4,file=evofile,form='unformatted',status='unknown')
	endif

	evofile=savefile
	return
1	stop
	end

	SUBROUTINE OPEN_COLOR_FILES(LUN,EVOFILE)

!	Open files to write colors and other model properties

	character evofile*(*)
	include 'stelib.dec'
        common /vazdekis/ jvaz,xmu

        ! Check if model is for Vazdekis et al. (1996) IMF
	if (jvaz == 0) then
		i = index(evofile,'_v')
		if (i > 0 .and. evofile(i+3:i+3) == 'p' .and. evofile(i+6:i+6) == '_') jvaz =1
	endif

!	File to write spectral indices computed from fits (1-13)
!	write (6,*) 'lun,lun+6 =',lun,lun+6,' ',evofile(1:largo(evofile))
	call chaext(evofile,'6lsindx_ffn',nm)
       	close (lun+6)
       	open (lun+6,file=evofile,status='unknown')
       	call file_header(lun+6,evofile,0)
	write (lun+6,'(a)') '#    (1)       (2)      (3)      (4)      (5)      (6)      (7)      (8)      (9)      (10)     (11)     (12)     (13)     (14)     (15)     (16)     (17)     (18)     (19)     (20)     (21)     (22)'
	write (lun+6,'(a)') '#'
	write (lun+6,'(a)') '# Index_No.      1:       2:       3:       4:       5:       6:       7:       8:       9:      10:      11:      12:      13:      14:      15:      16:      17:      18:      19:      20:      21:'
	write (lun+6,'(a)') '# log-age      CN_1     CN_2   Ca4227    G4300   Fe4383   Ca4455   Fe4531   Fe4668    Hbeta   Fe5015     Mg_1     Mg_2     Mg-b   Fe5270   Fe5335   Fe5406   Fe5709   Fe5782     Na-D    TiO_1    TiO_2'
	write (lun+6,'(a)') '#   (yr)      (mag)    (mag)     (A)      (A)      (A)      (A)      (A)      (A)      (A)      (A)     (mag)    (mag)     (A)      (A)      (A)      (A)      (A)      (A)      (A)     (mag)    (mag)'

!	File to write spectral indices measured from sed (1-13)
!	write (6,*) 'lun,lun+8 =',lun,lun+8,' ',evofile(1:largo(evofile))
	call chaext(evofile,'6lsindx_sed',nm)
       	close (lun+8)
       	open (lun+8,file=evofile,status='unknown')
       	call file_header(lun+8,evofile,0)
	write (lun+8,'(a)') '#    (1)       (2)      (3)      (4)      (5)      (6)      (7)      (8)      (9)      (10)     (11)     (12)     (13)     (14)     (15)     (16)     (17)     (18)     (19)     (20)     (21)     (22)'
	write (lun+8,'(a)') '#'
	write (lun+8,'(a)') '# Index_No.      1:       2:       3:       4:       5:       6:       7:       8:       9:      10:      11:      12:      13:      14:      15:      16:      17:      18:      19:      20:      21:'
	write (lun+8,'(a)') '# log-age      CN_1     CN_2   Ca4227    G4300   Fe4383   Ca4455   Fe4531   Fe4668    Hbeta   Fe5015     Mg_1     Mg_2     Mg-b   Fe5270   Fe5335   Fe5406   Fe5709   Fe5782     Na-D    TiO_1    TiO_2'
	write (lun+8,'(a)') '#   (yr)      (mag)    (mag)     (A)      (A)      (A)      (A)      (A)      (A)      (A)      (A)     (mag)    (mag)     (A)      (A)      (A)      (A)      (A)      (A)      (A)     (mag)    (mag)'

!	File to write spectral indices measured from sed transformed to lick system(1-13)
!	write (6,*) 'lun,lun+11 =',lun,lun+11,' ',evofile(1:largo(evofile))
	call chaext(evofile,'6lsindx_sed_lick_system',nm)
	close (lun+11)
	open (lun+11,file=evofile,status='unknown')
	call file_header(lun+11,evofile,0)
	write (lun+11,'(a)') '#    (1)      (2)     (3)     (4)     (5)     (6)     (7)     (8)     (9)     (10)    (11)    (12)    (13)    (14)    (15)    (16)    (17)    (18)    (19)    (20)    (21)    (22)'
	write (lun+11,'(a)') '#'
	write (lun+11,'(a)') '# Index_No.     1:      2:      3:      4:      5:      6:      7:      8:      9:     10:     11:     12:     13:     14:     15:     16:     17:     18:     19:     20:     21:'
	write (lun+11,'(a)') '# log-age     CN_1    CN_2  Ca4227   G4300  Fe4383  Ca4455  Fe4531  Fe4668   Hbeta  Fe5015    Mg_1    Mg_2    Mg-b  Fe5270  Fe5335  Fe5406  Fe5709  Fe5782    Na-D   TiO_1   TiO_2'
	write (lun+11,'(a)') '#   (yr)     (mag)   (mag)    (A)     (A)     (A)     (A)     (A)     (A)     (A)     (A)    (mag)   (mag)    (A)     (A)     (A)     (A)     (A)     (A)     (A)    (mag)   (mag)'

!	File to write spectral indices computed from fits (14-26)
!	write (6,*) 'lun,lun+7 =',lun,lun+7,' ',evofile(1:largo(evofile))
	call chaext(evofile,'7lsindx_ffn',nm)
       	close (lun+7)
       	open (lun+7,file=evofile,status='unknown')
       	call file_header(lun+7,evofile,0)
	write (lun+7,'(a)') '#    (1)          (2)        (3)        (4)        (5)        (6)'
	write (lun+7,'(a)') '#'
	write (lun+7,'(a)') '# Index_No.      WO-1:      WO-2:      WO-3:      WO-4:      GC-1:'
	write (lun+7,'(a)') '# log-age      HdeltaA    HgammaA    HdeltaF    HgammaF    D(4000)'
	write (lun+7,'(a)') '#   (yr)          (A)        (A)        (A)        (A)         . '

!	File to write spectral indices measured from sed (14-26)
!	write (6,*) 'lun,lun+9 =',lun,lun+9,' ',evofile(1:largo(evofile))
	call chaext(evofile,'7lsindx_sed',nm)
	close (lun+9)
	open (lun+9,file=evofile,status='unknown')
	call file_header(lun+9,evofile,0)
	write (lun+9,'(a)') '#    (1)        (2)      (3)      (4)      (5)      (6)      (7)       (8)        (9)        (10)      (11)      (12)      (13)      (14)      (15)'
	write (lun+9,'(a)') '#'
	write (lun+9,'(a)') '# Index_No.    WO-1:    WO-2:    WO-3:    WO-4:    GC-1:    4000A    DTT-Ca1    DTT-Ca2    DTT-Ca3   DTT-MgI     DM-04     DM-04     DM-04       BH:'
	write (lun+9,'(a)') '# log-age    HdeltaA  HgammaA  HdeltaF  HgammaF  D(4000)    B4_VN   CaII8498   CaII8542   CaII8662   MgI8807   H8_3889   H9_3835  H10_3798     BH-HK'
	write (lun+9,'(a)') '#   (yr)        (A)      (A)      (A)      (A)       .        .         (A)        (A)        (A)       (A)       (A)       (A)       (A)        (A)'

!	File to write spectral indices measured from sed transformed to lick system(14-26)
!	write (6,*) 'lun,lun+12 =',lun,lun+12,' ',evofile(1:largo(evofile))
	call chaext(evofile,'7lsindx_sed_lick_system',nm)
	close (lun+12)
	open (lun+12,file=evofile,status='unknown')
	call file_header(lun+12,evofile,0)
	write (lun+12,'(a)') '#    (1)        (2)      (3)      (4)      (5)      (6)      (7)       (8)        (9)        (10)      (11)      (12)      (13)      (14)      (15)'
	write (lun+12,'(a)') '#'
	write (lun+12,'(a)') '# Index_No.    WO-1:    WO-2:    WO-3:    WO-4:    GC-1:    4000A    DTT-Ca1    DTT-Ca2    DTT-Ca3   DTT-MgI     DM-04     DM-04     DM-04       BH:'
	write (lun+12,'(a)') '# log-age    HdeltaA  HgammaA  HdeltaF  HgammaF  D(4000)    B4_VN   CaII8498   CaII8542   CaII8662   MgI8807   H8_3889   H9_3835  H10_3798     BH-HK'
	write (lun+12,'(a)') '#   (yr)        (A)      (A)      (A)      (A)       .        .         (A)        (A)        (A)       (A)       (A)       (A)       (A)        (A)'

!	File to write Fanelli et al. UV spectral indices
!	write (6,*) 'lun,lun+87 =',lun,lun+87,' ',evofile(1:largo(evofile))
	call chaext (evofile,'8lsindx_sed',nm)
	close (lun+87)
	open (lun+87,file=evofile,status='unknown')
	call file_header(lun+87,evofile,0)
	write (lun+87,'(a)') '#'
	write (lun+87,'(a)') '#'
	write (lun+87,'(a)') '#    (1)        (2)        (3)        (4)        (5)        (6)        (7)        (8)        (9)        (10)       (11)       (12)       (13)       (14)       (15)       (16)       (17)       (18)       (19)       (20)       (21)       (22)'
	write (lun+87,'(a)') '#'
	write (lun+87,'(a)') '# log_t(yr)    BL1302      SiIV      BL1425     Fe1453    CIV1548a   CIV1548c   CIV1548e    BL1617     BL1664     BL1719     BL1853    FeII2402    BL2538    FeII2609     MgII       MgI       Mgwide      FeI       BL3096     CIVabs     HeIIems'

!	File to write fluxes used to compute spectral indices
!	write (6,*) 'lun,lun+10 =',lun,lun+10,' ',evofile(1:largo(evofile))
	call chaext(evofile,'9lsindx_sed_fluxes',nm)
       	close (lun+10)
       	open (lun+10,file=evofile,form='formatted',status='unknown')
       	call file_header(lun+10,evofile,0)
	write (lun+10,'(a)') '# log-age   Nx  Im   Flux_Blue   Flux_Red    Flux_Ctrl   Flux_Line       Index'

!	Open file .3color to write colors and other model properties
!	write (6,*) 'lun,lun+3 =',lun,lun+3,' ',evofile(1:largo(evofile))
	call chaext(evofile,'3color',nm)
	close (lun+3)
       	open (lun+3,file=evofile,status='unknown')
       	call file_header(lun+3,evofile,0)
	!rite (lun+3,'(a)') '#    (1)        (2)      (3)       (4)        (5)         (6)       (7)      (8)       (9)         (10)      (11)        (12)        (13)        (14)        (15)        (16)        (17)        (18)'
	!rite (lun+3,'(a)') '#'
       	!rite (lun+3,'(a)') '#log-age-yr    B4000     B4_VN   B4_SDSS      B912        NLy      NHeI      NHeII  NHeII/NLy       Mbol    Bol_Flux    SNR/yr/Lo    N(BH)       N(NS)     PNBR/yr/Lo    N(WD)    M(Remnants)   Lx(Lo)'
	!rite (lun+3,'(a)') '#    (1)        (2)      (3)       (4)        (5)         (6)       (7)      (8)       (9)       (10)        (11)        (12)        (13)        (14)        (15)        (16)        (17)        (18)        (19)        (20)'
	!rite (lun+3,'(a)') '#'
       	!rite (lun+3,'(a)') '#log-age-yr    B4000     B4_VN   B4_SDSS      B912        NLy      NHeI      NHeII  NHeII/NLy    SNe/yr     PISNe/yr    RegSNe/yr   IaSNe/yr  FailedSNe/yr   N(BH)       N(NS)      PNBR/yr      N(WD)     M(Remnants)   Lx(Lo)'
	write (lun+3,'(a)') '#    (1)       (2)     (3)      (4)      (5)        (6)      (7)      (8)      (9)      (10)        (11)        (12)        (13)        (14)        (15)        (16)        (17)        (18)        (19)'
	write (lun+3,'(a)') '#'
	write (lun+3,'(a)') '#log-age-yr   B4000    B4_VN  B4_SDSS    B912       NLy      NHeI    NHeII  NHeII/NLy   SNe/yr     PISNe/yr    RegSNe/yr   IaSNe/yr  FailedSNe/yr  PNBR/yr      N(BH)       N(NS)       N(WD)       Lx(Lo)'

!	Open file .4color to write colors and other model properties
!	write (6,*) 'lun,lun+4 =',lun,lun+4,' ',evofile(1:largo(evofile))
	call chaext(evofile,'4color',nm)
	close (lun+4)
	open (lun+4,file=evofile,status='unknown')
       	call file_header(lun+4,evofile,0)
	write (lun+4,'(a)') '#    (1)       (2)       (3)       (4)        (5)       (6)          (7)          (8)          (9)         (10)         (11)           (12)        (13)        (14)        (15)        (16)        (17)'
	write (lun+4,'(a)') '#                                                                                                                      M*_tot='
	write (lun+4,'(a)') '#log-age-yr    Mbol      Bmag      Vmag      Kmag      M*_liv     M_remnants    M_ret_gas    M_galaxy      SFR/yr    M*_liv+M_rem   M*_tot/Lb   M*_tot/Lv   M*_tot/Lk   M*_liv/Lb   M*_liv/Lv   M*_liv/Lk'

!	Open file .5color to write colors and other model properties
!	write (6,*) 'lun,lun+5 =',lun,lun+5,' ',evofile(1:largo(evofile))
	call chaext(evofile,'5color',nm)
	close (lun+5)
	open (lun+5,file=evofile,status='unknown')
       	call file_header(lun+5,evofile,0)
	write (lun+5,'(a)') '#    (1)       (2)         (3)          (4)          (5)          (6)          (7)          (8)'
	write (lun+5,'(a)') '#                                                                          Total Mass   Total Dust'
	write (lun+5,'(a)') '#log-age-yr    Mbol    b(t)*''s/yr   B(t)/yr/Lo  Turnoff_mass   BPMS/BMS     Loss Rate   Prod. Rate'

!	Open file .1VEGAmag to write Vega magnitudes
!	write (6,*) 'lun,lun+1 =',lun,lun+1,' ',evofile(1:largo(evofile))
	call chaext(evofile,'1VEGAmag',nm)
	close (lun+1)
	open (lun+1,file=evofile,status='unknown')
	call file_header(lun+1,evofile,0)
	write (lun+1,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)'
	write (lun+1,'(a)') '#'
	write (lun+1,'(a)') '#log-age-yr    Mbol       U         B2        B3        V         Rc        Ic        R         I         J         K         L        PalJ      PalH      PalK       K'''

!	Open file .2VEGAmag to write Vega magnitudes
!	write (6,*) 'lun,lun+2 =',lun,lun+2,' ',evofile(1:largo(evofile))
	call chaext(evofile,'2VEGAmag',nm)
	close (lun+2)
	open (lun+2,file=evofile,status='unknown')
	call file_header(lun+2,evofile,0)
	write (lun+2,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)'
	write (lun+2,'(a)') '#'
	write (lun+2,'(a)') '#log-age-yr    Mbol      Vmag     2MASSJ    2MASSH   2MASSKs   IRAC3.5   IRAC4.5   IRAC5.7   IRAC7.9    IRAS12    IRAS25    IRAS60   IRAS100    MIPS24    MIPS70   MIPS160'

!	Open files to write AB mag (SDSS, CFHT MegaCam, and GALEX)
!	write (6,*) 'lun,lun+13 =',lun,lun+13,' ',evofile(1:largo(evofile))
	call chaext(evofile,'1ABmag',nm)
	close (lun+13)
	open (lun+13,file=evofile,status='unknown')
	call file_header(lun+13,evofile,0)
	write (lun+13,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)          (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)         (17)      (18)         (19)      (20)      (21)         (22)         (23)'
	write (lun+13,'(a)') '#             <---------------- SDSS AB mag --------------->       <---------------------- CFHT MegaCam 1st and 3rd generation filters AB mag -------------------->       2Mass    CFHT Ks      <---------------- GALEX AB mag and flux --------------->'
	write (lun+13,'(a)') '#log-age-yr     u         g         r         i         z            u1        u3        g1        g3        r1        r3        i2        i3        z1        z3           H         Ks          FUV       NUV      F(FUV)       F(NUV)       F(1500A)'

!	Open file to write HST WFC and ACS bands in ABmag
!	write (6,*) 'lun,lun+88 =',lun,lun+88,' ',evofile(1:largo(evofile))
	call chaext(evofile,'2ABmag',nm)
	close (lun+88)
	open (lun+88,file=evofile,status='unknown')
	call file_header(lun+88,evofile,0)
	write (lun+88,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)      (18)      (19)      (20)      (21)      (22)      (23)      (24)      (25)      (26)      (27)      (28)      (29)      (30)      (31)      (32)      (33)'
	write (lun+88,'(a)') '#             <------------------------------------------- HST WFC3 UVIS1 LEGUS ABmag ------------------------------------------->    <----------------------------------- HST WFC3 ABmag --------------------------------->    <------------------------------------------- HST ACS WFC ABmag ------------------------------------------>'
	write (lun+88,'(a)') '#log-age-yr    F225w     F275w     F336w     F438w     F547m     F555w     F606w     F625w     F656n     F657n     F658n     F814w     F110W     F125W     F160W     F225W     F336W    FR388N     F438W     F555W     F814W     F220w     F250w     F330w     F410w     F435w     F475w     F555w     F606w     F625w     F775w     F814w'

	return
	end

	function jjdef(name)

!	Returns the value of jdef corresponding to the IMF selected in name

!	Array declaration
	character name*(*),imf(-1:18)*6
        common /vazdekis/ jvaz,xmu
	data imf /'chab','chab','kroup','KROUP','salp','scalo','SCALO','mi_sc','ktg93','kenn', 'leque','ferra','piott','carig','thIMF','myIMF','vxpxx','Schab','Skroup','Ssalp'/

!	Check for Vazdekis IMF
	if (jvaz > 0) then
		write (imf(15)(2:),'(f4.2)') xmu
		imf(15)(3:3) = 'p'
	endif

!	Get jjdef
	do jjdef=0,18
	k=index(name,imf(jjdef)(1:largo(imf(jjdef))))
	if (k.gt.0) return
	enddo
	write (6,*) 'Unknown IMF coding in this file name'
	jjdef=15
	return
	end

	subroutine check2007(name,cb07)
	
!	Check if Padova tracks are selected ( = Padova 1994/2000/2007/2008)
!	Check if Weiss 2008 tracks are selected
!	Modifies file name accordingly

	logical cb07
	character*(*) name
        character  a94(7)*3,b94(7)*3
!	character  a00(6)*4,b00(6)*4
!	character  a20(6)*3,b20(6)*4
!	character  a07(8)*4,b07(8)*3
        character  a08(10)*4,b08(10)*3
        character  a08b(10)*4,b08b(10)*3
        character  a10(10)*4,b10(10)*3
        character  a13(17)*4,a14(17)*4,a15(17)*4,a16(17)*4,a17(17)*4,b13(17)*9
        character  w07(6)*2,v07(6)*6
        character  w08(12)*4,v08(12)*6
        character  w09(12)*4,v09(12)*6

!	CB2007 models
	data a94/'m92', 'm93', 'm94', 'm95', 'm96', 'm97', 'm98'/
	data b94/'m22', 'm32', 'm42', 'm52', 'm62', 'm72', 'm82'/

!	CB2009 models (calibrated TP-AGB)
	data a08b/'m100','m101','m102','m103','m104','m105','m106','m107','m108','m109'/
	data b08b/ 'c00', 'c02', 'c12', 'c22', 'c32', 'c42', 'c52', 'c62', 'c72', 'c82'/

!	CB2009 models (uncalibrated TP-AGB)
	data a08 /'m110','m111','m112','m113','m114','m115','m116','m117','m118','m119'/
	data b08 / 'u00', 'u02', 'u12', 'u22', 'u32', 'u42', 'u52', 'u62', 'u72', 'u82'/

!	CB2010 models
	data a10 /'m210','m211','m212','m213','m214','m215','m216','m217','m218','m219'/
	data b10 / 'o00', 'o02', 'o12', 'o22', 'o32', 'o42', 'o52', 'o62', 'o72', 'o82'/

!	CB2013 models
	data b13 /'z0000y230','z0001y249','z0002y249','z0005y249','z001y250','z002y252','z004y256','z006y259','z008y263','z010y267','z014y273','z017y279','z020y284','z030y302','z040y321','z050y339','z060y356'/
	data a13 /     'm300',     'm301',     'm302',     'm303',    'm304',    'm305',    'm306',    'm307',    'm308',    'm309',    'm310',    'm311',    'm312',    'm313',    'm314',    'm315',    'm316'/  ! SET 1
!	CB2019 models
	data a14 /     'm400',     'm401',     'm402',     'm403',    'm404',    'm405',    'm406',    'm407',    'm408',    'm409',    'm410',    'm411',    'm412',    'm413',    'm414',    'm415',    'm416'/  ! VW PAGB
	data a15 /     'm500',     'm501',     'm502',     'm503',    'm504',    'm505',    'm506',    'm507',    'm508',    'm509',    'm510',    'm511',    'm512',    'm513',    'm514',    'm515',    'm516'/  ! MB PAGB
!	CB2021 models
	data a16 /     'm600',     'm601',     'm602',     'm603',    'm604',    'm605',    'm606',    'm607',    'm608',    'm609',    'm610',    'm611',    'm612',    'm613',    'm614',    'm615',    'm616'/  ! VW PAGB
	data a17 /     'm700',     'm701',     'm702',     'm703',    'm704',    'm705',    'm706',    'm707',    'm708',    'm709',    'm710',    'm711',    'm712',    'm713',    'm714',    'm715',    'm716'/  ! MB PAGB

!	Coelho et al models
	data w07/    'M1',     'M2',     'M3',     'M4',     'M5',     'M6'/
	data v07/'p00p00', 'p00p04', 'm05p00', 'm05p04', 'p02p00', 'p02p04'/

!	Coelho et al models
	data w08/  'm181',   'm182',   'm183',   'm184',   'm185',   'm186',   'm187',   'm188',   'm189',   'm190',   'm191',   'm192'/
	data v08/'m10p00', 'm13p04', 'm10p04', 'm05p00', 'm08p04', 'm05p04', 'p00p00', 'm03p04', 'p02p00', 'm01p04', 'p00p04', 'p02p04'/

!	Coelho et al models
	data w09/  'm281',   'm282',   'm283',   'm284',   'm285',   'm286',   'm287',   'm288',   'm289',   'm290',   'm291',   'm292'/
	data v09/'m10p00', 'm13p04', 'm10p04', 'm05p00', 'm08p04', 'm05p04', 'p00p00', 'm03p04', 'p02p00', 'm01p04', 'p00p04', 'p02p04'/

!	Store in common model option
	common /model/ kopt

!	Look for different option numbers to correct file name
	kopt=3
	cb07=.false.
	do i=1,7
	j=index(name,a94(i))
	if (j.gt.0) then
		name(j:j+2) = b94(i)
		name(1:2) = 'cb'
		name(6:6) = '7'
		cb07=.true.
		kopt=7
		return
	endif
	enddo
	do i=1,10
	j=index(name,a08(i))
	if (j.gt.0) then
		name(j:j+2) = b08(i)
		name(1:2) = 'cb'
		name(6:6) = '9'
		name(j+3:)=name(j+4:)
		cb07=.true.
		return
	endif
	enddo
	do i=1,10
	j=index(name,a08b(i))
	if (j.gt.0) then
		name(j:j+2) = b08b(i)
		name(1:2) = 'cb'
		name(6:6) = '9'
		name(j+3:)=name(j+4:)
		cb07=.true.
		return
	endif
	enddo
	do i=1,10
	j=index(name,a10(i))
	if (j.gt.0) then
		name(j:j+2) = b10(i)
		name(1:2) = 'cb'
		name(5:6) = '10'
		name(j+3:)=name(j+4:)
		cb07=.true.
		kopt=0
		return
	endif
	enddo
	do i=1,15
	j=index(name,a13(i))
	if (j.gt.0) then
		js=index(name,'_ssp')
		name(11:) = b13(i) // name(j+4:js-1) // name(7:j-1) // 'ssp'
		name(1:10) = 'cb2013_s1_'
		call no_bl(name)
		cb07=.true.
		kopt=13
		return
	endif
	enddo
	do i=1,17
	j=index(name,a14(i))
	if (j.gt.0) then
		js=index(name,'_ssp')
		name(11:) = b13(i) // name(j+4:js-1) // name(7:j-1) // 'ssp'
		name(1:10) = 'cb2013_s2_'
		call no_bl(name)
!		write (6,*)
!		write (6,*) 'Old style file name: ',trim(name)
		name(1:10) = 'bc2019_   '		! Vassiliadis & Wood (1994) WD tracks with 2020 version of PoWR models
		call no_bl(name)
		k = index(name,'y')
		if (name(k+4:k+4) /= '_') then
			name(k:k+4) = '     '		! No need for Y abundance in file name, nor "n" character related to TP-AGB choice
		else
			name(k:k+3) = '     '		! No need for Y abundance in file name
		endif
		call no_bl(name)
		cb07=.true.
		kopt=13
		return
	endif
	enddo
	do i=1,17
	j=index(name,a15(i))
	if (j.gt.0) then
		js=index(name,'_ssp')
		name(11:) = b13(i) // name(j+4:js-1) // name(7:j-1) // 'ssp'
		name(1:10) = 'cb2013_s2_'
		call no_bl(name)
!		write (6,*)
!		write (6,*) 'Old style file name: ',trim(name)
		name(1:10) = 'cb2019_   '		! Miller-Bertolami WD tracks with 2018 version of PoWR models
		call no_bl(name)
		k = index(name,'y')
		if (name(k+4:k+4) /= '_') then
			name(k:k+4) = '     '		! No need for Y abundance in file name, nor "n" character related to TP-AGB choice
		else
			name(k:k+3) = '     '		! No need for Y abundance in file name
		endif
		call no_bl(name)
		cb07=.true.
		kopt=13
		return
	endif
	enddo
	do i=1,17
	j=index(name,a16(i))
	if (j.gt.0) then
		js=index(name,'_ssp')
		name(11:) = b13(i) // name(j+4:js-1) // name(7:j-1) // 'ssp'
		name(1:10) = 'bc2021_   '		! Vassiliadis & Wood (1994) WD tracks with 2020 version of PoWR models
		call no_bl(name)
		k = index(name,'y')
		if (name(k+4:k+4) /= '_') then
			name(k:k+4) = '     '		! No need for Y abundance in file name, nor "n" character related to TP-AGB choice
		else
			name(k:k+3) = '     '		! No need for Y abundance in file name
		endif
		call no_bl(name)
		cb07=.true.
		kopt=13
		return
	endif
	enddo
	do i=1,17
	j=index(name,a17(i))
	if (j.gt.0) then
		js=index(name,'_ssp')
		name(11:) = b13(i) // name(j+4:js-1) // name(7:j-1) // 'ssp'
		call no_bl(name)
		name(1:10) = 'cb2021_   '		! Miller-Bertolami WD tracks with 2018 version of PoWR models
		call no_bl(name)
		k = index(name,'y')
		if (name(k+4:k+4) /= '_') then
			name(k:k+4) = '     '		! No need for Y abundance in file name, nor "n" character related to TP-AGB choice
		else
			name(k:k+3) = '     '		! No need for Y abundance in file name
		endif
		call no_bl(name)
		cb07=.true.
		kopt=13
		return
	endif
	enddo
	do i=1,6
	j=index(name,w07(i))
	if (j.gt.0) then
                name(j+6:)=name(j+2:)
		name(j:j+5) = v07(i)
		name(1:2) = 'bc'
		name(6:6) = '7'
		cb07=.true.
		kopt=7
		return
	endif
	enddo
	do i=1,12
	j=index(name,w08(i))
	if (j.gt.0) then
		name(j+6:)=name(j+4:)
!		name(j+10:)=name(j+4:)
		name(j:j+5) = v08(i)
!		name(j:j+9) = v08(i)
		name(1:2) = 'bc'
		name(6:6) = '9'
		name = 'c' // name
		cb07=.true.
		return
	endif
	enddo
	do i=1,12
	j=index(name,w09(i))
	if (j.gt.0) then
		name(j+6:)=name(j+4:)
!		name(j+10:)=name(j+4:)
		name(j:j+5) = v09(i)
!		name(j:j+9) = v09(i)
		name(1:2) = 'bc'
		name(5:6) = '12'
		name = 'c' // name
		cb07=.true.
		return
	endif
	enddo
!	write (6,*) 'Unmodified name: ',name(1:largo(name))
	return
	end

	SUBROUTINE NAME2016(A)

!	Changes bc03 and cb07 files name from old to new style (2016)

!	Variables
	character(*) a,b*512

	if (a == 'bc2003_p99_') return
	if ( index(a,'bc2003') > 0 .or. index(a,'cb2007') > 0 ) then
		i = index(a,'2_')-2
		j = index(a,'_ssp')
		if (a(i:i+1) == 'm2') then
			b = a(:6) //  '_z0001'
		elseif (a(i:i+1) == 'm3') then
			b = a(:6) //  '_z0004'
		elseif (a(i:i+1) == 'm4') then
			b = a(:6) //  '_z004'
		elseif (a(i:i+1) == 'm5') then
			b = a(:6) //  '_z008'
		elseif (a(i:i+1) == 'm6') then
			b = a(:6) //  '_z020'
		elseif (a(i:i+1) == 'm7') then
			b = a(:6) //  '_z050'
		elseif (a(i:i+1) == 'm8') then
			b = a(:6) //  '_z100'
		else
			return
		endif
!		write (6,*)
!		write (6,*) 'Old style file name: ',trim(a)
		b = trim(b) // a(i+3:j) // a(8:i-1) // a(j+1:)
		a = b
	else
		return
	endif
	end

	! Open file .1color to write colors and other model properties
	! write (6,*) 'lun,lun+1 =',lun,lun+1,' ',evofile(1:largo(evofile))
	! call chaext(evofile,'1color',nm)
	! close (lun+1)
	! open (lun+1,file=evofile,status='unknown')
	! call file_header(lun+1,evofile,0)
	! write (lun+1,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)        (7)       (8)       (9)      (10)      (11)      (12)      (13)      (14)      (15)'
	! write (lun+1,'(a)') '#'
	! write (lun+1,'(a)') '#log-age-yr    Mbol      Umag      Bmag      Vmag      Kmag      14-V      17-V      22-V      27-V       U-J       J-F       F-N       U-B       B-V'

	! Open file .2color to write colors and other model properties
	! write (6,*) 'lun,lun+2 =',lun,lun+2,' ',evofile(1:largo(evofile))
	! call chaext(evofile,'2color',nm)
	! close (lun+2)
	! open (lun+2,file=evofile,status='unknown')
	! call file_header(lun+2,evofile,0)
	! write (lun+2,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)      (14)      (15)'
	! write (lun+2,'(a)') '#'
	! write (lun+2,'(a)') '#log-age-yr    Rmag     J2Mmag     Kmag      V-R       V-I       V-J       V-K        R-I       J-H       H-K      V-K''      V-Ks    (J-H)2M   (J-Ks)2M'
	! write (lun+2,'(a)') '#log-age-yr    Rmag     J2Mmag     Kmag      V-R       V-I       V-J       V-K        R-I       J-H       H-K      V-K''      V-J2M     V-H2M    V-Ks2M'

	! Open file to write IR colors (IRAC, IRAS, MIPS)
	! write (6,*) 'lun,lun+14 =',lun,lun+14,' ',evofile(1:largo(evofile))
	! call chaext(evofile,'9color',nm)
	! close (lun+14)
	! open (lun+14,file=evofile,status='unknown')
	! call file_header(lun+14,evofile,0)
	! write (lun+14,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)       (14)         (15)'
	! write (lun+14,'(a)') '#                                                                                                                                    Total Mass   Total Dust'
	! write (lun+14,'(a)') '#log-age-yr    Kmag     K-I3.5  I3.5-I4.5 I4.5-I5.7 I5.7-I7.9 I7.9-I12   I12-I25   I25-I60  I60-I100   M24-M70  M70-M160 I100-M160   Loss Rate    Prod. Rate'

	! 1ABmag old header
	! write (lun+13,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)          (7)       (8)       (9)       (10)      (11)      (12)      (13)         (14)      (15)      (16)         (17)         (18)            (19)      (20)      (21)      (22)      (23)'
	! write (lun+13,'(a)') '#             <---------------- SDSS AB mag --------------->       <---------------------- CFHT MegaCam AB mag --------------------->      <---------------- GALEX AB mag and flux --------------->       <---- CFHT New Megaprime Filters AB mag ----->'
	! write (lun+13,'(a)') '#log-age-yr     u         g         r         i         z            u         g         r         i         y         z         Ks          FUV       NUV      F(FUV)       F(NUV)       F(1500A)          u''        g''        r''        i''        z'''
	! write (lun+13,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)          (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)         (17)         (18)      (19)      (20)         (21)         (22)'
	! write (lun+13,'(a)') '#             <---------------- SDSS AB mag --------------->       <---------------------- CFHT MegaCam 1st and 3rd generation filters AB mag -------------------->      CFHT Ks      <---------------- GALEX AB mag and flux --------------->'
	! write (lun+13,'(a)') '#log-age-yr     u         g         r         i         z            u1        u3        g1        g3        r1        r3        i2        i3        z1        z3           Ks          FUV       NUV      F(FUV)       F(NUV)       F(1500A)'

	! Open file to write HST WFPC2 + Johnson colors (requested by RAGL)
	! write (6,*) 'lun,lun+87 =',lun,lun+87,' ',evofile(1:largo(evofile))
	! call chaext(evofile,'wfpc2_johnson_color',nm)
	! close (lun+87)
	! open (lun+87,file=evofile,status='unknown')
	! call file_header(lun+87,evofile,0)
	! write (lun+87,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)      (18)'
	! write (lun+87,'(a)') '#'
	! write (lun+87,'(a)') '#log-age-yr    Vmag    V-F300w   V-F300rl  V-F336w   V-F439w   V-F450w   V-F555w   V-F606w   V-F675w   V-F814w     V-U       V-B       V-R       V-I       V-J       V-K       V-L'

	! Open file to write UVIS1 LEGUS colors (Aida Wofford) 2014
	! write (6,*) 'lun,lun+88 =',lun,lun+88,' ',evofile(1:largo(evofile))
	! call chaext(evofile,'wfc3_uvis1_legus_ABmag',nm)
	! close (lun+88)
	! open (lun+88,file=evofile,status='unknown')
	! call file_header(lun+88,evofile,0)
	! write (lun+88,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)'
	! write (lun+88,'(a)') '#'
	! write (lun+88,'(a)') '#log-age-yr    F225w     F275w     F336w     F438w     F547m     F555w     F606w     F625w     F656n     F657n     F658n     F814w'
	! write (lun+88,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)'
	! write (lun+88,'(a)') '#'
	! write (lun+88,'(a)') '#log-age-yr    F225w     F275w     F336w     F438w     F547m     F555w     F606w     F625w     F656n     F657n     F658n     F814w,    F110W     F125W     F160W     F225W     F336W    FR388N     F438W     F555W     F814W,    F220w     F250w     F330w     F410w     F435w     F475w     F555w     F606w     F625w     F775w     F814w'

	! Open file to write WFC3 colors
	! write (6,*) 'lun,lun+84 =',lun,lun+84,' ',evofile(1:largo(evofile))
	! call chaext(evofile,'wfc3_ABmag',nm)
	! close (lun+84)
	! open (lun+84,file=evofile,status='unknown')
	! call file_header(lun+84,evofile,0)
	! write (lun+84,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)'
	! write (lun+84,'(a)') '#'
	! write (lun+84,'(a)') '#log-age-yr    F110W     F125W     F160W     F225W     F336W    FR388N     F438W     F555W     F814W'

	! Open file to write HST ACS colors
	! write (6,*) 'lun,lun+86 =',lun,lun+86,' ',evofile(1:largo(evofile))
	! call chaext(evofile,'acs_wfc_ABmag',nm)
	! close (lun+86)
	! open (lun+86,file=evofile,status='unknown')
	! call file_header(lun+86,evofile,0)
	! write (lun+86,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)'
	! write (lun+86,'(a)') '#'
	! write (lun+86,'(a)') '#log-age-yr    F220w     F250w     F330w     F410w     F435w     F475w     F555w     F606w     F625w     F775w     F814w'

	! Open file to write WFC3.UVIS1 colors (Rupali Chandar) 2011
	! write (6,*) 'lun,lun+85 =',lun,lun+85,' ',evofile(1:largo(evofile))
	! !all chaext(evofile,'wfc3_uvis1_color',nm)
	! call chaext(evofile,'wfc3_uvis1_ABmag',nm)
	! close (lun+85)
	! open (lun+85,file=evofile,status='unknown')
	! call file_header(lun+85,evofile,0)
	! !rite (lun+85,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)      (18)'
	! !rite (lun+85,'(a)') '#'
	! !rite (lun+85,'(a)') '#log-age-yr    Vmag      Kmag    V-F225w   V-F336w   V-F438w   V-F547m   V-F555w   V-F606w   V-F625w   V-F656n   V-F657n   V-F658n   V-F814w   log Nly     Mt/Lb     Mt/Lv     Mt/Lk'
	! write (lun+85,'(a)') '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)      (14)      (15)'
	! write (lun+85,'(a)') '#'
	! write (lun+85,'(a)') '#log-age-yr    F225w     F336w     F390w     F438w     F475w     F547m     F555w     F606w     F625w     F656n     F657n     F658n     F775w     F814w'
