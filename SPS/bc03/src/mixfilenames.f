	SUBROUTINE MIXFILENAMES(B,A,ioption)

	! Builds script to perform correct file assignments for program mix_stelib

	! Variables
	character*(*) a,b,c*512

	! Find first and third occurrence of character '_'
        c = trim(a) // '.ised'
	if (index(b,'BaSeL') > 0) then
		if (index(a,'BaSeL') == 0) then
			a = trim(a) // '_lr_BaSeL_ssp'
		endif
		i=index(a,'_lr')
		j=index(a(i+1:),'_')
		k=index(a(i+j+1:),'_')
		l=largo(a)
		c = a(:i) // 'lr_BaSeL' // a(i+j+k:l) // '.ised'
		open (120,file='mix_stelib.tmp')
		! Write script to use program mix_stelib
		! write (120,'(5a)') 'setenv name0  ',a(:i-1)
		! write (120,'(5a)') 'setenv name1  ',a(:i),'hrs_stelib' ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name2  ',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name3  ',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised'
		! write (120,'(1a)') '#'
		! write (120,'(5a)') 'setenv name4  ',a(:i),'lrs_hngsl'  ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name5  ',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name6  ',a(:i),'lr_xhngsl'  ,a(i+j+k:l),'.ised'
		! write (120,'(1a)') '#'
		! write (120,'(5a)') 'setenv name7  ',a(:i),'hrs_indous' ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name8  ',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name9  ',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised'
		! write (120,'(1a)') '#'
		! write (120,'(5a)') 'setenv name10 ',a(:i),'hrx_miles'  ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name11 ',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name12 ',a(:i),'hr_xmiless' ,a(i+j+k:l),'.ised'
		! write (120,'(1a)') '#'
		! write (120,'(5a)') 'setenv name13 ',a(:i),'hrx_miles'  ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name14 ',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised'
		! write (120,'(5a)') 'setenv name15 ',a(:i),'er_xmilesi' ,a(i+j+k:l),'.ised'
		! write (120,'(1a)') '#'
		write (120,'(5a)') 'export  name1="',a(:i),'hrs_stelib' ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export  name2="',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export  name3="',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised"'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'export  name4="',a(:i),'lrs_hngsl'  ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export  name5="',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export  name6="',a(:i),'lr_xhngsl'  ,a(i+j+k:l),'.ised"'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'export  name7="',a(:i),'hrs_indous' ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export  name8="',a(:i),'lr_BaSeL'   ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export  name9="',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised"'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'export name10="',a(:i),'hrx_miles'  ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export name11="',a(:i),'hr_stelib'  ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export name12="',a(:i),'hr_xmiless' ,a(i+j+k:l),'.ised"'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'export name13="',a(:i),'hrx_miles'  ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export name14="',a(:i),'hr_xindous' ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export name15="',a(:i),'hr_xmilesi' ,a(i+j+k:l),'.ised"'
		write (120,'(1a)') '#'
		write (120,'(5a)') 'export name16="',a(:i),'er_xmilesi' ,a(i+j+k:l),'.ised"'
		write (120,'(5a)') 'export name17="',a(:i),'er_xmilesi' ,a(i+j+k:l),'.all"'
		write (120,'(5a)') 'export name18="',a(:i),'ht_xmilesi' ,a(i+j+k:l),'.ised"'
		write (120,'(1a)') '#'
		write (120,'(a,i5,a)') 'export ioptn="',ioption,'"'
		write (120,'(1a)') '#'
		close (120)
	endif
	open (120,file='mix_stelib2.tmp')
	write (120,'(a)') 'export name0="' // trim(c) //'"'
	write (120,'(a,i5,a)') 'export ioptn="',ioption,'"'
	close (120)
	return
	end
