	SUBROUTINE DARK_GALAXY(LUN,T)

	! Includes age record in output files indicating galaxy is still dark

	! Variables
	data icall/0/

	! Write age record
	if (icall == 0) then
		icall = 1
		write (lun+1 ,107) 0.
		write (lun+2 ,107) 0.
		write (lun+3 ,107) 0.
		write (lun+4 ,107) 0.
		write (lun+5 ,107) 0.
		write (lun+6 ,107) 0.
		write (lun+7 ,107) 0.
		write (lun+8 ,107) 0.
		write (lun+9 ,107) 0.
		write (lun+10,107)(0.,i=1,35)
		write (lun+11,107) 0.
		write (lun+12,107) 0.
		write (lun+13,107) 0.
		write (lun+87,107) 0.
		write (lun+88,107) 0.
		write (2002,  107) 0.
	endif
	if (t <= 0.) return
	tl = alog10(t)
	write (lun+1 ,107) tl
	write (lun+2 ,107) tl
	write (lun+3 ,107) tl
	write (lun+4 ,107) tl
	write (lun+5 ,107) tl
	write (lun+6 ,107) tl
	write (lun+7 ,107) tl
	write (lun+8 ,107) tl
	write (lun+9 ,107) tl
	write (lun+10,107)(tl,i=1,35)
	write (lun+11,107) tl
	write (lun+12,107) tl
	write (lun+13,107) tl
	write (lun+87,107) tl
	write (lun+88,107) tl
	write (2002,  107) tl
107	format ('#',f9.6,'    Dark galaxy')
	return
	end

		! write (lun+1 ,101) 0.,bolmag,umag,bmag,vmag,kmag,(col(i),i=1,9)
		! write (lun+2 ,101) 0.,rmag,jmag,kmag,            (col(i),i=10,20)
		! write (lun+14,108) 0.,     kmag,                 (col(i),i=21,31)
		! write (lun+87,108) 0.,vmag,                      (col(i),i=65,80)
		! write (lun+84,108) 0.,vmag,kmag,                 (col(i),i=32,42),cly,blr,vlr,klr
		! write (lun+85,108) 0.,vmag,kmag,                 (col(i),i=43,53),cly,blr,vlr,klr
		! write (lun+86,108) 0.,vmag,kmag,                 (col(i),i=54,64),cly,blr,vlr,klr
		! write (lun+88,108) 0.,vmag,kmag,                 (col(i),i=81,92),cly,blr,vlr,klr
		! write (lun+17,171) 0.,mabo(ic),mbel(ic),mabo(ic)+mbel(ic),mcut,rmup(ic),rmlo(ic)    ! Esto escribe los datos requeridos por Alba (opcion kdist = 5)

	! write (lun+1 ,101) tl,bolmag,umag,bmag,vmag,kmag,(col(i),i=1,9)
	! write (lun+2 ,101) tl,rmag,jmag,kmag,            (col(i),i=10,20)
	! write (lun+14,108) tl,     kmag,                 (col(i),i=21,31)
	! write (lun+87,108) tl,vmag,                      (col(i),i=65,80)
	! write (lun+84,108) tl,vmag,kmag,                 (col(i),i=32,42),cly,blr,vlr,klr
	! write (lun+85,108) tl,vmag,kmag,                 (col(i),i=43,53),cly,blr,vlr,klr
	! write (lun+86,108) tl,vmag,kmag,                 (col(i),i=54,64),cly,blr,vlr,klr
	! write (lun+88,108) tl,vmag,kmag,                 (col(i),i=81,92),cly,blr,vlr,klr
	! write (lun+17,171) tl,mabo(ic),mbel(ic),mabo(ic)+mbel(ic),mcut,rmup(ic),rmlo(ic)    ! Esto escribe los datos requeridos por Alba (opcion kdist = 5)
	! 101 format (f10.6,14f10.4,1pe13.4)
	! 171 format (f10.6,1p6e12.3)
