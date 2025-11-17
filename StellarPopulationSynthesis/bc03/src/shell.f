
	CHARACTER*32 FUNCTION SHELL()

	! Returns current shell

	! Variables
	implicit none
	character*32 csh,zsh

	! Get shell name from environment
        call getenv('SHELL',zsh)
        call getenv('shell',csh)
        shell = trim(csh) // trim(zsh)
	return
	end
