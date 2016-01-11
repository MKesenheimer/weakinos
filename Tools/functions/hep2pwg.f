cccccccccccccccccc transform hep-momenta to powheg-style
	subroutine phep2powheg(phep,ppwg)
	implicit none
	real*8 phep(1:4),ppwg(0:3)
	ppwg(1)=phep(1)
	ppwg(2)=phep(2)
	ppwg(3)=phep(3)
	ppwg(0)=phep(4)
	end