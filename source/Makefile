# Makefile for $(prog)
include Makeinit

prog=mkl_omp

a.out:	$(prog).o
	$(f90) $(mkldist) $(optf90) $(prog).o $(libf77)

$(prog).o:	$(prog).f90
	$(f90) $(optdir)$(dir) $(mkldist) -c $(mkldist) $(optf90) $(prog).f90


clean:
	rm *.o *.mod a.out
