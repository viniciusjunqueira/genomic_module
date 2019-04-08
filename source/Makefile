# Makefile for $(prog)
include Makeinit
#include Makeinit

prog=genomic_module

dir=../libs

a.out:	$(prog).o $(dir)/model.o $(dir)/blupsubs.a $(dir)/sparsem.a $(dir)/denseop.a $(dir)/genomic.a $(dir)/stat.a
	$(f90)  $(optf90) $(mkl) $(prog).o $(dir)/denseop.a $(dir)/sparsem.a $(dir)/model.o $(dir)/blupsubs.a $(dir)/genomic.a $(dir)/stat.a
	cp a.out ../bin/$(prog)
	cp a.out ../examples/ex3/$(prog)

$(prog).o:	$(prog).f90 $(dir)/denseop.a $(dir)/model.o $(dir)/sparsem.a $(dir)/blupsubs.a $(dir)/genomic.a $(dir)/stat.a
	$(f90) $(optdir)$(dir) $(mkl) -c $(optf90) $(prog).f90


clean:
	rm *.o *.mod a.out *.r
	#rm a.out *.
