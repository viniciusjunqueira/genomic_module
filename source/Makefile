# Makefile for $(prog)
include Makeinit

prog=eigen

#dir= ../../libs
dir= /Users/vinicius/OneDrive/Cursos/Computer_programming_Animal_Breeding/Vinicius_Junqueira/guest28/course2016/libs

a.out:	$(prog).o $(dir)/model.o $(dir)/blupsubs.a $(dir)/sparsem.a \
	$(dir)/denseop.a
	$(f90)  $(optf90) $(prog).o $(dir)/model.o $(dir)/blupsubs.a \
	$(dir)/sparsem.a $(dir)/denseop.a \
	$(dir)/iounf.o  $(libf77)

	mv a.out $(prog)
	cp $(prog) ..

$(prog).o:	$(prog).f90 $(dir)/model.o $(dir)/sparsem.a $(dir)/denseop.a
	$(f90) $(optdir)$(dir) -c $(optf90) $(prog).f90


clean:
	rm *.o *.mod
