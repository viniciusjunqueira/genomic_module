# Makes intergen

include Makeinit

all: denseop.a sparsem.a model.o blupsubs.a stat.a iounf.o 

denseop.a: lapack90r.o kind.o denseop.o  
	ar cr  denseop.a kind.o denseop.o lapack90r.o 

sparsem.a: kind.o hash.o sparse.o sparse2.o fspak90.o fspak.o fspaksub.o \
	pcg.o sparssub.o second.o 
	ar cr sparsem.a kind.o hash.o sparse.o sparse2.o fspak90.o \
           fspak.o fspaksub.o pcg.o sparssub.o second.o
	
blupsubs.a: blups1.o blups2.o ginv.o 
	ar cr blupsubs.a  blups1.o blups2.o ginv.o 
	
stat.a:	ranlib.o prob.o gibbs.o dcdflib.o
	ar cr stat.a prob.o  gibbs.o ranlib.o dcdflib.o

genomic.a: genomic.o genomicsub.o 
	ar cr genomic.a genomic.o genomicsub.o

genomicDist.a:	genomicDist.o genomicsub.o
	ar cr genomicDist.a genomicDist.o genomicsub.o

kind.o:	kind.f90
	$(f90) -c  $(optf90) kind.f90

model.o:	model.f90 
	$(f90) -c  $(optf90) model.f90

blups1.o:	blups1.f90 
	$(f90) -c  $(optf90) blups1.f90

blups2.o:	blups2.f 
	$(f90) -c  $(optf77) blups2.f

hash.o:	hash.f90
	$(f90) -c $(optf90)  hash.f90

sparse.o:	sparse.f90 
	$(f90) -c  $(optf90) sparse.f90

sparse2.o:	sparse2.f 
	$(f90) -c  $(optf77) sparse2.f

fspak90.o:	fspak90.f90 
	$(f90) -c  $(optf90) fspak90.f90

fspak.o:	fspak.f
	$(f90) -c  $(optf77) fspak.f

fspaksub.o:	fspaksub.f
	$(f90) -c  $(optf77) fspaksub.f

ginv.o:	ginv.f kind.f90
	$(f90) -c  $(optf77) ginv.f

sparssub.o:	sparssub.f 
	$(f90) $(optf77) -c sparssub.f

second.o:	second.f
	$(f90) -c  $(optf77) second.f

denseop.o:	denseop.f90
	$(f90) -c  $(optf90) denseop.f90

lapack90r.o:	lapack90r.f90 
	$(f90) -c  $(optf90) lapack90r.f90

pcg.o:	pcg.f90 
	$(f90) -c  $(optf90) pcg.f90

prob.o:	prob.f90
	$(f90) -c  $(optf90) prob.f90

problib.o:	problib.f90
	$(f90) -c  $(optf90) problib.f90
	
gibbs.o:	gibbs.f90
	$(f90) -c  $(optf90) gibbs.f90

ranlib.o:	ranlib.f
	$(f90) -c  $(optf77) ranlib.f

dcdflib.o:	dcdflib.f90
	$(f90) -c  $(optf90) dcdflib.f90

iounf.o:	iounf.f90
	$(f90) -c  $(optf90) iounf.f90

genomic.o:	genomic.f90
	$(f90) -c $(optf90) -fpp  genomic.f90 

genomicDist.o:	 genomic.f90
	 $(f90) -c $(optf90) -fpp -D _MYCHECK=1 genomic.f90 -o genomicDist.o

genomicsub.o:	genomicsub.F90
	$(f90) -c $(optf90)  genomicsub.F90

clean:
	rm *.o *.a *.mod #a.out
