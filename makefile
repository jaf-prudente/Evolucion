FLAGS = -O2 -ffree-form -w 

OBJS = arrays.o main.o initial.o initialDirac0.o initialDirac1.o reader.o evolve.o sources.o \
sourcesDirac.o metric.o metricDirac.o potential.o energy.o energyDirac.o \
save0Ddata.o save1Ddata.o 


qstar : dir link
test : dir_test link_test

dir :
	@ mkdir -p exe

dir_test :
	@ mkdir -p test

link : $(OBJS)
	gfortran $(FLAGS) -o xQStar $(OBJS) 
	@ mv xQStar exe

link_test : $(OBJS)
	gfortran $(FLAGS) -o xQStar $(OBJS)
	@ mv xQStar test


arrays.o : arrays.f90
	gfortran $(FLAGS) -c arrays.f90

main.o : main.f90
	gfortran $(FLAGS) -c main.f90

initial.o : initial.f90
	gfortran $(FLAGS) -c initial.f90

initialDirac0.o : initialDirac0.f90
	gfortran $(FLAGS) -c initialDirac0.f90

initialDirac1.o : initialDirac1.f90
	gfortran $(FLAGS) -c initialDirac1.f90

reader.o : reader.f90
	gfortran $(FLAGS) -c reader.f90

evolve.o : evolve.f90
	gfortran $(FLAGS) -c evolve.f90

sources.o : sources.f90
	gfortran $(FLAGS) -c sources.f90

sourcesDirac.o : sourcesDirac.f90
	gfortran $(FLAGS) -c sourcesDirac.f90

metric.o : metric.f90
	gfortran $(FLAGS) -c metric.f90

metricDirac.o : metricDirac.f90
	gfortran $(FLAGS) -c metricDirac.f90

potential.o : potential.f90
	gfortran $(FLAGS) -c potential.f90

energy.o : energy.f90
	gfortran $(FLAGS) -c energy.f90

energyDirac.o : energyDirac.f90
	gfortran $(FLAGS) -c energyDirac.f90

save0Ddata.o : save0Ddata.f90
	gfortran $(FLAGS) -c save0Ddata.f90

save1Ddata.o : save1Ddata.f90
	gfortran $(FLAGS) -c save1Ddata.f90

clean :
#	/bin/rm -r xxx *.o *.mod *.vo *.inc
	/bin/rm -r exe/*.rl exe/*.tl *.o *.mod *.vo *.inc





