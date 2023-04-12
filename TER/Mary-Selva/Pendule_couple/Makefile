FC=gfortran
FFLAGS=-O0 -g -Wall -ffpe-trap=invalid,overflow -fcheck=all -pedantic

EXEC=run
all:$(EXEC)

$(EXEC): pendule_couple_EE.o
	$(FC) -o $@ $^ $(FFLAGS)

%.o: %.f90
	$(FC) -c $< $(FFLAGS)

clean:
	rm -f *.o *.mod $(EXEC)
