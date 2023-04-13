FFLAGS = -Wall -fcheck=all -ffpe-trap=invalid -pedantic
all : exec
exec : constantes_modrk.o edo_modrk.o schema_modrk.o res_edork.o
	gfortran $(FFLAGS) -o exec constantes_modrk.o edo_modrk.o schema_modrk.o res_edork.o


constantes_modrk.o : constantes_modrk.f90
	gfortran $(FFLAGS) -c constantes_modrk.f90

edo_modrk.o : constantes_modrk.o edo_modrk.f90
	gfortran $(FFLAGS) -c edo_modrk.f90

schema_modrk.o : edo_modrk.o constantes_modrk.o schema_modrk.f90
	gfortran $(FFLAGS) -c  schema_modrk.f90

res_edork.o : constantes_modrk.o edo_modrk.o schema_modrk.o res_edork.f90
	gfortran $(FFLAGS) -c  res_edork.f90

clean :
	rm -f *.o *.mod exec
