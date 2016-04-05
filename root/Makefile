# compiler 
F90 = ifort 
FFLAGS = -O2 -i4 -r8 -warn -CB -traceback -vec-report1 -nogen-interfaces

SRC = solve.f90 main.f90
OBJ = $(SRC:.f90=.o)
MOD = $(SRC:.f90=.mod)

EXE = root.x

# suffix rules 
%.o: %.f90
	$(F90) $(FFLAGS) -c $<

# make rules 
root: $(OBJ)
	$(F90) $(FFLAGS) -o $(EXE) $(OBJ)

clean: 
	-rm -f $(OBJ) $(MOD)
