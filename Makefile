#Makefile

FC = ifort
LD = ifort 

FCFLAGS = -openmp
LDFLAGS = -openmp

OBJS = Fbase_Mod.o Reading_Mod.o Clementi.o Clementi_DM.o DM_Mod.o ECD_Mod.o EMD_Mod.o Bs_Mod.o CP_BS_Mod.o Wigner_Mod.o Moyal_Mod.o Parity_Mod.o D3_Mod.o main.o

Cry_DM.exe: $(OBJS) 
	$(LD) $(LDFLAGS) $^ -o $@ 

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@
clean:
	rm -f *.exe *.o *.mod *.dat
