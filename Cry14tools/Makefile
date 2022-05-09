#Makefile

FC = ifort
LD = ifort 


OBJS = cryapi_inp.o

Cry_out.exe: $(OBJS) 
	$(LD) $^ -o $@ 

%.o: %.f
	$(FC) -c $< -o $@
clean:
	rm -f *.exe *.o *.mod job-openmp.* *.dat

