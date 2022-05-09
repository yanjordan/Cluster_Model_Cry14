PROGRAM main

        USE D3_Mod

	USE EMD_Mod
        
        CHARACTER(len=150) :: out_file, d3_file
        
        CHARACTER(len=1) :: ur
        
        LOGICAL :: alive
        
        INTEGER :: narg
        
        REAL :: start, finish

	integer :: h_time, min_time, s_time
	
        CALL CPU_TIME(start)
        
        WRITE(*,*) 'This is a program to calculate the Density Matrix for Crystal14.'
        WRITE(*,*) 'We read the out put file of Cry_out program'
        WRITE(*,*) 'The command is like: xxx.exe Cry_out.out d3_file (U or R optional, U if no input)'
        WRITE(*,*) '***************************************************************'
        WRITE(*,*) ' '
        
        narg=iargc()
        
        if(narg==2) THEN
                CALL getarg(1, out_file)
                CALL getarg(2, d3_file)
                
                INQUIRE(file=trim(out_file), exist=alive)
                if(.NOT.alive) THEN
                        WRITE(*,*) 'Hey man, Can not find the output file?' 
                        STOP
                ENDIF
                
                INQUIRE(file=trim(d3_file), exist=alive)
                if(.NOT.alive) THEN
                        WRITE(*,*) 'Hey man, Can not find the d3 file?' 
                        STOP
                ENDIF
                
                CALL Read_File(trim(out_file), 'U')
               !CALL WRITE_para('U')
		
	        CALL normalisation()
	        STOP
	        CALL cal_dm(trim(d3_file))

        ELSEIF(narg==3) THEN
                CALL getarg(1, out_file)
                CALL getarg(2, d3_file)
                CALL getarg(3, ur)
                
                INQUIRE(file=trim(out_file), exist=alive)
                if(.NOT.alive) THEN
                        WRITE(*,*) 'Hey man, Can not find the output file?' 
                        STOP
                ENDIF
                
                INQUIRE(file=trim(d3_file), exist=alive)
                if(.NOT.alive) THEN
                        WRITE(*,*) 'Hey man, Can not find the d3 file?' 
                        STOP
                ENDIF
                
                if (ur=='U') THEN
                        CALL Read_File(trim(out_file), 'U')
                ELSEIF(ur=='R') THEN
                        CALL Read_File(trim(out_file), 'R')
                ELSE
                        WRITE(*,*) 'illegal input: ', ur
                        WRITE(*,*) 'Only two choices: U(Unrestricted) or R(Restricted)'
                        STOP
                ENDIF
                
                !CALL WRITE_para('U')
	

	        CALL normalisation()
		!call test_B0()

	        CALL cal_dm(trim(d3_file))
	
		!call test_p()
		!call test_r()
        ELSE	        
                WRITE(*,*) 'Wrong input, just outfile + d3 file + U or R optional'
                STOP
                
        ENDIF


        CALL CPU_TIME(finish)

	h_time=floor((finish-start)/3600)

	min_time=floor((finish-start)/60)-h_time*60d0

	s_time=floor(finish-start)-h_time*3600d0-min_time*60d0
        
        WRITE(*,'("CPU Time = ", I4," h", I4," min", I4," s" )') h_time, min_time, s_time

ENDPROGRAM main
