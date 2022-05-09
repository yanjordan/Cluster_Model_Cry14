MODULE D3_Mod
	!$ use OMP_LIB
        USE DM_Mod
        USE ECD_Mod
        USE Clementi
	USE Clementi_DM
	USE CP_BS_Mod
	use Wigner_Mod
	use Moyal_Mod
        use Parity_Mod
	
CONTAINS
     
!read d3 file
!the format can be read in the Xtal14 user munual 
        SUBROUTINE cal_dm(d3_file)
                IMPLICIT NONE
                CHARACTER(len=*), INTENT(in)  :: d3_file
                
                CHARACTER(len=20)  :: c
                
                OPEN(60, file=trim(d3_file), action='read')
                
                READ(60, '(A)') c
               
                SELECT CASE (trim(c))
                        CASE ('DENSMAT')
                                CALL DENSMAT(trim(c))
                        CASE ('SDENSMAT')
                                CALL DENSMAT(trim(c))
                        CASE ('CDENSMAT')
                                CALL DENSMAT(trim(c))
                        CASE ('ECHG')
                                CALL ECHG(trim(c))
                        CASE ('SECHG')
                                CALL ECHG(trim(c))
                        CASE ('CECHG')
                                CALL ECHG(trim(c))
			CASE ('BIDIERD')
                                CALL CP(trim(c))
			CASE ('SBIDIERD')
                                CALL CP(trim(c))
                        CASE ('CBIDIERD')
                                CALL CP(trim(c))
			CASE ('WIGNER')
                                CALL WIGNER(trim(c))
			CASE ('SWIGNER')
                                CALL WIGNER(trim(c))
                        CASE ('CWIGNER')
                                CALL WIGNER(trim(c))
			CASE ('MOYAL')
                                CALL MOYAL(trim(c))
			CASE ('SMOYAL')
                                CALL MOYAL(trim(c))
                        CASE ('CMOYAL')
                                CALL MOYAL(trim(c))
                        CASE ('PARITY')
                                CALL PARITY(trim(c))
                        CASE ('SPARITY')
                                CALL PARITY(trim(c))
                        CASE ('CPARITY')
                                CALL PARITY(trim(c))
                        CASE ('USER')
                                CALL USER(trim(c))
                ENDSELECT

        ENDSUBROUTINE 

!routine PARITY function(separation orbital) 
        SUBROUTINE PARITY(key_word)
                IMPLICIT NONE
                
                CHARACTER(len=*), INTENT(in) :: key_word
                
                INTEGER ::  num_row, num_col, num_orbital, i, j, ihferm, num_orbital1
                
                INTEGER, DIMENSION(:), ALLOCATABLE :: list_orbital, list_orbital1
                
                DOUBLEPRECISION, DIMENSION(:,:), ALLOCATABLE :: Den,  Den_cle, Spin

                INTEGER, DIMENSION(3,4)::  atom_cell

                DOUBLEPRECISION, DIMENSION(3) :: A,B,C, v_dx, v_dy, grid_r
                
                DOUBLEPRECISION, DIMENSION(4) :: margins
                
                DOUBLEPRECISION :: lamda, dx, dy, cosxy
                
                CHARACTER(20) :: c_str
                
                LOGICAL :: b_log
                
!continue read the d3 file  with the format Xtal14  

                READ(60, *) c_str
                
                READ(60, *) num_row

                READ(60, *) c_str
                
                if (trim(c_str)=='ATOMS') THEN
                        READ(60,*) atom_cell(1,:)
                        READ(60,*) atom_cell(2,:)
                        READ(60,*) atom_cell(3,:)
                        
                        A(:)=atom_coord_list(atom_cell(1,1),1:3)+&
                                atom_cell(1,2)*lattice(1,1:3)+&
                                atom_cell(1,3)*lattice(2,1:3)+&
                                atom_cell(1,4)*lattice(3,1:3)
                        B(:)=atom_coord_list(atom_cell(2,1),1:3)+&
                                atom_cell(2,2)*lattice(1,1:3)+&
                                atom_cell(2,3)*lattice(2,1:3)+&
                                atom_cell(2,4)*lattice(3,1:3)
                        C(:)=atom_coord_list(atom_cell(3,1),1:3)+&
                                atom_cell(3,2)*lattice(1,1:3)+&
                                atom_cell(3,3)*lattice(2,1:3)+&
                                atom_cell(3,4)*lattice(3,1:3)
                        
                ELSEIF (trim(c_str)=='COORDINA') THEN
                        READ(60,*) A(:)
                        READ(60,*) B(:)
                        READ(60,*) C(:)                        
                ENDIF 

                b_log=.FALSE.
                do WHILE(.NOT. b_log)
                        READ(60,*) c_str
                        b_log=(trim(c_str)=='MARGINS')
                ENDDO
                
                READ(60,*) margins(1:4)
                
                margins=margins/0.52917706
                
! key world for separate orbital                 
                if (trim(key_word)=='SPARITY' ) THEN
                        b_log= .FALSE.
	                do WHILE(.NOT. b_log)
	                        READ(60,*) c_str
	                        b_log=(trim(c_str)=='ORBITAL')
	                ENDDO             
   
                        READ(60, *) num_orbital
                        
                        ALLOCATE(list_orbital(num_orbital))
                        
                        do i=1, num_orbital
                                READ(60, *) list_orbital(i)
                        ENDDO
                        
			WRITE(*,*) num_orbital, num_total_AOh
                        write(*,*) list_orbital
                elseif (trim(key_word)=='CPARITY' ) THEN
                        b_log= .FALSE.
	                do WHILE(.NOT. b_log)
	                        READ(60,*) c_str
	                        b_log=(trim(c_str)=='ORBITAL1')
	                ENDDO             
   
                        READ(60, *) num_orbital
                        
                        ALLOCATE(list_orbital(num_orbital))
                        
                        do i=1, num_orbital
                                READ(60, *) list_orbital(i)
                        ENDDO
                        
			WRITE(*,*) num_orbital, num_total_AOh
                        write(*,*) list_orbital
                        
                        b_log= .FALSE.
	                do WHILE(.NOT. b_log)
	                        READ(60,*) c_str
	                        b_log=(trim(c_str)=='ORBITAL2')
	                ENDDO             
   
                        READ(60, *) num_orbital1
                        
                        ALLOCATE(list_orbital1(num_orbital1))
                        
                        do i=1, num_orbital1
                                READ(60, *) list_orbital1(i)
                        ENDDO
                        
			WRITE(*,*) num_orbital1, num_total_AOh
			write(*,*) list_orbital1
                ENDIF                
                CLOSE(60)
                
                lamda=DOT_PRODUCT((A-B),(C-B))/dis(B,C)
                v_dy=(C-B)/dis(B,C)
 
                if (lamda >=0) THEN
                        A=A-lamda*v_dy
                else
                        B=B+lamda*v_dy
                ENDIF               
                
                v_dx=(A-B)/dis(B,A)
                
                A=A-margins(1)*v_dy+margins(3)*v_dx
                B=B-margins(1)*v_dy-margins(4)*v_dx
                C=C+margins(2)*v_dy-margins(4)*v_dx

                num_col=NINT(dis(C,B)/dis(B,A)*num_row)
                dy=dis(C,B)/(num_col-1d0)
                dx=dis(A,B)/(num_row-1d0) 

                cosxy=DOT_PRODUCT((A-B),(C-B))/(dis(B,C)*dis(A,B))

                ALLOCATE(Den(num_row, num_col), Spin(num_row, num_col))
                
                Spin=0.d0
                write(*,*) 'begin the PARITY computation'
	        if (trim(key_word)=='PARITY' ) THEN 
                        call Parity_cluster(B, num_row, num_col, dx*v_dx, dy*v_dy)
	        ELSEIF (trim(key_word)=='SPARITY' ) THEN
                        call SParity_cluster(B, num_row, num_col, dx*v_dx, dy*v_dy, num_orbital, list_orbital)
                        deallocate(list_orbital)
                ELSEIF (trim(key_word)=='CPARITY' ) THEN
                        call CParity_cluster(B, num_row, num_col, dx*v_dx, dy*v_dy, num_orbital, list_orbital, num_orbital1, list_orbital1)
                        deallocate(list_orbital, list_orbital1)
	        ENDIF            
        ENDSUBROUTINE

!routine Moyal function	
		subroutine MOYAL(key_word)
			IMPLICIT NONE
                        CHARACTER(len=*), INTENT(in) :: key_word
			doubleprecision :: d_s
			integer :: n_s, n_hkl, i, num_orbital, num_orbital1
			doubleprecision, dimension(3) :: dir_s
			integer, dimension(:,:), allocatable :: hkl_list
			
			INTEGER, DIMENSION(:), ALLOCATABLE :: list_orbital, list_orbital1
			
                        character(20) :: c_str
                        
                        LOGICAL :: b_log
			
			READ(60, *) n_s, n_hkl, d_s
			
			allocate(hkl_list(n_hkl,3))
			
			READ(60, *) c_str
			
			if(trim(c_str)=='DIR') then
				READ(60, *) dir_s(1:3)
			else
				write(*,*) 'Wrong Format. DIR information??'
			endif
			
			READ(60, *) c_str
			if(trim(c_str)=='LIST') then
				do i=1, n_hkl
					READ(60, *) hkl_list(i,1:3)
				enddo
			endif
			
			if (trim(key_word)=='MOYAL') then
                                CLOSE(60)
				call Moyal_Cluster(n_s, n_hkl, d_s, dir_s, hkl_list)
			elseif (trim(key_word)=='SMOYAL') then
				b_log= .FALSE.
                                do WHILE(.NOT. b_log)
                                        READ(60,*) c_str
                                        b_log=(trim(c_str)=='ORBITAL')
                                ENDDO             
        
                                READ(60, *) num_orbital
                                
                                ALLOCATE(list_orbital(num_orbital))
                                
                                do i=1, num_orbital
                                        READ(60, *) list_orbital(i)
                                ENDDO
                                CLOSE(60)

                                WRITE(*,*) num_orbital, num_total_AOh
                                write(*,*) list_orbital

				call SMoyal_Cluster(n_s, n_hkl, d_s, dir_s, hkl_list, num_orbital, list_orbital)
				
                                deallocate(list_orbital)
                        elseif (trim(key_word)=='CMOYAL') then
                                b_log= .FALSE.
                                do WHILE(.NOT. b_log)
                                        READ(60,*) c_str
                                        b_log=(trim(c_str)=='ORBITAL1')
                                ENDDO             
        
                                READ(60, *) num_orbital
                                
                                ALLOCATE(list_orbital(num_orbital))
                                
                                do i=1, num_orbital
                                        READ(60, *) list_orbital(i)
                                ENDDO

                                WRITE(*,*) num_orbital, num_total_AOh
                                write(*,*) list_orbital

                                b_log= .FALSE.
                                do WHILE(.NOT. b_log)
                                        READ(60,*) c_str
                                        b_log=(trim(c_str)=='ORBITAL2')
                                ENDDO             
        
                                READ(60, *) num_orbital1
                                
                                ALLOCATE(list_orbital1(num_orbital1))
                                
                                do i=1, num_orbital1
                                        READ(60, *) list_orbital1(i)
                                ENDDO
                                CLOSE(60)

                                WRITE(*,*) num_orbital1, num_total_AOh
                                write(*,*) list_orbital1
                                
				call CMoyal_Cluster(n_s, n_hkl, d_s, dir_s, hkl_list, num_orbital, list_orbital, num_orbital1, list_orbital1)
				
                                deallocate(list_orbital, list_orbital1)
			endif
			
			deallocate(hkl_list)
		endsubroutine		
		
!routine Wigner function	
		subroutine WIGNER(key_word)
			IMPLICIT NONE
                        CHARACTER(len=*), INTENT(in) :: key_word
			doubleprecision :: d_row, d_col
			integer :: n_row, n_col, num_orbital, i, num_orbital1
			doubleprecision, dimension(6) :: xp1, xp2
			
			INTEGER, DIMENSION(:), ALLOCATABLE :: list_orbital, list_orbital1
			
                        character(20) :: c_str
                        
                        LOGICAL :: b_log
			
			READ(60, *) n_row, n_col, d_row, d_col
					
			READ(60, *) c_str
			
			if(trim(c_str)=='DIR') then
				READ(60, *) xp1(1:6)
				READ(60, *) xp2(1:6)
			else
				write(*,*) 'Wrong Format. DIR information??'
			endif
			write(*,*) 'begin the Wigner computation'
			if (trim(key_word)=='WIGNER') then
                                CLOSE(60)                                
				call Wigner_Cluster(n_row, n_col, d_row, d_col, xp1, xp2)
			elseif (trim(key_word)=='SWIGNER') then
				b_log= .FALSE.
                                do WHILE(.NOT. b_log)
                                        READ(60,*) c_str
                                        b_log=(trim(c_str)=='ORBITAL')
                                ENDDO             
        
                                READ(60, *) num_orbital
                                
                                ALLOCATE(list_orbital(num_orbital))
                                
                                do i=1, num_orbital
                                        READ(60, *) list_orbital(i)
                                ENDDO
                                CLOSE(60)

                                WRITE(*,*) num_orbital, num_total_AOh
                                write(*,*) list_orbital

				call SWigner_Cluster(n_row, n_col, d_row, d_col, xp1, xp2, num_orbital, list_orbital)
				
                                deallocate(list_orbital)
                        elseif (trim(key_word)=='CWIGNER') then
                                b_log= .FALSE.
                                do WHILE(.NOT. b_log)
                                        READ(60,*) c_str
                                        b_log=(trim(c_str)=='ORBITAL1')
                                ENDDO             
        
                                READ(60, *) num_orbital
                                
                                ALLOCATE(list_orbital(num_orbital))
                                
                                do i=1, num_orbital
                                        READ(60, *) list_orbital(i)
                                ENDDO

                                WRITE(*,*) num_orbital, num_total_AOh
                                write(*,*) list_orbital

                                b_log= .FALSE.
                                do WHILE(.NOT. b_log)
                                        READ(60,*) c_str
                                        b_log=(trim(c_str)=='ORBITAL2')
                                ENDDO             
        
                                READ(60, *) num_orbital1
                                
                                ALLOCATE(list_orbital1(num_orbital1))
                                
                                do i=1, num_orbital1
                                        READ(60, *) list_orbital1(i)
                                ENDDO
                                CLOSE(60)

                                WRITE(*,*) num_orbital1, num_total_AOh
                                write(*,*) list_orbital1

                                call CWigner_Cluster(n_row, n_col, d_row, d_col, xp1, xp2, num_orbital, list_orbital, num_orbital1, list_orbital1)
                                deallocate(list_orbital, list_orbital1)
			endif
		
		endsubroutine

!routine Compton profile
		SUBROUTINE CP(key_word)
			IMPLICIT NONE
                
                        CHARACTER(len=*), INTENT(in) :: key_word
                
                        INTEGER ::  Ndir, Npu, Imodo, Icaso, Npoip, i, num_orbital, num_orbital1
				
			INTEGER, DIMENSION(:), ALLOCATABLE :: list_orbital, list_orbital1
                
                        DOUBLEPRECISION :: step, stepc
                
                        DOUBLEPRECISION, DIMENSION(:,:), ALLOCATABLE :: dirs
				
			integer, DIMENSION(:,:), ALLOCATABLE :: atomes
                
                        CHARACTER(20) :: c_str
				
				LOGICAL :: b_log
                				
				READ(60, *) Ndir, Npu, step, Imodo, Icaso
				
				allocate(dirs(Ndir,3))
				if (Imodo/=0d0) then
					allocate(atomes(Ndir,4))
				endif
				
				READ(60, *) c_str
				
				if(trim(c_str)=='PROF') then
					READ(60, *) Npoip, stepc
					READ(60, *) c_str
					if(trim(c_str)=='DIR') then
						if (Imodo==0d0) then
							do i=1, Ndir
								READ(60, *) dirs(i,:)
							enddo
						else
							do i=1, Ndir
								READ(60, *) atomes(i,:)
								if (atomes(i,1)>num_atom_cell) then
									write(*,*) 'Wrong definition of atom label'
									stop
								else
									dirs(i,:)=atom_coord_list(atomes(i,1),:)+matmul(atomes(i,2:4), lattice)
								endif
							enddo						
                                                endif
                                                
						if (key_word=='SBIDIERD') then
							b_log= .FALSE.
							do WHILE(.NOT. b_log)
									READ(60,*) c_str
									b_log=(trim(c_str)=='ORBITAL')
							ENDDO             
		   
								READ(60, *) num_orbital
								
								ALLOCATE(list_orbital(num_orbital))
								
								do i=1, num_orbital
										READ(60, *) list_orbital(i)
								ENDDO
							close(60)	
							WRITE(*,*) num_orbital, num_total_AOh
							write(*,*) list_orbital	
							
							if (Icaso==1d0) then
								call SBIDIERD_Cluster(Ndir, Npu, step, Npoip, stepc, dirs, num_orbital, list_orbital)
							else
								write(*,*) 'The ICASO should only be 1:'
								write(*,*) '	1 represente the cluster computation'
								write(*,*) '	no other choice for SBIDIERD'
								stop
                                                        endif
                                                elseif (key_word=='CBIDIERD') then
							b_log= .FALSE.
							do WHILE(.NOT. b_log)
									READ(60,*) c_str
									b_log=(trim(c_str)=='ORBITAL1')
							ENDDO             
		   
								READ(60, *) num_orbital
								
								ALLOCATE(list_orbital(num_orbital))
								
								do i=1, num_orbital
										READ(60, *) list_orbital(i)
								ENDDO
								
							WRITE(*,*) num_orbital, num_total_AOh
                                                        write(*,*) list_orbital	
                                                        
                                                        b_log= .FALSE.
							do WHILE(.NOT. b_log)
									READ(60,*) c_str
									b_log=(trim(c_str)=='ORBITAL2')
							ENDDO             
		   
								READ(60, *) num_orbital1
								
								ALLOCATE(list_orbital1(num_orbital1))
								
								do i=1, num_orbital1
										READ(60, *) list_orbital1(i)
								ENDDO
							close(60)	
							WRITE(*,*) num_orbital1, num_total_AOh
							write(*,*) list_orbital1	
							
							if (Icaso==1d0) then
								call CBIDIERD_Cluster(Ndir, Npu, step, Npoip, stepc, dirs, num_orbital, list_orbital, num_orbital1, list_orbital1)
							else
								write(*,*) 'The ICASO should only be 1:'
								write(*,*) '	1 represente the cluster computation'
								write(*,*) '	no other choice for SBIDIERD'
								stop
							endif
						else 
							if (Icaso==1d0) then
								call BIDIERD_Cluster(Ndir, Npu, step, Npoip, stepc, dirs)
							elseif(Icaso==2d0) then
								call BIDIERD(Ndir, Npu, step, Npoip, stepc, dirs)
							else
								write(*,*) 'The ICASO should only be 1 or 2:'
								write(*,*) '	1 represente the cluster computation'
								write(*,*) '	2 represente the center cell computation'
								stop
							endif
							close(60)
						endif
					else
						write(*,*) 'format wrong of d3 file to define the DIR' 
						stop					
					endif
				else
					write(*,*) 'format wrong of d3 file to define the PROF' 
					stop
				endif
				
		ENDSUBROUTINE

!routine density matrix (separation orbital)        
        SUBROUTINE DENSMAT(key_word)
                IMPLICIT NONE
                
                DOUBLEPRECISION, DIMENSION(:,:), ALLOCATABLE :: knots_coord, grid, DM_t, DM_s, DM_Cle
                
                DOUBLEPRECISION, DIMENSION(3) :: displace
                
                INTEGER, DIMENSION(:), ALLOCATABLE :: atom_num, list_orbital, list_orbital1
                
                INTEGER, DIMENSION(:,:), ALLOCATABLE ::  atom_cell
                
                DOUBLEPRECISION :: lpr,dl
                
                INTEGER :: nkn, npu, iplot, imodo, i, j, IHFERM, num_orbital, num_orbital1
                
                CHARACTER(len=*), INTENT(in) :: key_word
                
                LOGICAL :: b_log
                
                CHARACTER(len=20) :: c_str
                
!continue read the file                
                
                READ(60, *) nkn, npu, iplot, imodo, lpr
                
!                write(*, *) nkn, npu, iplot, imodo, lpr
                
                ALLOCATE(knots_coord(nkn,3))
                
                if (imodo==0) THEN
                        do i=1, nkn
                                READ(60, *) knots_coord(i,1:3)
                        ENDDO
                ELSE
                        ALLOCATE(atom_num(nkn), atom_cell(nkn,3))
                        READ(60,*) displace(1:3)
                        
                        do i=1, nkn
                                READ(60, *) atom_num(i), atom_cell(i,1:3)
                                knots_coord(i,1:3)=atom_coord_list(atom_num(i),1:3)+atom_cell(i,1)*lattice(1,1:3)&
                                                                +atom_cell(i,2)*lattice(2,1:3)&
                                                                +atom_cell(i,3)*lattice(3,1:3)
                        ENDDO
						

			knots_coord(1,1:3)=knots_coord(1,1:3)+lpr*(knots_coord(1,1:3)-knots_coord(2,1:3))
			knots_coord(nkn,1:3)=knots_coord(nkn,1:3)+lpr*(knots_coord(nkn,1:3)-knots_coord(nkn-1,1:3))
	
                        
                        if (.NOT. all(displace(1:3)==(/0d0,0d0,0d0/))) THEN
                                do i=1,nkn
                                        knots_coord(i, 1:3)=knots_coord(i, 1:3)+displace(1:3)
                                ENDDO
                        ENDIF
                        
                ENDIF
                
                if (trim(key_word)=='SDENSMAT' ) THEN
                        b_log= .FALSE.
	                do WHILE(.NOT. b_log)
	                        READ(60,*) c_str
	                        b_log=(trim(c_str)=='ORBITAL')
	                ENDDO             
   
                        READ(60, *) num_orbital
                        
                        ALLOCATE(list_orbital(num_orbital))
                        
                        do i=1, num_orbital
                                READ(60, *) list_orbital(i)
                        ENDDO
                        
                        WRITE(*,*) num_orbital, num_total_AOh
                        write(*,*) list_orbital
                elseif (trim(key_word)=='CDENSMAT' ) THEN
                        b_log= .FALSE.
                        do WHILE(.NOT. b_log)
                                READ(60,*) c_str
                                b_log=(trim(c_str)=='ORBITAL1')
                        ENDDO             
        
                        READ(60, *) num_orbital
                        
                        ALLOCATE(list_orbital(num_orbital))
                        
                        do i=1, num_orbital
                                READ(60, *) list_orbital(i)
                        ENDDO
                        
                        WRITE(*,*) num_orbital, num_total_AOh
                        write(*,*) list_orbital

                        b_log= .FALSE.
                        do WHILE(.NOT. b_log)
                                READ(60,*) c_str
                                b_log=(trim(c_str)=='ORBITAL2')
                        ENDDO             
        
                        READ(60, *) num_orbital1
                        
                        ALLOCATE(list_orbital1(num_orbital1))
                        
                        do i=1, num_orbital1
                                READ(60, *) list_orbital1(i)
                        ENDDO
                        
                        WRITE(*,*) num_orbital1, num_total_AOh
                        write(*,*) list_orbital1
                ENDIF
                
                CLOSE(60)
                
!get the grid                
                ALLOCATE(grid(npu,3), DM_t(npu,npu), DM_Cle(npu,npu), DM_s(npu,npu))
                
                DM_s=0.d0
				DM_Cle=0.d0
                
                IHFERM=0d0
                
                CALL grid_route(nkn, npu, knots_coord, grid, dl)
                
! get the DM value

                if (trim(key_word)=='DENSMAT' ) THEN
	                call DM(npu, grid, 'R', DM_t)
!$OMP PARALLEL
!$OMP DO SCHEDULE(RUNTIME) 					
					do i=1, npu
						do j=1,i
							DM_Cle(i,j)=Cle_DM(grid(i,:),grid(j,:))
							DM_Cle(j,i)=DM_Cle(i,j)
						enddo
					enddo
!$OMP END DO
!$OMP END PARALLEL				
	                if (uorr=='U' ) THEN
	                        IHFERM=1d0
	                        call DM(npu, grid, 'U',DM_s)
	                ENDIF                
                ELSEIF (trim(key_word)=='SDENSMAT' ) THEN
	                call SDM(npu, grid, 'R' , num_orbital, list_orbital,DM_t)
			
	                if (uorr=='U' ) THEN
	                        IHFERM=1d0
	                        call SDM(npu, grid, 'U', num_orbital, list_orbital, DM_s)
                        ENDIF
                ELSEIF (trim(key_word)=='CDENSMAT' ) THEN
                        if (uorr=='R' ) THEN
                                call CDM(npu, grid, num_orbital, list_orbital, num_orbital1, list_orbital1, DM_t)
                        elseif (uorr=='U' ) THEN
                                IHFERM=1d0
                                call CDM(npu, grid, num_orbital, list_orbital, num_orbital1, list_orbital1, DM_t, DM_s)
                        ENDIF
                ENDIF
                

!print the result with the format of fort.25 (see the Xtal14 user manual)

                !change the unit to angstrom
                dl=dl*0.529177208
                atom_coord_list=atom_coord_list*0.529177208
                lattice=lattice*0.529177208
                
                OPEN(80, file='DM_CRY.data', status='replace')
                
                WRITE(80,'(A3, I1, A4, 2I5, 3E12.5)') '-%-', IHFERM, 'MAPN', npu, npu, dl, dl, 0.d0
                WRITE(80,'(1P, 6E12.5)') 0d0, (npu-1d0)*dl, 0.d0, 0.d0, 0.d0, 0.d0
                WRITE(80,'(1P, 3E12.5, 4X, 2I4)') (npu-1d0)*dl, 0.d0, 0.d0, num_atom_cell, 3
                WRITE(80,'(1P, 6E12.5)') ((DM_t(i,j), i=1, npu),j=1,npu)
                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(80,'(A3, I1, A4, 2I5, 3E12.5)') '-%-', IHFERM, 'MAPN', npu, npu, dl, dl, 0.d0
                WRITE(80,'(1P, 6E12.5)') 0d0, (npu-1d0)*dl, 0.d0, 0.d0, 0.d0, 0.d0
                WRITE(80,'(1P, 3E12.5, 4X, 2I4)') (npu-1d0)*dl, 0.d0, 0.d0, num_atom_cell, 3
                WRITE(80,'(1P, 6E12.5)') ((DM_s(i,j), i=1,npu),j=1,npu)
                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO               
                
                CLOSE(80)
	if (trim(key_word)=='DENSMAT' ) THEN	
		OPEN(90, file='DM_CRY_PATO.data', status='replace')
                
                WRITE(90,'(A3, I1, A4, 2I5, 3E12.5)') '-%-', IHFERM, 'MAPN', npu, npu, dl, dl, 0.d0
                WRITE(90,'(1P, 6E12.5)') 0.d0, (npu-1d0)*dl, 0.d0, 0.d0, 0.d0, 0.d0
                WRITE(90,'(1P, 3E12.5, 4X, 2I4)') (npu-1d0)*dl, 0.d0, 0.d0, num_atom_cell, 3
                WRITE(90,'(1P, 6E12.5)') ((DM_Cle(i,j), i=1, npu),j=1,npu)
                
                do i=1, num_atom_cell
                        WRITE(90, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(90,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(90,'(A3, I1, A4, 2I5, 3E12.5)') '-%-', IHFERM, 'MAPN', npu, npu, dl, dl, 0d0
                WRITE(90,'(1P, 6E12.5)') 0d0, (npu-1d0)*dl, 0d0, 0d0, 0d0, 0d0
                WRITE(90,'(1P, 3E12.5, 4X, 2I4)') (npu-1d0)*dl, 0d0, 0d0, num_atom_cell, 3
                WRITE(90,'(1P, 6E12.5)') ((0.d0, i=1,npu),j=1,npu)
                
                do i=1, num_atom_cell
                        WRITE(90, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(90,'(3E20.12)') lattice(i,:)
                ENDDO               
                
                CLOSE(90)

		OPEN(100, file='DM_CRY_Defo.data', status='replace')
                
                WRITE(100,'(A3, I1, A4, 2I5, 3E12.5)') '-%-', IHFERM, 'MAPN', npu, npu, dl, dl, 0.d0
                WRITE(100,'(1P, 6E12.5)') 0.d0, (npu-1d0)*dl, 0.d0, 0.d0, 0.d0, 0.d0
                WRITE(100,'(1P, 3E12.5, 4X, 2I4)') (npu-1d0)*dl, 0.d0, 0.d0, num_atom_cell, 3
                WRITE(100,'(1P, 6E12.5)') ((DM_t(i,j)-DM_Cle(i,j), i=1, npu),j=1,npu)
                
                do i=1, num_atom_cell
                        WRITE(100, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(100,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(100,'(A3, I1, A4, 2I5, 3E12.5)') '-%-', IHFERM, 'MAPN', npu, npu, dl, dl, 0d0
                WRITE(100,'(1P, 6E12.5)') 0d0, (npu-1d0)*dl, 0d0, 0d0, 0d0, 0d0
                WRITE(100,'(1P, 3E12.5, 4X, 2I4)') (npu-1d0)*dl, 0d0, 0d0, num_atom_cell, 3
                WRITE(100,'(1P, 6E12.5)') ((DM_s(i,j), i=1,npu),j=1,npu)
                
                do i=1, num_atom_cell
                        WRITE(100, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(100,'(3E20.12)') lattice(i,:)
                ENDDO               
                
                CLOSE(100)
            endif
			
			DEALLOCATE(knots_coord, grid, DM_t, DM_Cle, DM_s)
				
			if (imodo/=0) THEN
				DEALLOCATE(atom_num, atom_cell)
			endif

			if (trim(key_word)=='SDENSMAT' ) THEN
                                DEALLOCATE(list_orbital)
                        elseif (trim(key_word)=='CDENSMAT' ) THEN
				DEALLOCATE(list_orbital, list_orbital1)
			endif

        ENDSUBROUTINE
 
!electron density (separation orbital)  
        SUBROUTINE ECHG(key_word)
                IMPLICIT NONE
                
                CHARACTER(len=*), INTENT(in) :: key_word
                
                INTEGER ::  num_row, num_col, num_orbital, i, j, ihferm, num_orbital1
                
                INTEGER, DIMENSION(:), ALLOCATABLE :: list_orbital, list_orbital1
                
                DOUBLEPRECISION, DIMENSION(:,:), ALLOCATABLE :: Den,  Den_cle, Spin

                INTEGER, DIMENSION(3,4)::  atom_cell

                DOUBLEPRECISION, DIMENSION(3) :: A,B,C, v_dx, v_dy, grid_r
                
                DOUBLEPRECISION, DIMENSION(4) :: margins
                
                DOUBLEPRECISION :: lamda, dx, dy, cosxy
                
                CHARACTER(20) :: c_str
                
                LOGICAL :: b_log
                
!continue read the d3 file  with the format Xtal14  

                READ(60, *) c_str
                
                READ(60, *) num_row

                READ(60, *) c_str
                
                if (trim(c_str)=='ATOMS') THEN
                        READ(60,*) atom_cell(1,:)
                        READ(60,*) atom_cell(2,:)
                        READ(60,*) atom_cell(3,:)
                        
                        A(:)=atom_coord_list(atom_cell(1,1),1:3)+&
                                atom_cell(1,2)*lattice(1,1:3)+&
                                atom_cell(1,3)*lattice(2,1:3)+&
                                atom_cell(1,4)*lattice(3,1:3)
                        B(:)=atom_coord_list(atom_cell(2,1),1:3)+&
                                atom_cell(2,2)*lattice(1,1:3)+&
                                atom_cell(2,3)*lattice(2,1:3)+&
                                atom_cell(2,4)*lattice(3,1:3)
                        C(:)=atom_coord_list(atom_cell(3,1),1:3)+&
                                atom_cell(3,2)*lattice(1,1:3)+&
                                atom_cell(3,3)*lattice(2,1:3)+&
                                atom_cell(3,4)*lattice(3,1:3)
                        
                ELSEIF (trim(c_str)=='COORDINA') THEN
                        READ(60,*) A(:)
                        READ(60,*) B(:)
                        READ(60,*) C(:)                        
                ENDIF 

                b_log=.FALSE.
                do WHILE(.NOT. b_log)
                        READ(60,*) c_str
                        b_log=(trim(c_str)=='MARGINS')
                ENDDO
                
                READ(60,*) margins(1:4)
                
                margins=margins/0.52917706
                
! key world for separate orbital                 
                if (trim(key_word)=='SECHG' ) THEN
                        b_log= .FALSE.
	                do WHILE(.NOT. b_log)
	                        READ(60,*) c_str
	                        b_log=(trim(c_str)=='ORBITAL')
	                ENDDO             
   
                        READ(60, *) num_orbital
                        
                        ALLOCATE(list_orbital(num_orbital))
                        
                        do i=1, num_orbital
                                READ(60, *) list_orbital(i)
                        ENDDO
                        
			WRITE(*,*) num_orbital, num_total_AOh
                        write(*,*) list_orbital
                elseif (trim(key_word)=='CECHG' ) THEN
                        b_log= .FALSE.
                        do WHILE(.NOT. b_log)
                                READ(60,*) c_str
                                b_log=(trim(c_str)=='ORBITAL1')
                        ENDDO             
        
                        READ(60, *) num_orbital
                        
                        ALLOCATE(list_orbital(num_orbital))
                        
                        do i=1, num_orbital
                                READ(60, *) list_orbital(i)
                        ENDDO
                        
                        WRITE(*,*) num_orbital, num_total_AOh
                        write(*,*) list_orbital

                        b_log= .FALSE.
                        do WHILE(.NOT. b_log)
                                READ(60,*) c_str
                                b_log=(trim(c_str)=='ORBITAL2')
                        ENDDO             
        
                        READ(60, *) num_orbital1
                        
                        ALLOCATE(list_orbital1(num_orbital1))
                        
                        do i=1, num_orbital1
                                READ(60, *) list_orbital1(i)
                        ENDDO
                        
                        WRITE(*,*) num_orbital1, num_total_AOh
                        write(*,*) list_orbital1
                ENDIF                
                CLOSE(60)
                
                lamda=DOT_PRODUCT((A-B),(C-B))/dis(B,C)
                v_dy=(C-B)/dis(B,C)
 
                if (lamda >=0) THEN
                        A=A-lamda*v_dy
                else
                        B=B+lamda*v_dy
                ENDIF               
                
                v_dx=(A-B)/dis(B,A)
                
                A=A-margins(1)*v_dy+margins(3)*v_dx
                B=B-margins(1)*v_dy-margins(4)*v_dx
                C=C+margins(2)*v_dy-margins(4)*v_dx

                num_col=NINT(dis(C,B)/dis(B,A)*num_row)
                dy=dis(C,B)/(num_col-1d0)
                dx=dis(A,B)/(num_row-1d0) 

                cosxy=DOT_PRODUCT((A-B),(C-B))/(dis(B,C)*dis(A,B))

                ALLOCATE(Den(num_row, num_col), Den_cle(num_row, num_col), Spin(num_row, num_col))
                
                Spin=0.d0
                write(*,*) 'begin the ECHG computation'
                if (trim(key_word)=='ECHG' ) THEN 
                        call ECHG_cluster(B, num_row, num_col, dx*v_dx, dy*v_dy)
	        ELSEIF (trim(key_word)=='SECHG' ) THEN
                        call SECHG_cluster(B, num_row, num_col, dx*v_dx, dy*v_dy, num_orbital, list_orbital)
                        deallocate(list_orbital)
                ELSEIF (trim(key_word)=='CECHG' ) THEN
                        call CECHG_cluster(B, num_row, num_col, dx*v_dx, dy*v_dy, num_orbital, list_orbital, num_orbital1, list_orbital1)
                        deallocate(list_orbital, list_orbital1)
	        ENDIF  
        ENDSUBROUTINE


!get the grid of the route of DM by the coordinate
        SUBROUTINE grid_route(nkn, npu, knots_coord, grid_rout, dl)
                IMPLICIT NONE
        
                INTEGER, INTENT(in) :: nkn, npu
                DOUBLEPRECISION, DIMENSION(nkn,3) :: knots_coord
                
                DOUBLEPRECISION, DIMENSION(npu, 3) :: grid_rout
                
                DOUBLEPRECISION :: dl
                
                INTEGER :: i,j, num_seg, pos_num
                DOUBLEPRECISION :: longeur, res
                DOUBLEPRECISION, DIMENSION(3) :: vec
                DOUBLEPRECISION, DIMENSION(nkn-1) :: l_seg
                
                longeur=0d0
                l_seg=0d0
                
                do i=1,nkn-1
                        l_seg(i)=dis(knots_coord(i,1:3), knots_coord(i+1,1:3))
                        longeur=longeur+l_seg(i)
                ENDDO
                
!                WRITE(*,*) longeur
                grid_rout=0d0
                
                dl=longeur/(npu-1)
                
                grid_rout(1,1:3)=knots_coord(1,1:3)

                pos_num=1d0
                
                res=0d0
!                vec(1:3)=(knots_coord(2,1:3)-knots_coord(1,1:3))/l_seg(1)
                
                do i=1, nkn-1
                        num_seg=(l_seg(i)-res)/dl
                        vec(1:3)=(knots_coord(i+1,1:3)-knots_coord(i,1:3))/l_seg(i)
                        
                        do j=1,num_seg
                                grid_rout(pos_num+j,1:3)=knots_coord(i,1:3)+(dl*j+res)*vec(1:3)
                        ENDDO
                        pos_num = pos_num+ num_seg
                        if (i<nkn-1) THEN
	                        res=dl-(l_seg(i)-res-num_seg*dl)
	                        vec(1:3)=(knots_coord(i+2,1:3)-knots_coord(i+1,1:3))/l_seg(i+1)
	                        pos_num =  pos_num+ 1
	                        grid_rout(pos_num,1:3)=knots_coord(i+1,1:3)+res*vec(1:3)
                        ENDIF

                ENDDO

        ENDSUBROUTINE        
   
!User definition function for some test
        subroutine USER(key_word)
                IMPLICIT NONE

                CHARACTER(len=*), INTENT(in) :: key_word
                
                INTEGER ::  num_row, num_col, num_orbital, i, j, ihferm
                
                INTEGER, DIMENSION(:), ALLOCATABLE :: list_orbital
                
                DOUBLEPRECISION, DIMENSION(:,:), ALLOCATABLE :: Den,  Den_cle, Spin

                INTEGER, DIMENSION(3,4)::  atom_cell

                DOUBLEPRECISION, DIMENSION(3) :: A,B,C, v_dx, v_dy, grid_r
                
                DOUBLEPRECISION, DIMENSION(4) :: margins
                
                DOUBLEPRECISION :: lamda, dx, dy, cosxy
                
                read(60,*) num_row, num_col
                read(60,*) dx, dy
                read(60,*) A(:)
                read(60,*) B(:)
                read(60,*) C(:)
		A=A/0.529177208
		B=B/0.529177208
		C=C/0.529177208
		dx=dx/0.529177208
		dy=dy/0.529177208

                v_dx=(A-B)/dis(B,A)
                v_dy=(C-B)/dis(B,C)

                call ECHG_cluster(B, num_row, num_col, dx*v_dx, dy*v_dy)


                close(60)



        endsubroutine    


ENDMODULE D3_Mod