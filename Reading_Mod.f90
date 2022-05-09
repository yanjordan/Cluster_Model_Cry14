MODULE Reading_Mod

        USE Fbase_Mod
        
        
        IMPLICIT NONE
        
	Type(GTOs), dimension(:), allocatable :: cart_gauss
		                    
        DOUBLEPRECISION, DIMENSION(:,:), allocatable :: mat
        
        DOUBLEPRECISION, DIMENSION(:,:,:) , ALLOCATABLE :: Pop_tot, Pop_spin

	DOUBLEPRECISION, DIMENSION(:,:,:) , ALLOCATABLE :: Pop_tot_car, Pop_spin_car
                
        DOUBLEPRECISION, DIMENSION(3,3) :: lattice
		
		DOUBLEPRECISION, DIMENSION(3,3) :: lattice_rec
		
	DOUBLEPRECISION, DIMENSION(:,:), ALLOCATABLE :: coord_frac
        
        INTEGER :: num_cell
                
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: cell_dir
        
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: mat_cell
        
        INTEGER, DIMENSION(200) :: atome_num_list
        
        CHARACTER(len=2) , DIMENSION(200) :: atome_ele_list
        
        DOUBLEPRECISION, DIMENSION(200,3) :: atom_coord_list
        
        CHARACTER(len=1) :: uorr
        
!type_line: 0--atome, 1--type shell, 2--type cofficient                
        INTEGER :: num_total_AOc, num_total_AOh, num_atom_cell
        
      
CONTAINS




        SUBROUTINE Read_File(filename, type_cal)
                IMPLICIT NONE
                1080  FORMAT(I4,1X,A,3F7.3)
		1090  FORMAT(31X,I4,A)
		1100  FORMAT(26X,I4,'-',I4,A)
		1120  FORMAT(36X,1P,E14.7,3E10.3)
		!1120  FORMAT(40X,1P,4E10.3)
		1071  FORMAT(1X,79('*'))
		1200  FORMAT(I4,3X,1P,10E16.8)
		1082  FORMAT(I4,2X,I4,1X,A,1P,3E20.12)
                
		DOUBLEPRECISION, DIMENSION(:,:,:), allocatable :: coeff
		INTEGER, DIMENSION(:), allocatable :: num_GTOs
		DOUBLEPRECISION, DIMENSION(:,:), allocatable :: mat_temps
				
                CHARACTER(len=*), INTENT(in) :: filename
                
                CHARACTER(len=1),  INTENT(in) :: type_cal
                
                DOUBLEPRECISION, DIMENSION(3) :: atom_coord
				
				DOUBLEPRECISION, DIMENSION(3) :: cross_vec_temps
                
                INTEGER, DIMENSION(3) :: vec_dir
                
                INTEGER :: i, j, k,l,m, num_line, num_space,num_atom_AO, num_GTO, num_atom, num_temp, num1,num2
				
                INTEGER :: num_coeff,  order_AO1,order_AO2, type_line, num_atom_AOh
                
                DOUBLEPRECISION :: exp_c, s_c, p_c, dfg_c
                
                CHARACTER(len=3) :: sym_ele, str_ele, type_AO
                
                CHARACTER(len=80) :: c1, symbol  

		logical :: b              
                
                b= .FALSE.

        	allocate(coeff(1000,10,10), num_GTOs(1000), mat_temps(1000,1000))

		uorr=type_cal
        
	        mat_temps=0.d0
	        coeff=0.d0
	        num_total_AOh=0d0
	        num_total_AOc=0d0
	        Pop_tot=0.d0
	        Pop_spin=0.d0
	        num_GTOs=0d0
	        num_atom_cell=0d0
                
                WRITE(symbol,1071)
                
		OPEN(35, file='G_Vector',action='read')
		Read(35, *) num_cell

		allocate(cell_dir(num_cell,3), mat_cell(num_cell,num_cell))

		do i=1, num_cell
			Read(35, *) cell_dir(i,:)
		enddo
		close(35)

                do i=1, num_cell
                        do j=1,num_cell
                                vec_dir(:)=cell_dir(j,:)-cell_dir(i,:)
                                mat_cell(i,j)=num_row(vec_dir, cell_dir(1:num_cell,:), num_cell)
                        ENDDO
                ENDDO
                                
!                do i=1,num_cell
!                        WRITE(*,*) mat_cell(i,:)
!                        WRITE(*,*) cell_dir(i,:)
!                ENDDO
                
!                STOP
                
                OPEN(40,file=trim(filename),action='read')
                
                do WHILE(.NOT. b)
                        READ(40,*) c1
                        b=(trim(c1) =='DIRECT')
                ENDDO
                
                READ(40,*) c1
                
                do i=1,3
                        READ(40,*) c1,lattice(i,1:3)
                ENDDO
						
				lattice_rec(1,:)=cross_3d(lattice(2,1:3), lattice(3,1:3))*&
					2*pi/(dot_product(lattice(1,1:3), cross_3d(lattice(2,1:3), lattice(3,1:3))))
				lattice_rec(2,:)=cross_3d(lattice(3,1:3), lattice(1,1:3))*&
					2*pi/(dot_product(lattice(1,1:3), cross_3d(lattice(2,1:3), lattice(3,1:3))))				
				lattice_rec(3,:)=cross_3d(lattice(1,1:3), lattice(2,1:3))*&
					2*pi/(dot_product(lattice(1,1:3), cross_3d(lattice(2,1:3), lattice(3,1:3))))
                

		b= .FALSE.
                do WHILE(.NOT. b)
                        READ(40,*) c1
                        b=(trim(c1) =='ATOM')
                ENDDO
                
                do i=1,3
                        READ(40,'(A)') c1
                ENDDO
		
		b= .FALSE.
                do WHILE(.NOT. b)
                        READ(40,'(A)') c1
                        if (c1==symbol) THEN
				b=.TRUE.
			else
				READ(c1, '(I4)') num_temp
				READ(c1, 1082) num1, num2, atome_ele_list(num_temp), &
						atom_coord_list(num_temp,1:3)
			endif
                ENDDO




                
                b= .FALSE.
                do WHILE(.NOT. b)
                        READ(40,*) c1
                        b=(trim(c1) =='LOCAL')
                ENDDO
                
                do i=1,3
                        READ(40,'(A)') c1
                ENDDO
                
                num_atom_AO=0
                
               
                num_coeff=0
                type_AO='  '
                type_line=3

                
                b= .FALSE.
                do WHILE(.NOT. b)
	                READ(40,'(A)') c1
	                num_space=LEN_TRIM(c1)-LEN_TRIM(ADJUSTL(c1))
	                
	                if (c1==symbol) THEN
	                        if (type_line==2 ) THEN
		                        if (TRIM(ADJUSTL(type_AO))=='S') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                
		                                num_atom_AO=num_atom_AO+1
		                                num_total_AOc=num_total_AOc+1
		                                
		                                num_atom_AOh=num_atom_AOh+1
		                                num_total_AOh =num_total_AOh+1
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='SP') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+2, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+2, :,:)
		                                coeff(num_total_AOc+3, 1,3:5)=(/0,1,0/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/0,0,1/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+2)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+3)=1d0
		                                mat_temps(num_total_AOh+4, num_total_AOc+4)=1d0
		                                	                                
		                                
		                                num_atom_AO=num_atom_AO+4
		                                num_total_AOc=num_total_AOc+4
		                                
		                                num_atom_AOh=num_atom_AOh+4
		                                num_total_AOh =num_total_AOh+4
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+2, 1,3:5)=(/0,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/0,0,1/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+2)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+3)=1d0
		                                	                                
		                                
		                                num_atom_AO=num_atom_AO+3
		                                num_total_AOc=num_total_AOc+3
		                                
		                                num_atom_AOh=num_atom_AOh+3
		                                num_total_AOh =num_total_AOh+3
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+5, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+6, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+2, 1,3:5)=(/1,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/1,0,1/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/0,2,0/)
		                                coeff(num_total_AOc+5, 1,3:5)=(/0,1,1/)
		                                coeff(num_total_AOc+6, 1,3:5)=(/0,0,2/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)= -0.5
		                                mat_temps(num_total_AOh+1, num_total_AOc+4)= -0.5
		                                mat_temps(num_total_AOh+1, num_total_AOc+6)= 1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+3)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+5)=1d0
		                                mat_temps(num_total_AOh+4, num_total_AOc+1)= 0.5*sqrt(3d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+4)= -0.5*sqrt(3d0)
		                                mat_temps(num_total_AOh+5, num_total_AOc+2)=1d0
		                                
		                                num_atom_AO=num_atom_AO+6
		                                num_total_AOc=num_total_AOc+6
		                                
		                                num_atom_AOh=num_atom_AOh+5
		                                num_total_AOh =num_total_AOh+5
		                                
		                                
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+5, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+6, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+7, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+8, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+9, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+10, :,:)=coeff(num_total_AOc+1, :,:)
		                                
		                                coeff(num_total_AOc+2, 1,3:5)=(/2,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/2,0,1/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/1,2,0/)
		                                coeff(num_total_AOc+5, 1,3:5)=(/1,1,1/)
		                                coeff(num_total_AOc+6, 1,3:5)=(/1,0,2/)
		                                coeff(num_total_AOc+7, 1,3:5)=(/0,3,0/)
		                                coeff(num_total_AOc+8, 1,3:5)=(/0,2,1/)
		                                coeff(num_total_AOc+9, 1,3:5)=(/0,1,2/)
		                                coeff(num_total_AOc+10, 1,3:5)=(/0,0,3/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+3)= -3d0/10d0*sqrt(5d0)
		                                mat_temps(num_total_AOh+1, num_total_AOc+8)= -3d0/10d0*sqrt(5d0)
		                                mat_temps(num_total_AOh+1, num_total_AOc+10)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+1)= -1d0/4d0*sqrt(6d0)
		                                mat_temps(num_total_AOh+2, num_total_AOc+4)= -1d0/20d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+2, num_total_AOc+6)= 1d0/5d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+2)=-1d0/20d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+7)=-1d0/4d0*sqrt(6d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+9)=1d0/5d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+3)= 1d0/2d0*sqrt(3d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+8)= 1d0/2d0*sqrt(3d0)
		                                mat_temps(num_total_AOh+5, num_total_AOc+5)= 1d0
		                                mat_temps(num_total_AOh+6, num_total_AOc+1)= sqrt(10d0)/4d0
		                                mat_temps(num_total_AOh+6, num_total_AOc+4)= -3d0/4d0*sqrt(2d0)
		                                mat_temps(num_total_AOh+7, num_total_AOc+2)= 3d0/4d0*sqrt(2d0)
		                                mat_temps(num_total_AOh+7, num_total_AOc+7)= -1d0/4d0*sqrt(10d0)
		                                
		                                num_atom_AO=num_atom_AO+10
		                                num_total_AOc=num_total_AOc+10
		                                
		                                num_atom_AOh=num_atom_AOh+7
		                                num_total_AOh =num_total_AOh+7
		                        
		                        ENDIF

! num_GTOs 	                        
	                                if(TRIM(ADJUSTL(type_AO))=='S') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='SP') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                        num_GTOs(num_total_AOc-4)=num_GTO
	                                        num_GTOs(num_total_AOc-5)=num_GTO
	                                        
	                                
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                        num_GTOs(num_total_AOc-4)=num_GTO
	                                        num_GTOs(num_total_AOc-5)=num_GTO
	                                        num_GTOs(num_total_AOc-6)=num_GTO
	                                        num_GTOs(num_total_AOc-7)=num_GTO
	                                        num_GTOs(num_total_AOc-8)=num_GTO
	                                        num_GTOs(num_total_AOc-9)=num_GTO
	                                       
	                                
	                                ENDIF	                        
	                        
	                        ENDIF
                                b= .TRUE.

!	                        WRITE(*,*) 'It is the end'
	                ELSEIF (num_space<4) THEN
	                        if (type_line==2) THEN
		                        if (TRIM(ADJUSTL(type_AO))=='S') THEN
		                        

		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                
		                                num_atom_AO=num_atom_AO+1
		                                num_total_AOc=num_total_AOc+1
		                                
		                                num_atom_AOh=num_atom_AOh+1
		                                num_total_AOh =num_total_AOh+1
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='SP') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+2, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+2, :,:)
		                                coeff(num_total_AOc+3, 1,3:5)=(/0,1,0/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/0,0,1/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+2)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+3)=1d0
		                                mat_temps(num_total_AOh+4, num_total_AOc+4)=1d0
		                                		                                
		                                
		                                num_atom_AO=num_atom_AO+4
		                                num_total_AOc=num_total_AOc+4
		                                
		                                num_atom_AOh=num_atom_AOh+4
		                                num_total_AOh =num_total_AOh+4
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+2, 1,3:5)=(/0,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/0,0,1/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+2)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+3)=1d0

		                                
		                                num_atom_AO=num_atom_AO+3
		                                num_total_AOc=num_total_AOc+3
		                                
		                                num_atom_AOh=num_atom_AOh+3
		                                num_total_AOh =num_total_AOh+3
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+5, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+6, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+2, 1,3:5)=(/1,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/1,0,1/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/0,2,0/)
		                                coeff(num_total_AOc+5, 1,3:5)=(/0,1,1/)
		                                coeff(num_total_AOc+6, 1,3:5)=(/0,0,2/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)= -0.5
		                                mat_temps(num_total_AOh+1, num_total_AOc+4)= -0.5
		                                mat_temps(num_total_AOh+1, num_total_AOc+6)= 1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+3)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+5)=1d0
		                                mat_temps(num_total_AOh+4, num_total_AOc+1)= 0.5*sqrt(3d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+4)= -0.5*sqrt(3d0)
		                                mat_temps(num_total_AOh+5, num_total_AOc+2)=1d0
		                                
		                                num_atom_AO=num_atom_AO+6
		                                num_total_AOc=num_total_AOc+6
		                                
		                                num_atom_AOh=num_atom_AOh+5
		                                num_total_AOh =num_total_AOh+5
		                                
		                                
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+5, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+6, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+7, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+8, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+9, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+10, :,:)=coeff(num_total_AOc+1, :,:)
		                                
		                                coeff(num_total_AOc+2, 1,3:5)=(/2,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/2,0,1/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/1,2,0/)
		                                coeff(num_total_AOc+5, 1,3:5)=(/1,1,1/)
		                                coeff(num_total_AOc+6, 1,3:5)=(/1,0,2/)
		                                coeff(num_total_AOc+7, 1,3:5)=(/0,3,0/)
		                                coeff(num_total_AOc+8, 1,3:5)=(/0,2,1/)
		                                coeff(num_total_AOc+9, 1,3:5)=(/0,1,2/)
		                                coeff(num_total_AOc+10, 1,3:5)=(/0,0,3/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+3)= -3d0/10d0*sqrt(5d0)
		                                mat_temps(num_total_AOh+1, num_total_AOc+8)= -3d0/10d0*sqrt(5d0)
		                                mat_temps(num_total_AOh+1, num_total_AOc+10)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+1)= -1d0/4d0*sqrt(6d0)
		                                mat_temps(num_total_AOh+2, num_total_AOc+4)= -1d0/20d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+2, num_total_AOc+6)= 1d0/5d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+2)=-1d0/20d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+7)=-1d0/4d0*sqrt(6d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+9)=1d0/5d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+3)= 1d0/2d0*sqrt(3d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+8)= 1d0/2d0*sqrt(3d0)
		                                mat_temps(num_total_AOh+5, num_total_AOc+5)= 1d0
		                                mat_temps(num_total_AOh+6, num_total_AOc+1)= sqrt(10d0)/4d0
		                                mat_temps(num_total_AOh+6, num_total_AOc+4)= -3d0/4d0*sqrt(2d0)
		                                mat_temps(num_total_AOh+7, num_total_AOc+2)= 3d0/4d0*sqrt(2d0)
		                                mat_temps(num_total_AOh+7, num_total_AOc+7)= -1d0/4d0*sqrt(10d0)
		                                
		                                num_atom_AO=num_atom_AO+10
		                                num_total_AOc=num_total_AOc+10
		                                
		                                num_atom_AOh=num_atom_AOh+7
		                                num_total_AOh =num_total_AOh+7	                        
	                                ENDIF
	                                
! num_GTOs 	                        
	                                if(TRIM(ADJUSTL(type_AO))=='S') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='SP') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                        num_GTOs(num_total_AOc-4)=num_GTO
	                                        num_GTOs(num_total_AOc-5)=num_GTO
	                                        
	                                
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                        num_GTOs(num_total_AOc-4)=num_GTO
	                                        num_GTOs(num_total_AOc-5)=num_GTO
	                                        num_GTOs(num_total_AOc-6)=num_GTO
	                                        num_GTOs(num_total_AOc-7)=num_GTO
	                                        num_GTOs(num_total_AOc-8)=num_GTO
	                                        num_GTOs(num_total_AOc-9)=num_GTO
	                                       
	                                
	                                ENDIF	                        
	                        
	                        
	                        
	                        
	                        ENDIF



	                        
                                READ(c1, 1080) num_atom, str_ele
                                
                                num_atom_cell=num_atom_cell+1d0
                                !READ(c1, 1080) atome_num_list(num_atom_cell), atome_ele_list(num_atom_cell),&
                                !                         atom_coord_list(num_atom_cell,1:3)
                                
!                                WRITE(*,*) str_ele, sym_ele
                                
                                if (trim(ADJUSTL(str_ele))==trim(ADJUSTL(sym_ele))) THEN

	                                coeff(num_total_AOc+1:num_total_AOc+num_atom_AO,:,:)= & 
                                                coeff(num_total_AOc - num_atom_AO+1 : num_total_AOc,:,:)
                                                
                                        num_GTOs(num_total_AOc+1 :num_total_AOc+num_atom_AO ) = &
                                                num_GTOs(num_total_AOc - num_atom_AO+1 : num_total_AOc)
                                        
                                        
                                        
	                                do i= 1, num_atom_AO
						coeff(num_total_AOc+i, 1, 6:8)=atom_coord_list(num_atom, 1:3)
	                                        !READ(c1, 1080) num_atom, sym_ele, coeff(num_total_AOc+i, 1, 6:8)
!	                                        coeff(num_total_AOc+i, 1, 6:8)=coeff(num_total_AOc+i, 1, 6:8)*0.52917720859
	                                ENDDO
!	                                WRITE(*,*) coeff(num_total_AOc+1+1, 1, 6:8)

                                        mat_temps(num_total_AOh+1 :num_total_AOh+num_atom_AOh, &
                                                num_total_AOc+1 :num_total_AOc+num_atom_AO) =&
                                        mat_temps(num_total_AOh- num_atom_AOh+1 :num_total_AOh, &
                                                num_total_AOc-num_atom_AO+1 :num_total_AOc)        
	                                
                                        num_total_AOc = num_total_AOc + num_atom_AO
                                        num_total_AOh = num_total_AOh + num_atom_AOh
                                        
!                                        WRITE(*,*) num_total_AOc, num_atom_AO
!                                        WRITE(*,*) num_total_AOh, num_atom_AOh
                                ELSE
!change to Angstrom                                
                                        READ(c1, 1080) num_atom, sym_ele, atom_coord(:)
					atom_coord(:)=atom_coord_list(num_atom, 1:3)
!                                        atom_coord(:)=atom_coord(:)*0.52917720859
!                                        WRITE(*,*) atom_coord(:)
                                ENDIF
                                
!                                WRITE(*,*) num_total_AOc, num_atom_AO
!                                WRITE(*,*) num_total_AOh, num_atom_AOh   
                                
                                
	                        type_line=0

!	                        WRITE(*,*) 'atome and coordonne'
!	                ELSEIF (num_space>4 .AND. num_space <30) THEN
!	                        WRITE(*,*) 'type SP, P, D, F'
!	                        READ(c1, 1100) order_AO1, order_AO2, type_AO
!	                        type_line=1
	                
	                ELSEIF (num_space>4 .AND. num_space<35) THEN
	                
	                        
	                        
	                        if (type_line ==2) THEN
	                        
	                               
	                                
		                        if (TRIM(ADJUSTL(type_AO))=='S') THEN
		                        

		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                
		                                num_atom_AO=num_atom_AO+1
		                                num_total_AOc=num_total_AOc+1
		                                
		                                num_atom_AOh=num_atom_AOh+1
		                                num_total_AOh =num_total_AOh+1
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='SP') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+2, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+2, :,:)
		                                coeff(num_total_AOc+3, 1,3:5)=(/0,1,0/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/0,0,1/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+2)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+3)=1d0
		                                mat_temps(num_total_AOh+4, num_total_AOc+4)=1d0
		                                
		                                
		                                num_atom_AO=num_atom_AO+4
		                                num_total_AOc=num_total_AOc+4
		                                
		                                num_atom_AOh=num_atom_AOh+4
		                                num_total_AOh =num_total_AOh+4
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+2, 1,3:5)=(/0,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/0,0,1/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+2)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+3)=1d0

		                                
		                                num_atom_AO=num_atom_AO+3
		                                num_total_AOc=num_total_AOc+3
		                                
		                                num_atom_AOh=num_atom_AOh+3
		                                num_total_AOh =num_total_AOh+3
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+5, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+6, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+2, 1,3:5)=(/1,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/1,0,1/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/0,2,0/)
		                                coeff(num_total_AOc+5, 1,3:5)=(/0,1,1/)
		                                coeff(num_total_AOc+6, 1,3:5)=(/0,0,2/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+1)= -0.5
		                                mat_temps(num_total_AOh+1, num_total_AOc+4)= -0.5
		                                mat_temps(num_total_AOh+1, num_total_AOc+6)= 1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+3)=1d0
		                                mat_temps(num_total_AOh+3, num_total_AOc+5)=1d0
		                                mat_temps(num_total_AOh+4, num_total_AOc+1)= 0.5*sqrt(3d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+4)= -0.5*sqrt(3d0)
		                                mat_temps(num_total_AOh+5, num_total_AOc+2)=1d0
		                                
		                                num_atom_AO=num_atom_AO+6
		                                num_total_AOc=num_total_AOc+6
		                                
		                                num_atom_AOh=num_atom_AOh+5
		                                num_total_AOh =num_total_AOh+5
		                                
		                                
		                        
		                        ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
		                                coeff(num_total_AOc+1, 1, 6:8)=atom_coord(:)
		                                coeff(num_total_AOc+2, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+3, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+4, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+5, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+6, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+7, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+8, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+9, :,:)=coeff(num_total_AOc+1, :,:)
		                                coeff(num_total_AOc+10, :,:)=coeff(num_total_AOc+1, :,:)
		                                
		                                coeff(num_total_AOc+2, 1,3:5)=(/2,1,0/)
		                                coeff(num_total_AOc+3, 1,3:5)=(/2,0,1/)
		                                coeff(num_total_AOc+4, 1,3:5)=(/1,2,0/)
		                                coeff(num_total_AOc+5, 1,3:5)=(/1,1,1/)
		                                coeff(num_total_AOc+6, 1,3:5)=(/1,0,2/)
		                                coeff(num_total_AOc+7, 1,3:5)=(/0,3,0/)
		                                coeff(num_total_AOc+8, 1,3:5)=(/0,2,1/)
		                                coeff(num_total_AOc+9, 1,3:5)=(/0,1,2/)
		                                coeff(num_total_AOc+10, 1,3:5)=(/0,0,3/)
		                                
		                                mat_temps(num_total_AOh+1, num_total_AOc+3)= -3d0/10d0*sqrt(5d0)
		                                mat_temps(num_total_AOh+1, num_total_AOc+8)= -3d0/10d0*sqrt(5d0)
		                                mat_temps(num_total_AOh+1, num_total_AOc+10)=1d0
		                                mat_temps(num_total_AOh+2, num_total_AOc+1)= -1d0/4d0*sqrt(6d0)
		                                mat_temps(num_total_AOh+2, num_total_AOc+4)= -1d0/20d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+2, num_total_AOc+6)= 1d0/5d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+2)=-1d0/20d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+7)=-1d0/4d0*sqrt(6d0)
		                                mat_temps(num_total_AOh+3, num_total_AOc+9)=1d0/5d0*sqrt(30d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+3)= 1d0/2d0*sqrt(3d0)
		                                mat_temps(num_total_AOh+4, num_total_AOc+8)= 1d0/2d0*sqrt(3d0)
		                                mat_temps(num_total_AOh+5, num_total_AOc+5)= 1d0
		                                mat_temps(num_total_AOh+6, num_total_AOc+1)= sqrt(10d0)/4d0
		                                mat_temps(num_total_AOh+6, num_total_AOc+4)= -3d0/4d0*sqrt(2d0)
		                                mat_temps(num_total_AOh+7, num_total_AOc+2)= 3d0/4d0*sqrt(2d0)
		                                mat_temps(num_total_AOh+7, num_total_AOc+7)= -1d0/4d0*sqrt(10d0)
		                                
		                                num_atom_AO=num_atom_AO+10
		                                num_total_AOc=num_total_AOc+10
		                                
		                                num_atom_AOh=num_atom_AOh+7
		                                num_total_AOh =num_total_AOh+7
		                        
		                        ENDIF
                                ELSEIF (type_line ==0) THEN
                                        num_atom_AO=0
                                        num_atom_AOh=0
        
                                ENDIF

                                if (type_line /=0 ) THEN
	                                if(TRIM(ADJUSTL(type_AO))=='S') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='SP') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                        num_GTOs(num_total_AOc-4)=num_GTO
	                                        num_GTOs(num_total_AOc-5)=num_GTO
	                                        
	                                
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
	                                        num_GTOs(num_total_AOc)=num_GTO
	                                        num_GTOs(num_total_AOc-1)=num_GTO
	                                        num_GTOs(num_total_AOc-2)=num_GTO
	                                        num_GTOs(num_total_AOc-3)=num_GTO
	                                        num_GTOs(num_total_AOc-4)=num_GTO
	                                        num_GTOs(num_total_AOc-5)=num_GTO
	                                        num_GTOs(num_total_AOc-6)=num_GTO
	                                        num_GTOs(num_total_AOc-7)=num_GTO
	                                        num_GTOs(num_total_AOc-8)=num_GTO
	                                        num_GTOs(num_total_AOc-9)=num_GTO
	                                       
	                                
	                                ENDIF
                                ENDIF 
                                
	                        
	                        if (num_space>30 ) THEN
		                        
		                        READ(c1, 1090) order_AO2, type_AO
		                        coeff(num_total_AOc+1, 1, 3:5)=(/0d0,0d0,0d0/)
		                        
!		                        WRITE(*,*) type_AO
		                        
                                ELSE
                                        
	                                READ(c1, 1100) order_AO1, order_AO2, type_AO
	                                
!	                                WRITE(*,*) type_AO
	                                
	                                if (TRIM(ADJUSTL(type_AO))=='SP') THEN
	                                        coeff(num_total_AOc+1, 1, 3:5)=(/0d0, 0d0, 0d0/)
	                                        coeff(num_total_AOc+2, 1, 3:5)=(/1d0, 0d0, 0d0/)
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
	                                        coeff(num_total_AOc+1, 1, 3:5)=(/1d0, 0d0, 0d0/)
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
	                                        coeff(num_total_AOc+1, 1, 3:5)=(/2d0, 0d0, 0d0/)
	                                ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
	                                        coeff(num_total_AOc+1, 1, 3:5)=(/3d0,0d0,0d0/)

	                                ENDIF
	                              
	                                
!                                WRITE(*,*) num_total_AOc, num_atom_AO
!                                WRITE(*,*) num_total_AOh, num_atom_AOh   
     
                                ENDIF
                                 
                                type_line=1
	                
	                ELSEIF (num_space>=36) THEN
	                        if (type_line==1) THEN
	                                num_GTO=1

	                        ELSEIF (type_line==2) THEN
	                                num_GTO=num_GTO+1
	                                
	                        ENDIF
	                        
!	                        WRITE(*,*) num_GTO
	                        
	                        READ(c1, 1120) exp_c , s_c,  p_c, dfg_c
	                        
	                        if (TRIM(ADJUSTL(type_AO))=='S') THEN
	                                coeff(num_total_AOc+1, num_GTO, 1)= exp_c
	                                coeff(num_total_AOc+1, num_GTO, 2)= s_c
	                        ELSEIF(TRIM(ADJUSTL(type_AO))=='P') THEN
	                                coeff(num_total_AOc+1, num_GTO, 1)= exp_c
	                                coeff(num_total_AOc+1, num_GTO, 2)= p_c
	                        ELSEIF(TRIM(ADJUSTL(type_AO))=='SP') THEN
	                                coeff(num_total_AOc+1, num_GTO, 1)= exp_c
	                                coeff(num_total_AOc+1, num_GTO, 2)= s_c
	                                coeff(num_total_AOc+2, num_GTO, 1)= exp_c
	                                coeff(num_total_AOc+2, num_GTO, 2)= p_c
	                        ELSEIF(TRIM(ADJUSTL(type_AO))=='D') THEN
	                                coeff(num_total_AOc+1, num_GTO, 1)= exp_c
	                                coeff(num_total_AOc+1, num_GTO, 2)= dfg_c
	                        ELSEIF(TRIM(ADJUSTL(type_AO))=='F') THEN
	                                coeff(num_total_AOc+1, num_GTO, 1)= exp_c
	                                coeff(num_total_AOc+1, num_GTO, 2)= dfg_c
	                                
                                ENDIF
                                
!	                        WRITE(*,*) 'coefficient'
	                        
	                        type_line=2
	                
	                ENDIF
	                

                ENDDO
                
                ALLOCATE (Pop_tot(num_cell,num_total_AOh, num_total_AOh))
		
                if (type_cal=='U') THEN
			ALLOCATE(Pop_spin(num_cell, num_total_AOh, num_total_AOh))
                ENDIF                
                
                Pop_tot=0d0
                Pop_spin=0d0
            do m= 1,num_cell
! density matrix  alpha+beta                
                b= .FALSE.
                do WHILE( .NOT.b)
                        READ(40,'(A)') c1
                        b=(TRIM(ADJUSTL(c1))=='ALPHA+BETA ELECTRONS')
                
                ENDDO
                
                k= num_total_AOh/10
                
                if (k>=1) THEN
                        do i=1, k
                                do l=1,3
			                READ(40,'(A)') c1
!			                WRITE(*,*) c1
		                ENDDO
		                
                                do j=1, num_total_AOh
                                        READ(40, 1200) num_line, Pop_tot(m, j, (i-1)*10+1: i*10)
                                ENDDO
                        ENDDO
                ENDIF
                
                
                if (MOD(num_total_AOh, 10)>0) THEN
                        do l=1,2
		                READ(40,'(A)') c1
!		                WRITE(*,*) c1
	                ENDDO
                
                        do j=1, num_total_AOh
                                READ(40, 1200) num_line, Pop_tot(m, j, k*10+1 : num_total_AOh)
                        ENDDO
                ENDIF
                

! density matrix  alpha-beta                  
                if (type_cal=='U') THEN
	                b= .FALSE.
	                do WHILE( .NOT.b)
	                        READ(40,'(A)') c1
	                        b=(TRIM(ADJUSTL(c1))=='ALPHA-BETA ELECTRONS')
	                        
	                ENDDO
!	                WRITE(*,'(A)') c1
!	                WRITE(*,*) 'Right'
	                k= num_total_AOh/10
	                
	                if (k>=1) THEN
	                        do i=1, k
	                                do l=1,3
				                READ(40,'(A)') c1
			                ENDDO
			                
	                                do j=1, num_total_AOh
	                                        READ(40, 1200) num_line, Pop_spin(m, j, (i-1)*10+1: i*10)
	                                ENDDO
	                        ENDDO
	                ENDIF
	                
	                
	                if (MOD(num_total_AOh, 10)>0) THEN
	                        do l=1,2
			                READ(40,'(A)') c1
		                ENDDO
	                
	                        do j=1, num_total_AOh
	                                READ(40, 1200) num_line, Pop_spin(m, j, k*10+1 : num_total_AOh)
	                        ENDDO
	                ENDIF
                ENDIF
            ENDDO

        
                CLOSE(40)
				
            allocate(cart_gauss(num_total_AOc), mat(num_total_AOh, num_total_AOc))
			do i=1, num_total_AOc
				cart_gauss(i)%num_gauss=num_GTOs(i)
				cart_gauss(i)%l_gauss(1)=NINT(coeff(i,1,3))
				cart_gauss(i)%l_gauss(2)=NINT(coeff(i,1,4))
				cart_gauss(i)%l_gauss(3)=NINT(coeff(i,1,5))
				cart_gauss(i)%coord_atom=coeff(i,1,6:8)
				do j=1,num_GTOs(i)
					cart_gauss(i)%alpha(j)=coeff(i,j,1)
					cart_gauss(i)%coefficient(j)=coeff(i,j,2)
				enddo
			enddo
			mat(:,:)=mat_temps(1:num_total_AOh, 1:num_total_AOc)
						
			deallocate(coeff,num_GTOs, mat_temps)  
        ENDSUBROUTINE
        
        
        
        SUBROUTINE normalisation()
                IMPLICIT NONE
                
                INTEGER :: i, j, k, m,n, l
                DOUBLEPRECISION, DIMENSION(3,3) :: inv_lattice

                ALLOCATE (coord_frac(num_total_AOh,3))
		coord_frac=0.d0
				
				call inv_mat3(transpose(lattice), inv_lattice)
				k=1d0
                do i= 1, num_total_AOc
				
!re normalisation V3
			l=cart_gauss(i)%l_gauss(1)+cart_gauss(i)%l_gauss(2)+cart_gauss(i)%l_gauss(3)
			do j= 1, cart_gauss(i)%num_gauss
				cart_gauss(i)%coefficient(j)= cart_gauss(i)%coefficient(j)*(4* cart_gauss(i)%alpha(j))**(-l/2.d0-3.d0/4.d0)*&
					2**(1.d0/2.d0-l/2.d0)*pi**(3.d0/4.d0-5.d0/8.d0)*&
					sqrt(fact(2*l)/fact(l))
				cart_gauss(i)%coefficient(j) = cart_gauss(i)%coefficient(j)* &
					A(cart_gauss(i)%alpha(j), cart_gauss(i)%l_gauss(1), cart_gauss(i)%l_gauss(2), cart_gauss(i)%l_gauss(3))
			ENDDO

! list de coord fraction pour l'impulsion calcul
			if (l==2) then
				if (cart_gauss(i)%l_gauss(3)/=2d0) then
					coord_frac(k,:)=matmul(inv_lattice,cart_gauss(i)%coord_atom)
					k=k+1
				endif
			elseif(l==3) then
				if (cart_gauss(i)%l_gauss(1)>0d0 .OR. cart_gauss(i)%l_gauss(2)==3d0) then
					coord_frac(k,:)=matmul(inv_lattice,cart_gauss(i)%coord_atom)
					k=k+1
				endif	
			else
				coord_frac(k,:)=matmul(inv_lattice,cart_gauss(i)%coord_atom)
				k=k+1			
			endif

            ENDDO
        
        ENDSUBROUTINE normalisation
        
        SUBROUTINE WRITE_para(type_cal)
                IMPLICIT NONE
                CHARACTER(len=5) :: Num_col
                
                CHARACTER(len=1),  INTENT(in) :: type_cal
                
                INTEGER :: i, j
                
                
                OPEN(30, file='basis.dat', status='replace')
                
                do  i=1, num_total_AOc
                        WRITE(30,'(I5)') i
                        do j=1, cart_gauss(i)%num_gauss
                                WRITE(30, '(2E10.3)') cart_gauss(i)%alpha(j), cart_gauss(i)%coefficient(j)
                        ENDDO
						WRITE(30, '(3I4)') cart_gauss(i)%l_gauss
						WRITE(30, '(3E20.12)') cart_gauss(i)%coord_atom
                        WRITE(30,*) ' '
                ENDDO
                CLOSE(30)
                
                OPEN(20, file='mat_temps.dat', status='replace')
                write( Num_col,'(i5)') num_total_AOc
                do  i=1, num_total_AOh
                        WRITE(20, '('//trim(Num_col)//'E12.4)') mat(i, :)
                
                ENDDO
                CLOSE(20)
                
                OPEN(10, file='Pop_tot.dat', status='replace')
                write( Num_col,'(i5)') num_total_AOh
                do  i=1, num_total_AOh
                        WRITE(10, '('//trim(Num_col)//'E12.3)') Pop_tot(1,i,:)
                
                ENDDO
                CLOSE(10)
                

                if (type_cal=='U' ) THEN
               
	                OPEN(5, file='Pop_spin.dat', status='replace')
	                write( Num_col,'(i5)') num_total_AOh
	                do  i=1, num_total_AOh
	                        WRITE(5, '('//trim(Num_col)//'E12.3)') Pop_spin(1,i,:)
	                
	                ENDDO
	                CLOSE(5)
                ENDIF
              
        ENDSUBROUTINE




ENDMODULE Reading_Mod
