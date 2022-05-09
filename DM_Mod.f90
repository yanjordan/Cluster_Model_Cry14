MODULE DM_Mod
        USE Reading_Mod
        !$ use OMP_LIB
        
CONTAINS

        SUBROUTINE DM_r(num_p, r, cell1, cell2, Pop_mat, DM)
                IMPLICIT NONE
                
                INTEGER, intent(in) :: cell1, cell2, num_p
                
                DOUBLEPRECISION, DIMENSION(num_p, 3), intent(in) :: r
                
                DOUBLEPRECISION, DIMENSION(3) :: dep1, dep2
                
                DOUBLEPRECISION, DIMENSION(num_p, num_p) :: DM
                
                DOUBLEPRECISION, DIMENSION(num_total_AOh, num_total_AOh), intent(in) :: Pop_mat
                
                DOUBLEPRECISION, DIMENSION(num_total_AOc) :: AOc
                
                DOUBLEPRECISION, DIMENSION(num_total_AOh, num_p) :: AOh, AOh1
                                
                INTEGER :: i
                               
                dep1=cell_dir(cell1,1)*lattice(1,:)+cell_dir(cell1,2)*lattice(2,:)+cell_dir(cell1,3)*lattice(3,:)
                dep2=cell_dir(cell2,1)*lattice(1,:)+cell_dir(cell2,2)*lattice(2,:)+cell_dir(cell2,3)*lattice(3,:)

!vector of wave function
                do i=1,num_p
                      call AOc_total(r(i,:), num_total_AOc, cart_gauss, AOc, dep1)
                      if (num_total_AOc==num_total_AOh) THEN
	                     AOh(:,i)=AOc(:)
	              ELSE
	                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh(:,i))
	              ENDIF

                      call AOc_total(r(i,:), num_total_AOc, cart_gauss, AOc, dep2)
                      if (num_total_AOc==num_total_AOh) THEN
	                     AOh1(:,i)=AOc(:)
	              ELSE
	                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh1(:,i))
	              ENDIF
                         
                ENDDO                
                DM=0.d0
                
                DM=MATMUL(TRANSPOSE(AOh), MATMUL(Pop_mat, AOh1))
                                
        ENDSUBROUTINE


        SUBROUTINE  DM(num_p, r, tors, DM_Cry)
                IMPLICIT NONE

                INTEGER, intent(in) :: num_p
                
                CHARACTER(len=1), intent(in) :: tors
                
                DOUBLEPRECISION, DIMENSION(num_p, num_p) :: DM_Cry, DM_temps
                               
                DOUBLEPRECISION, DIMENSION(num_p,3), intent(in) :: r
                
                INTEGER :: i,j
                
                
                DM_Cry=0.d0
				DM_temps=0.d0
				
                if(tors=='R') THEN
	                do i=1, num_cell
	                        do j=1, num_cell
	                                if (Mat_cell(i,j)/=0d0) THEN
                                                call DM_r(num_p, r, i, j, Pop_tot(Mat_cell(i,j),:,:), DM_temps)
	                                        DM_Cry=DM_Cry+DM_temps
	                                ENDIF
	                        ENDDO
	                ENDDO
	                
                ELSEIF(tors=='U') THEN
	                do i=1, num_cell
	                        do j=1, num_cell
	                                if (Mat_cell(i,j)/=0d0) THEN
                                                call DM_r(num_p, r, i, j, Pop_spin(Mat_cell(i,j),:,:), DM_temps)
	                                        DM_Cry=DM_Cry+DM_temps									
	                                ENDIF
	                        ENDDO
	                ENDDO
                ENDIF

        ENDSUBROUTINE

        
! separate the orbital contribution
        SUBROUTINE SDM_r(num_p, r, cell1, cell2, num_orbital, list_orbital, Pop_mat, SDM)
                IMPLICIT NONE
                
                INTEGER, intent(in) :: cell1, cell2, num_p, num_orbital
                
                DOUBLEPRECISION, DIMENSION(num_p, 3), intent(in) :: r
                
                INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
                
                DOUBLEPRECISION, DIMENSION(3) :: dep1, dep2
                
                DOUBLEPRECISION, DIMENSION(num_p, num_p) :: SDM
                
                DOUBLEPRECISION, DIMENSION(num_total_AOh, num_total_AOh), intent(in) :: Pop_mat

				DOUBLEPRECISION, DIMENSION(num_orbital, num_orbital):: Pop_mat_s
                
                DOUBLEPRECISION, DIMENSION(num_total_AOc) :: AOc
                
                DOUBLEPRECISION, DIMENSION(num_total_AOh, num_p) :: AOh, AOh1

				DOUBLEPRECISION, DIMENSION(num_orbital,num_p) :: SAOh, SAOh1
                                
                INTEGER :: i

                
                dep1=cell_dir(cell1,1)*lattice(1,:)+cell_dir(cell1,2)*lattice(2,:)+cell_dir(cell1,3)*lattice(3,:)
                dep2=cell_dir(cell2,1)*lattice(1,:)+cell_dir(cell2,2)*lattice(2,:)+cell_dir(cell2,3)*lattice(3,:)

!vector of wave function
                do i=1,num_p
                      call AOc_total(r(i,:), num_total_AOc, cart_gauss, AOc, dep1)
                      if (num_total_AOc==num_total_AOh) THEN
	                     AOh(:,i)=AOc(:)
	              ELSE
	                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh(:,i))
	              ENDIF

                      call AOc_total(r(i,:), num_total_AOc, cart_gauss, AOc, dep2)
                      if (num_total_AOc==num_total_AOh) THEN
	                     AOh1(:,i)=AOc(:)
	              ELSE
	                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh1(:,i))
	              ENDIF
                         
                ENDDO 

		Pop_mat_s=Pop_mat(list_orbital(:),list_orbital(:))          

                SAOh=AOh(list_orbital(:),:)
                SAOh1=AOh1(list_orbital(:),:)

                SDM=0.d0
                
                SDM=MATMUL(TRANSPOSE(SAOh), MATMUL(Pop_mat_s, SAOh1))
                                
        ENDSUBROUTINE
        

        SUBROUTINE  SDM(num_p, r, tors, num_orbital, list_orbital, SDM_Cry)
                IMPLICIT NONE

                INTEGER, intent(in) :: num_p, num_orbital
                
                CHARACTER(len=1), intent(in) :: tors
                
                INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
                
                DOUBLEPRECISION, DIMENSION(num_p, num_p) :: SDM_Cry, SDM_temps
                
                DOUBLEPRECISION, DIMENSION(num_p,3) :: r
                
                INTEGER :: i,j
		                
                SDM_Cry=0.d0
				SDM_temps=0.d0
		
		
                if(tors=='R') THEN
	                do i=1, num_cell
	                        do j=1, num_cell
	                                if (Mat_cell(i,j)/=0d0) THEN
                                                call SDM_r(num_p, r, i, j, num_orbital, list_orbital, Pop_tot(Mat_cell(i,j),:,:), SDM_temps)
	                                        SDM_Cry=SDM_Cry+SDM_temps
	                                ENDIF
	                        ENDDO
	                ENDDO
                ELSEIF(tors=='U') THEN
                        do i=1, num_cell
	                        do j=1, num_cell
	                                if (Mat_cell(i,j)/=0d0) THEN
                                                call SDM_r(num_p, r, i, j, num_orbital, list_orbital, Pop_spin(Mat_cell(i,j),:,:), SDM_temps)
	                                        SDM_Cry=SDM_Cry+SDM_temps									
	                                ENDIF
	                        ENDDO
	                ENDDO
                ENDIF

        ENDSUBROUTINE

        SUBROUTINE CDM_r(num_p, r, cell1, cell2, CTDM, CSDM)
                IMPLICIT NONE
                
                INTEGER, intent(in) :: cell1, cell2, num_p
                
                DOUBLEPRECISION, DIMENSION(num_p, 3), intent(in) :: r
                
                DOUBLEPRECISION, DIMENSION(3) :: dep1, dep2
                
                DOUBLEPRECISION, DIMENSION(num_p, num_p) :: CTDM

                DOUBLEPRECISION, DIMENSION(num_p, num_p), optional :: CSDM

                DOUBLEPRECISION, DIMENSION(num_total_AOc, num_p) :: AOc, AOc1
                                
                INTEGER :: i,j

                dep1=cell_dir(cell1,1)*lattice(1,:)+cell_dir(cell1,2)*lattice(2,:)+cell_dir(cell1,3)*lattice(3,:)
                dep2=cell_dir(cell2,1)*lattice(1,:)+cell_dir(cell2,2)*lattice(2,:)+cell_dir(cell2,3)*lattice(3,:)

!vector of wave function
                do i=1,num_p
                      call AOc_total(r(i,:), num_total_AOc, cart_gauss, AOc(:,i), dep1)

                      call AOc_total(r(i,:), num_total_AOc, cart_gauss, AOc1(:,i), dep2)                         
                ENDDO 

                CTDM=0.d0        
                CTDM=MATMUL(TRANSPOSE(AOc), MATMUL(Pop_tot_car(mat_cell(cell1, cell2),:,:), AOc1))  

                if (uorr=='U') then
                        CSDM=0.d0
                        CSDM=MATMUL(TRANSPOSE(AOc), MATMUL(Pop_spin_car(mat_cell(cell1, cell2),:,:), AOc1))    
                endif       
        ENDSUBROUTINE
        

        SUBROUTINE  CDM(num_p, r, num_orbital, list_orbital, num_orbital1, list_orbital1, CTDM, CSDM)
                IMPLICIT NONE

                INTEGER, intent(in) :: num_p, num_orbital, num_orbital1
                
                INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital

                INTEGER, DIMENSION(num_orbital1), INTENT(in) :: list_orbital1
                
                DOUBLEPRECISION, DIMENSION(num_p, num_p) :: CTDM, CTDM_temps, CSDM_temps

                DOUBLEPRECISION, DIMENSION(num_p, num_p), optional :: CSDM
                
                DOUBLEPRECISION, DIMENSION(num_p,3) :: r

                doubleprecision, dimension(num_total_AOh, num_total_AOh) :: pop_temps
                
                INTEGER :: i,j
		                
                allocate(Pop_tot_car(num_cell,num_total_AOc, num_total_AOc))
	
		do i=1, num_cell
			pop_temps=0.d0
			pop_temps(list_orbital(:),list_orbital1(:))=Pop_tot(i,list_orbital(:),list_orbital1(:))
                        pop_temps(list_orbital1(:),list_orbital(:))=Pop_tot(i,list_orbital1(:),list_orbital(:))			
                        Pop_tot_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
		enddo
		
		if (uorr=='U') THEN
			ALLOCATE(Pop_spin_car(num_cell, num_total_AOc, num_total_AOc))
			
			do i=1, num_cell
				pop_temps=0.d0
				pop_temps(list_orbital(:),list_orbital1(:))=Pop_spin(i,list_orbital(:),list_orbital1(:))
                                pop_temps(list_orbital1(:),list_orbital(:))=Pop_spin(i,list_orbital1(:),list_orbital(:))
				Pop_spin_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
			enddo			
		ENDIF
                !do i=1, 436
                !        write(*,'(436F8.4)') pop_temps(i,:)
                !enddo


                if (uorr=='R') THEN
                        CTDM=0.d0
                        do i=1,num_cell
                                do j=1,num_cell
                                        if (Mat_cell(i, j)/=0d0) THEN
                                                call CDM_r(num_p, r, i, j, CTDM_temps)
                                                CTDM=CTDM+CTDM_temps
                                        endif
                                enddo
                        enddo
                        deallocate(Pop_tot_car)
                elseif(uorr=='U') then
                        CTDM=0.d0
                        CSDM=0.d0
                        do i=1,num_cell
                                do j=1,num_cell
                                        if (Mat_cell(i, j)/=0d0) THEN
                                                call CDM_r(num_p, r, i, j, CTDM_temps, CSDM_temps)
                                                CTDM=CTDM+CTDM_temps
                                                CSDM=CSDM+CSDM_temps
                                        endif
                                enddo
                        enddo
                        deallocate(Pop_tot_car,Pop_spin_car)
                endif



        ENDSUBROUTINE

ENDMODULE