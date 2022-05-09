MODULE ECD_Mod
        USE Reading_Mod
        USE Clementi
	!$ use OMP_LIB
        
CONTAINS

! ECHG function
        SUBROUTINE ECHG_Cluster(ori, num_row, num_col, v_dx, v_dy)
                IMPLICIT NONE
                doubleprecision, dimension(3), intent(in) :: ori, v_dx, v_dy
		integer, intent(in) :: num_row, num_col

		doubleprecision, dimension(num_row, num_col) :: Den, Spin, Den_cle
                doubleprecision, dimension(3) :: grid_r, A, B, C

		doubleprecision :: den_temps, spin_temps, dx, dy, den_ij, spin_ij

		integer :: i, j, k, l, IHFERM
		
		IHFERM=0d0
		allocate(Pop_tot_car(num_cell,num_total_AOh, num_total_AOh))
                Pop_tot_car=Pop_tot

		if (uorr=='U') THEN
			ALLOCATE(Pop_spin_car(num_cell, num_total_AOh, num_total_AOh))
			IHFERM=1d0
			Pop_spin_car=Pop_spin		
		ENDIF                

		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(grid_r, den_temps)
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: den_ij)		
		        do i=1, num_row
				do j=1, num_col
					grid_r=ori+(i-1)*v_dx+(j-1)*v_dy
                                        Den_cle(i,j)=Cle_DEN(grid_r)
                                        den_ij=0.d0
					do k=1,num_cell
						do l=1, num_cell
							if (Mat_cell(k, l)/=0d0) THEN
								call ECHG_restricted(grid_r, k, l, den_temps)
                                                                den_ij=den_ij+den_temps
							ENDIF
						enddo
					enddo
                                        Den(i,j)=den_ij
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
		elseif (uorr=='U') then
!$OMP PARALLEL PRIVATE(grid_r, den_temps, spin_temps)
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: den_ij,spin_ij)				
			do i=1, num_row
				do j=1, num_col
					grid_r=ori+(i-1)*v_dx+(j-1)*v_dy
                                        Den_cle(i,j)=Cle_DEN(grid_r)
                                        den_ij=0.d0
                                        spin_ij=0.d0
					do k=1,num_cell
						do l=1, num_cell
							if (Mat_cell(k, l)/=0d0) THEN
								call ECHG_unrestricted(grid_r, k, l, den_temps, spin_temps)
                                                                den_ij=den_ij+den_temps
                                                                spin_ij=spin_ij+spin_temps								
							ENDIF
						enddo
					enddo
                                        Den(i,j)=den_ij
					Spin(i,j)=spin_ij
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
		ENDIF

                deallocate(Pop_tot_car)
		if (uorr=='U') THEN
			deallocate(Pop_spin_car)
		ENDIF
		write(*,*) 'end the ECHG computation' 

! write the data to the file
                B=ori*0.529177208
		A=(ori+(num_row-1)*v_dx)*0.529177208
		C=(ori+(num_col-1)*v_dy)*0.529177208

                dx=(sqrt(v_dx(1)**2+v_dx(2)**2+v_dx(3)**2))*0.529177208
                dy=(sqrt(v_dy(1)**2+v_dy(2)**2+v_dy(3)**2))*0.529177208
                atom_coord_list=atom_coord_list*0.529177208
                lattice=lattice*0.529177208
                
                OPEN(70, file='ECHG.data', status='replace')
                WRITE(70,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(70,'(1P,6E12.5)') A(:), B(:)
                WRITE(70, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(70, '(1P,6E12.5)') ((Den(i,j),i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(70, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(70,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(70,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(70,'(1P,6E12.5)') A(:), B(:)
                WRITE(70, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(70, '(1P,6E12.5)') ((Spin(i,j),i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(70, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(70,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(70)

                OPEN(80, file='ECHG_PATO.data', status='replace')
                WRITE(80,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(80,'(1P,6E12.5)') A(:), B(:)
                WRITE(80, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(80, '(1P,6E12.5)') ((Den_cle(i,j),i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(80,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(80,'(1P,6E12.5)') A(:), B(:)
                WRITE(80, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(80, '(1P,6E12.5)') ( (0.d0 , i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(80)
        ENDSUBROUTINE


! SECHG function
        SUBROUTINE SECHG_Cluster(ori, num_row, num_col, v_dx, v_dy, num_orbital, list_orbital)
                IMPLICIT NONE
                doubleprecision, dimension(3), intent(in) :: ori, v_dx, v_dy
		integer, intent(in) :: num_row, num_col, num_orbital
                integer, dimension(num_orbital), intent(in) :: list_orbital

		doubleprecision, dimension(num_row, num_col) :: Den, Spin 
                doubleprecision, dimension(3) :: grid_r, A, B, C

		doubleprecision :: den_temps, spin_temps, dx, dy, den_ij, spin_ij

		integer :: i, j, k, l, IHFERM
		
		IHFERM=0d0
		allocate(Pop_tot_car(num_cell,num_total_AOh, num_total_AOh))
                Pop_tot_car=0.d0
                do i=1, num_cell
                        Pop_tot_car(i,list_orbital(:),list_orbital(:))=Pop_tot(i,list_orbital(:),list_orbital(:))
                enddo

		if (uorr=='U') THEN
			ALLOCATE(Pop_spin_car(num_cell, num_total_AOh, num_total_AOh))
			IHFERM=1d0
                        Pop_spin_car=0.d0
                        do i=1,num_cell
			        Pop_spin_car(i,list_orbital(:),list_orbital(:))=Pop_spin(i,list_orbital(:),list_orbital(:))
                        enddo	
		ENDIF                

		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(grid_r, den_temps)
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: den_ij)			
		        do i=1, num_row
				do j=1, num_col
					grid_r=ori+(i-1)*v_dx+(j-1)*v_dy
                                        den_ij=0.d0
					do k=1,num_cell
						do l=1, num_cell
							if (Mat_cell(k, l)/=0d0) THEN
								call ECHG_restricted(grid_r, k, l, den_temps)
								den_ij=den_ij+den_temps
							ENDIF
						enddo
					enddo
                                        Den(i,j)=den_ij
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
		elseif (uorr=='U') then
!$OMP PARALLEL PRIVATE(grid_r, den_temps, spin_temps)
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: den_ij,spin_ij)			
			do i=1, num_row
				do j=1, num_col
					grid_r=ori+(i-1)*v_dx+(j-1)*v_dy
                                        den_ij=0.d0
                                        spin_ij=0.d0
					do k=1,num_cell
						do l=1, num_cell
							if (Mat_cell(k, l)/=0d0) THEN
								call ECHG_unrestricted(grid_r, k, l, den_temps, spin_temps)
								den_ij=den_ij+den_temps
                                                                spin_ij=spin_ij+spin_temps
							ENDIF
						enddo
					enddo
                                        Den(i,j)=den_ij
					Spin(i,j)=spin_ij
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
		ENDIF

                deallocate(Pop_tot_car)
		if (uorr=='U') THEN
			deallocate(Pop_spin_car)
		ENDIF
		write(*,*) 'end the ECHG computation' 

! write the data to the file
                B=ori*0.529177208
		A=(ori+(num_row-1)*v_dx)*0.529177208
		C=(ori+(num_col-1)*v_dy)*0.529177208

                dx=(sqrt(v_dx(1)**2+v_dx(2)**2+v_dx(3)**2))*0.529177208
                dy=(sqrt(v_dy(1)**2+v_dy(2)**2+v_dy(3)**2))*0.529177208
                atom_coord_list=atom_coord_list*0.529177208
                lattice=lattice*0.529177208
                
                OPEN(70, file='ECHG.data', status='replace')
                WRITE(70,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(70,'(1P,6E12.5)') A(:), B(:)
                WRITE(70, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(70, '(1P,6E12.5)') ((Den(i,j),i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(70, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(70,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(70,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(70,'(1P,6E12.5)') A(:), B(:)
                WRITE(70, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(70, '(1P,6E12.5)') ((Spin(i,j),i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(70, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(70,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(70)

                OPEN(80, file='ECHG_PATO.data', status='replace')
                WRITE(80,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(80,'(1P,6E12.5)') A(:), B(:)
                WRITE(80, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(80, '(1P,6E12.5)') ((0.d0,i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(80,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(80,'(1P,6E12.5)') A(:), B(:)
                WRITE(80, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(80, '(1P,6E12.5)') ( (0.d0 , i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(80)
        ENDSUBROUTINE

! CECHG function
        SUBROUTINE CECHG_Cluster(ori, num_row, num_col, v_dx, v_dy, num_orbital, list_orbital, num_orbital1, list_orbital1)
                IMPLICIT NONE
                doubleprecision, dimension(3), intent(in) :: ori, v_dx, v_dy
		integer, intent(in) :: num_row, num_col, num_orbital, num_orbital1
                integer, dimension(num_orbital), intent(in) :: list_orbital
                integer, dimension(num_orbital1), intent(in) :: list_orbital1

		doubleprecision, dimension(num_row, num_col) :: Den, Spin 
                doubleprecision, dimension(3) :: grid_r, A, B, C

		doubleprecision :: den_temps, spin_temps, dx, dy, den_ij, spin_ij

		integer :: i, j, k, l, IHFERM
		
		IHFERM=0d0
		allocate(Pop_tot_car(num_cell,num_total_AOh, num_total_AOh))
                Pop_tot_car=0.d0
                do i=1, num_cell
                        Pop_tot_car(i,list_orbital1(:),list_orbital(:))=Pop_tot(i,list_orbital1(:),list_orbital(:))
                        Pop_tot_car(i,list_orbital(:),list_orbital1(:))=Pop_tot(i,list_orbital(:),list_orbital1(:))
                enddo

		if (uorr=='U') THEN
			ALLOCATE(Pop_spin_car(num_cell, num_total_AOh, num_total_AOh))
			IHFERM=1d0
                        Pop_spin_car=0.d0
                        do i=1,num_cell
                                Pop_spin_car(i,list_orbital1(:),list_orbital(:))=Pop_spin(i,list_orbital1(:),list_orbital(:))
                                Pop_spin_car(i,list_orbital(:),list_orbital1(:))=Pop_spin(i,list_orbital(:),list_orbital1(:))
                        enddo	
		ENDIF                

		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(grid_r, den_temps)
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: den_ij)			
		        do i=1, num_row
				do j=1, num_col
					grid_r=ori+(i-1)*v_dx+(j-1)*v_dy
                                        den_ij=0.d0
					do k=1,num_cell
						do l=1, num_cell
							if (Mat_cell(k, l)/=0d0) THEN
								call ECHG_restricted(grid_r, k, l, den_temps)
								den_ij=den_ij+den_temps
							ENDIF
						enddo
					enddo
                                        Den(i,j)=den_ij
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
		elseif (uorr=='U') then
!$OMP PARALLEL PRIVATE(grid_r, den_temps, spin_temps)
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: den_ij,spin_ij)			
			do i=1, num_row
				do j=1, num_col
					grid_r=ori+(i-1)*v_dx+(j-1)*v_dy
                                        den_ij=0.d0
                                        spin_ij=0.d0
					do k=1,num_cell
						do l=1, num_cell
							if (Mat_cell(k, l)/=0d0) THEN
								call ECHG_unrestricted(grid_r, k, l, den_temps, spin_temps)
								den_ij=den_ij+den_temps
                                                                spin_ij=spin_ij+spin_temps
							ENDIF
						enddo
					enddo
                                        Den(i,j)=den_ij
					Spin(i,j)=spin_ij
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
		ENDIF

                deallocate(Pop_tot_car)
		if (uorr=='U') THEN
			deallocate(Pop_spin_car)
		ENDIF
		write(*,*) 'end the ECHG computation' 

! write the data to the file
                B=ori*0.529177208
		A=(ori+(num_row-1)*v_dx)*0.529177208
		C=(ori+(num_col-1)*v_dy)*0.529177208

                dx=(sqrt(v_dx(1)**2+v_dx(2)**2+v_dx(3)**2))*0.529177208
                dy=(sqrt(v_dy(1)**2+v_dy(2)**2+v_dy(3)**2))*0.529177208
                atom_coord_list=atom_coord_list*0.529177208
                lattice=lattice*0.529177208
                
                OPEN(70, file='ECHG.data', status='replace')
                WRITE(70,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(70,'(1P,6E12.5)') A(:), B(:)
                WRITE(70, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(70, '(1P,6E12.5)') ((Den(i,j),i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(70, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(70,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(70,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(70,'(1P,6E12.5)') A(:), B(:)
                WRITE(70, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(70, '(1P,6E12.5)') ((Spin(i,j),i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(70, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(70,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(70)

                OPEN(80, file='ECHG_PATO.data', status='replace')
                WRITE(80,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(80,'(1P,6E12.5)') A(:), B(:)
                WRITE(80, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(80, '(1P,6E12.5)') ((0.d0,i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO

                WRITE(80,'(A3,I1,A4,2I5,3E12.5)') '-%-', ihferm, 'MAPN', num_row, num_col, dx, dy, 0.d0
                WRITE(80,'(1P,6E12.5)') A(:), B(:)
                WRITE(80, '(1P,3E12.5, 4X,2I4)') C(:), num_atom_cell, 3
                
                WRITE(80, '(1P,6E12.5)') ( (0.d0 , i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(80)
        ENDSUBROUTINE

!ECHG_restricted
        SUBROUTINE ECHG_restricted(r, cell1, cell2, den_r)
                IMPLICIT NONE
                doubleprecision, dimension(3), intent(in) :: r
                integer, intent(in) :: cell1, cell2

                doubleprecision :: den_r

                DOUBLEPRECISION, DIMENSION(3) :: dep1, dep2

                DOUBLEPRECISION, DIMENSION(num_total_AOc) :: AOc
                
                DOUBLEPRECISION, DIMENSION(num_total_AOh,1) :: AOh, AOh1
                                
                DOUBLEPRECISION, DIMENSION(1,1) :: Dens

                dep1=cell_dir(cell1,1)*lattice(1,:)+cell_dir(cell1,2)*lattice(2,:)+cell_dir(cell1,3)*lattice(3,:)
                dep2=cell_dir(cell2,1)*lattice(1,:)+cell_dir(cell2,2)*lattice(2,:)+cell_dir(cell2,3)*lattice(3,:)

! wave function cell1                
                call AOc_total(r, num_total_AOc, cart_gauss, AOc, dep1)
                if (num_total_AOc==num_total_AOh) THEN
                     AOh(:,1)=AOc(:)
                ELSE
                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh(:,1))
                ENDIF                
                
! wave function cell2               
                call AOc_total(r, num_total_AOc, cart_gauss, AOc, dep2)
                if (num_total_AOc==num_total_AOh) THEN
                     AOh1(:,1)=AOc(:)
                ELSE
                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh1(:,1))
                ENDIF 

                Dens=MATMUL(TRANSPOSE(AOh), MATMUL(Pop_tot_car(Mat_cell(cell1,cell2),:,:), AOh1))

                den_r=Dens(1,1)

        ENDSUBROUTINE

!ECHG_unrestricted
        SUBROUTINE ECHG_unrestricted(r, cell1, cell2, den_r, spin_r)
                IMPLICIT NONE
                doubleprecision, dimension(3), intent(in) :: r
                integer, intent(in) :: cell1, cell2

                doubleprecision :: den_r, spin_r

                DOUBLEPRECISION, DIMENSION(3) :: dep1, dep2

                DOUBLEPRECISION, DIMENSION(num_total_AOc) :: AOc
                
                DOUBLEPRECISION, DIMENSION(num_total_AOh,1) :: AOh, AOh1
                                
                DOUBLEPRECISION, DIMENSION(1,1) :: Dens, Spins

                dep1=cell_dir(cell1,1)*lattice(1,:)+cell_dir(cell1,2)*lattice(2,:)+cell_dir(cell1,3)*lattice(3,:)
                dep2=cell_dir(cell2,1)*lattice(1,:)+cell_dir(cell2,2)*lattice(2,:)+cell_dir(cell2,3)*lattice(3,:)

! wave function cell1                
                call AOc_total(r, num_total_AOc, cart_gauss, AOc, dep1)
                if (num_total_AOc==num_total_AOh) THEN
                     AOh(:,1)=AOc(:)
                ELSE
                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh(:,1))
                ENDIF                
                
! wave function cell2               
                call AOc_total(r, num_total_AOc, cart_gauss, AOc, dep2)
                if (num_total_AOc==num_total_AOh) THEN
                     AOh1(:,1)=AOc(:)
                ELSE
                     call AOh_total(num_total_AOc, num_total_AOh, mat, AOc, AOh1(:,1))
                ENDIF 

                Dens=MATMUL(TRANSPOSE(AOh), MATMUL(Pop_tot_car(Mat_cell(cell1,cell2),:,:), AOh1))
                Spins=MATMUL(TRANSPOSE(AOh), MATMUL(Pop_spin_car(Mat_cell(cell1,cell2),:,:), AOh1))

                den_r=Dens(1,1)
                spin_r=Spins(1,1)

        ENDSUBROUTINE
ENDMODULE