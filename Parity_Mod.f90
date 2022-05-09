MODULE Parity_Mod
        USE Reading_Mod
		!$ use OMP_LIB

CONTAINS
!Parity cluster
		subroutine Parity_cluster(ori, num_row, num_col, v_dx, v_dy)
			implicit none
			doubleprecision, dimension(3), intent(in) :: ori, v_dx, v_dy
			integer, intent(in) :: num_row, num_col

			doubleprecision, dimension(num_row, num_col) :: Den, Spin 
			doubleprecision, dimension(3) :: grid_r, A, B, C

			doubleprecision :: den_temps, spin_temps, dx, dy, den_ij, spin_ij


			integer :: i, j, k, l, IHFERM
		
			IHFERM=0d0

			allocate(Pop_tot_car(num_cell,num_total_AOc, num_total_AOc))

			do i=1, num_cell
				Pop_tot_car(i,:,:)=matmul(transpose(mat), matmul(Pop_tot(i,:,:), mat))
			enddo

			if (uorr=='U') THEN
				ALLOCATE(Pop_spin_car(num_cell, num_total_AOc, num_total_AOc))
				IHFERM=1d0
				
				do i=1, num_cell
					Pop_spin_car(i,:,:)=matmul(transpose(mat), matmul(Pop_spin(i,:,:), mat))
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
									call parity_restricted(grid_r, k, l, den_temps)
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
									call parity_unrestricted(grid_r, k, l, den_temps, spin_temps)
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
			write(*,*) 'end the Parity computation'                
             
                
! write the data to the file
                B=ori*0.529177208
				A=(ori+(num_row-1)*v_dx)*0.529177208
				C=(ori+(num_col-1)*v_dy)*0.529177208

                dx=(sqrt(v_dx(1)**2+v_dx(2)**2+v_dx(3)**2))*0.529177208
                dy=(sqrt(v_dy(1)**2+v_dy(2)**2+v_dy(3)**2))*0.529177208
                atom_coord_list=atom_coord_list*0.529177208
                lattice=lattice*0.529177208
                
                OPEN(70, file='PARITY.data', status='replace')
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
! used for drawing by the ECHG program
		OPEN(80, file='PARITY_PATO.data', status='replace')
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
                
                WRITE(80, '(1P,6E12.5)') ((0.d0,i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(80)
		endsubroutine
!SParity cluster
		subroutine SParity_cluster(ori, num_row, num_col, v_dx, v_dy, num_orbital, list_orbital)
			implicit none
			doubleprecision, dimension(3), intent(in) :: ori, v_dx, v_dy
			integer, intent(in) :: num_row, num_col, num_orbital
			integer, dimension(num_orbital), intent(in) :: list_orbital

			doubleprecision, dimension(num_row, num_col) :: Den, Spin 
			doubleprecision, dimension(3) :: grid_r, A, B, C

			doubleprecision :: den_temps, spin_temps, dx, dy, den_ij, spin_ij

			doubleprecision, dimension(num_total_AOh, num_total_AOh) :: pop_temps


			integer :: i, j, k, l, IHFERM
		
			IHFERM=0d0

			allocate(Pop_tot_car(num_cell,num_total_AOc, num_total_AOc))

			do i=1, num_cell
				pop_temps=0.d0
				pop_temps(list_orbital(:),list_orbital(:))=Pop_tot(i,list_orbital(:),list_orbital(:))
				Pop_tot_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
			enddo

			if (uorr=='U') THEN
				ALLOCATE(Pop_spin_car(num_cell, num_total_AOc, num_total_AOc))
				IHFERM=1d0
				
				do i=1, num_cell
					pop_temps=0.d0
					pop_temps(list_orbital(:),list_orbital(:))=Pop_spin(i,list_orbital(:),list_orbital(:))
					Pop_spin_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
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
									call parity_restricted(grid_r, k, l, den_temps)
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
									call parity_unrestricted(grid_r, k, l, den_temps, spin_temps)
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

			write(*,*) 'end the Parity computation'                
             
                
! write the data to the file
                B=ori*0.529177208
				A=(ori+(num_row-1)*v_dx)*0.529177208
				C=(ori+(num_col-1)*v_dy)*0.529177208

                dx=(sqrt(v_dx(1)**2+v_dx(2)**2+v_dx(3)**2))*0.529177208
                dy=(sqrt(v_dy(1)**2+v_dy(2)**2+v_dy(3)**2))*0.529177208
                atom_coord_list=atom_coord_list*0.529177208
                lattice=lattice*0.529177208
                
                OPEN(70, file='PARITY.data', status='replace')
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

! used for drawing by the ECHG program
				OPEN(80, file='PARITY_PATO.data', status='replace')
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
                
                WRITE(80, '(1P,6E12.5)') ((0.d0,i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(80)
		endsubroutine

!CParity cluster
		subroutine CParity_cluster(ori, num_row, num_col, v_dx, v_dy, num_orbital, list_orbital, num_orbital1, list_orbital1)
			implicit none
			doubleprecision, dimension(3), intent(in) :: ori, v_dx, v_dy
			integer, intent(in) :: num_row, num_col, num_orbital, num_orbital1
			integer, dimension(num_orbital), intent(in) :: list_orbital
			integer, dimension(num_orbital1), intent(in) :: list_orbital1

			doubleprecision, dimension(num_row, num_col) :: Den, Spin 
			doubleprecision, dimension(3) :: grid_r, A, B, C

			doubleprecision :: den_temps, spin_temps, dx, dy, den_ij, spin_ij

			doubleprecision, dimension(num_total_AOh, num_total_AOh) :: pop_temps


			integer :: i, j, k, l, IHFERM
		
			IHFERM=0d0

			allocate(Pop_tot_car(num_cell,num_total_AOc, num_total_AOc))

			do i=1, num_cell
				pop_temps=0.d0
				pop_temps(list_orbital(:),list_orbital1(:))=Pop_tot(i,list_orbital(:),list_orbital1(:))
				pop_temps(list_orbital1(:),list_orbital(:))=Pop_tot(i,list_orbital1(:),list_orbital(:))
				Pop_tot_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
			enddo

			if (uorr=='U') THEN
				ALLOCATE(Pop_spin_car(num_cell, num_total_AOc, num_total_AOc))
				IHFERM=1d0
				
				do i=1, num_cell
					pop_temps=0.d0
					pop_temps(list_orbital(:),list_orbital1(:))=Pop_spin(i,list_orbital(:),list_orbital1(:))
					pop_temps(list_orbital1(:),list_orbital(:))=Pop_spin(i,list_orbital1(:),list_orbital(:))
					Pop_spin_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
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
									call parity_restricted(grid_r, k, l, den_temps)
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
									call parity_unrestricted(grid_r, k, l, den_temps, spin_temps)
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

			write(*,*) 'end the Parity computation'                
             
                
! write the data to the file
                B=ori*0.529177208
				A=(ori+(num_row-1)*v_dx)*0.529177208
				C=(ori+(num_col-1)*v_dy)*0.529177208

                dx=(sqrt(v_dx(1)**2+v_dx(2)**2+v_dx(3)**2))*0.529177208
                dy=(sqrt(v_dy(1)**2+v_dy(2)**2+v_dy(3)**2))*0.529177208
                atom_coord_list=atom_coord_list*0.529177208
                lattice=lattice*0.529177208
                
                OPEN(70, file='PARITY.data', status='replace')
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

! used for drawing by the ECHG program
				OPEN(80, file='PARITY_PATO.data', status='replace')
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
                
                WRITE(80, '(1P,6E12.5)') ((0.d0,i=num_row,1,-1), j=1,num_col)
	                
                do i=1, num_atom_cell
                        WRITE(80, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)
                ENDDO
                
                do i=1,3
                        WRITE(80,'(3E20.12)') lattice(i,:)
                ENDDO
                
                CLOSE(80)
		endsubroutine


!parity function for the contribution between center cell1 and cell2 (restricted case)
	subroutine parity_restricted(x, cell1, cell2, parity_total)
		implicit none
		doubleprecision, dimension(3), intent(in) :: x
		integer, intent(in) :: cell1, cell2
		doubleprecision :: parity_total
		integer :: i, j
		
		parity_total=0.d0
		
		do i=1, num_total_AOc
			do j=1, num_total_AOc
				parity_total=parity_total+parity_nm_3D(x, cell1, cell2, i, j)*Pop_tot_car(Mat_cell(cell1, cell2),i,j)
			enddo
		enddo	
		
	endsubroutine

!parity function for the contribution between center cell1 and cell2 (unrestricted case)
	subroutine parity_unrestricted(x, cell1, cell2, parity_total, parity_spin)
		implicit none
		doubleprecision, dimension(3), intent(in) :: x
		integer, intent(in) :: cell1, cell2
		doubleprecision :: parity_total, parity_spin, parity_nm_temps
		integer :: i, j
		
		parity_total=0.d0
		parity_spin=0.d0
		
		do i=1, num_total_AOc
			do j=1, num_total_AOc
				parity_nm_temps=parity_nm_3D(x, cell1, cell2, i, j)
				parity_total=parity_total+parity_nm_temps*Pop_tot_car(Mat_cell(cell1, cell2),i,j)
				parity_spin=parity_spin+parity_nm_temps*Pop_spin_car(Mat_cell(cell1, cell2),i,j)
			enddo
		enddo			
	endsubroutine

! parity_nm(r) 3D (cell1 cell2) contribution 
	DOUBLEPRECISION FUNCTION parity_nm_3D(x, cell1, cell2, num_bs_n, num_bs_m)
		implicit none
		doubleprecision, dimension(3), intent(in) :: x
		integer, intent(in) :: num_bs_n, num_bs_m, cell1, cell2
		doubleprecision, dimension(3) :: dep1, dep2
		integer :: i,j
		
		parity_nm_3D=(0.d0, 0.d0)
		
		dep1=cell_dir(cell1,1)*lattice(1,:)+cell_dir(cell1,2)*lattice(2,:)+cell_dir(cell1,3)*lattice(3,:)
		dep2=cell_dir(cell2,1)*lattice(1,:)+cell_dir(cell2,2)*lattice(2,:)+cell_dir(cell2,3)*lattice(3,:)
		
		do i=1, cart_gauss(num_bs_n)%num_gauss
			do j=1, cart_gauss(num_bs_m)%num_gauss
				parity_nm_3D=parity_nm_3D+cart_gauss(num_bs_n)%coefficient(i)* cart_gauss(num_bs_m)%coefficient(j) * &
				parity_nm(x(1), cart_gauss(num_bs_n)%l_gauss(1), cart_gauss(num_bs_m)%l_gauss(1), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(1)+dep1(1), cart_gauss(num_bs_m)%coord_atom(1)+dep2(1))&
				* &
				parity_nm(x(2), cart_gauss(num_bs_n)%l_gauss(2), cart_gauss(num_bs_m)%l_gauss(2), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(2)+dep1(2), cart_gauss(num_bs_m)%coord_atom(2)+dep2(2))&
				* &
				parity_nm(x(3), cart_gauss(num_bs_n)%l_gauss(3), cart_gauss(num_bs_m)%l_gauss(3), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(3)+dep1(3), cart_gauss(num_bs_m)%coord_atom(3)+dep2(3))
			enddo
		enddo
	ENDFUNCTION

!Winger_nm 1D like for s, y, z
    DOUBLEPRECISION FUNCTION parity_nm( x, l_n, l_m, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            
            INTEGER, INTENT(in) :: l_n, l_m
            
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
            if (l_n==0d0) THEN
                if (l_m==0d0) THEN
                        parity_nm=parity_00(x, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        parity_nm=parity_01(x, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        parity_nm=parity_02(x, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        parity_nm=parity_03(x, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function parity_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                        
            elseif(l_n==1d0) THEN
                if (l_m==0d0) THEN
                        parity_nm=parity_10(x, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        parity_nm=parity_11(x, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        parity_nm=parity_12(x, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        parity_nm=parity_13(x, alpha_n, alpha_m, X_n, X_m)  
                else
                        WRITE(*,*) 'function parity_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            elseif(l_n==2d0) THEN
                if (l_m==0d0) THEN
                        parity_nm=parity_20(x, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        parity_nm=parity_21(x, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        parity_nm=parity_22(x, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        parity_nm=parity_23(x, alpha_n, alpha_m, X_n, X_m)  
                else
                        WRITE(*,*) 'function parity_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            elseif(l_n==3d0) THEN
                if (l_m==0d0) THEN
                        parity_nm=parity_30(x, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        parity_nm=parity_31(x, alpha_n, alpha_m, X_n, X_m)
                elseif(l_m==2d0) THEN
                        parity_nm=parity_32(x, alpha_n, alpha_m, X_n, X_m)
                elseif(l_m==3d0) THEN
                        parity_nm=parity_33(x, alpha_n, alpha_m, X_n, X_m) 
                else
                        WRITE(*,*) 'function parity_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            else
                    WRITE(*,*) 'function parity_nm wrong!!!'
                    WRITE(*,*) 'l_n = ', l_n
                    WRITE(*,*) 'the l should between 0 and 3 '
                    STOP
            ENDIF
    ENDFUNCTION



!Parity_nmij for gaussians
   DOUBLEPRECISION FUNCTION parity_00(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_00=(Exp(C_Parity)*Sqrt(Pi))/Sqrt(A_Parity)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_01(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_01=(Exp(C_Parity)*Sqrt(Pi)*(B_Parity + Y_Parity))/Sqrt(A_Parity)

    ENDFUNCTION

   DOUBLEPRECISION FUNCTION parity_02(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_02=(Exp(C_Parity)*Sqrt(Pi)*(1 + 2*A_Parity*(B_Parity + Y_Parity)**2))/(2.*A_Parity**1.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_03(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_03=(Exp(C_Parity)*Sqrt(Pi)*(B_Parity + Y_Parity)*(3 + 2*A_Parity*(B_Parity + Y_Parity)**2))/(2.*A_Parity**1.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_10(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_10=(Exp(C_Parity)*Sqrt(Pi)*(-B_Parity + X_Parity))/Sqrt(A_Parity)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_11(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_11=-(Exp(C_Parity)*Sqrt(Pi)*(1 + 2*A_Parity*(B_Parity - X_Parity)*(B_Parity + Y_Parity)))/(2.*A_Parity**1.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_12(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_12=(Exp(C_Parity)*Sqrt(Pi)*(X_Parity + B_Parity*(-3 + 2*A_Parity*B_Parity*(-B_Parity + X_Parity))&
						- 2*Y_Parity + 4*A_Parity*B_Parity*(-B_Parity + X_Parity)*Y_Parity + &
						2*A_Parity*(-B_Parity + X_Parity)*Y_Parity**2))/(2.*A_Parity**1.5)


    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_13(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_13=-(Exp(C_Parity)*Sqrt(Pi)*(3 + 4*A_Parity**2*(B_Parity - X_Parity)*(B_Parity + Y_Parity)**3 + &
						6*A_Parity*(B_Parity + Y_Parity)*(2*B_Parity - X_Parity + Y_Parity)))/(4.*A_Parity**2.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_20(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_20=(Exp(C_Parity)*Sqrt(Pi)*(1 + 2*A_Parity*(B_Parity - X_Parity)**2))/(2.*A_Parity**1.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_21(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_21=(Exp(C_Parity)*Sqrt(Pi)*(B_Parity*(3 + 2*A_Parity*(B_Parity - X_Parity)**2) - 2*X_Parity + &
						Y_Parity + 2*A_Parity*(B_Parity - X_Parity)**2*Y_Parity))/(2.*A_Parity**1.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_22(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_22=(Exp(C_Parity)*Sqrt(Pi)*(3 + 4*A_Parity**2*(B_Parity - X_Parity)**2*(B_Parity + Y_Parity)**2 + &
						2*A_Parity*(6*B_Parity**2 - 6*B_Parity*X_Parity + X_Parity**2 + 6*B_Parity*Y_Parity - &
						4*X_Parity*Y_Parity + Y_Parity**2)))/(4.*A_Parity**2.5)


    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_23(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_23=(Exp(C_Parity)*Sqrt(Pi)*(-6*X_Parity + B_Parity*(15 + 4*A_Parity**2*B_Parity**2*(B_Parity - &
					X_Parity)**2 + A_Parity*(20*B_Parity**2 - 24*B_Parity*X_Parity + 6*X_Parity**2)) + &
					9*Y_Parity + 6*A_Parity*(2*A_Parity*B_Parity**4 - 6*B_Parity*X_Parity - &
					4*A_Parity*B_Parity**3*X_Parity + X_Parity**2 + 2*B_Parity**2*(3 + &
					A_Parity*X_Parity**2))*Y_Parity + 6*A_Parity*(B_Parity*(3 + 2*A_Parity*(B_Parity - &
					X_Parity)**2) - 2*X_Parity)*Y_Parity**2 + &
					2*A_Parity*(1 + 2*A_Parity*(B_Parity - X_Parity)**2)*Y_Parity**3))/(4.*A_Parity**2.5)


    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_30(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_30=-(Exp(C_Parity)*Sqrt(Pi)*(3 + 2*A_Parity*(B_Parity - X_Parity)**2)*(B_Parity - X_Parity))/(2.*A_Parity**1.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_31(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_31=-(Exp(C_Parity)*Sqrt(Pi)*(3 + 4*A_Parity**2*(B_Parity - X_Parity)**3*(B_Parity + Y_Parity) + &
						6*A_Parity*(B_Parity - X_Parity)*(2*B_Parity - X_Parity + Y_Parity)))/(4.*A_Parity**2.5)

    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_32(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_32=-(Exp(C_Parity)*Sqrt(Pi)*(15*B_Parity + 20*A_Parity*B_Parity**3 + 4*A_Parity**2*B_Parity**5 - 9*X_Parity - &
						36*A_Parity*B_Parity**2*X_Parity - 12*A_Parity**2*B_Parity**4*X_Parity + &
						18*A_Parity*B_Parity*X_Parity**2 + 12*A_Parity**2*B_Parity**3*X_Parity**2 - &
						2*A_Parity*X_Parity**3 - 4*A_Parity**2*B_Parity**2*X_Parity**3 + 6*Y_Parity + &
						4*A_Parity*(2*B_Parity*(3 + A_Parity*(B_Parity - X_Parity)**2) - 3*X_Parity)*(B_Parity - &
						X_Parity)*Y_Parity + 2*A_Parity*(3 + 2*A_Parity*(B_Parity - X_Parity)**2)*(B_Parity - &
						X_Parity)*Y_Parity**2))/(4.*A_Parity**2.5)


    ENDFUNCTION
	
   DOUBLEPRECISION FUNCTION parity_33(x, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Parity, Y_Parity, A_Parity, B_Parity, C_Parity
			
			X_Parity=x-X_n
			Y_Parity=x+X_m
			A_Parity=alpha_n+alpha_m
			B_Parity=(alpha_n*X_Parity-alpha_m*Y_Parity)/(alpha_n+alpha_m)
			C_Parity=-alpha_n* X_Parity**2 -alpha_m* Y_Parity**2 +(alpha_n*X_Parity-alpha_m*Y_Parity)**2 /(alpha_n+alpha_m)
			
			parity_33=-(Exp(C_Parity)*Sqrt(Pi)*(15 + 8*A_Parity**3*(B_Parity - X_Parity)**3*(B_Parity + Y_Parity)**3 + &
						18*A_Parity*(5*B_Parity**2 - 5*B_Parity*X_Parity + X_Parity**2 + 5*B_Parity*Y_Parity - &
						3*X_Parity*Y_Parity + Y_Parity**2) + 12*A_Parity**2*(B_Parity - X_Parity)*(B_Parity + &
						Y_Parity)*(5*B_Parity**2 - 5*B_Parity*X_Parity + X_Parity**2 + 5*B_Parity*Y_Parity - &
						3*X_Parity*Y_Parity + Y_Parity**2)))/(8.*A_Parity**3.5)

    ENDFUNCTION

END MODULE