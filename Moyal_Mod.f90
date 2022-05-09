MODULE Moyal_Mod
        USE Reading_Mod
		!$ use OMP_LIB

CONTAINS

!Moyal(s,k)
	subroutine Moyal_Cluster(npu_row, npu_col, d_row, dir_s, hkl_list)
		implicit none
		integer, intent(in) :: npu_row, npu_col
		doubleprecision, intent(in) :: d_row
		doubleprecision, dimension(3), intent(in) :: dir_s
		integer, dimension(npu_col,3), intent(in) :: hkl_list

		integer :: i, j, k, IHFERM
		doubleprecision, dimension(3) :: u_dir_s, v_s, v_k
		!double complex, dimension(:,:), allocatable :: moyal_total, moyal_spin
		doubleprecision, dimension(:,:), allocatable :: moyal_total_real, moyal_spin_real
		doubleprecision, dimension(:,:), allocatable :: moyal_total_imag, moyal_spin_imag
		double complex ::moyal_ij_t, moyal_ij_s, moyal_t_temps, moyal_s_temps
		
		CHARACTER(len=3) :: Num_col
		
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
		
		u_dir_s=uni_vector(dir_s)
		!u_sk2=uni_vector(s_k2)
		
		!allocate(moyal_total(npu_row, npu_col))
		allocate(moyal_total_real(npu_row, npu_col), moyal_total_imag(npu_row, npu_col))
		
		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(v_s, v_k, moyal_t_temps, moyal_ij_t)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					!v_sk=d_row*(i-1)*u_sk1+d_col*(j-1)*u_sk2
					v_s=d_row*(i-1)*u_dir_s
					v_k=matmul(hkl_list(j,1:3),lattice_rec(:,:))
					moyal_ij_t=(0.d0, 0.d0)
					do k=1, num_cell
						call Moyal_restricted(v_s(1:3), v_k(1:3), moyal_t_temps, k)
						moyal_ij_t=moyal_ij_t+moyal_t_temps
						!moyal_total(i,j)=moyal_total(i,j)+moyal_t_temps
					enddo
					moyal_total_real(i,j)=dreal(moyal_ij_t)
					moyal_total_imag(i,j)=dimag(moyal_ij_t)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
			open(40, file='Moyal_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Moyal_Imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
				
		else
			!allocate(moyal_spin(npu_row, npu_col))
			allocate(moyal_spin_real(npu_row, npu_col), moyal_spin_imag(npu_row, npu_col))
		
!$OMP PARALLEL PRIVATE(v_s, v_k, moyal_t_temps, moyal_ij_t, moyal_s_temps, moyal_ij_s)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					!v_sk=d_row*(i-1)*u_sk1+d_col*(j-1)*u_sk2
					v_s=d_row*(i-1)*u_dir_s
					v_k=matmul(hkl_list(j,1:3),lattice_rec(:,:))
					moyal_ij_t=(0.d0, 0.d0)
					moyal_ij_s=(0.d0, 0.d0)
					do k=1, num_cell
						call Moyal_unrestricted(v_s(1:3), v_k(1:3), moyal_t_temps, moyal_s_temps, k)
						moyal_ij_t=moyal_ij_t+moyal_t_temps
						moyal_ij_s=moyal_ij_s+moyal_s_temps
						!moyal_total(i,j)=moyal_total(i,j)+moyal_t_temps
						!moyal_spin(i,j)=moyal_spin(i,j)+moyal_s_temps
					enddo
					moyal_total_real(i,j)=dreal(moyal_ij_t)
					moyal_total_imag(i,j)=dimag(moyal_ij_t)
					moyal_spin_real(i,j)=dreal(moyal_ij_s)
					moyal_spin_imag(i,j)=dimag(moyal_ij_s)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL

			open(40, file='Moyal_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_spin_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Moyal_Imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_spin_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
			deallocate(moyal_spin_real, moyal_spin_imag, Pop_spin_car)
			!deallocate(moyal_spin, Pop_spin_car)
		endif
		deallocate(moyal_total_real, moyal_total_imag, Pop_tot_car)
		!deallocate(moyal_total, Pop_tot_car)	
	endsubroutine
	
!SMoyal(s,k)
	subroutine SMoyal_Cluster(npu_row, npu_col, d_row, dir_s, hkl_list, num_orbital, list_orbital)
		implicit none
		integer, intent(in) :: npu_row, npu_col, num_orbital
		INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
		doubleprecision, intent(in) :: d_row
		doubleprecision, dimension(3), intent(in) :: dir_s
		integer, dimension(npu_col,3), intent(in) :: hkl_list

		integer :: i, j, k, IHFERM
		doubleprecision, dimension(3) :: u_dir_s, v_s, v_k
		!double complex, dimension(:,:), allocatable :: moyal_total, moyal_spin
		doubleprecision, dimension(:,:), allocatable :: moyal_total_real, moyal_spin_real
		doubleprecision, dimension(:,:), allocatable :: moyal_total_imag, moyal_spin_imag
		double complex ::moyal_ij_t, moyal_ij_s, moyal_t_temps, moyal_s_temps
		
		doubleprecision, dimension(num_total_AOh, num_total_AOh) :: pop_temps
		
		CHARACTER(len=3) :: Num_col
		
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
		
		u_dir_s=uni_vector(dir_s)
		!u_sk2=uni_vector(s_k2)
		
		!allocate(moyal_total(npu_row, npu_col))
		allocate(moyal_total_real(npu_row, npu_col), moyal_total_imag(npu_row, npu_col))
		
		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(v_s, v_k, moyal_t_temps, moyal_ij_t)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					!v_sk=d_row*(i-1)*u_sk1+d_col*(j-1)*u_sk2
					v_s=d_row*(i-1)*u_dir_s
					v_k=matmul(hkl_list(j,1:3),lattice_rec(:,:))
					moyal_ij_t=(0.d0, 0.d0)
					do k=1, num_cell
						call Moyal_restricted(v_s(1:3), v_k(1:3), moyal_t_temps, k)
						moyal_ij_t=moyal_ij_t+moyal_t_temps
						!moyal_total(i,j)=moyal_total(i,j)+moyal_t_temps
					enddo
					moyal_total_real(i,j)=dreal(moyal_ij_t)
					moyal_total_imag(i,j)=dimag(moyal_ij_t)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
			open(40, file='Moyal_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Moyal_Imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
				
		else
			!allocate(moyal_spin(npu_row, npu_col))
			allocate(moyal_spin_real(npu_row, npu_col), moyal_spin_imag(npu_row, npu_col))
		
!$OMP PARALLEL PRIVATE(v_s, v_k, moyal_t_temps, moyal_ij_t, moyal_s_temps, moyal_ij_s)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					!v_sk=d_row*(i-1)*u_sk1+d_col*(j-1)*u_sk2
					v_s=d_row*(i-1)*u_dir_s
					v_k=matmul(hkl_list(j,1:3),lattice_rec(:,:))
					moyal_ij_t=(0.d0, 0.d0)
					moyal_ij_s=(0.d0, 0.d0)
					do k=1, num_cell
						call Moyal_unrestricted(v_s(1:3), v_k(1:3), moyal_t_temps, moyal_s_temps, k)
						moyal_ij_t=moyal_ij_t+moyal_t_temps
						moyal_ij_s=moyal_ij_s+moyal_s_temps
						!moyal_total(i,j)=moyal_total(i,j)+moyal_t_temps
						!moyal_spin(i,j)=moyal_spin(i,j)+moyal_s_temps
					enddo
					moyal_total_real(i,j)=dreal(moyal_ij_t)
					moyal_total_imag(i,j)=dimag(moyal_ij_t)
					moyal_spin_real(i,j)=dreal(moyal_ij_s)
					moyal_spin_imag(i,j)=dimag(moyal_ij_s)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL

			open(40, file='Moyal_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_spin_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Moyal_Imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_spin_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
			deallocate(moyal_spin_real, moyal_spin_imag, Pop_spin_car)
			!deallocate(moyal_spin, Pop_spin_car)
		endif
		deallocate(moyal_total_real, moyal_total_imag, Pop_tot_car)
		!deallocate(moyal_total, Pop_tot_car)	
	endsubroutine	

	
!CMoyal(s,k)
	subroutine CMoyal_Cluster(npu_row, npu_col, d_row, dir_s, hkl_list, num_orbital, list_orbital, num_orbital1, list_orbital1)
		implicit none
		integer, intent(in) :: npu_row, npu_col, num_orbital, num_orbital1
		INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
		INTEGER, DIMENSION(num_orbital1), INTENT(in) :: list_orbital1
		doubleprecision, intent(in) :: d_row
		doubleprecision, dimension(3), intent(in) :: dir_s
		integer, dimension(npu_col,3), intent(in) :: hkl_list

		integer :: i, j, k, IHFERM
		doubleprecision, dimension(3) :: u_dir_s, v_s, v_k
		!double complex, dimension(:,:), allocatable :: moyal_total, moyal_spin
		doubleprecision, dimension(:,:), allocatable :: moyal_total_real, moyal_spin_real
		doubleprecision, dimension(:,:), allocatable :: moyal_total_imag, moyal_spin_imag
		double complex ::moyal_ij_t, moyal_ij_s, moyal_t_temps, moyal_s_temps
		
		doubleprecision, dimension(num_total_AOh, num_total_AOh) :: pop_temps
		
		CHARACTER(len=3) :: Num_col
		
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
		
		u_dir_s=uni_vector(dir_s)
		!u_sk2=uni_vector(s_k2)
		
		!allocate(moyal_total(npu_row, npu_col))
		allocate(moyal_total_real(npu_row, npu_col), moyal_total_imag(npu_row, npu_col))
		
		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(v_s, v_k, moyal_t_temps, moyal_ij_t)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					!v_sk=d_row*(i-1)*u_sk1+d_col*(j-1)*u_sk2
					v_s=d_row*(i-1)*u_dir_s
					v_k=matmul(hkl_list(j,1:3),lattice_rec(:,:))
					moyal_ij_t=(0.d0, 0.d0)
					do k=1, num_cell
						call Moyal_restricted(v_s(1:3), v_k(1:3), moyal_t_temps, k)
						moyal_ij_t=moyal_ij_t+moyal_t_temps
						!moyal_total(i,j)=moyal_total(i,j)+moyal_t_temps
					enddo
					moyal_total_real(i,j)=dreal(moyal_ij_t)
					moyal_total_imag(i,j)=dimag(moyal_ij_t)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
			open(40, file='Moyal_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Moyal_Imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
				
		else
			!allocate(moyal_spin(npu_row, npu_col))
			allocate(moyal_spin_real(npu_row, npu_col), moyal_spin_imag(npu_row, npu_col))
		
!$OMP PARALLEL PRIVATE(v_s, v_k, moyal_t_temps, moyal_ij_t, moyal_s_temps, moyal_ij_s)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					!v_sk=d_row*(i-1)*u_sk1+d_col*(j-1)*u_sk2
					v_s=d_row*(i-1)*u_dir_s
					v_k=matmul(hkl_list(j,1:3),lattice_rec(:,:))
					moyal_ij_t=(0.d0, 0.d0)
					moyal_ij_s=(0.d0, 0.d0)
					do k=1, num_cell
						call Moyal_unrestricted(v_s(1:3), v_k(1:3), moyal_t_temps, moyal_s_temps, k)
						moyal_ij_t=moyal_ij_t+moyal_t_temps
						moyal_ij_s=moyal_ij_s+moyal_s_temps
						!moyal_total(i,j)=moyal_total(i,j)+moyal_t_temps
						!moyal_spin(i,j)=moyal_spin(i,j)+moyal_s_temps
					enddo
					moyal_total_real(i,j)=dreal(moyal_ij_t)
					moyal_total_imag(i,j)=dimag(moyal_ij_t)
					moyal_spin_real(i,j)=dreal(moyal_ij_s)
					moyal_spin_imag(i,j)=dimag(moyal_ij_s)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL

			open(40, file='Moyal_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(40,'(1P, 3E12.5)') u_dir_s
			WRITE(40,'(1P, 6E12.5)') ((moyal_spin_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Moyal_Imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, num_atom_cell
			WRITE(50,'(1P, 3E12.5)') u_dir_s
			WRITE(50,'(1P, 6E12.5)') ((moyal_spin_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
			deallocate(moyal_spin_real, moyal_spin_imag, Pop_spin_car)
			!deallocate(moyal_spin, Pop_spin_car)
		endif
		deallocate(moyal_total_real, moyal_total_imag, Pop_tot_car)
		!deallocate(moyal_total, Pop_tot_car)	
	endsubroutine
	

!Moyal^{g_num}(s,k) for the contribution between center cell and g_num cell (restricted case)
	subroutine Moyal_restricted(s, k, moyal, g_num)
		implicit none
		doubleprecision, dimension(3), intent(in) :: s, k
		integer, intent(in) :: g_num
		double complex :: moyal
		integer :: i, j
		
		moyal=(0.d0, 0.d0)
		
		do i=1, num_total_AOc
			do j=1, num_total_AOc
				moyal=moyal+moyal_nm_3D(s, k, i, j, g_num)*Pop_tot_car(g_num,i,j)
			enddo
		enddo	
		
	endsubroutine

!Moyal^{g_num}(s,k) for the contribution between center cell and g_num cell (unrestricted case)
	subroutine Moyal_unrestricted(s, k, moyal_total, moyal_spin, g_num)
		implicit none
		doubleprecision, dimension(3), intent(in) :: s, k
		integer, intent(in) :: g_num
		double complex :: moyal_total, moyal_spin, moyal_nm_temps
		integer :: i, j
		
		moyal_total=(0.d0, 0.d0)
		moyal_spin=(0.d0, 0.d0)
		
		do i=1, num_total_AOc
			do j=1, num_total_AOc
				moyal_nm_temps=moyal_nm_3D(s, k, i, j, g_num)
				moyal_total=moyal_total+moyal_nm_temps*Pop_tot_car(g_num,i,j)
				moyal_spin=moyal_spin+moyal_nm_temps*Pop_spin_car(g_num,i,j)
			enddo
		enddo			
	endsubroutine
	
	
! Moyal_nm(r) 3D
	DOUBLE COMPLEX FUNCTION moyal_nm_3D(s, k, num_bs_n, num_bs_m, g_num)
		implicit none
		doubleprecision, dimension(3), intent(in) :: s, k
		integer, intent(in) :: num_bs_n, num_bs_m, g_num
		doubleprecision, dimension(3) :: dep
		integer :: i,j
		
		moyal_nm_3D=(0.d0, 0.d0)
		
		dep=cell_dir(g_num,1)*lattice(1,:)+cell_dir(g_num,2)*lattice(2,:)+cell_dir(g_num,3)*lattice(3,:)
		
		do i=1, cart_gauss(num_bs_n)%num_gauss
			do j=1, cart_gauss(num_bs_m)%num_gauss
				moyal_nm_3D=moyal_nm_3D+cart_gauss(num_bs_n)%coefficient(i)* cart_gauss(num_bs_m)%coefficient(j) * &
				moyal_nm(s(1), k(1), cart_gauss(num_bs_n)%l_gauss(1), cart_gauss(num_bs_m)%l_gauss(1), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(1), cart_gauss(num_bs_m)%coord_atom(1)+dep(1))&
				* &
				moyal_nm(s(2), k(2), cart_gauss(num_bs_n)%l_gauss(2), cart_gauss(num_bs_m)%l_gauss(2), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(2), cart_gauss(num_bs_m)%coord_atom(2)+dep(2))&
				* &
				moyal_nm(s(3), k(3), cart_gauss(num_bs_n)%l_gauss(3), cart_gauss(num_bs_m)%l_gauss(3), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(3), cart_gauss(num_bs_m)%coord_atom(3)+dep(3))
			enddo
		enddo
	ENDFUNCTION

!Winger_nm 1D like for s, y, z
    DOUBLE COMPLEX FUNCTION moyal_nm( s, k, l_n, l_m, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            
            INTEGER, INTENT(in) :: l_n, l_m
            
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
            
            if (l_n==0d0) THEN
                if (l_m==0d0) THEN
                            moyal_nm=moyal_00(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        moyal_nm=moyal_01(s, k, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        moyal_nm=moyal_02(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        moyal_nm=moyal_03(s, k, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function moyal_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                        
            elseif(l_n==1d0) THEN
                if (l_m==0d0) THEN
                            moyal_nm=moyal_10(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        moyal_nm=moyal_11(s, k, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        moyal_nm=moyal_12(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        moyal_nm=moyal_13(s, k, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function moyal_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            elseif(l_n==2d0) THEN
                if (l_m==0d0) THEN
                            moyal_nm=moyal_20(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        moyal_nm=moyal_21(s, k, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        moyal_nm=moyal_22(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        moyal_nm=moyal_23(s, k, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function moyal_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            elseif(l_n==3d0) THEN
                if (l_m==0d0) THEN
                            moyal_nm=moyal_30(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        moyal_nm=moyal_31(s, k, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        moyal_nm=moyal_32(s, k, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        moyal_nm=moyal_33(s, k, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function moyal_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            else
                    WRITE(*,*) 'function moyal_nm wrong!!!'
                    WRITE(*,*) 'l_n = ', l_n
                    WRITE(*,*) 'the l should between 0 and 3 '
                    STOP
            ENDIF
    ENDFUNCTION

!Moyal_nmij for gaussians
   DOUBLE COMPLEX FUNCTION moyal_00(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
            
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
			
            moyal_00=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi))/Sqrt(A_Moyal)

    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION moyal_01(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_01=((0,0.5)*EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/&
					(4.*A_Moyal))*Sqrt(Pi)*(k + (0,2)*A_Moyal*(B_Moyal + Y_Moyal)))/A_Moyal**1.5

    ENDFUNCTION
    

   DOUBLE COMPLEX FUNCTION moyal_02(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_02=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/&
					(4.*A_Moyal))*Sqrt(Pi)*(-k**2 + 4*A_Moyal**2*(B_Moyal + Y_Moyal)**2 +&
					A_Moyal*(2 - (0,4)*k*(B_Moyal + Y_Moyal))))/(4.*A_Moyal**2.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_03(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_03=-(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/&
				(4.*A_Moyal))*Sqrt(Pi)*((0,-1)*k + 2*A_Moyal*(B_Moyal + Y_Moyal))*&
				(-k**2 + 4*A_Moyal**2*(B_Moyal + Y_Moyal)**2 + &
				A_Moyal*(6 - (0,4)*k*(B_Moyal + Y_Moyal))))/(8.*A_Moyal**3.5)
    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION moyal_10(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_10=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/&
					(4.*A_Moyal))*Sqrt(Pi)*(-2*A_Moyal*B_Moyal + &
					(0,1)*k + 2*A_Moyal*X_Moyal))/(2.*A_Moyal**1.5)
    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION moyal_11(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_11=-(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/&
					(4.*A_Moyal))*Sqrt(Pi)*(k**2 - 4*A_Moyal**2*(B_Moyal - X_Moyal)*(B_Moyal + Y_Moyal) + &
					(0,2)*A_Moyal*((0,1) + k*(2*B_Moyal - X_Moyal + Y_Moyal))))/(4.*A_Moyal**2.5)
    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION moyal_12(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_12=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/&
				(4.*A_Moyal))*Sqrt(Pi)*((0,-1)*k**3 - 8*A_Moyal**3*(B_Moyal - X_Moyal)*&
				(B_Moyal + Y_Moyal)**2 + 2*A_Moyal*k*((0,3) + 3*B_Moyal*k - k*X_Moyal + 2*k*Y_Moyal) +&
				(0,4)*A_Moyal**2*(3*B_Moyal**2*k - (0,1)*(X_Moyal - 2*Y_Moyal) + &
				k*Y_Moyal*(-2*X_Moyal + Y_Moyal) + B_Moyal*((0,3) - 2*k*X_Moyal + 4*k*Y_Moyal))))/(8.*A_Moyal**3.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_13(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_13=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					(k**4 + 16*A_Moyal**4*(B_Moyal - X_Moyal)*(B_Moyal + Y_Moyal)**3 + &
					(0,2)*A_Moyal*k**2*((0,6) + 4*B_Moyal*k - k*X_Moyal + 3*k*Y_Moyal) - &
					(0,8)*A_Moyal**3*(B_Moyal + Y_Moyal)*(4*B_Moyal**2*k - 3*X_Moyal*((0,1) +&
					k*Y_Moyal) + Y_Moyal*((0,3) + k*Y_Moyal) + B_Moyal*((0,6) - &
					3*k*X_Moyal + 5*k*Y_Moyal)) - 12*A_Moyal**2*(-1 + k*(2*B_Moyal**2*k - &
					X_Moyal*((0,1) + k*Y_Moyal) + Y_Moyal*((0,3) + k*Y_Moyal) + B_Moyal*((0,4) - &
					k*X_Moyal + 3*k*Y_Moyal)))))/(16.*A_Moyal**4.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_20(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_20=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					(-k**2 + A_Moyal*(2 - (0,4)*k*(B_Moyal - X_Moyal)) + &
					4*A_Moyal**2*(B_Moyal - X_Moyal)**2))/(4.*A_Moyal**2.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_21(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_21=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					((0,-1)*k**3 - 8*A_Moyal**3*(B_Moyal - X_Moyal)**2*(B_Moyal + Y_Moyal) + &
					2*A_Moyal*k*((0,3) + k*(3*B_Moyal - 2*X_Moyal + Y_Moyal)) + &
					(0,4)*A_Moyal**2*(3*B_Moyal**2*k + X_Moyal*((0,-2) + k*(X_Moyal - 2*Y_Moyal)) + &
					(0,1)*Y_Moyal + B_Moyal*((0,3) - 4*k*X_Moyal + 2*k*Y_Moyal))))/(8.*A_Moyal**3.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_22(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_22=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					(k**4 + 16*A_Moyal**4*(B_Moyal - X_Moyal)**2*(B_Moyal + Y_Moyal)**2 + &
					(0,4)*A_Moyal*k**2*((0,3) + 2*B_Moyal*k + k*(-X_Moyal + Y_Moyal)) - &
					4*A_Moyal**2*(-3 + 6*B_Moyal**2*k**2 - (0,6)*k*(X_Moyal - Y_Moyal) + &
					k**2*(X_Moyal**2 - 4*X_Moyal*Y_Moyal + Y_Moyal**2) + 6*B_Moyal*k*((0,2) + &
					k*(-X_Moyal + Y_Moyal))) + 8*A_Moyal**3*((0,-4)*B_Moyal**3*k + &
					B_Moyal**2*(6 + (0,6)*k*(X_Moyal - Y_Moyal)) + Y_Moyal**2 + &
					X_Moyal**2*(1 - (0,2)*k*Y_Moyal) + (0,2)*X_Moyal*Y_Moyal*((0,2) + k*Y_Moyal) -&
					(0,2)*B_Moyal*(k*X_Moyal**2 + Y_Moyal*((0,3) + k*Y_Moyal) - &
					X_Moyal*((0,3) + 4*k*Y_Moyal)))))/(16.*A_Moyal**4.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_23(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_23=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					((0,1)*k**5 - 32*A_Moyal**5*(B_Moyal - X_Moyal)**2*(B_Moyal + Y_Moyal)**3 - &
					2*A_Moyal*k**3*((0,10) + 5*B_Moyal*k - 2*k*X_Moyal + 3*k*Y_Moyal) - &
					(0,4)*A_Moyal**2*k*(-15 + 10*B_Moyal**2*k**2 - (0,6)*k*(2*X_Moyal - 3*Y_Moyal) + &
					2*B_Moyal*k*((0,15) - 4*k*X_Moyal + 6*k*Y_Moyal) + k**2*(X_Moyal**2 - &
					6*X_Moyal*Y_Moyal + 3*Y_Moyal**2)) + (0,16)*A_Moyal**4*(B_Moyal + &
					Y_Moyal)*(5*B_Moyal**3*k + (0,1)*Y_Moyal**2 + 3*X_Moyal**2*((0,1) + k*Y_Moyal) - &
					2*X_Moyal*Y_Moyal*((0,3) + k*Y_Moyal) + B_Moyal**2*((0,10) - 8*k*X_Moyal + 7*k*Y_Moyal) +&
					B_Moyal*(3*k*X_Moyal**2 + 2*Y_Moyal*((0,4) + k*Y_Moyal) - 2*X_Moyal*((0,6) + 5*k*Y_Moyal))) + &
					8*A_Moyal**3*(10*B_Moyal**3*k**2 - 6*B_Moyal**2*k*((0,-5) + 2*k*X_Moyal - 3*k*Y_Moyal) + &
					3*k*X_Moyal**2*((0,1) + k*Y_Moyal) - 6*X_Moyal*(-1 + (0,3)*k*Y_Moyal + k**2*Y_Moyal**2) + &
					Y_Moyal*(-9 + (0,9)*k*Y_Moyal + k**2*Y_Moyal**2) + 3*B_Moyal*(-5 - (0,4)*k*(2*X_Moyal - 3*Y_Moyal) +&
					k**2*(X_Moyal**2 - 6*X_Moyal*Y_Moyal + 3*Y_Moyal**2)))))/(32.*A_Moyal**5.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_30(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_30=-(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					((0,-1)*k + 2*A_Moyal*(B_Moyal - X_Moyal))*(-k**2 + A_Moyal*(6 - (0,4)*k*(B_Moyal - X_Moyal)) +&
					4*A_Moyal**2*(B_Moyal - X_Moyal)**2))/(8.*A_Moyal**3.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_31(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_31= (EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					(k**4 + 16*A_Moyal**4*(B_Moyal - X_Moyal)**3*(B_Moyal + Y_Moyal) + &
					(0,2)*A_Moyal*k**2*((0,6) + k*(4*B_Moyal - 3*X_Moyal + Y_Moyal)) - &
					(0,8)*A_Moyal**3*(B_Moyal - X_Moyal)*(4*B_Moyal**2*k + X_Moyal*((0,-3) +&
					k*(X_Moyal - 3*Y_Moyal)) + (0,3)*Y_Moyal + B_Moyal*((0,6) - 5*k*X_Moyal + &
					3*k*Y_Moyal)) - 12*A_Moyal**2*(-1 + k*(2*B_Moyal**2*k + X_Moyal*((0,-3) + &
					k*(X_Moyal - Y_Moyal)) + (0,1)*Y_Moyal + B_Moyal*((0,4) + &
					k*(-3*X_Moyal + Y_Moyal))))))/(16.*A_Moyal**4.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_32(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_32=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					((0,1)*k**5 - 32*A_Moyal**5*(B_Moyal - X_Moyal)**3*(B_Moyal + Y_Moyal)**2 - &
					2*A_Moyal*k**3*((0,10) + 5*B_Moyal*k - 3*k*X_Moyal + 2*k*Y_Moyal) - &
					(0,4)*A_Moyal**2*k*(-15 + 10*B_Moyal**2*k**2 - (0,6)*k*(3*X_Moyal - &
					2*Y_Moyal) + 2*B_Moyal*k*((0,15) - 6*k*X_Moyal + 4*k*Y_Moyal) + &
					k**2*(3*X_Moyal**2 - 6*X_Moyal*Y_Moyal + Y_Moyal**2)) + &
					(0,16)*A_Moyal**4*(B_Moyal - X_Moyal)*(5*B_Moyal**3*k + (0,3)*Y_Moyal**2 - &
					3*X_Moyal*Y_Moyal*((0,2) + k*Y_Moyal) + X_Moyal**2*((0,1) + 2*k*Y_Moyal) + &
					B_Moyal**2*((0,10) - 7*k*X_Moyal + 8*k*Y_Moyal) + B_Moyal*(2*k*X_Moyal**2 + &
					3*Y_Moyal*((0,4) + k*Y_Moyal) - 2*X_Moyal*((0,4) + 5*k*Y_Moyal))) + &
					8*A_Moyal**3*(10*B_Moyal**3*k**2 - k**2*X_Moyal**3 - 6*B_Moyal**2*k*((0,-5) + &
					3*k*X_Moyal - 2*k*Y_Moyal) + (0,3)*Y_Moyal*((0,2) + k*Y_Moyal) + &
					3*k*X_Moyal**2*((0,3) + 2*k*Y_Moyal) - 3*X_Moyal*(-3 + (0,6)*k*Y_Moyal + k**2*Y_Moyal**2) + &
					3*B_Moyal*(-5 - (0,4)*k*(3*X_Moyal - 2*Y_Moyal) + k**2*(3*X_Moyal**2 - &
					6*X_Moyal*Y_Moyal + Y_Moyal**2)))))/(32.*A_Moyal**5.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION moyal_33(s, k, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: s, k, alpha_n, alpha_m, X_n, X_m
			
			DOUBLEPRECISION :: X_Moyal, Y_Moyal, A_Moyal, B_Moyal, C_Moyal
			
			X_Moyal=s/2-X_n
			Y_Moyal=s/2+X_m
			A_Moyal=alpha_n+alpha_m
			B_Moyal=(alpha_n*X_Moyal-alpha_m*Y_Moyal)/(alpha_n+alpha_m)
			C_Moyal=-alpha_n* X_Moyal**2 -alpha_m* Y_Moyal**2 +(alpha_n*X_Moyal-alpha_m*Y_Moyal)**2 /(alpha_n+alpha_m)
            
            moyal_33=(EXP(C_Moyal - (k*((0,4)*A_Moyal*B_Moyal + k))/(4.*A_Moyal))*Sqrt(Pi)*&
					(15*(1 + (3*(2*A_Moyal*B_Moyal - (0,1)*k)**2)/(2.*A_Moyal) + &
					(2*A_Moyal*B_Moyal - (0,1)*k)**4/(4.*A_Moyal**2) + (2*A_Moyal*B_Moyal - &
					(0,1)*k)**6/(120.*A_Moyal**3)) - 45*(1 + (2*A_Moyal*B_Moyal - &
					(0,1)*k)**2/(3.*A_Moyal) + (2*A_Moyal*B_Moyal - (0,1)*k)**4/(60.*A_Moyal**2))*&
					(2*A_Moyal*B_Moyal - (0,1)*k)*(X_Moyal - Y_Moyal) - 12*A_Moyal**2*(2*A_Moyal*B_Moyal - &
					(0,1)*k)*X_Moyal**2*(X_Moyal - Y_Moyal)*Y_Moyal**2 - 8*A_Moyal**3*X_Moyal**3*Y_Moyal**3 + &
					18*A_Moyal*(1 + (2*A_Moyal*B_Moyal - (0,1)*k)**2/A_Moyal + (2*A_Moyal*B_Moyal - &
					(0,1)*k)**4/(12.*A_Moyal**2))*(X_Moyal**2 - 3*X_Moyal*Y_Moyal + Y_Moyal**2) - &
					6*A_Moyal*(4*A_Moyal**2*B_Moyal**2 - k**2 + A_Moyal*(2 - &
					(0,4)*B_Moyal*k))*X_Moyal*Y_Moyal*(X_Moyal**2 - 3*X_Moyal*Y_Moyal + Y_Moyal**2) - &
					6*A_Moyal*(1 + (2*A_Moyal*B_Moyal - (0,1)*k)**2/(6.*A_Moyal))*(2*A_Moyal*B_Moyal - &
					(0,1)*k)*(X_Moyal**3 - 9*X_Moyal**2*Y_Moyal + &
					9*X_Moyal*Y_Moyal**2 - Y_Moyal**3)))/(8.*A_Moyal**3.5)

    ENDFUNCTION   	

END MODULE