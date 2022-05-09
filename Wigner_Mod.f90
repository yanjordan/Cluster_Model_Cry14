MODULE Wigner_Mod
        USE Reading_Mod
		!$ use OMP_LIB

CONTAINS

!Wigner(x,p)
	subroutine Wigner_Cluster(npu_row, npu_col, d_row, d_col, x_p1, x_p2)
		implicit none
		integer, intent(in) :: npu_row, npu_col
		doubleprecision, intent(in) :: d_row, d_col
		doubleprecision, dimension(6), intent(in) :: x_p1, x_p2

		integer :: i, j, k, IHFERM
		doubleprecision, dimension(6) :: u_xp1, u_xp2, v_xp
		!double complex, dimension(:,:), allocatable :: wigner_total, wigner_spin
		doubleprecision, dimension(:,:), allocatable :: wigner_total_real, wigner_spin_real
		doubleprecision, dimension(:,:), allocatable :: wigner_total_imag, wigner_spin_imag
		double complex ::wigner_ij_t, wigner_ij_s, wigner_t_temps, wigner_s_temps
		
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
			
		u_xp1=uni_vector(x_p1)
		u_xp2=uni_vector(x_p2)
		
		!allocate(wigner_total(npu_row, npu_col))
		allocate(wigner_total_real(npu_row, npu_col), wigner_total_imag(npu_row, npu_col))
		
		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(v_xp, wigner_t_temps, wigner_ij_t)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					v_xp=d_row*(i-1)*u_xp1+d_col*(j-1)*u_xp2
					wigner_ij_t=(0.d0, 0.d0)
					do k=1, num_cell
						call Wigner_restricted(v_xp(1:3), v_xp(4:6), wigner_t_temps, k)
						wigner_ij_t=wigner_ij_t+wigner_t_temps
						!wigner_total(i,j)=wigner_total(i,j)+wigner_t_temps
					enddo
					wigner_total_real(i,j)=dreal(wigner_ij_t)
					wigner_total_imag(i,j)=dimag(wigner_ij_t)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
			open(40, file='Wigner_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Wigner_imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
				
		else
	
			!allocate(wigner_spin(npu_row, npu_col))
			allocate(wigner_spin_real(npu_row, npu_col), wigner_spin_imag(npu_row, npu_col))
		
!$OMP PARALLEL PRIVATE(v_xp, wigner_t_temps, wigner_ij_t, wigner_s_temps, wigner_ij_s)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					v_xp=d_row*(i-1)*u_xp1+d_col*(j-1)*u_xp2
					wigner_ij_t=(0.d0, 0.d0)
					wigner_ij_s=(0.d0, 0.d0)
					do k=1, num_cell
						call Wigner_unrestricted(v_xp(1:3), v_xp(4:6), wigner_t_temps, wigner_s_temps, k)
						wigner_ij_t=wigner_ij_t+wigner_t_temps
						wigner_ij_s=wigner_ij_s+wigner_s_temps
						!wigner_total(i,j)=wigner_total(i,j)+wigner_t_temps
						!wigner_spin(i,j)=wigner_spin(i,j)+wigner_s_temps
					enddo
					wigner_total_real(i,j)=dreal(wigner_ij_t)
					wigner_total_imag(i,j)=dimag(wigner_ij_t)
					wigner_spin_real(i,j)=dreal(wigner_ij_s)
					wigner_spin_imag(i,j)=dimag(wigner_ij_s)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL

			open(40, file='Wigner_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_spin_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Wigner_imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_spin_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
			deallocate(wigner_spin_real, wigner_spin_imag, Pop_spin_car)
			!deallocate(wigner_spin, Pop_spin_car)
		endif
		deallocate(wigner_total_real, wigner_total_imag, Pop_tot_car)
		!deallocate(wigner_total, Pop_tot_car)	
		write(*,*) 'END the Wigner computation'
	endsubroutine

	
!SWigner(x,p)
	subroutine SWigner_Cluster(npu_row, npu_col, d_row, d_col, x_p1, x_p2, num_orbital, list_orbital)
		implicit none
		integer, intent(in) :: npu_row, npu_col, num_orbital
		INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
		doubleprecision, intent(in) :: d_row, d_col
		doubleprecision, dimension(6), intent(in) :: x_p1, x_p2

		integer :: i, j, k, IHFERM
		doubleprecision, dimension(6) :: u_xp1, u_xp2, v_xp
		!double complex, dimension(:,:), allocatable :: wigner_total, wigner_spin
		doubleprecision, dimension(:,:), allocatable :: wigner_total_real, wigner_spin_real
		doubleprecision, dimension(:,:), allocatable :: wigner_total_imag, wigner_spin_imag
		double complex ::wigner_ij_t, wigner_ij_s, wigner_t_temps, wigner_s_temps
		
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
		
		u_xp1=uni_vector(x_p1)
		u_xp2=uni_vector(x_p2)
		
		!allocate(wigner_total(npu_row, npu_col))
		allocate(wigner_total_real(npu_row, npu_col), wigner_total_imag(npu_row, npu_col))
		
		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(v_xp, wigner_t_temps, wigner_ij_t)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					v_xp=d_row*(i-1)*u_xp1+d_col*(j-1)*u_xp2
					wigner_ij_t=(0.d0, 0.d0)
					do k=1, num_cell
						call Wigner_restricted(v_xp(1:3), v_xp(4:6), wigner_t_temps, k)
						wigner_ij_t=wigner_ij_t+wigner_t_temps
						!wigner_total(i,j)=wigner_total(i,j)+wigner_t_temps
					enddo
					wigner_total_real(i,j)=dreal(wigner_ij_t)
					wigner_total_imag(i,j)=dimag(wigner_ij_t)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
			open(40, file='Wigner_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Wigner_imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
				
		else
			!allocate(wigner_spin(npu_row, npu_col))
			allocate(wigner_spin_real(npu_row, npu_col), wigner_spin_imag(npu_row, npu_col))
		
!$OMP PARALLEL PRIVATE(v_xp, wigner_t_temps, wigner_ij_t, wigner_s_temps, wigner_ij_s)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					v_xp=d_row*(i-1)*u_xp1+d_col*(j-1)*u_xp2
					wigner_ij_t=(0.d0, 0.d0)
					wigner_ij_s=(0.d0, 0.d0)
					do k=1, num_cell
						call Wigner_unrestricted(v_xp(1:3), v_xp(4:6), wigner_t_temps, wigner_s_temps, k)
						wigner_ij_t=wigner_ij_t+wigner_t_temps
						wigner_ij_s=wigner_ij_s+wigner_s_temps
						!wigner_total(i,j)=wigner_total(i,j)+wigner_t_temps
						!wigner_spin(i,j)=wigner_spin(i,j)+wigner_s_temps
					enddo
					wigner_total_real(i,j)=dreal(wigner_ij_t)
					wigner_total_imag(i,j)=dimag(wigner_ij_t)
					wigner_spin_real(i,j)=dreal(wigner_ij_s)
					wigner_spin_imag(i,j)=dimag(wigner_ij_s)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL

			open(40, file='Wigner_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_spin_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Wigner_imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_spin_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
			deallocate(wigner_spin_real, wigner_spin_imag, Pop_spin_car)
			!deallocate(wigner_spin, Pop_spin_car)
		endif
		deallocate(wigner_total_real, wigner_total_imag, Pop_tot_car)
		!deallocate(wigner_total, Pop_tot_car)	
	endsubroutine

!cross term
!CWigner(x,p)
	subroutine CWigner_Cluster(npu_row, npu_col, d_row, d_col, x_p1, x_p2, num_orbital, list_orbital, num_orbital1, list_orbital1)
		implicit none
		integer, intent(in) :: npu_row, npu_col, num_orbital, num_orbital1
		INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
		INTEGER, DIMENSION(num_orbital1), INTENT(in) :: list_orbital1
		doubleprecision, intent(in) :: d_row, d_col
		doubleprecision, dimension(6), intent(in) :: x_p1, x_p2

		integer :: i, j, k, IHFERM
		doubleprecision, dimension(6) :: u_xp1, u_xp2, v_xp
		!double complex, dimension(:,:), allocatable :: wigner_total, wigner_spin
		doubleprecision, dimension(:,:), allocatable :: wigner_total_real, wigner_spin_real
		doubleprecision, dimension(:,:), allocatable :: wigner_total_imag, wigner_spin_imag
		double complex ::wigner_ij_t, wigner_ij_s, wigner_t_temps, wigner_s_temps
		
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
		
		u_xp1=uni_vector(x_p1)
		u_xp2=uni_vector(x_p2)
		
		!allocate(wigner_total(npu_row, npu_col))
		allocate(wigner_total_real(npu_row, npu_col), wigner_total_imag(npu_row, npu_col))
		
		if (uorr=='R') then
!$OMP PARALLEL PRIVATE(v_xp, wigner_t_temps, wigner_ij_t)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					v_xp=d_row*(i-1)*u_xp1+d_col*(j-1)*u_xp2
					wigner_ij_t=(0.d0, 0.d0)
					do k=1, num_cell
						call Wigner_restricted(v_xp(1:3), v_xp(4:6), wigner_t_temps, k)
						wigner_ij_t=wigner_ij_t+wigner_t_temps
						!wigner_total(i,j)=wigner_total(i,j)+wigner_t_temps
					enddo
					wigner_total_real(i,j)=dreal(wigner_ij_t)
					wigner_total_imag(i,j)=dimag(wigner_ij_t)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL
			open(40, file='Wigner_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Wigner_imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((0.d0, i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
				
		else
			!allocate(wigner_spin(npu_row, npu_col))
			allocate(wigner_spin_real(npu_row, npu_col), wigner_spin_imag(npu_row, npu_col))
		
!$OMP PARALLEL PRIVATE(v_xp, wigner_t_temps, wigner_ij_t, wigner_s_temps, wigner_ij_s)
!$OMP DO SCHEDULE(RUNTIME)			
			do i=1, npu_row
				do j=1, npu_col
					v_xp=d_row*(i-1)*u_xp1+d_col*(j-1)*u_xp2
					wigner_ij_t=(0.d0, 0.d0)
					wigner_ij_s=(0.d0, 0.d0)
					do k=1, num_cell
						call Wigner_unrestricted(v_xp(1:3), v_xp(4:6), wigner_t_temps, wigner_s_temps, k)
						wigner_ij_t=wigner_ij_t+wigner_t_temps
						wigner_ij_s=wigner_ij_s+wigner_s_temps
						!wigner_total(i,j)=wigner_total(i,j)+wigner_t_temps
						!wigner_spin(i,j)=wigner_spin(i,j)+wigner_s_temps
					enddo
					wigner_total_real(i,j)=dreal(wigner_ij_t)
					wigner_total_imag(i,j)=dimag(wigner_ij_t)
					wigner_spin_real(i,j)=dreal(wigner_ij_s)
					wigner_spin_imag(i,j)=dimag(wigner_ij_s)
				enddo
			enddo
!$OMP END DO
!$OMP END PARALLEL

			open(40, file='Wigner_Real.dat', status='replace')
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_total_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(40,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(40,'(1P, 6E12.5)') u_xp1
			WRITE(40,'(1P, 6E12.5)') u_xp2
			WRITE(40,'(1P, 6E12.5)') ((wigner_spin_real(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(40, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(40,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(40)
			
			open(50, file='Wigner_imag.dat', status='replace')
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_total_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			
			write(50,'(A3, I1, A4, 2I5, 2E12.5, I4)') '-%-', IHFERM, 'MAPN', npu_row, npu_col, d_row, d_col, num_atom_cell
			WRITE(50,'(1P, 6E12.5)') u_xp1
			WRITE(50,'(1P, 6E12.5)') u_xp2
			WRITE(50,'(1P, 6E12.5)') ((wigner_spin_imag(i,j), i=1, npu_row),j=1,npu_col)
			
			do i=1, num_atom_cell
                WRITE(50, '(I4, 1X, A, 1P,3E20.12)') a_num(atome_ele_list(i)), atome_ele_list(i), atom_coord_list(i,1:3)*0.529177208
            ENDDO
			
			do i=1,3
                WRITE(50,'(3E20.12)') lattice(i,:)*0.529177208
            ENDDO
			close(50)
			deallocate(wigner_spin_real, wigner_spin_imag, Pop_spin_car)
			!deallocate(wigner_spin, Pop_spin_car)
		endif
		deallocate(wigner_total_real, wigner_total_imag, Pop_tot_car)
		!deallocate(wigner_total, Pop_tot_car)	
	endsubroutine



!Wigner^{g_num}(x,p) for the contribution between center cell and g_num cell (restricted case)
	subroutine Wigner_restricted(x, p, wigner, g_num)
		implicit none
		doubleprecision, dimension(3), intent(in) :: x, p
		integer, intent(in) :: g_num
		double complex :: wigner
		integer :: i, j
		
		wigner=(0.d0, 0.d0)
		
		do i=1, num_total_AOc
			do j=1, num_total_AOc
				wigner=wigner+wigner_nm_3D(x, p, i, j, g_num)*Pop_tot_car(g_num,i,j)
			enddo
		enddo	
		
	endsubroutine

!Wigner^{g_num}(x,p) for the contribution between center cell and g_num cell (unrestricted case)
	subroutine Wigner_unrestricted(x, p, wigner_total, wigner_spin, g_num)
		implicit none
		doubleprecision, dimension(3), intent(in) :: x, p
		integer, intent(in) :: g_num
		double complex :: wigner_total, wigner_spin, wigner_nm_temps
		integer :: i, j
		
		wigner_total=(0.d0, 0.d0)
		wigner_spin=(0.d0, 0.d0)
		
		do i=1, num_total_AOc
			do j=1, num_total_AOc
				wigner_nm_temps=wigner_nm_3D(x, p, i, j, g_num)
				wigner_total=wigner_total+wigner_nm_temps*Pop_tot_car(g_num,i,j)
				wigner_spin=wigner_spin+wigner_nm_temps*Pop_spin_car(g_num,i,j)
			enddo
		enddo			
	endsubroutine
	
	
! Wigner_nm(r) 3D
	DOUBLE COMPLEX FUNCTION wigner_nm_3D(x, p, num_bs_n, num_bs_m, g_num)
		implicit none
		doubleprecision, dimension(3), intent(in) :: x, p
		integer, intent(in) :: num_bs_n, num_bs_m, g_num
		doubleprecision, dimension(3) :: dep
		integer :: i,j
		
		wigner_nm_3D=(0.d0, 0.d0)
		
		dep=cell_dir(g_num,1)*lattice(1,:)+cell_dir(g_num,2)*lattice(2,:)+cell_dir(g_num,3)*lattice(3,:)
		
		do i=1, cart_gauss(num_bs_n)%num_gauss
			do j=1, cart_gauss(num_bs_m)%num_gauss
				wigner_nm_3D=wigner_nm_3D+cart_gauss(num_bs_n)%coefficient(i)* cart_gauss(num_bs_m)%coefficient(j) * &
				wigner_nm(x(1), p(1), cart_gauss(num_bs_n)%l_gauss(1), cart_gauss(num_bs_m)%l_gauss(1), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(1), cart_gauss(num_bs_m)%coord_atom(1)+dep(1))&
				* &
				wigner_nm(x(2), p(2), cart_gauss(num_bs_n)%l_gauss(2), cart_gauss(num_bs_m)%l_gauss(2), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(2), cart_gauss(num_bs_m)%coord_atom(2)+dep(2))&
				* &
				wigner_nm(x(3), p(3), cart_gauss(num_bs_n)%l_gauss(3), cart_gauss(num_bs_m)%l_gauss(3), &
				cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
				cart_gauss(num_bs_n)%coord_atom(3), cart_gauss(num_bs_m)%coord_atom(3)+dep(3))
			enddo
		enddo
	ENDFUNCTION

!Winger_nm 1D like for x, y, z
    DOUBLE COMPLEX FUNCTION wigner_nm( x, p, l_n, l_m, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            
            INTEGER, INTENT(in) :: l_n, l_m
            
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            if (l_n==0d0) THEN
                if (l_m==0d0) THEN
                            wigner_nm=wigner_00(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        wigner_nm=wigner_01(x, p, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        wigner_nm=wigner_02(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        wigner_nm=wigner_03(x, p, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function wigner_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                        
            elseif(l_n==1d0) THEN
                if (l_m==0d0) THEN
                            wigner_nm=wigner_10(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        wigner_nm=wigner_11(x, p, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        wigner_nm=wigner_12(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        wigner_nm=wigner_13(x, p, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function wigner_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            elseif(l_n==2d0) THEN
                if (l_m==0d0) THEN
                            wigner_nm=wigner_20(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        wigner_nm=wigner_21(x, p, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        wigner_nm=wigner_22(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        wigner_nm=wigner_23(x, p, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function wigner_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            elseif(l_n==3d0) THEN
                if (l_m==0d0) THEN
                            wigner_nm=wigner_30(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==1d0) THEN
                        wigner_nm=wigner_31(x, p, alpha_n, alpha_m, X_n, X_m)  
                elseif(l_m==2d0) THEN
                        wigner_nm=wigner_32(x, p, alpha_n, alpha_m, X_n, X_m) 
                elseif(l_m==3d0) THEN
                        wigner_nm=wigner_33(x, p, alpha_n, alpha_m, X_n, X_m)   
                else
                        WRITE(*,*) 'function wigner_nm wrong!!!'
                        WRITE(*,*) 'l_m = ', l_m
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF                 
            else
                    WRITE(*,*) 'function wigner_nm wrong!!!'
                    WRITE(*,*) 'l_n = ', l_n
                    WRITE(*,*) 'the l should between 0 and 3 '
                    STOP
            ENDIF
    ENDFUNCTION

!Wigner_nmij for gaussians
   DOUBLE COMPLEX FUNCTION wigner_00(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_00=(2*Sqrt(Pi))*(EXP(-(p**2 + (-2*x + X_n + X_m)**2*alpha_n*alpha_m &
					+ (0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/Sqrt(alpha_n + alpha_m))

    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION wigner_01(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_01=(2*Sqrt(Pi)*((0,-1)*p + (2*X_n - X_n - X_m)*alpha_n))*&
					(EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m + &
					(0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + alpha_m)**1.5)

    ENDFUNCTION
    

   DOUBLE COMPLEX FUNCTION wigner_02(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_02=(Sqrt(Pi)*(-2*p**2 + alpha_n - (0,4)*p*(2*X_n - X_n - X_m)&
					*alpha_n + 2*(-2*X_n + X_n + X_m)**2*alpha_n**2 + alpha_m))*&
					(EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m + &
					(0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + alpha_m)**2.5)

    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_03(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_03=(Sqrt(Pi)*((0,1)*p + (-2*X_n + X_n + X_m)*alpha_n)*&
					(2*p**2 - 3*alpha_n + (0,4)*p*(2*X_n - X_n - X_m)*alpha_n&
					- 2*(-2*X_n + X_n + X_m)**2*alpha_n**2 - 3*alpha_m))*&
					(EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m + &
					(0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + alpha_m)**3.5)
    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION wigner_10(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_10=(2*Sqrt(Pi)*((0,1)*p + (2*X_n - X_n - X_m)*alpha_m))*&
					(EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m + &
					(0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + alpha_m)**1.5)
    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION wigner_11(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_11=(Sqrt(Pi)*(2*p**2 + (0,2)*p*(2*X_n - X_n - X_m)*&
					(alpha_n - alpha_m) - alpha_m + alpha_n*(-1 + &
					8*X_n**2*alpha_m + 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m +&
					2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*alpha_m)))*(EXP(-(p**2 &
					+ (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m + (0,2)*p*&
					(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + X_m*alpha_m))/&
					(alpha_n + alpha_m))/(alpha_n + alpha_m)**2.5)
    ENDFUNCTION


   DOUBLE COMPLEX FUNCTION wigner_12(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_12=(Sqrt(Pi)*((0,-2)*p**3 + 2*p**2*(2*X_n - X_n - X_m)*&
					(2*alpha_n - alpha_m) + (2*X_n - X_n - X_m)*&
					(-(alpha_n*alpha_m) + alpha_m**2 + 2*alpha_n**2*(-1 &
					+ 4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m &
					+ X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m)) + &
					(0,1)*p*(2*(-2*X_n + X_n + X_m)**2*alpha_n**2 + &
					3*alpha_m + alpha_n*(3 - 16*X_n**2*alpha_m - &
					4*X_n**2*alpha_m - 8*X_n*X_m*alpha_m - 4*X_m**2*alpha_m &
					+ 16*X_n*(X_n + X_m)*alpha_m))))*&
					(EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m + &
					(0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + alpha_m)**3.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_13(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_13=(Sqrt(Pi)*(-4*p**4 - 3*alpha_n**2 - (0,4)*p**3*&
					(2*X_n - X_n - X_m)*(3*alpha_n - alpha_m) - 3*alpha_m**2 &
					+ 2*(-2*X_n + X_n + X_m)**2*alpha_n**3*(-3 + 8*X_n**2*alpha_m &
					+ 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - &
					8*X_n*(X_n + X_m)*alpha_m) + 6*alpha_n*alpha_m*(-1 + &
					4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m + &
					X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) + &
					(0,2)*p*(2*X_n - X_n - X_m)*(2*(-2*X_n + X_n + X_m)**2&
					*alpha_n**3 + 6*alpha_n*alpha_m - 3*alpha_m**2 - &
					3*alpha_n**2*(-3 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + &
					4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*&
					alpha_m)) + 12*p**2*((-2*X_n + X_n + X_m)**2*alpha_n**2 + &
					alpha_m - alpha_n*(-1 + 4*X_n**2*alpha_m + X_n**2*alpha_m &
					+ 2*X_n*X_m*alpha_m + X_m**2*alpha_m - 4*X_n*(X_n + &
					X_m)*alpha_m))))*(0.5*EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*&
					alpha_n*alpha_m + (0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n &
					- alpha_m) + X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + &
					alpha_m)**4.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_20(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_20=(Sqrt(Pi)*(-2*p**2 + alpha_n + (0,4)*p*(2*X_n - X_n - &
					X_m)*alpha_m + alpha_m*(1 + 8*X_n**2*alpha_m + &
					2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m &
					- 8*X_n*(X_n + X_m)*alpha_m)))*(EXP(-(p**2 + (-2*X_n + X_n &
					+ X_m)**2*alpha_n*alpha_m + (0,2)*p*(-(X_n*alpha_n) + &
					X_n*(alpha_n - alpha_m) + X_m*alpha_m))/(alpha_n + &
					alpha_m))/(alpha_n + alpha_m)**2.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_21(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_21=(Sqrt(Pi)*((0,2)*p**3 - 2*p**2*(2*X_n - X_n - X_m)*&
					(alpha_n - 2*alpha_m) + (2*X_n - X_n - X_m)*(alpha_n**2 &
					- 2*alpha_m**2 + alpha_n*alpha_m*(-1 + 8*X_n**2*alpha_m &
					+ 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m &
					- 8*X_n*(X_n + X_m)*alpha_m)) + (0,1)*p*(alpha_n*(-3 + &
					16*X_n**2*alpha_m + 4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m &
					+ 4*X_m**2*alpha_m - 16*X_n*(X_n + X_m)*alpha_m) - &
					alpha_m*(3 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + &
					4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + &
					X_m)*alpha_m))))*(EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*&
					alpha_n*alpha_m + (0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n &
					- alpha_m) + X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + &
					alpha_m)**3.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_22(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_22=(Sqrt(Pi)*(4*p**4 + 2*(-2*X_n + X_n + X_m)**2*alpha_n**3 &
					+ (0,8)*p**3*(2*X_n - X_n - X_m)*(alpha_n - alpha_m) + &
					alpha_m**2*(3 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + &
					4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + &
					X_m)*alpha_m) - 6*alpha_n*alpha_m*(-1 + 4*X_n**2*alpha_m &
					+ X_n**2*alpha_m + 2*X_n*X_m*alpha_m + X_m**2*alpha_m - &
					4*X_n*(X_n + X_m)*alpha_m) + alpha_n**2*(3 - 6*X_m**2*alpha_m &
					+ 64*X_n**4*alpha_m**2 + 4*X_n**4*alpha_m**2 + &
					16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*alpha_m**2 - &
					128*X_n**3*(X_n + X_m)*alpha_m**2 + 4*X_n*X_m*alpha_m*(-3 + &
					4*X_m**2*alpha_m) + 6*X_n**2*alpha_m*(-1 + 4*X_m**2*alpha_m)&
					- 8*X_n*(X_n + X_m)*alpha_m*(-3 + 4*X_n**2*alpha_m + &
					8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m) + 24*X_n**2*alpha_m*(-1 &
					+ 4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m)) + &
					(0,4)*p*(2*X_n - X_n - X_m)*(3*alpha_m**2 - 2*(-2*X_n + X_n + &
					X_m)**2*alpha_n*alpha_m**2 + alpha_n**2*(-3 + 8*X_n**2*alpha_m &
					+ 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - &
					8*X_n*(X_n + X_m)*alpha_m)) - 4*p**2*((-2*X_n + X_n + X_m)**2*&
					alpha_n**2 + alpha_m*(3 + 4*X_n**2*alpha_m + X_n**2*alpha_m + &
					2*X_n*X_m*alpha_m + X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) &
					+ alpha_n*(3 - 16*X_n**2*alpha_m - 4*X_n**2*alpha_m - &
					8*X_n*X_m*alpha_m - 4*X_m**2*alpha_m + 16*X_n*(X_n + X_m)*&
					alpha_m))))*(0.5*EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*&
					alpha_m + (0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + alpha_m)**4.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_23(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_23=(Sqrt(Pi)*((0,-4)*p**5 + 4*p**4*(2*X_n - X_n - X_m)*&
					(3*alpha_n - 2*alpha_m) - 4*p**2*(2*X_n - X_n - X_m)*&
					((-2*X_n + X_n + X_m)**2*alpha_n**3 - 6*alpha_m**2 - &
					3*alpha_n**2*(-3 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + &
					4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*alpha_m)&
					+ 3*alpha_n*alpha_m*(1 + 4*X_n**2*alpha_m + X_n**2*alpha_m +&
					2*X_n*X_m*alpha_m + X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m)) +& 
					(0,4)*p**3*(3*(-2*X_n + X_n + X_m)**2*alpha_n**2 + alpha_m*(5 + &
					4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m + &
					X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) + alpha_n*(5 - &
					24*X_n**2*alpha_m - 6*X_n**2*alpha_m - 12*X_n*X_m*alpha_m - &
					6*X_m**2*alpha_m + 24*X_n*(X_n + X_m)*alpha_m)) + &
					(0,1)*p*(2*(-2*X_n + X_n + X_m)**2*alpha_n**3*(-9 + &
					16*X_n**2*alpha_m + 4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m + &
					4*X_m**2*alpha_m - 16*X_n*(X_n + X_m)*alpha_m) - 3*alpha_m**2*(5 &
					+ 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + &
					2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*alpha_m) + 30*alpha_n*&
					alpha_m*(-1 + 4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*&
					alpha_m + X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) - &
					3*alpha_n**2*(5 - 6*X_m**2*alpha_m + 64*X_n**4*alpha_m**2 +&
					4*X_n**4*alpha_m**2 + 16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*&
					alpha_m**2 - 128*X_n**3*(X_n + X_m)*alpha_m**2 + 4*X_n*X_m*&
					alpha_m*(-3 + 4*X_m**2*alpha_m) + 6*X_n**2*alpha_m*(-1 + &
					4*X_m**2*alpha_m) - 8*X_n*(X_n + X_m)*alpha_m*(-3 + &
					4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m) + &
					24*X_n**2*alpha_m*(-1 + 4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m &
					+ 4*X_m**2*alpha_m))) + (2*X_n - X_n - X_m)*(2*(-2*X_n + X_n &
					+ X_m)**2*alpha_n**4 - 6*alpha_m**3 + 3*alpha_n*alpha_m**2*(-1 &
					+ 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + &
					2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*alpha_m) - 6*alpha_n**2*&
					alpha_m*(-2 + 4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*&
					alpha_m + X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) + &
					alpha_n**3*(9 - 10*X_m**2*alpha_m + 64*X_n**4*alpha_m**2 + &
					4*X_n**4*alpha_m**2 + 16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*&
					alpha_m**2 - 128*X_n**3*(X_n + X_m)*alpha_m**2 + 4*X_n*X_m*&
					alpha_m*(-5 + 4*X_m**2*alpha_m) - 8*X_n*(X_n + X_m)*alpha_m*(-5 &
					+ 4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m) + &
					2*X_n**2*alpha_m*(-5 + 12*X_m**2*alpha_m) + 8*X_n**2*alpha_m*(-5 &
					+ 12*X_n**2*alpha_m + 24*X_n*X_m*alpha_m + 12*X_m**2*alpha_m)))))*&
					(0.5*EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m + (0,2)*p*&
					(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + X_m*alpha_m))/&
					(alpha_n + alpha_m))/(alpha_n + alpha_m)**5.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_30(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_30=(Sqrt(Pi)*((0,-1)*p + (-2*X_n + X_n + X_m)*alpha_m)*&
					(2*p**2 - 3*alpha_n - (0,4)*p*(2*X_n - X_n - X_m)*alpha_m &
					- alpha_m*(3 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + &
					4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + &
					X_m)*alpha_m)))*(EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*&
					alpha_n*alpha_m + (0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n &
					- alpha_m) + X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n &
					+ alpha_m)**3.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_31(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_31=(Sqrt(Pi)*(-4*p**4 - (0,4)*p**3*(2*X_n - X_n - X_m)*&
					(alpha_n - 3*alpha_m) + 3*alpha_n**2*(-1 + 8*X_n**2*alpha_m &
					+ 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - &
					8*X_n*(X_n + X_m)*alpha_m) - 3*alpha_m**2*(1 + &
					8*X_n**2*alpha_m + 2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + &
					2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*alpha_m) + &
					2*alpha_n*alpha_m*(-3 + 32*X_n**4*alpha_m**2 + &
					2*X_n**4*alpha_m**2 + 8*X_n**3*X_m*alpha_m**2 + &
					12*X_n**2*X_m**2*alpha_m**2 + 8*X_n*X_m**3*alpha_m**2 +&
					2*X_m**4*alpha_m**2 - 64*X_n**3*(X_n + X_m)*alpha_m**2 + &
					48*X_n**2*(X_n + X_m)**2*alpha_m**2 - 16*X_n*(X_n + X_m)**3*&
					alpha_m**2) + (0,2)*p*(2*X_n - X_n - X_m)*(3*alpha_n**2 - &
					alpha_m**2*(9 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + &
					4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + &
					X_m)*alpha_m) + 6*alpha_n*alpha_m*(-1 + 4*X_n**2*alpha_m +&
					X_n**2*alpha_m + 2*X_n*X_m*alpha_m + X_m**2*alpha_m - &
					4*X_n*(X_n + X_m)*alpha_m)) - 12*p**2*(alpha_n*(-1 + &
					4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m + &
					X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) - alpha_m*(1 +&
					4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m + &
					X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m))))*&
					(0.5*EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*alpha_m +&
					(0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n - alpha_m) + &
					X_m*alpha_m))/(alpha_n + alpha_m))/(alpha_n + alpha_m)**4.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_32(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_32=(Sqrt(Pi)*((0,4)*p**5 - 4*p**4*(2*X_n - X_n - X_m)*&
					(2*alpha_n - 3*alpha_m) - 4*p**2*(2*X_n - X_n - X_m)*&
					(-3*alpha_n*alpha_m*(-1 + 8*X_n**2*alpha_m + &
					2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m &
					- 8*X_n*(X_n + X_m)*alpha_m) + 3*alpha_n**2*(-2 + &
					4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m +&
					X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) + &
					alpha_m**2*(9 + 4*X_n**2*alpha_m + X_n**2*alpha_m + &
					2*X_n*X_m*alpha_m + X_m**2*alpha_m - 4*X_n*(X_n + X_m)*&
					alpha_m)) - (0,4)*p**3*((-2*X_n + X_n + X_m)**2*alpha_n**2 &
					+ alpha_m*(5 + 12*X_n**2*alpha_m + 3*X_n**2*alpha_m + &
					6*X_n*X_m*alpha_m + 3*X_m**2*alpha_m - 12*X_n*(X_n + X_m)*&
					alpha_m) + alpha_n*(5 - 24*X_n**2*alpha_m - 6*X_n**2*alpha_m &
					- 12*X_n*X_m*alpha_m - 6*X_m**2*alpha_m + 24*X_n*(X_n + X_m)*&
					alpha_m)) + (2*X_n - X_n - X_m)*(-2*alpha_n*alpha_m**2*&
					(-6 + 20*X_n**2*alpha_m + 5*X_n**2*alpha_m + &
					10*X_n*X_m*alpha_m + 5*X_m**2*alpha_m - 20*X_n*(X_n + X_m)*&
					alpha_m) + alpha_m**3*(9 + 8*X_n**2*alpha_m + &
					2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - &
					8*X_n*(X_n + X_m)*alpha_m) + 6*alpha_n**3*(-1 + &
					4*X_n**2*alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m +&
					X_m**2*alpha_m - 4*X_n*(X_n + X_m)*alpha_m) + &
					alpha_n**2*alpha_m*(-3 - 6*X_m**2*alpha_m + &
					64*X_n**4*alpha_m**2 + 4*X_n**4*alpha_m**2 + &
					16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*alpha_m**2 - &
					128*X_n**3*(X_n + X_m)*alpha_m**2 + 4*X_n*X_m*alpha_m*(-3 +&
					4*X_m**2*alpha_m) + 6*X_n**2*alpha_m*(-1 + 4*X_m**2*alpha_m) &
					- 8*X_n*(X_n + X_m)*alpha_m*(-3 + 4*X_n**2*alpha_m + &
					8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m) + 24*X_n**2*alpha_m*(-1&
					+ 4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m))) + &
					(0,1)*p*(6*(-2*X_n + X_n + X_m)**2*alpha_n**3 + &
					3*alpha_m**2*(5 + 24*X_n**2*alpha_m + 6*X_n**2*alpha_m + &
					12*X_n*X_m*alpha_m + 6*X_m**2*alpha_m - 24*X_n*(X_n + &
					X_m)*alpha_m) - 2*alpha_n*alpha_m*(-15 + 9*X_m**2*alpha_m + &
					64*X_n**4*alpha_m**2 + 4*X_n**4*alpha_m**2 + &
					16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*alpha_m**2 - &
					128*X_n**3*(X_n + X_m)*alpha_m**2 + 3*X_n**2*alpha_m*(3 +&
					8*X_m**2*alpha_m) + 2*X_n*X_m*alpha_m*(9 + 8*X_m**2*alpha_m)&
					+ 12*X_n**2*alpha_m*(3 + 8*X_n**2*alpha_m + &
					16*X_n*X_m*alpha_m + 8*X_m**2*alpha_m) - 4*X_n*(X_n + X_m)*&
					alpha_m*(9 + 8*X_n**2*alpha_m + 16*X_n*X_m*alpha_m + &
					8*X_m**2*alpha_m)) + 3*alpha_n**2*(5 - 10*X_m**2*alpha_m + &
					64*X_n**4*alpha_m**2 + 4*X_n**4*alpha_m**2 + 16*X_n**3*X_m*&
					alpha_m**2 + 4*X_m**4*alpha_m**2 - 128*X_n**3*(X_n + X_m)*&
					alpha_m**2 + 4*X_n*X_m*alpha_m*(-5 + 4*X_m**2*alpha_m) - &
					8*X_n*(X_n + X_m)*alpha_m*(-5 + 4*X_n**2*alpha_m + &
					8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m) + 2*X_n**2*alpha_m*&
					(-5 + 12*X_m**2*alpha_m) + 8*X_n**2*alpha_m*(-5 + &
					12*X_n**2*alpha_m + 24*X_n*X_m*alpha_m + &
					12*X_m**2*alpha_m)))))*(0.5*EXP(-(p**2 + (-2*X_n + X_n + &
					X_m)**2*alpha_n*alpha_m + (0,2)*p*(-(X_n*alpha_n) + &
					X_n*(alpha_n - alpha_m) + X_m*alpha_m))/(alpha_n + &
					alpha_m))/(alpha_n + alpha_m)**5.5)
    ENDFUNCTION

   DOUBLE COMPLEX FUNCTION wigner_33(x, p, alpha_n, alpha_m, X_n, X_m)
            IMPLICIT NONE
            DOUBLEPRECISION, INTENT(in) :: x, p, alpha_n, alpha_m, X_n, X_m
            
            wigner_33=(Sqrt(Pi)*(8*p**6 + (0,24)*p**5*(2*X_n - X_n - X_m)*&
					(alpha_n - alpha_m) - 3*alpha_m**3*(5 + 24*X_n**2*alpha_m +&
					6*X_n**2*alpha_m + 12*X_n*X_m*alpha_m + 6*X_m**2*alpha_m - &
					24*X_n*(X_n + X_m)*alpha_m) + 6*(-2*X_n + X_n + X_m)**2*&
					alpha_n**4*(-3 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + &
					4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*&
					alpha_m) - 3*alpha_n**2*alpha_m*(15 - 24*X_m**2*alpha_m + &
					128*X_n**4*alpha_m**2 + 8*X_n**4*alpha_m**2 + &
					32*X_n**3*X_m*alpha_m**2 + 8*X_m**4*alpha_m**2 - &
					256*X_n**3*(X_n + X_m)*alpha_m**2 + 16*X_n*X_m*alpha_m*&
					(-3 + 2*X_m**2*alpha_m) + 24*X_n**2*alpha_m*(-1 + &
					2*X_m**2*alpha_m) - 32*X_n*(X_n + X_m)*alpha_m*(-3 + &
					2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m) +&
					96*X_n**2*alpha_m*(-1 + 2*X_n**2*alpha_m + 4*X_n*X_m*&
					alpha_m + 2*X_m**2*alpha_m)) + 3*alpha_n*alpha_m**2*(-15 &
					+ 6*X_m**2*alpha_m + 64*X_n**4*alpha_m**2 + 4*X_n**4*&
					alpha_m**2 + 16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*&
					alpha_m**2 - 128*X_n**3*(X_n + X_m)*alpha_m**2 + &
					6*X_n**2*alpha_m*(1 + 4*X_m**2*alpha_m) + 4*X_n*X_m*&
					alpha_m*(3 + 4*X_m**2*alpha_m) + 24*X_n**2*alpha_m*(1 &
					+ 4*X_n**2*alpha_m + 8*X_n*X_m*alpha_m + 4*X_m**2*&
					alpha_m) - 8*X_n*(X_n + X_m)*alpha_m*(3 + 4*X_n**2*&
					alpha_m + 8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m)) - &
					12*p**4*(2*(-2*X_n + X_n + X_m)**2*alpha_n**2 + alpha_m*&
					(5 + 8*X_n**2*alpha_m + 2*X_n**2*alpha_m + 4*X_n*X_m*&
					alpha_m + 2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*alpha_m) &
					+ alpha_n*(5 - 24*X_n**2*alpha_m - 6*X_n**2*alpha_m - &
					12*X_n*X_m*alpha_m - 6*X_m**2*alpha_m + 24*X_n*(X_n + &
					X_m)*alpha_m)) - (0,8)*p**3*(2*X_n - X_n - X_m)*&
					(alpha_n - alpha_m)*((-2*X_n + X_n + X_m)**2*alpha_n**2 &
					+ alpha_m*(15 + 4*X_n**2*alpha_m + X_n**2*alpha_m + &
					2*X_n*X_m*alpha_m + X_m**2*alpha_m - 4*X_n*(X_n + X_m)*&
					alpha_m) + alpha_n*(15 - 32*X_n**2*alpha_m - 8*X_n**2*&
					alpha_m - 16*X_n*X_m*alpha_m - 8*X_m**2*alpha_m + &
					32*X_n*(X_n + X_m)*alpha_m)) - 6*p**2*(-3*alpha_m**2*(5 &
					+ 16*X_n**2*alpha_m + 4*X_n**2*alpha_m + 8*X_n*X_m*&
					alpha_m + 4*X_m**2*alpha_m - 16*X_n*(X_n + X_m)*alpha_m)&
					+ 4*(-2*X_n + X_n + X_m)**2*alpha_n**3*(-3 + 4*X_n**2*&
					alpha_m + X_n**2*alpha_m + 2*X_n*X_m*alpha_m + X_m**2*&
					alpha_m - 4*X_n*(X_n + X_m)*alpha_m) + 2*alpha_n*alpha_m*&
					(-15 + 12*X_m**2*alpha_m + 32*X_n**4*alpha_m**2 + &
					2*X_n**4*alpha_m**2 + 8*X_n**3*X_m*alpha_m**2 + &
					2*X_m**4*alpha_m**2 - 64*X_n**3*(X_n + X_m)*alpha_m**2 &
					+ 12*X_n**2*alpha_m*(1 + X_m**2*alpha_m) + 8*X_n*X_m*&
					alpha_m*(3 + X_m**2*alpha_m) + 48*X_n**2*alpha_m*(1 + &
					X_n**2*alpha_m + 2*X_n*X_m*alpha_m + X_m**2*alpha_m) - &
					16*X_n*(X_n + X_m)*alpha_m*(3 + X_n**2*alpha_m + &
					2*X_n*X_m*alpha_m + X_m**2*alpha_m)) - 3*alpha_n**2*&
					(5 - 8*X_m**2*alpha_m + 64*X_n**4*alpha_m**2 + &
					4*X_n**4*alpha_m**2 + 16*X_n**3*X_m*alpha_m**2 + &
					4*X_m**4*alpha_m**2 - 128*X_n**3*(X_n + X_m)*&
					alpha_m**2 + 16*X_n*X_m*alpha_m*(-1 + X_m**2*alpha_m) &
					- 32*X_n*(X_n + X_m)*alpha_m*(-1 + X_n**2*alpha_m + &
					2*X_n*X_m*alpha_m + X_m**2*alpha_m) + 8*X_n**2*&
					alpha_m*(-1 + 3*X_m**2*alpha_m) + 32*X_n**2*alpha_m*&
					(-1 + 3*X_n**2*alpha_m + 6*X_n*X_m*alpha_m + &
					3*X_m**2*alpha_m))) + alpha_n**3*(-15 + 18*X_m**2*&
					alpha_m - 24*X_m**4*alpha_m**2 + 512*X_n**6*alpha_m**3 &
					+ 8*X_n**6*alpha_m**3 + 48*X_n**5*X_m*alpha_m**3 + &
					8*X_m**6*alpha_m**3 - 1536*X_n**5*(X_n + X_m)*&
					alpha_m**3 + 32*X_n**3*X_m*alpha_m**2*(-3 + &
					5*X_m**2*alpha_m) + 24*X_n**4*alpha_m**2*(-1 + &
					5*X_m**2*alpha_m) - 256*X_n**3*(X_n + X_m)*&
					alpha_m**2*(-3 + 5*X_n**2*alpha_m + 10*X_n*X_m*&
					alpha_m + 5*X_m**2*alpha_m) + 384*X_n**4*alpha_m**2*&
					(-1 + 5*X_n**2*alpha_m + 10*X_n*X_m*alpha_m + &
					5*X_m**2*alpha_m) + 12*X_n*X_m*alpha_m*(3 - &
					8*X_m**2*alpha_m + 4*X_m**4*alpha_m**2) + 6*X_n**2*&
					alpha_m*(3 - 24*X_m**2*alpha_m + 20*X_m**4*alpha_m**2) &
					- 24*X_n*(X_n + X_m)*alpha_m*(3 - 8*X_m**2*alpha_m + &
					4*X_n**4*alpha_m**2 + 16*X_n**3*X_m*alpha_m**2 + &
					4*X_m**4*alpha_m**2 + 16*X_n*X_m*alpha_m*(-1 + &
					X_m**2*alpha_m) + 8*X_n**2*alpha_m*(-1 + &
					3*X_m**2*alpha_m)) + 24*X_n**2*alpha_m*(3 - &
					24*X_m**2*alpha_m + 20*X_n**4*alpha_m**2 + &
					80*X_n**3*X_m*alpha_m**2 + 20*X_m**4*alpha_m**2 +&
					16*X_n*X_m*alpha_m*(-3 + 5*X_m**2*alpha_m) + &
					24*X_n**2*alpha_m*(-1 + 5*X_m**2*alpha_m))) + &
					(0,6)*p*(2*X_n - X_n - X_m)*(alpha_n - alpha_m)*&
					(2*(-2*X_n + X_n + X_m)**2*alpha_n**3 - &
					2*alpha_n*alpha_m*(-15 + 28*X_n**2*alpha_m + &
					7*X_n**2*alpha_m + 14*X_n*X_m*alpha_m + &
					7*X_m**2*alpha_m - 28*X_n*(X_n + X_m)*alpha_m) +&
					alpha_m**2*(15 + 8*X_n**2*alpha_m + &
					2*X_n**2*alpha_m + 4*X_n*X_m*alpha_m + &
					2*X_m**2*alpha_m - 8*X_n*(X_n + X_m)*alpha_m) + &
					alpha_n**2*(15 - 14*X_m**2*alpha_m + &
					64*X_n**4*alpha_m**2 + 4*X_n**4*alpha_m**2 + &
					16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*alpha_m**2&
					- 128*X_n**3*(X_n + X_m)*alpha_m**2 + &
					4*X_n*X_m*alpha_m*(-7 + 4*X_m**2*alpha_m) - &
					8*X_n*(X_n + X_m)*alpha_m*(-7 + 4*X_n**2*alpha_m &
					+ 8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m) + &
					2*X_n**2*alpha_m*(-7 + 12*X_m**2*alpha_m) + &
					8*X_n**2*alpha_m*(-7 + 12*X_n**2*alpha_m + &
					24*X_n*X_m*alpha_m + 12*X_m**2*alpha_m)))))*&
					(0.25*EXP(-(p**2 + (-2*X_n + X_n + X_m)**2*alpha_n*&
					alpha_m + (0,2)*p*(-(X_n*alpha_n) + X_n*(alpha_n &
					- alpha_m) + X_m*alpha_m))/(alpha_n + alpha_m))/&
					(alpha_n + alpha_m)**6.5)
    ENDFUNCTION   	

END MODULE