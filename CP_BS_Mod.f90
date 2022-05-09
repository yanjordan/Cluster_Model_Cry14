MODULE CP_BS_Mod
	USE Bs_Mod
	!$ use OMP_LIB
	
CONTAINS

	!only the center cell
	subroutine BIDIERD(ndir, npu, step, npoip, stepc, dirs)
		IMPLICIT NONE
		integer, intent(in) :: ndir, npu, npoip
		doubleprecision, intent(in) :: step, stepc
		doubleprecision, dimension(ndir, 3) :: dirs
		
		integer :: i,j
		doubleprecision, dimension(:,:), allocatable :: Bdir_t, Bdir_s
		doubleprecision, dimension(:,:), allocatable :: CPdir_t, CPdir_s
		doubleprecision, dimension(3) :: v, s
		CHARACTER(len=2) :: Num_col
		
		allocate(Bdir_t(npu,ndir), CPdir_t(npoip,ndir), Pop_tot_car(1,num_total_AOc, num_total_AOc))
		Pop_tot_car(1,:,:)=matmul(transpose(mat), matmul(Pop_tot(1,:,:), mat))
		if (uorr=='U') THEN
			ALLOCATE(Pop_spin_car(1, num_total_AOc, num_total_AOc))		
			Pop_spin_car(1,:,:)=matmul(transpose(mat), matmul(Pop_spin(1,:,:), mat))
        ENDIF

		Bdir_t=0.d0
		CPdir_t=0.d0
		write( Num_col,'(i2)') ndir
		
		if (uorr=='R') then
			do i=1, ndir
				v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
!$OMP PARALLEL PRIVATE(j,s)
!$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(/0.d0,0.d0,0.d0/)+(j-1)*v*step
					call Bs_restricted(s, Bdir_t(j,i))
				enddo
!$OMP END DO
!$OMP END PARALLEL	
!normalisation B(0) --> electron number
				Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))		
			enddo
			
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)
			
			open(45, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(45)
			
		else
			allocate(Bdir_s(npu,ndir), CPdir_s(npoip,ndir))
			CPdir_s=0.d0
			Bdir_s=0.d0
			do i=1, ndir
				v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
!$OMP PARALLEL PRIVATE(j,s)
!$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(/0.d0,0.d0,0.d0/)+(j-1)*v*step
					call Bs_unrestricted(s, Bdir_t(j,i), Bdir_s(j,i))
				enddo
!$OMP END DO
!$OMP END PARALLEL	
!normalisation B(0) --> electron number
				Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))
				Bdir_s(:,i)=Bdir_s(:,i)*NINT(Bdir_s(1,i))/(Bdir_s(1,i))				
			enddo
			
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_s, CPdir_s)
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)

			open(41, file='B_spin.dat', status='replace')
			do i=1,npu
				write(41, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_s(i,:)
			enddo
			close(41)
			
			open(42, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(42, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(42)	
			
			open(45, file='CP_spin.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_s(i,:)
			enddo
			close(45)	
		deallocate(Bdir_s, CPdir_s, Pop_spin_car)
		endif	
		deallocate(Bdir_t, CPdir_t, Pop_tot_car)
	endsubroutine
	
	!cluster contribution
	subroutine BIDIERD_Cluster(ndir, npu, step, npoip, stepc, dirs)
		IMPLICIT NONE
		integer, intent(in) :: ndir, npu, npoip
		doubleprecision, intent(in) :: step, stepc
		doubleprecision, dimension(ndir, 3) :: dirs
		
		integer :: i,j,k,m
		doubleprecision, dimension(:,:), allocatable :: Bdir_t, Bdir_s
		doubleprecision, dimension(:,:), allocatable :: CPdir_t, CPdir_s
		doubleprecision, dimension(3) :: v, s
		
		doubleprecision :: bs_ij_t, bs_ij_s
		
		CHARACTER(len=2) :: Num_col
		
		allocate(Bdir_t(npu,ndir), CPdir_t(npoip,ndir), Pop_tot_car(num_cell,num_total_AOc, num_total_AOc))
	
		do i=1, num_cell
			Pop_tot_car(i,:,:)=matmul(transpose(mat), matmul(Pop_tot(i,:,:), mat))
		enddo
		
		if (uorr=='U') THEN
			ALLOCATE(Pop_spin_car(num_cell, num_total_AOc, num_total_AOc))
			
			do i=1, num_cell
				Pop_spin_car(i,:,:)=matmul(transpose(mat), matmul(Pop_spin(i,:,:), mat))
			enddo	
			
        ENDIF
		
		Bdir_t=0.d0
		CPdir_t=0.d0
		write( Num_col,'(i2)') ndir
		
		if (uorr=='R') then
			do i=1, ndir
				!v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
				v=uni_vector(dirs(i,:))
!$OMP PARALLEL PRIVATE(j,s,bs_ij_t)
!$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(j-1)*v*step
					do k=1, num_cell
						call Bs_restricted(s, bs_ij_t, k)
						Bdir_t(j,i)=Bdir_t(j,i)+bs_ij_t
					enddo
				enddo
!$OMP END DO
!$OMP END PARALLEL	
!normalisation B(0) --> electron number
				Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))		
			enddo
			
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)
			
			open(45, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(45)
			
		else
			allocate(Bdir_s(npu,ndir), CPdir_s(npoip,ndir))
			CPdir_s=0.d0
			Bdir_s=0.d0
			do i=1, ndir
				!v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
				v=uni_vector(dirs(i,:))
 !$OMP PARALLEL PRIVATE(j,s,bs_ij_t,bs_ij_s)
 !$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(j-1)*v*step
					do k=1, num_cell
						call Bs_unrestricted(s, bs_ij_t, bs_ij_s, k)
						Bdir_t(j,i)=Bdir_t(j,i)+bs_ij_t
						Bdir_s(j,i)=Bdir_s(j,i)+bs_ij_s
					enddo
				enddo
 !$OMP END DO
 !$OMP END PARALLEL	
 !normalisation B(0) --> electron number
				Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))
				Bdir_s(:,i)=Bdir_s(:,i)*NINT(Bdir_s(1,i))/(Bdir_s(1,i))	
			enddo
						
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_s, CPdir_s)
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)

			open(41, file='B_spin.dat', status='replace')
			do i=1,npu
				write(41, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_s(i,:)
			enddo
			close(41)
			
			open(42, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(42, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(42)	
			
			open(45, file='CP_spin.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_s(i,:)
			enddo
			close(45)	
			deallocate(Bdir_s, CPdir_s, Pop_spin_car)
		endif	
		deallocate(Bdir_t, CPdir_t, Pop_tot_car)
	endsubroutine	

!cluster contribution (separation orbitals)
	subroutine SBIDIERD_Cluster(ndir, npu, step, npoip, stepc, dirs, num_orbital, list_orbital)
		IMPLICIT NONE
		integer, intent(in) :: ndir, npu, npoip, num_orbital
		INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
		doubleprecision, intent(in) :: step, stepc
		doubleprecision, dimension(ndir, 3) :: dirs
		
		integer :: i,j,k,m
		doubleprecision, dimension(:,:), allocatable :: Bdir_t, Bdir_s
		doubleprecision, dimension(:,:), allocatable :: CPdir_t, CPdir_s
		
		doubleprecision, dimension(num_total_AOh, num_total_AOh) :: pop_temps
		
		doubleprecision, dimension(3) :: v, s
		
		doubleprecision :: bs_ij_t, bs_ij_s
		
		CHARACTER(len=2) :: Num_col
		
		allocate(Bdir_t(npu,ndir), CPdir_t(npoip,ndir), Pop_tot_car(num_cell,num_total_AOc, num_total_AOc))
	
		do i=1, num_cell
			pop_temps=0.d0
			pop_temps(list_orbital(:),list_orbital(:))=Pop_tot(i,list_orbital(:),list_orbital(:))
			Pop_tot_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
		enddo
		
		if (uorr=='U') THEN
			ALLOCATE(Pop_spin_car(num_cell, num_total_AOc, num_total_AOc))
			
			do i=1, num_cell
				pop_temps=0.d0
				pop_temps(list_orbital(:),list_orbital(:))=Pop_spin(i,list_orbital(:),list_orbital(:))
				Pop_spin_car(i,:,:)=matmul(transpose(mat), matmul(pop_temps, mat))
			enddo	
	
        	ENDIF
		
		Bdir_t=0.d0
		CPdir_t=0.d0
		write( Num_col,'(i2)') ndir
		
		if (uorr=='R') then
			do i=1, ndir
				!v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
				v=uni_vector(dirs(i,:))
!$OMP PARALLEL PRIVATE(j,s,bs_ij_t)
!$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(j-1)*v*step
					do k=1, num_cell
						call Bs_restricted(s, bs_ij_t, k)
						Bdir_t(j,i)=Bdir_t(j,i)+bs_ij_t
					enddo
				enddo
!$OMP END DO
!$OMP END PARALLEL	
!normalisation B(0) --> electron number
				!Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))		
			enddo
			
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)
			
			open(45, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(45)
			
		else
			allocate(Bdir_s(npu,ndir), CPdir_s(npoip,ndir))
			CPdir_s=0.d0
			Bdir_s=0.d0
			do i=1, ndir
				!v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
				v=uni_vector(dirs(i,:))
 !$OMP PARALLEL PRIVATE(j,s,bs_ij_t,bs_ij_s)
 !$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(j-1)*v*step
					do k=1, num_cell
						call Bs_unrestricted(s, bs_ij_t, bs_ij_s, k)
						Bdir_t(j,i)=Bdir_t(j,i)+bs_ij_t
						Bdir_s(j,i)=Bdir_s(j,i)+bs_ij_s
					enddo
				enddo
 !$OMP END DO
 !$OMP END PARALLEL	
 !normalisation B(0) --> electron number
				!Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))
				!Bdir_s(:,i)=Bdir_s(:,i)*NINT(Bdir_s(1,i))/(Bdir_s(1,i))	
			enddo
						
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_s, CPdir_s)
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)

			open(41, file='B_spin.dat', status='replace')
			do i=1,npu
				write(41, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_s(i,:)
			enddo
			close(41)
			
			open(42, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(42, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(42)	
			
			open(45, file='CP_spin.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_s(i,:)
			enddo
			close(45)	
			deallocate(Bdir_s, CPdir_s, Pop_spin_car)
		endif	
		deallocate(Bdir_t, CPdir_t, Pop_tot_car)
	endsubroutine	

! cross terim
!cluster contribution (separation orbitals)
	subroutine CBIDIERD_Cluster(ndir, npu, step, npoip, stepc, dirs, num_orbital, list_orbital, num_orbital1, list_orbital1)
		IMPLICIT NONE
		integer, intent(in) :: ndir, npu, npoip, num_orbital, num_orbital1
		INTEGER, DIMENSION(num_orbital), INTENT(in) :: list_orbital
		INTEGER, DIMENSION(num_orbital1), INTENT(in) :: list_orbital1
		doubleprecision, intent(in) :: step, stepc
		doubleprecision, dimension(ndir, 3) :: dirs
		
		integer :: i,j,k,m
		doubleprecision, dimension(:,:), allocatable :: Bdir_t, Bdir_s
		doubleprecision, dimension(:,:), allocatable :: CPdir_t, CPdir_s
		
		doubleprecision, dimension(num_total_AOh, num_total_AOh) :: pop_temps
		
		doubleprecision, dimension(3) :: v, s
		
		doubleprecision :: bs_ij_t, bs_ij_s
		
		CHARACTER(len=2) :: Num_col
		
		allocate(Bdir_t(npu,ndir), CPdir_t(npoip,ndir), Pop_tot_car(num_cell,num_total_AOc, num_total_AOc))
	
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
		
		Bdir_t=0.d0
		CPdir_t=0.d0
		write( Num_col,'(i2)') ndir
		
		if (uorr=='R') then
			do i=1, ndir
				!v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
				v=uni_vector(dirs(i,:))
!$OMP PARALLEL PRIVATE(j,s,bs_ij_t)
!$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(j-1)*v*step
					do k=1, num_cell
						call Bs_restricted(s, bs_ij_t, k)
						Bdir_t(j,i)=Bdir_t(j,i)+bs_ij_t
					enddo
				enddo
!$OMP END DO
!$OMP END PARALLEL	
!normalisation B(0) --> electron number
				!Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))		
			enddo
			
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)
			
			open(45, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(45)
			
		else
			allocate(Bdir_s(npu,ndir), CPdir_s(npoip,ndir))
			CPdir_s=0.d0
			Bdir_s=0.d0
			do i=1, ndir
				!v=dirs(i,:)/dis((/0.d0,0.d0,0.d0/), dirs(i,:))
				v=uni_vector(dirs(i,:))
 !$OMP PARALLEL PRIVATE(j,s,bs_ij_t,bs_ij_s)
 !$OMP DO SCHEDULE(RUNTIME)				
				do j=1, npu
					s=(j-1)*v*step
					do k=1, num_cell
						call Bs_unrestricted(s, bs_ij_t, bs_ij_s, k)
						Bdir_t(j,i)=Bdir_t(j,i)+bs_ij_t
						Bdir_s(j,i)=Bdir_s(j,i)+bs_ij_s
					enddo
				enddo
 !$OMP END DO
 !$OMP END PARALLEL	
 !normalisation B(0) --> electron number
				!Bdir_t(:,i)=Bdir_t(:,i)*NINT(Bdir_t(1,i))/(Bdir_t(1,i))
				!Bdir_s(:,i)=Bdir_s(:,i)*NINT(Bdir_s(1,i))/(Bdir_s(1,i))	
			enddo
						
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_t, CPdir_t)
			call Bs2Cp(ndir, npu, step, npoip, stepc, Bdir_s, CPdir_s)
			
			open(40, file='B_total.dat', status='replace')
			do i=1,npu
				write(40, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_t(i,:)
			enddo
			close(40)

			open(41, file='B_spin.dat', status='replace')
			do i=1,npu
				write(41, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*step, Bdir_s(i,:)
			enddo
			close(41)
			
			open(42, file='CP_total.dat', status='replace')
			do i=1,npoip
				write(42, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_t(i,:)
			enddo
			close(42)	
			
			open(45, file='CP_spin.dat', status='replace')
			do i=1,npoip
				write(45, '(F6.2, '//trim(Num_col)//'F20.12)') (i-1)*stepc, CPdir_s(i,:)
			enddo
			close(45)	
			deallocate(Bdir_s, CPdir_s, Pop_spin_car)
		endif	
		deallocate(Bdir_t, CPdir_t, Pop_tot_car)
	endsubroutine
	
	subroutine Bs2Cp(ndir, npu, step, npoip, stepc, Bdir, CPdir)
		IMPLICIT NONE
		integer, intent(in) :: ndir, npu, npoip
		doubleprecision, intent(in) :: step, stepc
		doubleprecision, dimension(npu,ndir), intent(in) :: Bdir
		doubleprecision, dimension(npoip,ndir), intent(out) :: CPdir
		
		integer :: i,j, k
		
		CPdir=0.d0
		
		do i=1, ndir
			do j=1, npoip
				CPdir(j,i) = CPdir(j,i)+ Bdir(1,i)*0.5d0
				do k=2, npu
					CPdir(j,i) = CPdir(j,i)+ Bdir(k,i)*cos((j-1)*stepc*step*(k-1))
				enddo
			enddo
		enddo
		
		CPdir=CPdir*step/pi		
	endsubroutine
	
	subroutine Bs_restricted(x, Bs_total, g_num)
		IMPLICIT NONE
		doubleprecision, dimension(3), intent(in) :: x
		integer, intent(in), optional :: g_num
		doubleprecision :: Bs_total
		integer :: i, j
		
		Bs_total=0.d0
		if (present(g_num)) then
			do i=1, num_total_AOc
				do j=1, num_total_AOc
					Bs_total=Bs_total+Bs_nm_3D(x,i,j, g_num)*Pop_tot_car(g_num,i,j)
				enddo
			enddo			
		else
			do i=1, num_total_AOc
				do j=1, num_total_AOc
					Bs_total=Bs_total+Bs_nm_3D(x,i,j)*Pop_tot_car(1d0,i,j)
				enddo
			enddo			
			
		endif

	endsubroutine
	
	subroutine Bs_unrestricted(x, Bs_total, Bs_spin, g_num)
		IMPLICIT NONE
		doubleprecision, dimension(3), intent(in) :: x
		integer, intent(in), optional :: g_num
		doubleprecision :: Bs_total, Bs_spin, Bs_nm_temps
		integer :: i, j
		
		Bs_total=0.d0
		Bs_spin=0.d0
		
		if (present(g_num)) then
			do i=1, num_total_AOc
				do j=1, num_total_AOc
					Bs_nm_temps=Bs_nm_3D(x,i,j,g_num)
					Bs_total=Bs_total+Bs_nm_temps*Pop_tot_car(g_num,i,j)
					Bs_spin=Bs_spin+Bs_nm_temps*Pop_spin_car(g_num,i,j)
				enddo
			enddo			
		else
			do i=1, num_total_AOc
				do j=1, num_total_AOc
					Bs_nm_temps=Bs_nm_3D(x,i,j,g_num)
					Bs_total=Bs_total+Bs_nm_temps*Pop_tot_car(1d0,i,j)
					Bs_spin=Bs_spin+Bs_nm_temps*Pop_spin_car(1d0,i,j)
				enddo
			enddo
		endif

	endsubroutine
	
ENDMODULE