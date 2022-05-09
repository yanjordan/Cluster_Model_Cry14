MODULE EMD_Mod
    USE Reading_Mod
	!$ use OMP_LIB
	
CONTAINS

!electron density in the point p
        SUBROUTINE Density_p(p, g_num, Pop_mat, den_p)
                IMPLICIT NONE

                DOUBLEPRECISION, DIMENSION(3), INTENT(in) :: p
                
                INTEGER, INTENT(in) :: g_num
		
		integer :: i, j
                
                DOUBLEPRECISION, DIMENSION(num_total_AOh, num_total_AOh), intent(in) :: Pop_mat
				
				DOUBLE COMPLEX, DIMENSION(num_total_AOh, num_total_AOh) :: inter_term

				DOUBLE COMPLEX :: den_p
                
                DOUBLEPRECISION, DIMENSION(3) :: g_vec
                
                DOUBLE COMPLEX, DIMENSION(num_total_AOc) :: AOc
                
                DOUBLE COMPLEX, DIMENSION(num_total_AOh,1) :: AOh
                                
                DOUBLE COMPLEX, DIMENSION(1,1) :: Dens

				g_vec=matmul(cell_dir(g_num,:),lattice)
! wave function cell1                
                call AOc_p_total(p, num_total_AOc, cart_gauss, AOc)
                if (num_total_AOc==num_total_AOh) THEN
                     AOh(:,1)=AOc(:)
                ELSE
                     call AOh_p_total(num_total_AOc, num_total_AOh, mat, AOc, AOh(:,1))
                ENDIF                     

				inter_term=(0.d0,0.d0)
				
				do i=1, num_total_AOh
					do j=1, num_total_AOh
						inter_term(i,j)=exp(-Img*dot_product(p, g_vec))
					enddo
				enddo
				
                Dens=MATMUL(TRANSPOSE(AOh), MATMUL(Pop_mat*inter_term, AOh))

                den_p=Dens(1,1)

        ENDSUBROUTINE	
	
       SUBROUTINE Density_p_cluster(p, tors, denp_cluster_fin)
                IMPLICIT NONE
                               
                CHARACTER(len=1), INTENT(in) :: tors
                                
                DOUBLEPRECISION, DIMENSION(3), intent(in) :: p
				
				DOUBLE COMPLEX :: denp_cluster, denp_temps
				
				DOUBLEPRECISION :: denp_cluster_fin
                
                INTEGER :: i,j

                
                denp_cluster=(0.d0,0.d0)
				
                if(tors=='R') THEN
	                do i=1, num_cell
						call Density_p(p, i, Pop_tot(i,:,:),denp_temps)
						denp_cluster=denp_cluster+denp_temps
	                ENDDO
				
                ELSEIF(tors=='U') THEN			
	                do i=1, num_cell
						call Density_p(p, i, Pop_spin(i,:,:),denp_temps)
						denp_cluster=denp_cluster+denp_temps
	                ENDDO				
                ENDIF
				denp_cluster_fin=abs(denp_cluster)
        ENDSUBROUTINE


		SUBROUTINE test_p()
			IMPLICIT NONE
			DOUBLEPRECISION, DIMENSION(3) :: p
			DOUBLEPRECISION :: dpx, dpy, dpz, n_ele, n_ele_temps
			integer :: i, j, k
			
			n_ele=0.d0
			p=0.d0
			n_ele_temps=0.d0
			dpx=0.01d0
			!dpy=0.01d0
			!dpz=0.01d0
!$OMP PARALLEL PRIVATE(p, n_ele_temps)
!$OMP DO SCHEDULE(RUNTIME)	REDUCTION(+: n_ele)				
			do i=0,1000
				!do j=0,300
					!do k=0,300
						p= i*dpx*(/1.0,0.0,0.0/)!+ j*dpy*(/0.0,1.0,0.0/) +k*dpz*(/0.0,0.0,1.0/)
						call Density_p_cluster(p, 'U', n_ele_temps)
						n_ele=n_ele+n_ele_temps*4*pi*i*i*dpx*dpx
					!enddo
				!enddo	
			enddo
!$OMP END DO
!$OMP END PARALLEL	
			!n_ele=n_ele*dpx*dpy*dpz*8
			write(*,*) n_ele*dpx
		ENDSUBROUTINE
		
ENDMODULE	