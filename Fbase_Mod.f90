MODULE Fbase_Mod

        IMPLICIT NONE
        DOUBLEPRECISION, PARAMETER :: pi=3.141592654
        
        COMPLEX, PARAMETER :: Img=(0.d0, 1.d0)
		
		Type GTOs
			Integer :: num_gauss
			doubleprecision, dimension(10) :: alpha
			doubleprecision, dimension(10) :: coefficient
			INTEGER, dimension(3) :: l_gauss
			doubleprecision, dimension(3) :: coord_atom
		end type GTOs

	CHARACTER(len=2), DIMENSION(54), PARAMETER :: element=(/'H ','HE','LI',&
		'BE','B ','C ','N ','O ','F ','NE','NA','MG','AL','SI','P ','S ','CL',&
		'AR','K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',&
		'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR','NB','MO','TC',&
		'RU','RH','PD','AG','CD','IN','SN','SB','TE','I ','XE'/)

CONTAINS
        
        !fonction de base g
        DOUBLEPRECISION FUNCTION g(r, l1, l2, l3, alpha, c, r0)
				IMPLICIT NONE
                DOUBLEPRECISION, DIMENSION(3), intent(in) :: r, r0
                INTEGER, intent(in) :: l1, l2, l3
                DOUBLEPRECISION, intent(in) :: alpha, c

				g =  c*(r(1) - r0(1))**l1 * (r(2) - r0(2))**l2 * (r(3)-r0(3))**l3
				g = g * exp(-alpha*dot_product(r-r0, r-r0))
        ENDFUNCTION g

		!derivation n times of exp(-beta p^2) beta=1/(4 alpha)
        DOUBLEPRECISION FUNCTION F(p, alpha, n)
                IMPLICIT NONE
                DOUBLEPRECISION, intent(in) :: p, alpha
                INTEGER, intent(in) :: n
                DOUBLEPRECISION :: somme, beta, x
                INTEGER :: k, j
                
				beta = 1/(4.d0*alpha)
				x = p/1.d0
				j = floor(n/2.d0)
				somme = 0.d0
				do k = 0, j
					somme = somme + x**(n-2*k)/((-4*beta)**k * fact(k) * fact(n-2*k) * 1.d0)
				end do

				F = 2**n * fact(n) * (-beta)**n * exp(-beta*x**2) * somme
        ENDFUNCTION F
		
		! 1D  transforme fourier of  g
        DOUBLE COMPLEX FUNCTION TF(p,l,alpha,r0)
				IMPLICIT NONE
                DOUBLEPRECISION, intent(in) :: p, r0, alpha
                INTEGER, intent(in) :: l
                
                TF=Img**l * SQRT(pi/alpha) * F(p,alpha,l) * EXP(-Img*(p*r0)/1.d0)
                TF=TF * 1/SQRT(2*pi*1.d0)
        ENDFUNCTION TF		
		
        !3D transforme fourier of g
        DOUBLE COMPLEX FUNCTION TF_g(p, l1, l2, l3, alpha, c, r0)
				IMPLICIT NONE
                DOUBLEPRECISION, DIMENSION(3), intent(in) :: p, r0
				
                INTEGER, intent(in) :: l1, l2, l3
                DOUBLEPRECISION, intent(in) :: alpha, c

				TF_g = c * TF(p(1), l1, alpha, r0(1))
				TF_g = TF_g * TF(p(2), l2, alpha, r0(2))
				TF_g = TF_g * TF(p(3), l3, alpha, r0(3))
        ENDFUNCTION TF_g		
		

        ! AO postion cartesian
        SUBROUTINE AOc_total(r, num_AOc, cart_gauss, phi, dep)
				IMPLICIT NONE
                INTEGER, INTENT(in) :: num_AOc
                
                DOUBLEPRECISION, DIMENSION(3), intent(in) :: r
		
				DOUBLEPRECISION, DIMENSION(3), intent(in), optional :: dep
                
                DOUBLEPRECISION, DIMENSION(num_AOc) :: phi
                
                Type(GTOs), dimension(num_AOc), INTENT(in) :: cart_gauss
                
                
                INTEGER :: i, j
                
                phi=0.d0
		
		if (present(dep)) then
        	        do i = 1, num_AOc
	                        do j= 1 , cart_gauss(i)%num_gauss
                        	        phi(i)= phi(i) + g(r, cart_gauss(i)%l_gauss(1), &
									cart_gauss(i)%l_gauss(2), cart_gauss(i)%l_gauss(3),&
									cart_gauss(i)%alpha(j), cart_gauss(i)%coefficient(j), &
									(cart_gauss(i)%coord_atom+dep)) 
        	                ENDDO
	                ENDDO
		else                
        	        do i = 1, num_AOc
	                        do j= 1 , cart_gauss(i)%num_gauss
                        	        phi(i)= phi(i) + g(r, cart_gauss(i)%l_gauss(1),&
									cart_gauss(i)%l_gauss(2), cart_gauss(i)%l_gauss(3), &
									cart_gauss(i)%alpha(j), cart_gauss(i)%coefficient(j), &
									cart_gauss(i)%coord_atom) 
        	                ENDDO
	                ENDDO
		endif
        ENDSUBROUTINE AOc_total

        ! AO momentum cartesian
        SUBROUTINE AOc_p_total(p, num_AOc, cart_gauss, phi_p)
				IMPLICIT NONE
                INTEGER, INTENT(in) :: num_AOc
                
                DOUBLEPRECISION, DIMENSION(3), intent(in) :: p
				
                
                DOUBLE COMPLEX, DIMENSION(num_AOc) :: phi_p
                
                Type(GTOs), dimension(num_AOc), INTENT(in) :: cart_gauss
                                
                INTEGER :: i, j
                
                phi_p=(0.d0,0.d0)

                do i = 1, num_AOc
                        do j= 1 , cart_gauss(i)%num_gauss
                                phi_p(i)= phi_p(i) +  TF_g(p, cart_gauss(i)%l_gauss(1), &
								cart_gauss(i)%l_gauss(2), cart_gauss(i)%l_gauss(3),&
								cart_gauss(i)%alpha(j), cart_gauss(i)%coefficient(j), &
								cart_gauss(i)%coord_atom)
                        ENDDO
                ENDDO
        ENDSUBROUTINE AOc_p_total

         ! AO_ r pure
        SUBROUTINE AOh_total(num_AOc, num_AOh, mat, AOc, phi)
				IMPLICIT NONE
                INTEGER, INTENT(in) :: num_AOc, num_AOh
                
                DOUBLEPRECISION, DIMENSION(num_AOh) :: phi
                
                DOUBLEPRECISION, DIMENSION(num_AOc), INTENT(in) :: AOc
                
                DOUBLEPRECISION, DIMENSION(num_AOh, num_AOc), INTENT(in) :: mat
                                
                phi=0.d0
                
			phi=MATMUL(mat(:,:),AOc(:) )

        ENDSUBROUTINE AOh_total

         ! AO_p pure
        SUBROUTINE AOh_p_total(num_AOc, num_AOh, mat, AOc, phi)
		IMPLICIT NONE
                INTEGER, INTENT(in) :: num_AOc, num_AOh
                
                DOUBLE COMPLEX, DIMENSION(num_AOh) :: phi
                
                DOUBLE COMPLEX, DIMENSION(num_AOc), INTENT(in) :: AOc
                
                DOUBLEPRECISION, DIMENSION(num_AOh, num_AOc), INTENT(in) :: mat
                                
                phi=(0.d0,0.d0)
                
		phi=MATMUL(mat(:,:),AOc(:) )

        ENDSUBROUTINE AOh_p_total
        
!inv 3X3 matrix		
		SUBROUTINE inv_mat3(mat_in, mat_out)
			IMPLICIT NONE
			DOUBLEPRECISION, DIMENSION(3,3), INTENT(in) :: mat_in
			DOUBLEPRECISION, DIMENSION(3,3), INTENT(out) :: mat_out
			
			mat_out(1,1)=(mat_in(2,2)*mat_in(3,3) - mat_in(2,3)*mat_in(3,2))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) - &
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))
			
			mat_out(1,2)= -(mat_in(1,2)*mat_in(3,3) - mat_in(1,3)*mat_in(3,2))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))
			
			mat_out(1,3)=(mat_in(1,2)*mat_in(2,3) - mat_in(1,3)*mat_in(2,2))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))

			mat_out(2,1)=-(mat_in(2,1)*mat_in(3,3) - mat_in(2,3)*mat_in(3,1))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))
			
			mat_out(2,2)=(mat_in(1,1)*mat_in(3,3) - mat_in(1,3)*mat_in(3,1))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))
			
			mat_out(2,3)= -(mat_in(1,1)*mat_in(2,3) - mat_in(1,3)*mat_in(2,1))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))
			
			mat_out(3,1)=(mat_in(2,1)*mat_in(3,2) - mat_in(2,2)*mat_in(3,1))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))
			
			mat_out(3,2)= -(mat_in(1,1)*mat_in(3,2) - mat_in(1,2)*mat_in(3,1))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))

			mat_out(3,3)=(mat_in(1,1)*mat_in(2,2) - mat_in(1,2)*mat_in(2,1))&
			/(mat_in(1,1)*mat_in(2,2)*mat_in(3,3) - mat_in(1,1)*mat_in(2,3)&
			*mat_in(3,2) - mat_in(1,2)*mat_in(2,1)*mat_in(3,3) + mat_in(1,2)&
			*mat_in(2,3)*mat_in(3,1) + mat_in(1,3)*mat_in(2,1)*mat_in(3,2) -&
			mat_in(1,3)*mat_in(2,2)*mat_in(3,1))

		ENDSUBROUTINE
		
        !factoriell
        REAL FUNCTION fact(n)
				IMPLICIT NONE
                INTEGER, intent(in) :: n
                INTEGER :: i
                
                fact=1d0
                if(n>=2) THEN
                        do i=2,n
                                fact=fact*i
                        ENDDO
                ENDIF
        ENDFUNCTION fact      


        !facteur de normalisation de fonction base 
        DOUBLEPRECISION FUNCTION A(alpha, l1, l2, l3)
				IMPLICIT NONE
                DOUBLEPRECISION, intent(in) :: alpha
                INTEGER :: l1, l2, l3
                                
                A=(2.0d0*alpha/pi)**(3.0d0/4.0d0)
                A = A*sqrt((4.0d0*alpha)**(l1+l2+l3))
!		A=2**(-5d0/4d0)
                A = A*sqrt(2**(l1+l2+l3) * fact(l1)*fact(l2)*fact(l3)*1.d0)
                A = A/sqrt(fact(2*l1)*fact(2*l2)*fact(2*l3)*1.d0)
        ENDFUNCTION A  
       
       !find the row number of vector in a matrix
        INTEGER FUNCTION num_row(x, Mat, num_cell)
                IMPLICIT NONE
                
                INTEGER, intent(in) :: num_cell
                
                INTEGER, DIMENSION(3) :: x
                
                INTEGER, DIMENSION(num_cell,3) :: Mat
                
                INTEGER :: i
                
                num_row=0d0

                do i=1,num_cell
                        if (ALL(x(:) == Mat(i,:))) THEN
                                num_row=i
                                EXIT
                        ENDIF
                ENDDO

        ENDFUNCTION

        ! distance between two points DOUBLEPRECISION
        DOUBLEPRECISION FUNCTION dis(r, r1)
                IMPLICIT NONE
                DOUBLEPRECISION, DIMENSION(3), INTENT(in) :: r, r1
                
                dis=sqrt( (r(1)-r1(1))**2 +(r(2)-r1(2))**2 +(r(3)-r1(3))**2)
                
        ENDFUNCTION 
		
		function uni_vector(v) result(u)
			implicit none
			doubleprecision, dimension(:), intent(in) :: v
			
			doubleprecision, dimension(size(v)) :: u
			doubleprecision :: norm_v
			integer :: dim_v, i
			
			dim_v=size(v)
			norm_v=0.d0
			do i=1,dim_v
				norm_v=norm_v+v(i)**2
			enddo
			
			u=v/(sqrt(norm_v))
		endfunction
        
        !find the atome number
        INTEGER FUNCTION a_num(a_ele)
                IMPLICIT NONE
                
                CHARACTER(len=*) :: a_ele
                
                INTEGER :: i
                
                a_num=0d0
                
                do i=1,54
                        if (TRIM(a_ele)==TRIM(element(i))) THEN
                                a_num=i
                                EXIT
                        ENDIF
                ENDDO
                
        ENDFUNCTION
		
		!c= cross_product(a,b)
		function cross_3d(a, b) result(c)
			implicit none
			doubleprecision, dimension(3), intent(in) :: a, b
			doubleprecision, dimension(3) :: c
			
			c=0.d0
			
			c(1)=a(2)*b(3)-a(3)*b(2)
			c(2)=a(3)*b(1)-a(1)*b(3)
			c(3)=a(1)*b(2)-a(2)*b(1)
			
		endfunction


ENDMODULE Fbase_Mod