MODULE Bs_Mod
        USE Reading_Mod

CONTAINS

       DOUBLEPRECISION FUNCTION Bs_00(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_00=Sqrt(Pi)/(EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
					(alpha_n + alpha_m))*Sqrt(alpha_n + alpha_m))

        ENDFUNCTION


       DOUBLEPRECISION FUNCTION Bs_01(x,  alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_01=(Sqrt(Pi)*(x + X_n - X_m)*alpha_n)/&
				(EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**1.5)

        ENDFUNCTION
        

       DOUBLEPRECISION FUNCTION Bs_02(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_02=(Sqrt(Pi)*(alpha_n + 2*(x + X_n - X_m)**2*alpha_n**2 + alpha_m))/&
				(2.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**2.5)

        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_03(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_03=(Sqrt(Pi)*(x + X_n - X_m)*alpha_n*(alpha_n*&
				(3 + 2*(x + X_n - X_m)**2*alpha_n) + 3*alpha_m))/&
				(2.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**3.5)
        ENDFUNCTION


       DOUBLEPRECISION FUNCTION Bs_10(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_10=-((Sqrt(Pi)*(x + X_n - X_m)*alpha_m)/&
				(EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**1.5))
        ENDFUNCTION


       DOUBLEPRECISION FUNCTION Bs_11(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_11=(Sqrt(Pi)*(alpha_n + alpha_m - &
				2*(x + X_n - X_m)**2*alpha_n*alpha_m))/&
				(2.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**2.5)
        ENDFUNCTION


       DOUBLEPRECISION FUNCTION Bs_12(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_12=-(Sqrt(Pi)*(x + X_n - X_m)*(-(alpha_n*alpha_m) +&
				alpha_m**2 + 2*alpha_n**2*(-1 + (x + X_n - X_m)**2*alpha_m)))/&
				(2.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**3.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_13(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_13=(Sqrt(Pi)*(3*alpha_n**2 + 3*alpha_m**2 + &
				2*(x + X_n - X_m)**2*alpha_n**3*(3 - 2*(x + X_n - X_m)**2&
				*alpha_m) + 6*alpha_n*alpha_m*(1 - (x + X_n - X_m)**2*alpha_m)))/&
				(4.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**4.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_20(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_20=(Sqrt(Pi)*(alpha_n + alpha_m + 2*(x + X_n - X_m)**2*alpha_m**2))/&
				(2.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**2.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_21(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_21=(Sqrt(Pi)*(x + X_n - X_m)*(alpha_n**2 - 2*alpha_m**2 + &
				alpha_n*alpha_m*(-1 + 2*(x + X_n - X_m)**2*alpha_m)))/&
				(2.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**3.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_22(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_22=(Sqrt(Pi)*(2*(x + X_n - X_m)**2*alpha_n**3 + &
				6*alpha_n*alpha_m*(1 - (x + X_n - X_m)**2*alpha_m) + &
				alpha_m**2*(3 + 2*(x + X_n - X_m)**2*alpha_m) + &
				alpha_n**2*(3 - 6*(x + X_n - X_m)**2*alpha_m + &
				4*(x + X_n - X_m)**4*alpha_m**2)))/&
				(4.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**4.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_23(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_23=(Sqrt(Pi)*(x + X_n - X_m)*(alpha_n**3*(9 + 2*(x + X_n - &
				X_m)**2*alpha_n) - 2*alpha_n**2*(-6 + 5*(x + X_n - X_m)**2*alpha_n)&
				*alpha_m + alpha_n*(-3 - 6*(x + X_n - X_m)**2*alpha_n + &
				4*(x + X_n - X_m)**4*alpha_n**2)*alpha_m**2 + &
				6*(-1 + (x + X_n - X_m)**2*alpha_n)*alpha_m**3))/&
				(4.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**5.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_30(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_30=(Sqrt(Pi)*(x + X_n - X_m)*alpha_m*(-3*alpha_n + &
				alpha_m*(-3 - 2*(x + X_n - X_m)**2*alpha_m)))/&
				(2.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**3.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_31(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_31=(Sqrt(Pi)*(6*alpha_n*alpha_m - 4*(x + X_n - X_m)**4*alpha_n&
				*alpha_m**3 + 3*alpha_n**2*(1 - 2*(x + X_n - X_m)**2*alpha_m) + &
				3*alpha_m**2*(1 + 2*(x + X_n - X_m)**2*alpha_m)))/&
				(4.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**4.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_32(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_32=-(Sqrt(Pi)*(x + X_n - X_m)*(2*alpha_n*alpha_m**2*(6 - &
				5*(x + X_n - X_m)**2*alpha_m) + 6*alpha_n**3*(-1 + &
				(x + X_n - X_m)**2*alpha_m) + alpha_m**3*(9 + 2*(x + X_n - X_m)**2*alpha_m)&
				+ alpha_n**2*alpha_m*(-3 - 6*(x + X_n - X_m)**2*alpha_m + &
				4*(x + X_n - X_m)**4*alpha_m**2)))/(4.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**5.5)
        ENDFUNCTION

       DOUBLEPRECISION FUNCTION Bs_33(x, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                Bs_33=-(Sqrt(Pi)*(6*(x + X_n - X_m)**2*alpha_n**4*(-3 + 2*x**2*alpha_m &
				+ 2*X_n**2*alpha_m + 4*x*(X_n - X_m)*alpha_m - 4*X_n*X_m*alpha_m + &
				2*X_m**2*alpha_m) - 3*alpha_m**3*(5 + 6*x**2*alpha_m + 6*X_n**2*alpha_m +&
				12*x*(X_n - X_m)*alpha_m - 12*X_n*X_m*alpha_m + 6*X_m**2*alpha_m) - &
				3*alpha_n**2*alpha_m*(15 - 24*X_m**2*alpha_m + 8*x**4*alpha_m**2 + &
				8*X_n**4*alpha_m**2 + 32*x**3*(X_n - X_m)*alpha_m**2 - &
				32*X_n**3*X_m*alpha_m**2 + 8*X_m**4*alpha_m**2 + &
				16*X_n*X_m*alpha_m*(3 - 2*X_m**2*alpha_m) + &
				24*X_n**2*alpha_m*(-1 + 2*X_m**2*alpha_m) + &
				16*x*(X_n - X_m)*alpha_m*(-3 + 2*X_n**2*alpha_m - 4*X_n*X_m*alpha_m + &
				2*X_m**2*alpha_m) + 24*x**2*alpha_m*(-1 + 2*X_n**2*alpha_m - &
				4*X_n*X_m*alpha_m + 2*X_m**2*alpha_m)) + 3*alpha_n*alpha_m**2*(-15 + &
				6*X_m**2*alpha_m + 4*x**4*alpha_m**2 + 4*X_n**4*alpha_m**2 + &
				16*x**3*(X_n - X_m)*alpha_m**2 - 16*X_n**3*X_m*alpha_m**2 + &
				4*X_m**4*alpha_m**2 + 6*X_n**2*alpha_m*(1 + 4*X_m**2*alpha_m) - &
				4*X_n*X_m*alpha_m*(3 + 4*X_m**2*alpha_m) + 6*x**2*alpha_m*(1 + &
				4*X_n**2*alpha_m - 8*X_n*X_m*alpha_m + 4*X_m**2*alpha_m) + &
				4*x*(X_n - X_m)*alpha_m*(3 + 4*X_n**2*alpha_m - 8*X_n*X_m*alpha_m + &
				4*X_m**2*alpha_m)) + alpha_n**3*(-15 + 18*X_m**2*alpha_m - &
				24*X_m**4*alpha_m**2 + 8*x**6*alpha_m**3 + 8*X_n**6*alpha_m**3 + &
				48*x**5*(X_n - X_m)*alpha_m**3 - 48*X_n**5*X_m*alpha_m**3 + &
				8*X_m**6*alpha_m**3 - 32*X_n**3*X_m*alpha_m**2*(-3 + 5*X_m**2*alpha_m) + &
				24*X_n**4*alpha_m**2*(-1 + 5*X_m**2*alpha_m) + 32*x**3*(X_n - X_m)*&
				alpha_m**2*(-3 + 5*X_n**2*alpha_m - 10*X_n*X_m*alpha_m + &
				5*X_m**2*alpha_m) + 24*x**4*alpha_m**2*(-1 + 5*X_n**2*alpha_m - &
				10*X_n*X_m*alpha_m + 5*X_m**2*alpha_m) - 12*X_n*X_m*alpha_m*(3 - &
				8*X_m**2*alpha_m + 4*X_m**4*alpha_m**2) + &
				6*X_n**2*alpha_m*(3 - 24*X_m**2*alpha_m + 20*X_m**4*alpha_m**2) + &
				12*x*(X_n - X_m)*alpha_m*(3 - 8*X_m**2*alpha_m + 4*X_n**4*alpha_m**2 - &
				16*X_n**3*X_m*alpha_m**2 + 4*X_m**4*alpha_m**2 - 16*X_n*X_m*alpha_m*(-1 + &
				X_m**2*alpha_m) + 8*X_n**2*alpha_m*(-1 + 3*X_m**2*alpha_m)) + &
				6*x**2*alpha_m*(3 - 24*X_m**2*alpha_m + 20*X_n**4*alpha_m**2 - &
				80*X_n**3*X_m*alpha_m**2 + 20*X_m**4*alpha_m**2 + 16*X_n*X_m*alpha_m*(3 - &
				5*X_m**2*alpha_m) + 24*X_n**2*alpha_m*(-1 + 5*X_m**2*alpha_m)))))/&
				(8.*EXP(((x + X_n - X_m)**2*alpha_n*alpha_m)/&
				(alpha_n + alpha_m))*(alpha_n + alpha_m)**6.5)
        ENDFUNCTION
       

        DOUBLEPRECISION FUNCTION Bs_nm( x, l_n, l_m, alpha_n, alpha_m, X_n, X_m)
                IMPLICIT NONE
                
                INTEGER, INTENT(in) :: l_n, l_m
                
                DOUBLEPRECISION, INTENT(in) :: x, alpha_n, alpha_m, X_n, X_m
                
                if (l_n==0d0) THEN
	                if (l_m==0d0) THEN
                                Bs_nm=Bs_00(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==1d0) THEN
	                        Bs_nm=Bs_01(x, alpha_n, alpha_m, X_n, X_m)  
	                elseif(l_m==2d0) THEN
	                        Bs_nm=Bs_02(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==3d0) THEN
	                        Bs_nm=Bs_03(x, alpha_n, alpha_m, X_n, X_m)   
	                else
	                        WRITE(*,*) 'function Bs_nm wrong!!!'
	                        WRITE(*,*) 'l_m = ', l_m
	                        WRITE(*,*) 'the l should between 0 and 3 '
	                        STOP
	                ENDIF                        
                elseif(l_n==1d0) THEN
	                if (l_m==0d0) THEN
                                Bs_nm=Bs_10(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==1d0) THEN
	                        Bs_nm=Bs_11(x, alpha_n, alpha_m, X_n, X_m)  
	                elseif(l_m==2d0) THEN
	                        Bs_nm=Bs_12(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==3d0) THEN
	                        Bs_nm=Bs_13(x, alpha_n, alpha_m, X_n, X_m)   
	                else
	                        WRITE(*,*) 'function Bs_nm wrong!!!'
	                        WRITE(*,*) 'l_m = ', l_m
	                        WRITE(*,*) 'the l should between 0 and 3 '
	                        STOP
	                ENDIF                 
                elseif(l_n==2d0) THEN
	                if (l_m==0d0) THEN
                                Bs_nm=Bs_20(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==1d0) THEN
	                        Bs_nm=Bs_21(x, alpha_n, alpha_m, X_n, X_m)  
	                elseif(l_m==2d0) THEN
	                        Bs_nm=Bs_22(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==3d0) THEN
	                        Bs_nm=Bs_23(x, alpha_n, alpha_m, X_n, X_m)   
	                else
	                        WRITE(*,*) 'function Bs_nm wrong!!!'
	                        WRITE(*,*) 'l_m = ', l_m
	                        WRITE(*,*) 'the l should between 0 and 3 '
	                        STOP
	                ENDIF                 
                elseif(l_n==3d0) THEN
	                if (l_m==0d0) THEN
                                Bs_nm=Bs_30(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==1d0) THEN
	                        Bs_nm=Bs_31(x, alpha_n, alpha_m, X_n, X_m)  
	                elseif(l_m==2d0) THEN
	                        Bs_nm=Bs_32(x, alpha_n, alpha_m, X_n, X_m) 
	                elseif(l_m==3d0) THEN
	                        Bs_nm=Bs_33(x, alpha_n, alpha_m, X_n, X_m)   
	                else
	                        WRITE(*,*) 'function Bs_nm wrong!!!'
	                        WRITE(*,*) 'l_m = ', l_m
	                        WRITE(*,*) 'the l should between 0 and 3 '
	                        STOP
	                ENDIF                 
                else
                        WRITE(*,*) 'function Bs_nm wrong!!!'
                        WRITE(*,*) 'l_n = ', l_n
                        WRITE(*,*) 'the l should between 0 and 3 '
                        STOP
                ENDIF
      
        ENDFUNCTION
			
		doubleprecision function Bs_nm_3D(x, num_bs_n, num_bs_m, g_num)
			IMPLICIT NONE
			doubleprecision, dimension(3), intent(in) :: x
			integer, intent(in) :: num_bs_n, num_bs_m
			integer, intent(in), optional :: g_num
			doubleprecision, dimension(3) :: dep
			integer :: i, j
			
			Bs_nm_3D=0.d0
			
			if (present(g_num)) then
				dep=cell_dir(g_num,1)*lattice(1,:)+cell_dir(g_num,2)*lattice(2,:)+cell_dir(g_num,3)*lattice(3,:)
				do i=1, cart_gauss(num_bs_n)%num_gauss
					do j=1, cart_gauss(num_bs_m)%num_gauss
						Bs_nm_3D=Bs_nm_3D + cart_gauss(num_bs_n)%coefficient(i)* &
						cart_gauss(num_bs_m)%coefficient(j) * &
						Bs_nm(x(1),  cart_gauss(num_bs_n)%l_gauss(1), cart_gauss(num_bs_m)%l_gauss(1), &
						cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
						cart_gauss(num_bs_n)%coord_atom(1), cart_gauss(num_bs_m)%coord_atom(1)+dep(1))&
						* &
						Bs_nm(x(2),  cart_gauss(num_bs_n)%l_gauss(2), cart_gauss(num_bs_m)%l_gauss(2), &
						cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
						cart_gauss(num_bs_n)%coord_atom(2), cart_gauss(num_bs_m)%coord_atom(2)+dep(2))&
						* &
						Bs_nm(x(3),  cart_gauss(num_bs_n)%l_gauss(3), cart_gauss(num_bs_m)%l_gauss(3), &
						cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
						cart_gauss(num_bs_n)%coord_atom(3), cart_gauss(num_bs_m)%coord_atom(3)+dep(3))
					enddo
				enddo			
			else
				do i=1, cart_gauss(num_bs_n)%num_gauss
					do j=1, cart_gauss(num_bs_m)%num_gauss
						Bs_nm_3D=Bs_nm_3D + cart_gauss(num_bs_n)%coefficient(i)* &
						cart_gauss(num_bs_m)%coefficient(j) * &
						Bs_nm(x(1),  cart_gauss(num_bs_n)%l_gauss(1), cart_gauss(num_bs_m)%l_gauss(1), &
						cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
						cart_gauss(num_bs_n)%coord_atom(1), cart_gauss(num_bs_m)%coord_atom(1))&
						* &
						Bs_nm(x(2),  cart_gauss(num_bs_n)%l_gauss(2), cart_gauss(num_bs_m)%l_gauss(2), &
						cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
						cart_gauss(num_bs_n)%coord_atom(2), cart_gauss(num_bs_m)%coord_atom(2))&
						* &
						Bs_nm(x(3),  cart_gauss(num_bs_n)%l_gauss(3), cart_gauss(num_bs_m)%l_gauss(3), &
						cart_gauss(num_bs_n)%alpha(i), cart_gauss(num_bs_m)%alpha(j),&
						cart_gauss(num_bs_n)%coord_atom(3), cart_gauss(num_bs_m)%coord_atom(3))
					enddo
				enddo
			endif
		endfunction


END MODULE