MODULE Clementi_DM
        USE Reading_Mod
        
CONTAINS


        DOUBLEPRECISION FUNCTION norm_cle_DM(alpha, n)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: alpha
                INTEGER, INTENT(in) :: n
                norm_cle_DM=0.d0
                
                norm_cle_DM=sqrt((2*alpha)**(2*n+1) /fact(2*n))

        
        ENDFUNCTION

        DOUBLEPRECISION FUNCTION DM_H(r,r1)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: r,r1
                DM_H=0.d0
                
!               DM_H=1d0*(norm_cle_DM(1d0,1)*1d0*exp(-1d0*r))**2
                DM_H=4*exp(-r-r1)
                
                
        ENDFUNCTION
        
        DOUBLEPRECISION FUNCTION DM_C(r,r1)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: r,r1
                DOUBLEPRECISION :: wf, wf1
                
                DOUBLEPRECISION, DIMENSION(6) :: alphas, sc1, sc2
                DOUBLEPRECISION, DIMENSION(4) :: alphap, pc1
                INTEGER, DIMENSION(6) :: sn
                INTEGER, DIMENSION(4) :: pn
                
                INTEGER :: i
                
                alphas=(/5.43599, 9.48256, 1.05749, 1.52427, 2.68435, 4.20096/)
                sc1=(/ .93262, .06931, .00083, -.00176, .00559, .00382/)
                sc2=(/-.20814, -.01071, .08099, .75045, .33549, -.14765/)
                sn=(/1, 1, 2, 2, 2, 2/)

                
                alphap=(/0.98073, 1.44361, 2.60051, 6.51003/)
                pc1=(/.28241, .54697, .23195, .01025/)
                pn=(/2, 2, 2, 2/)
                
                
                DM_C=0.d0
                
                
                wf=0.d0
				wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_C=DM_C+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_C=DM_C+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,4
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_C=DM_C+2*wf*wf1
        ENDFUNCTION
        
        DOUBLEPRECISION FUNCTION DM_N(r,r1)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: r,r1
                DOUBLEPRECISION :: wf,wf1
                
                DOUBLEPRECISION, DIMENSION(6) :: alphas, sc1, sc2
                DOUBLEPRECISION, DIMENSION(4) :: alphap, pc1
                INTEGER, DIMENSION(6) :: sn
                INTEGER, DIMENSION(4) :: pn
                
                INTEGER :: i
                
                alphas=(/ 6.45739, 11.17200, 1.36405, 1.89734, 3.25291, 5.08238/)
                sc1=(/  .93780, .05849, .00093, -.00170, .00574, .00957/)
                sc2=(/-.21677, -.00846, .17991, .67416, .31297, -.14497/)
                sn=(/1, 1, 2, 2, 2, 2/)

                
                alphap=(/1.16068, 1.70472, 3.03935, 7.17482/)
                pc1=(/.26639, .52319, .27353, .01292/)
                pn=(/2, 2, 2, 2/)
                
                
                DM_N=0.d0
                
                
                wf=0.d0
				wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_N=DM_N+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_N=DM_N+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,4
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_N=DM_N+3*wf*wf1
        ENDFUNCTION
        

        DOUBLEPRECISION FUNCTION DM_O(r,r1)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: r, r1
                DOUBLEPRECISION :: wf, wf1
                
                DOUBLEPRECISION, DIMENSION(6) :: alphas, sc1, sc2
                DOUBLEPRECISION, DIMENSION(4) :: alphap, pc1
                INTEGER, DIMENSION(6) :: sn
                INTEGER, DIMENSION(4) :: pn
                
                INTEGER :: i
                
                alphas=(/ 7.61413, 13.75740, 1.69824, 2.48022, 4.31196, 5.86596/)
                sc1=(/ .94516, .03391, -.00034, .00241, -.00486, .03681/)
                sc2=(/-.22157, -.00476, .34844, .60807, .25365, -.19183/)
                sn=(/1, 1, 2, 2, 2, 2/)

                
                alphap=(/1.14394, 1.81730, 3.44988, 7.56484/)
                pc1=(/.16922, .57974, .32352, .01660/)
                pn=(/2, 2, 2, 2/)
                
                
                DM_O=0.d0
                
                
                wf=0.d0
				wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_O=DM_O+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_O=DM_O+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,4
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_O=DM_O+4*wf*wf1
        ENDFUNCTION        
        
        DOUBLEPRECISION FUNCTION DM_S(r,r1)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: r,r1
                DOUBLEPRECISION :: wf,wf1
                
                DOUBLEPRECISION, DIMENSION(8) :: alphas, sc1, sc2, sc3, alphap, pc1, pc2
                INTEGER, DIMENSION(8) :: sn,pn

                INTEGER :: i
                
                alphas=(/ 16.00000, 18.12080,12.81680,9.62700, 6.41240, 3.77420, 2.70950, 1.67090/)
                sc1=(/.97218, .02807, .01574,-.00142, .00083,-.00011, .00001,-.00001/)
                sc2=(/-.26848,-.00207, .07475, .35394, .61978, .04624,-.00873, .00157/)
                sc3=(/ .07790, .00204,-.03170,-.09454,-.25461,-.07603, .65648, .51654/)
                sn=(/1,3,3,3,3,3,3,3/)

                
                alphap=(/8.00000,14.12590,11.57180, 7.91710, 5.60390, 2.89060, 1.62750, 0.86650/)
                pc1=(/  .61338,  .00712,  .02562,  .28388,  .18293,  .00529, -.00009,  .00016/)
                pc2=(/-.15546,-.00995, .01287,-.09992, .05174, .55282, .50896, .02556/)
                pn=(/2,4,4,4,4,4,4,4/)
                
                
                DM_S=0.d0
                
                
                wf=0.d0
				wf1=0.d0
                do i=1,8
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_S=DM_S+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,8
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_S=DM_S+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,8
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc3(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc3(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_S=DM_S+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,8
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_S=DM_S+6*wf*wf1
                
                 wf=0.d0
				 wf1=0.d0
                do i=1,8
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc2(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc2(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_S=DM_S+4*wf*wf1
        ENDFUNCTION

        DOUBLEPRECISION FUNCTION DM_Y(r,r1)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: r,r1
                DOUBLEPRECISION :: wf,wf1
                
                DOUBLEPRECISION, DIMENSION(13) :: alphas, sc1, sc2, sc3, sc4, sc5
                INTEGER, DIMENSION(13) :: sn
                
                DOUBLEPRECISION, DIMENSION(10) :: alphap, pc1, pc2, pc3
                INTEGER, DIMENSION(10) :: pn
                
                DOUBLEPRECISION, DIMENSION(8) :: alphad, dc1, dc2
                INTEGER, DIMENSION(8) :: dn

                INTEGER :: i
                
                alphas=(/60.664173,50.975223,40.935841,24.747326,18.046217,18.044325,10.201763, 7.100605, &
                        6.098655, 3.847410, 2.856269, 2.102583, 0.664503/)
                sc1=(/-0.0176001,-0.0316712,-0.9176934, 0.1149616,-0.2122013, 0.0675445,-0.0100817, &
                        0.0047046,-0.0011301, 0.0011457,-0.0013048, 0.0001949,-0.0000028/)
                sc2=(/ 0.0221354, 0.0359175, 0.2909009,-0.3201723,-0.2528916,-0.4305314,-0.1012115, &
                        0.0270250,-0.0066148, 0.0070465,-0.0079563, 0.0012401,-0.0000181/)
                sc3=(/ 0.0105274, 0.0175844, 0.1198556,-0.1108259,-0.1311709,-0.3614456,-0.7066031, &
                        1.6281982, 0.1456418, 0.0061994,-0.0025391, 0.0013244,-0.0000085/)
                sc4=(/-0.0041604,-0.0070581,-0.0441032, 0.0330659, 0.0629593, 0.1601787, 0.2040213,&
                        -0.7545000,-0.7567802,-0.6203207, 1.7680472, 0.5018628, 0.0001533/)
                sc5=(/ 0.001077, 0.001807, 0.011706,-0.009415,-0.015542,-0.042224,-0.058771, &
                        0.110009, 0.075629,-0.448223, 0.784950,-0.889130, 1.174469/)
                sn=(/ 2, 3, 1, 3, 2, 3, 2, 2, 3, 3, 2, 2, 2/)

                
                alphap=(/51.271036,44.352773,27.799532,19.823485,16.018534,11.183812, &
                        7.826405, 5.241793, 2.098300, 1.281011/)
                pc1=(/ 0.0152145, 0.0219494, 0.1311566,-0.9467697,-0.2390735,-0.0102678,&
                        -0.0015598, 0.0000634, 0.0000101,-0.0000002/)
                pc2=(/ 0.0170185, 0.0238146, 0.1213548,-0.5297530,-0.1814548, 0.1240941,&
                         0.5661652, 0.4784007, 0.0019550,-0.0003711/)
                pc3=(/-0.008785,-0.012265,-0.060639, 0.232584, 0.093393, 0.029127, &
                        0.042259,-0.885882, 1.214234, 0.086854/)
                pn=(/ 2, 3, 3, 2, 3, 3, 3, 2, 2, 2/)
                
                alphad=(/20.763265,11.356318, 8.671457, 5.239435, 3.167295, 1.989968, 1.280570, 0.850768/)
                dc1=(/ 0.0161944, 0.4085581, 0.3277136, 0.3385215, 0.0164784,-0.0010850, 0.0003789,-0.0000529/)
                dc2=(/-0.003205 ,-0.088972 ,-0.058387 ,-0.130652 , 0.137859 , 0.390520 , 0.418799 , 0.2159011 /)
                dn=(/ 3, 3, 4, 3, 3, 3, 3, 3/)
                
                
                DM_Y=0.d0
                
! s orbital                
                wf=0.d0
				wf1=0.d0
                do i=1,13
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Y=DM_Y+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,13
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Y=DM_Y+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,13
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc3(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc3(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Y=DM_Y+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,13
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc4(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc4(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Y=DM_Y+2*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,13
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc5(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc5(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Y=DM_Y+2*wf*wf1
! p orbital                
                wf=0.d0
				wf1=0.d0
                do i=1,10
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_Y=DM_Y+6*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,10
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc2(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc2(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_Y=DM_Y+6*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,10
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc3(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc3(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_Y=DM_Y+6*wf*wf1
! d orbital                
                wf=0.d0
				wf1=0.d0
                do i=1,8
                        wf=wf+norm_cle_DM(alphad(i),dn(i))*r**(dn(i)-1)*dc1(i)*exp(-alphad(i)*r)
						wf1=wf1+norm_cle_DM(alphad(i),dn(i))*r1**(dn(i)-1)*dc1(i)*exp(-alphad(i)*r1)
                ENDDO
                DM_Y=DM_Y+10*wf*wf1
                
                wf=0.d0
				wf1=0.d0
                do i=1,8
                        wf=wf+norm_cle_DM(alphad(i),dn(i))*r**(dn(i)-1)*dc2(i)*exp(-alphad(i)*r)
						wf1=wf1+norm_cle_DM(alphad(i),dn(i))*r1**(dn(i)-1)*dc2(i)*exp(-alphad(i)*r1)
                ENDDO
                DM_Y=DM_Y+1*wf*wf1
 
        ENDFUNCTION

        
        DOUBLEPRECISION FUNCTION DM_Ti(r,r1)
		IMPLICIT NONE
                DOUBLEPRECISION, INTENT(in) :: r,r1
                DOUBLEPRECISION :: wf,wf1
                
                DOUBLEPRECISION, DIMENSION(11) :: alphas, sc1, sc2, sc3, sc4
                INTEGER, DIMENSION(11) :: sn
                
                DOUBLEPRECISION, DIMENSION(6) :: alphap, pc1, pc2
                INTEGER, DIMENSION(6) :: pn
                
                DOUBLEPRECISION, DIMENSION(5) :: alphad, dc1
                INTEGER, DIMENSION(5) :: dn

                INTEGER :: i
                
                alphas=(/ 21.78320,34.57420,18.70590, 9.30981, 8.23718, 4.83162, 3.37930, 3.48770, 1.53912, 0.93846, 0.70026/)
                sc1=(/ .95029, .02149, .03563, .00444,-.00286, .00338,-.00487, .00267,-.00013, .00011,-.00006/)
                sc2=(/-.29231,-.00057,-.15767, .93686, .22976,-.00795, .02460,-.01355, .00082,-.00072, .00036/)
                sc3=(/ .10649,-.00084, .05797,-.35705,-.25133, .54343, .48850, .17860, .01333,-.00925, .00438/)
                sc4=(/-.02379,-.00029,-.01531, .08760, .05373,-.06232,-.33674, .12273, .56471, .49596, .02654/)
                sn=(/ 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4/)

                
                alphap=(/10.00560,16.79230, 8.32711, 4.73525, 2.92461, 1.34176/)
                pc1=(/ .69283, .05359, .29944, .01300,-.00116, .00017/)
                pc2=(/-.23557,-.01963,-.16304, .39692, .71972, .03091/)
                pn=(/2,2,3,3,3,3/)
                
                alphad=(/2.00850,9.10000,4.40390,3.42000,1.10460/)
                dc1=(/ .55179, .03278, .25284, .17865, .14505/)
                dn=(/3,3,3,3,3/)
                
                
                DM_Ti=0.d0
                
! s orbital                
                wf=0.d0
                wf1=0.d0
                do i=1,11
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc1(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Ti=DM_Ti+2*wf*wf1
                
                wf=0.d0
                wf1=0.d0
                do i=1,11
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc2(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Ti=DM_Ti+2*wf*wf1
                
                wf=0.d0
                wf1=0.d0
                do i=1,11
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc3(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc3(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Ti=DM_Ti+2*wf*wf1
                
                wf=0.d0
                wf1=0.d0
                do i=1,11
                        wf=wf+norm_cle_DM(alphas(i),sn(i))*r**(sn(i)-1)*sc4(i)*exp(-alphas(i)*r)
						wf1=wf1+norm_cle_DM(alphas(i),sn(i))*r1**(sn(i)-1)*sc4(i)*exp(-alphas(i)*r1)
                ENDDO
                
                DM_Ti=DM_Ti+2*wf*wf1
! p orbital                
                wf=0.d0
                wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc1(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_Ti=DM_Ti+6*wf*wf1
                
                wf=0.d0
                wf1=0.d0
                do i=1,6
                        wf=wf+norm_cle_DM(alphap(i),pn(i))*r**(pn(i)-1)*pc2(i)*exp(-alphap(i)*r)
						wf1=wf1+norm_cle_DM(alphap(i),pn(i))*r1**(pn(i)-1)*pc2(i)*exp(-alphap(i)*r1)
                ENDDO
                DM_Ti=DM_Ti+6*wf*wf1
! d orbital                
                wf=0.d0
                wf1=0.d0
                do i=1,5
                        wf=wf+norm_cle_DM(alphad(i),dn(i))*r**(dn(i)-1)*dc1(i)*exp(-alphad(i)*r)
						wf1=wf1+norm_cle_DM(alphad(i),dn(i))*r1**(dn(i)-1)*dc1(i)*exp(-alphad(i)*r1)
                ENDDO
                DM_Ti=DM_Ti+2*wf*wf1
        ENDFUNCTION


	DOUBLEPRECISION FUNCTION Cle_DM(r,r1)
	        IMPLICIT NONE
	        
		DOUBLEPRECISION, DIMENSION(3) :: r,r1
		DOUBLEPRECISION :: dis_r, dis_r1
                
                DOUBLEPRECISION,DIMENSION(1,3) :: deplace
                
		INTEGER :: i,j,k, m
		
		DOUBLEPRECISION, DIMENSION(num_atom_cell,3) :: atom_list
		
		DOUBLEPRECISION, DIMENSION(num_atom_cell,1) ::vec_dep
		
		vec_dep=1.d0

		Cle_DM=0.d0
		
                       do i= -1,1
                                do j= -1,1
                                        do k= -1,1
                                                deplace(1,:)=i*lattice(1,:) + j*lattice(2,:)+k*lattice(3,:)
                                                atom_list(:,:)=atom_coord_list(:,:) + MATMUL(vec_dep, deplace)
						do m=1,num_atom_cell
						        dis_r=dis(atom_list(m,:),r)!/0.52917706
								dis_r1=dis(atom_list(m,:),r1)
							if (trim(atome_ele_list(m))=='H') then
								Cle_DM=Cle_DM+DM_H(dis_r,dis_r1)/(4*pi)
							elseif (trim(atome_ele_list(m))=='C') then
								Cle_DM=Cle_DM+DM_C(dis_r,dis_r1)/(4*pi)
							elseif (trim(atome_ele_list(m))=='N') then
								Cle_DM=Cle_DM+DM_N(dis_r,dis_r1)/(4*pi)
							elseif (trim(atome_ele_list(m))=='O') then
								Cle_DM=Cle_DM+DM_O(dis_r,dis_r1)/(4*pi)
							elseif (trim(atome_ele_list(m))=='S') then
								Cle_DM=Cle_DM+DM_S(dis_r,dis_r1)/(4*pi)
							elseif (trim(atome_ele_list(m))=='TI') then
								Cle_DM=Cle_DM+DM_Ti(dis_r,dis_r1)/(4*pi)
				                        elseif (trim(atome_ele_list(m))=='Y') then
								Cle_DM=Cle_DM+DM_Y(dis_r,dis_r1)/(4*pi)
							endif
						enddo
                                        ENDDO
                                ENDDO
                        ENDDO
                

	ENDFUNCTION

ENDMODULE Clementi_DM
