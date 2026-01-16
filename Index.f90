PROGRAM Index   ! *** December-2025.
USE INITIAL_ShW_CLOUD 
USE DISORD
!***************************************************************
!*   This program calculates DOSES                             *
!*   for the A,B etc.UV region by KD & Discr. Ord. method      *
!***************************************************************

! REAL*8 :: dob_sum, dob_mean
! INTEGER :: dob_cnt, i

INTEGER :: i, k, left, right, nvalid
REAL*8  :: dob_sum, dob_mean, t

CHARACTER*50 FIA,FLE_AL,FIT,FIR,FIRS &
            ,F_DOB,F_OROG,F_ALBGLS,F_AEROSGLS,F_CLOUDGLS,F_DAY   
PARAMETER(IOUT=147) !***
PARAMETER (JMAX_F=54) 
PARAMETER (LEG_GA=1)   !  UVI at midday will be calculated only  -*** December-2025.
DIMENSION XLGA(LEG_GA),WLGA(LEG_GA),SOLAR_TIME(LEG_GA),CSSZA(LEG_GA)
ALLOCATABLE DOS_A(:),DOS_B(:),DOS_E(:),DOS_S(:) &
           ,DOB(:),OROG(:),ALBGLS(:),AEROSGLS(:),CLOUDGLS(:) &
		   ,grid_s(:),grid_d(:) ! *** grds -> shirota  dolgota 	 
! ***        OCO  altitude  albedo     aerosol *** !
COMMON/DOB/DOBSON_F,O3_F(JMAX_F) 
COMMON/SZA/ALC,O,ALBGLS1,RAY_WN
DATA XLGA,WLGA/1.0,1.0/  !  *** December-2025.
!* ---------------------------------------------------------------- *
 DtoR=PI/180.0
WRITE(*,*) ' =============== REMEMBER   0.0 =>  300.0 DU =============='
!*************  ALL FILEs ***********************!
OPEN(IOUT,FILE='UV-CONTROL_FILE')
READ(IOUT,100)F_day
READ(IOUT,100)F_DOB
READ(IOUT,100)F_OROG
READ(IOUT,100)F_AEROSGLS 
READ(IOUT,100)F_ALBGLS
READ(IOUT,100)F_CLOUDGLS
READ(IOUT,100)FIR
CLOSE(IOUT)
!************************************************!
OPEN(310514,FILE=F_DAY) ; READ(310514,*)DAY ; CLOSE(310514)
OPEN(310515,FILE=F_DOB)
OPEN(310516,FILE=F_OROG)
OPEN(310517,FILE=F_AEROSGLS)
OPEN(310518,FILE=F_ALBGLS)
OPEN(310519,FILE=F_CLOUDGLS)
OPEN(310514,FILE='./Input/grid') 
READ(310514,*)nshir,DSHIR,SHIR1,SHIR2
READ(310514,*)ndolg,DDOLG,DOLG1,DOLG2
allocate(grid_s(nshir),grid_d(ndolg))
do j=1,nshir ;  read(310514,*)grid_s(j) ;  end do
do j=1,ndolg ;  read(310514,*)grid_d(j) ; end do
close(310514)

ALLOCATE (DOS_A(NDOLG),DOS_B(NDOLG),DOS_E(NDOLG),DOS_S(NDOLG) &
,DOB(NDOLG),OROG(NDOLG),ALBGLS(NDOLG),AEROSGLS(NDOLG),CLOUDGLS(NDOLG))
 FIRS=TRIM(FIR)//'.UV'
 100    FORMAT(A50)
S_Insol=0.0 
! ############################################## !
 !### READ(IOUT,*)S_Insol ! Extraterrestrial Solar Flux
! ############################################## ! 
			  IOT=8432+1        
	   OPEN(IOT+22,FILE='./Output/E_'//TRIM(FIRS))    
    CALL ATM_PROF_READING
! ======================================== !
NCASE=LEG_GA

!***************************************************************!
!***************************************************************!
!***************************************************************!
! *** Loop over Latitude *** !

DOB=0. 
DO LAT=1,NSHIR
READ(310515,*)ATITUDE ;  READ(310515,*)DOB

! Count valid + compute mean (fallback)
dob_sum = 0.0D0
nvalid  = 0
DO i = 1, NDOLG
  IF (DOB(i) >= 1.0) THEN
    dob_sum = dob_sum + DOB(i)
    nvalid  = nvalid + 1
  END IF
END DO
IF (nvalid > 0) THEN
  dob_mean = dob_sum / nvalid
ELSE
  dob_mean = 300.0D0
END IF

IF (nvalid < 2) THEN
  ! Not enough points to interpolate: fill with mean
  DO i = 1, NDOLG
    IF (DOB(i) < 1.0) DOB(i) = dob_mean
  END DO
ELSE
  ! For each missing point, find nearest valid neighbors to the left and right (cyclic)
  DO i = 1, NDOLG
    IF (DOB(i) < 1.0) THEN

      left = i
      DO k = 1, NDOLG
        left = left - 1
        IF (left < 1) left = NDOLG
        IF (DOB(left) >= 1.0) EXIT
      END DO

      right = i
      DO k = 1, NDOLG
        right = right + 1
        IF (right > NDOLG) right = 1
        IF (DOB(right) >= 1.0) EXIT
      END DO

      IF (DOB(left) < 1.0 .OR. DOB(right) < 1.0) THEN
        ! Should not happen if nvalid>=2, but keep it safe
        DOB(i) = dob_mean
      ELSE
        ! Distance to neighbors along cyclic index
        IF (right > left) THEN
          t = REAL(i - left, 8) / REAL(right - left, 8)
        ELSE
          ! Wrap-around case: interval crosses NDOLG->1
          IF (i > left) THEN
            t = REAL(i - left, 8) / REAL((NDOLG - left) + right, 8)
          ELSE
            t = REAL((NDOLG - left) + i, 8) / REAL((NDOLG - left) + right, 8)
          END IF
        END IF

        DOB(i) = (1.0D0 - t) * DOB(left) + t * DOB(right)
      END IF

    END IF
  END DO
END IF
! -------------------------------------------------------------------------------


! ! --- Fill missing ozone (DOB < 1) using zonal mean at this latitude ---

! dob_sum = 0.0D0
! dob_cnt = 0
! DO i = 1, NDOLG
!   IF (DOB(i) >= 1.0) THEN
!     dob_sum = dob_sum + DOB(i)
!     dob_cnt = dob_cnt + 1
!   END IF
! END DO

! IF (dob_cnt > 0) THEN
!   dob_mean = dob_sum / dob_cnt
! ELSE
!   dob_mean = 300.0D0   ! fallback: if the entire latitude row is missing
! END IF

! DO i = 1, NDOLG
!   IF (DOB(i) < 1.0) DOB(i) = dob_mean
! END DO
! ! ----------------------------------------------------------------------

 if(LAT/10*10==LAT)write(*,*)ATITUDE 
READ(310516,*)ATITUDE2 ; READ(310516,*)OROG ; OROG=OROG/1000.0 ! km <== m
READ(310517,*)ATITUDE3 ; READ(310517,*)AEROSGLS
READ(310518,*)ATITUDE4 ; READ(310518,*)ALBGLS
READ(310519,*)ATITUDE5 ; READ(310519,*)CLOUDGLS
CRT=(ATITUDE2-ATITUDE)**2+(ATITUDE3-ATITUDE)**2+(ATITUDE4-ATITUDE)**2
IF(CRT>0.001)THEN
WRITE(*,*)'*** Latitude Grids NOT the SAME ***' ; STOP
END IF
CALL TIME_SZA(LEG_GA,DAY,ATITUDE,XLGA,SOLAR_TIME,CSSZA,DLD2)
 DOS_A=0. ; DOS_B=0. ; DOS_E=0. ; DOS_S=0. 
IF(DLD2<0.)GOTO  12345
  
  DO LON=1,NDOLG
! -------------------------------------------------- !
SUM_A=0. ; SUM_B=0. ; SUM_E=0. ; SUM_S=0.


! *** the Given Points *** !
DOBSON=DOB(LON)
ALBGLS1=ALBGLS(LON)
AEROS_FACTOR=1.-AEROSGLS(LON)
POINT_ALT=OROG(LON)
CLOUD_BALL=CLOUDGLS(LON)
CLOUD_FACTOR=1.

IF(CLOUD_BALL>50.)THEN
   IF(CLOUD_BALL>95.) CLOUD_FACTOR=0.5
   IF(CLOUD_BALL>85.0.AND.CLOUD_BALL<=95.)CLOUD_FACTOR=0.7
   IF(CLOUD_BALL>75.0.AND.CLOUD_BALL<=85.)CLOUD_FACTOR=0.8
   IF(CLOUD_BALL>65.0.AND.CLOUD_BALL<=75.)CLOUD_FACTOR=0.85
   IF(CLOUD_BALL<=65.)CLOUD_FACTOR=0.9
END IF
! *** WRITE(*,*)CLOUD_BALL,CLOUD_FACTOR
! ************************ ! 

! *********************** ATTENTION ********************* !
! IF(DOBSON<1.)DOBSON=300. ! no O3  ! *** ATTENTION!!! *** 
! *********************** ATTENTION ********************* !

! *** Ozone profife for the given poin *** !
DO J=1,JMAX_F
 RO_O3(J)=O3_F(J)*DOBSON/DOBSON_F
END DO

ROFINT_O3(JMAX)=0.
DO J=JMAX,2,-1
ROFINT_O3(J-1)=ROFINT_O3(J)+0.5*(RO_O3(J)+RO_O3(J-1))*(Z(J)-Z(J-1))
END DO
DO J=JMAX,2,-1
ROFINT_O3(J-1)=ALOG(ROFINT_O3(J-1))
END DO


! =================================== !
DO NCA=1,NCASE  ! Loop over SZA 
TIME=SOLAR_TIME(NCA)
ANGLE=ACOS(CSSZA(NCA))/DtoR

!### WRITE(*,*)NCA,TIME,ANGLE,DOBSON

! *** Profile for KD ***  !
ALC=CSSZA(NCA)
IF(ALC<0.0001)ALC=0.0001
ALCL=ALOG(ALC)   
ALDALC=-ALCL
DO J=1,JMAX
ALO=ROFINT_O3(J)+ALDALC
RO_INT_O3(J)=AMAX1(ALO,0.)
END DO

DO NV=3,3  !  *** Erithema ONLY !!! N_VALUES => Loop over KD for the given (9) values.
! ***  4 values is not used *** !
 CALL K_COEF_KD

! *** RAYLEIGH ***  !
BDA=10000./RAY_WN
      R=7.68246E-4/BDA**(3.55212+1.35579*BDA+0.11563/BDA)*288.15
     DO J=1,JMAX
	 S_MOL(J)=R*P1(J)/T1(J)
IF(S_MOL(J)*RABMA(J)>0.)THEN
WW00(J)=S_MOL(J)/(S_MOL(J)+RABMA(J))
ELSE
WW00(J)=0.
END IF
      END DO
TAUMA(JMAX)=0.
DO J=JMAX,2,-1
TAUMA(J-1)=TAUMA(J)+0.5*(S_MOL(J)+RABMA(J)+S_MOL(J-1) &
 +RABMA(J-1))*(Z(J)-Z(J-1))
END DO

! *** TAU-grid! Here the ALTITUDE is TAKEN INTO ACCOUNT  ****
IF(POINT_ALT<0.)POINT_ALT=0.
DO JA=2,JMAX
IF(Z(JA-1)<=POINT_ALT.AND.Z(JA)>=POINT_ALT)EXIT
END DO
C2=(POINT_ALT-Z(JA-1))/(Z(JA)-Z(JA-1)) ; C1=1.-C2
TAU1=C1*TAUMA(JA-1)+C2*TAUMA(JA) 
STEP=TAU1/(MP-1)
DO J=1,MP
XXX=TAU1-STEP*(J-1) ; IF(XXX<0.)XXX=0. ; X(J)=XXX !Homogen. TAU-grid.
   DO JJ=2,JMAX
   IF(TAUMA(JJ-1)>=XXX.AND.TAUMA(JJ)<=XXX)EXIT
   END DO 
   C1=(XXX-TAUMA(JJ))/(TAUMA(JJ-1)-TAUMA(JJ)) ; C2=1.-C1
ZZ(J)=C1*Z(JJ-1)+C2*Z(JJ)
W0(J)=C1*WW00(JJ-1)+C2*WW00(JJ)
END DO
	FLUXUP0=0. ; FLUXUPN=0. ; FLUXD0=0. ; FLUXdN =0.  
VCHAN=RAY_WN
 CALL DISCRORD 
FDO9(:,NV)=FLUXD0+FLUXdN
!### FUP9(:,NV)=FLUXUP0+FLUXUPN   ! Is't used!
       	END DO  ! Loop over "values" (UV-A, UV-B, Erithema etc.
!* ---------------------------------- *

!*			Writing			*
IF(S_Insol>0.0)THEN 
cs=S_INSOL  ! <=== ATTENTION! it's FACTOR 
FDO9=FDO9*cs
!### FUP9=FUP9*cs
END IF
       
     !*******  Writing ******
			  IOT=8432+1        
! *** 
ErDO=FDO9(1,3)
CMWT=1000.*AEROS_FACTOR*CLOUD_FACTOR ! Wt->mWt and Aerosol&Cloud Factor
2765 FORMAT(F8.3,F8.2,5F10.1)
  2         FORMAT(F7.2,'  km ')
SUM_E=ErDO*CMWT ! *** December-2025. 
END DO ! DO NCA=1,NCASE
59 CONTINUE
DOS_E(LON)=SUM_E/25.   !  Index UVI at midday *** December-2025. 
! ************************** !
    END DO ! NDOLG
12345 CONTINUE 
WRITE(IOT+22,*)ATITUDE
WRITE(IOT+22,5)DOS_E
5 FORMAT(60E10.3)

END DO ! NSHIR
CLOSE(IOT+20) ; CLOSE(IOT+21) ; CLOSE(IOT+22) ; CLOSE(IOT+23) 
CLOSE(310515)
!* ---------------------------------- *
       END
   SUBROUTINE ATM_PROF_READING
  USE INITIAL_ShW_CLOUD
   CHARACTER*1 A_ZO
PARAMETER (JMAX_F=54)
   DIMENSION Z_F(JMAX_F),P_F(JMAX_F),T_F(JMAX_F)
COMMON/DOB/DOBSON_F,O3_F(JMAX_F) 
DATA DOBSON_F,Z_F,P_F,T_F,O3_F/300.0, &
 0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0, &
 11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0, &
 24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0, &
 37.0,38.0,39.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0, &
0.1000E+01,0.9452E+00,0.8904E+00,0.8411E+00,0.7917E+00,0.7463E+00, &
0.7009E+00,0.6604E+00,0.6199E+00,0.5834E+00,0.5469E+00,0.4808E+00, &
0.4205E+00,0.3672E+00,0.3198E+00,0.2774E+00,0.2399E+00,0.2063E+00, &
0.1767E+00,0.1510E+00,0.1283E+00,0.1096E+00,0.9378E-01,0.8016E-01, &
0.6861E-01,0.5874E-01,0.5035E-01,0.4314E-01,0.3712E-01,0.3179E-01, &
0.2734E-01,0.2354E-01,0.2028E-01,0.1748E-01,0.1509E-01,0.1303E-01, &
0.1133E-01,0.9851E-02,0.8555E-02,0.7423E-02,0.6436E-02,0.5612E-02, &
0.4900E-02,0.4284E-02,0.3750E-02,0.3287E-02,0.1737E-02,0.9388E-03, &
0.5084E-03,0.2685E-03,0.1372E-03,0.6614E-04,0.2962E-04,0.1185E-04, &
294.2,292.0,289.7,287.5,285.2,282.2,279.2,276.2,273.2,270.2,267.2, &
261.2,254.7,248.2,241.7,235.3,228.8,222.3,215.8,215.7,215.7,215.7, &
215.7,216.8,217.9,219.2,220.4,221.6,222.8,223.9,225.1,226.2,227.6, &
229.3,231.4,233.7,235.7,237.9,240.2,242.6,245.2,247.6,250.1,252.5, &
255.0,257.5,269.9,275.7,269.3,257.1,240.1,218.1,196.1,174.1, &
0.6719E+17,0.6720E+17,0.6721E+17,0.6719E+17,0.6718E+17,0.6831E+17, &
0.6943E+17,0.7055E+17,0.7167E+17,0.7280E+17,0.7393E+17,0.7725E+17, &
0.8397E+17,0.8844E+17,0.9629E+17,0.1007E+18,0.1231E+18,0.1356E+18, &
0.1610E+18,0.2019E+18,0.1949E+18,0.1997E+18,0.1994E+18,0.2423E+18, &
0.3095E+18,0.3511E+18,0.3592E+18,0.3698E+18,0.3710E+18,0.3719E+18, &
0.3820E+18,0.3614E+18,0.3370E+18,0.3103E+18,0.2829E+18,0.2557E+18, &
0.2354E+18,0.2143E+18,0.1933E+18,0.1727E+18,0.1530E+18,0.1326E+18, &
0.1132E+18,0.9495E+17,0.7824E+17,0.6313E+17,0.1897E+17,0.6244E+16, &
0.2226E+16,0.8892E+15,0.2994E+15,0.7945E+14,0.1879E+14,0.8912E+13/
JMAX=JMAX_-1
DO J=1,JMAX
Z(J)=Z_F(J) ; P1(J)=P_F(J) ; T1(J)=T_F(J) 
END DO
     DO J=1,JMAX
	 PdT(J)=P1(J)/T1(J)
	 END DO
	 S_PdT=0.
	 DO J=JMAX,2,-1
	 A=(Z(J)-Z(J-1))/(ALOG(PdT(J-1)/PdT(J)))
	 DTU=PdT(J-1)*A*(1.-EXP(-(Z(J)-Z(J-1))/A))
     S_PdT(J-1)=S_PdT(J)+DTU
	 END DO
 END SUBROUTINE ATM_PROF_READING
! *** 29 May,2015 ***
! *** Time_SZA_CREATION *** !
SUBROUTINE TIME_SZA(NT,DAY,ATITUDE,XLGA,SOLAR_TIME,COSSZA,DLD2)
PARAMETER (Pi=3.14159265359D0,DtoR=PI/180.0, AX=23.45*DtoR)
                         !Earth Axis- Ecliptic
DIMENSION XLGA(NT),SOLAR_TIME(NT),COSSZA(NT)
SAX=SIN(AX)

! --- Initial Data --- !
ATITUDE=AMIN1(ATITUDE,89.9)
FI=ATITUDE*DtoR
SIN_FI=SIN(FI) ; COS_FI=COS(FI)

! ------------------------------------------- !
! *** Sklonenie DELT **** !
DLT=23.45*SIN(2.*PI*(DAY-81.)/365.)*DtoR
SIN_DELT=SIN(DLT) ; COS_DELT=COS(DLT)
ARGUM=-SIN_DELT*SIN_FI/COS_DELT/COS_FI
! Polar Night 
IF(ARGUM>=1.)THEN
DLD2=-1.
RETURN
END IF
IF(ARGUM>-1.) THEN
SUN_RISE_UG=ACOS(ARGUM)
SUN_RISE=12.-SUN_RISE_UG/(2.*PI/24.)
DLD2=12.-SUN_RISE
ELSE
! Polar Day
SUN_RISE_UG=Pi
SUN_RISE=0.
DLD2=12.
END IF
DO N=1,NT
SSOLAR_TIME=SUN_RISE+DLD2*XLGA(N)
HRA=PI/12.*(SSOLAR_TIME-12.)
! *** SZA *** !
COS_SZA=SIN_DELT*SIN_FI+COS_DELT*COS_FI*COS(HRA)
SOLAR_TIME(N)=SSOLAR_TIME-SUN_RISE
COSSZA(N)=COS_SZA
END DO
END
SUBROUTINE K_COEF_KD
 COMMON/SZA/COS_SZA,O,ALBED,RAY_WN 
       CALL UV_E_RABMA
RETURN
 END 

! UV-E, (RABMA (now for O3 only) )  (03 Apr.,2017).
     SUBROUTINE UV_E_RABMA     ! <-- Erithema
 USE INITIAL_ShW_CLOUD 
    	COMMON/SZA/COS_SZA,O,ALBED,RAY_WN
         O=10.22   ! <-- E
		RAY_WN=31702.0 * 1.013  !1.015
         JZM_1=JMAX-1
         RABMA(JMAX)=0.
         DO J=JZM_1,1,-1
	   WWW1=RO_INT_O3(J) ! this is LOG for ammaunt
	   RK_KD1=RK_E_03(WWW1)*RO_O3(J)
RABMA(J)=RK_KD1 ! + other components if necessary  
         END DO
		 END
FUNCTION RK_E_03(WWW)
! W-> LOG of ammount along the solar ray (LOG(...[mol/cm^2])
PARAMETER (NW=17,DW=0.25,W1=16.0,W2=20.0,DW2=DW*DW) 
DIMENSION X(NW),Y(NW)
DATA X/16.00,16.25,16.50,16.75,17.00,17.25,17.50,17.75,18.00,18.25, &
18.50,18.75,19.00,19.25,19.50,19.75,20.00/, &
Y/117.34,117.77,117.32,115.75,112.20,105.86,95.49,82.39,65.75,49.58, &
34.01,21.72,12.58,6.83,3.19,1.60,0.80/
DATA IBG/0/
IF(IBG==0)THEN ; Y=Y*1E-20 ; IBG=1 ;END IF
WW=WWW/2.3025825  ! Ln -> Log
! ------------------------------------------------------------- !
IF(WW<=W1)THEN ; RK_E_03=Y(1) ; RETURN ; END IF
IF(WW>=W2)THEN ; RK_E_03=Y(NW) ; RETURN ; END IF
! --- !
I=(WW-W1)/DW+2
IF(I+1>NW)I=I-1
   X2X=X(I)-WW ; X1X=WW-X(I-1) ; X3X=X(I+1)-WW
C1=0.5*X2X*X3X/DW2
C2=X1X*X3X/DW2
C3=-0.5*X2X*X1X/DW2
RK_E_03=Y(I-1)*C1+Y(I)*C2+Y(I+1)*C3
END
!*******************************************************************!
SUBROUTINE DISCRORD 
! *** 4 Discrete Ordinates - 11 Sept., 2017  *** !
USE DISORD
 COMMON/SZA/UM0,F0,ALBEDO,RAY_WN 
! *** Matrixes etc. ***
CALL Bn_creating(UM0,Bn,UMn,Gn,NM) 
Bn=0.25*F0*Bn/3.14159265359
PIPI1=0.
DO N=1,NM2
PIPI1=PIPI1+Gn(N)*UMn(N)   
END DO
PIPI=PIPI1/0.5*Pi 

 ! ********************* !
! --- Direct Solar Flux --- !
  DDXX=X(MP-1)
FLUXd0(MP)=F0*UM0
OSLAB=EXP(-DDXX/UM0)
DO J=MP,2,-1
FLUXd0(J-1)=FLUXd0(J)*OSLAB
END DO

! --- ' nonSCATERRED part from Surface ' --- ! 

DO N=1,NM2
Fu0(1,N)=ALBEDO*FLUXd0(1)/PIPI ! Correction for N_streem summation
END DO
DO N=1,NM2
OSLAB=EXP(-DDXX/UMn(N))
   DO J=2,MP
   Fu0(J,N)=Fu0(J-1,N)*OSLAB
   END DO
END DO
DO J=1,MP
   SUM=0.
   DO N=1,NM2
   SUM=SUM+Gn(N)*UMn(N)*Fu0(J,N)
   END DO
   FLUXup0(J)= 2.*PI*SUM                    
END DO

! *** -----  Y0  (Single Scattering approximation) *** ----- !
 ! --- W0 - Function ---- !
DD=ALBEDO*F0*EXP(-X(1)/UM0)/(32.*Pi)
  DO J=1,MP
    DO N=1,NM
    WA(J,N)=-W0(J)*Bn(N)/UMn(N)*EXP(-X(J)/UM0)
	XXS=X(1)-X(J) ; CALL E2E4(XXS,E2,E4)
WB(J,N)=-DD*W0(J)/UMn(N)*((3.-UMn(N)**2)*E2+(3.*UMn(N)**2)*E4)
	END DO
  END DO
W=WA+WB   !

! --- DOWNward --- !
DO N=NM_2,NM
 Y0(MP,N)=0. ! <--- Boundary Conditions
  GRAL=0.
 OSLUMN1=EXP(-DDXX/UMn(N))
   EP=1.
   G1=W(MP,N)              
  DO J=MP,2,-1
    EP=EP*OSLUMN1
    G2=W(J-1,N)*EP
GRAL=GRAL+0.5*(G1+G2)*DDXX
	G1=G2
Y0(J-1,N)=GRAL/EP 
  END DO
END DO

 
! --- UPward --- !
DO N=1,NM2
Y0(1,N)=0. ! <--- Boundary Conditions for UPWARD radiation
   GRAL=0.
  OSLUMN1=EXP(DDXX/UMn(N))
   EP=EXP(-X(1)/UMn(N))
   EP_S=EP
   G1=W(1,N)*EP              
       DO J=2,MP
       EP=EP*OSLUMN1
	   G2=W(J,N)*EP
	   GRAL=GRAL-0.5*(G1+G2)*DDXX
	   G1=G2
Y0(J,N)=GRAL/EP 
	   END DO
END DO
! *** Taking into account REFLECTED UPWARD RADIATION *** !
REFL=0
DO N=NM_2,NM
REFL=REFL-Gn(N)*UMn(N)*Y0(1,N)
END DO 
REFL=ALBEDO*REFL*2. 
DO N=1,NM2
 OSLUMN1=EXP(-DDXX/UMn(N))
   EP=1.
  DO J=1,MP
  Y0(J,N)=Y0(J,N)+REFL*EP !
  EP=EP*OSLUMN1
  END DO
END DO

! ----------------------------- !
! *** Fluxes *** !
   DO J=1,MP
FF00UP=0.
FF00DO=0.
DO N=1,NM2
FF00UP=FF00UP+Gn(N)*UMn(N)*Y0(J,N)
END DO
FLUXup1(J)=FF00UP*2.*Pi 
DO N=NM_2,NM
FF00DO=FF00DO-Gn(N)*UMn(N)*Y0(J,N)
END DO 
FLUXd1(J)=FF00DO*2.*Pi
   END DO
! ************ N-scattering **************** !
FLUXD0=FLUXD0+FLUXd1 ; FLUXup0=FLUXup0+FLUXup1
FLUXdN=0. ; FLUXupN=0. 
DO N_CRAT=1,NSCATEVENT  ! Loop over Scatterring Events 
 ! --- Wn - Function ---- !
  DO J=1,MP
    DO N=1,NM
    SUM=0.
          DO M=1,NM
		  SUM=SUM+Anm(N,M)*Y0(J,M)
		  END DO
	W(J,N)=-SUM*W0(J)/UMn(N)
	END DO
  END DO
! ----------------------- !
! --- DOWNward --- !
DO N=NM_2,NM
 Y0(MP,N)=0. ! <--- Boundary Conditions
  GRAL=0.
 OSLUMN1=EXP(-DDXX/UMn(N))
   EP=1.
   G1=W(MP,N)              
  DO J=MP,2,-1
    EP=EP*OSLUMN1
    G2=W(J-1,N)*EP
GRAL=GRAL+0.5*(G1+G2)*DDXX
	G1=G2
Y0(J-1,N)=GRAL/EP        
  END DO
END DO

! *** Fluxes DOWN *** !
   DO J=1,MP
FF00DO=0.
DO N=NM_2,NM
FF00DO=FF00DO-Gn(N)*UMn(N)*Y0(J,N)
END DO 
IF(J==1)YnRefl=ALBEDO*FF00DO *2. 
FLUXdN(J)=FLUXdN(J)+FF00DO*Pi2 
   END DO

! --- UPward --- !
DO N=1,NM2
Y0(1,N)=YnRefl  ! <--- Boundary Conditions for UPWARD radiation
   GRAL=0.
  OSLUMN1=EXP(DDXX/UMn(N))
   EP=EXP(-X(1)/UMn(N))
   EP_S=EP
   G1=W(1,N)*EP              
       DO J=2,MP
       EP=EP*OSLUMN1
	   G2=W(J,N)*EP
	   GRAL=GRAL-0.5*(G1+G2)*DDXX
	   G1=G2
Y0(J,N)=GRAL/EP 
! *** Taking into account REFLECTED UPWARD RADIATION *** !
Y0(J,N)=Y0(J,N)+YnRefl*EP_S/EP  
	   END DO
END DO
! *** Fluxes UP *** !
   DO J=1,MP
FF00UP=0.
DO N=1,NM2
FF00UP=FF00UP+Gn(N)*UMn(N)*Y0(J,N)
END DO
FLUXupN(J)=FLUXupN(J)+FF00UP*Pi2 
   END DO
END DO          ! END Loop over Scatterring Events 
END
SUBROUTINE E2E4(X,E2,E4)
! *** E2(X) and E2(X) Functions; 18 August, 2017 *** !
PARAMETER(NP=71,STEP0=0.0004,DECR=1.08,DECR1=DECR-1,DECRL=7.6961078E-2)
PARAMETER (NP2=91,STEP2=0.1)
DIMENSION XG(NP),Y_E2(NP),Y_E4(NP)
DIMENSION X2(NP2),Y2E2(NP2),Y2E4(NP2)
DATA XG/0.0,0.4E-03,0.832E-03,0.12986E-02,0.18024E-02,0.23466E-02, &
0.29344E-02,0.35691E-02,0.42547E-02,0.49950E-02, &
 0.57946E-02,0.66582E-02,0.75909E-02,0.85981E-02,0.96860E-02, &
 0.10861E-01,0.12130E-01,0.13500E-01,0.14980E-01,0.16579E-01, &
 0.18305E-01,0.20169E-01,0.22183E-01,0.24357E-01,0.26706E-01, &
 0.29242E-01,0.31982E-01,0.34940E-01,0.38136E-01,0.41586E-01, &
 0.45313E-01,0.49338E-01,0.53685E-01,0.58380E-01,0.63451E-01, &
 0.68927E-01,0.74841E-01,0.81228E-01,0.88126E-01,0.95576E-01, &
 0.10362E+00,0.11231E+00,0.12170E+00,0.13183E+00,0.14278E+00, &
 0.15460E+00,0.16737E+00,0.18116E+00,0.19605E+00,0.21214E+00, &
 0.22951E+00,0.24827E+00,0.26853E+00,0.29041E+00,0.31405E+00, &
 0.33957E+00,0.36713E+00,0.39691E+00,0.42906E+00,0.46378E+00, &
 0.50129E+00,0.54179E+00,0.58553E+00,0.63277E+00,0.68380E+00, &
 0.73890E+00,0.79841E+00,0.86268E+00,0.93210E+00,0.10071E+01, &
 0.10880E+01/
 DATA Y_E2/1.0,0.99661,0.99367,0.99075,0.98778,0.98474,0.98160, &
 0.97821,0.97495,0.97141,0.96760,0.96373,0.95966,0.95541,0.95094, &
 0.94615,0.94122,0.93605,0.93062,0.92482,0.91883,0.91244,0.90584, &
 0.89891,0.89154,0.88391,0.87581,0.86739,0.85854,0.84925, &
 0.83953,0.82946,0.81884,0.80775,0.79619,0.78413,0.77157,0.75841, &
 0.74482,0.73062,0.71598,0.70071,0.68491,0.66858,0.65173, &
 0.63435,0.61639,0.59793,0.57899,0.55961,0.53980,0.51952,0.49888, &
 0.47784,0.45662,0.43510,0.41335,0.39151,0.36963,0.34773, &
 0.32588,0.30423,0.28281,0.26175,0.24109,0.22092,0.20133,0.18242, &
 0.16426,0.14693,0.13052/
 DATA Y_E4/0.33333,0.33304,0.33283,0.33261,0.33237,0.33211,0.33183, &
 0.33152,0.33120,0.33085,0.33037,0.32996,0.32952,0.32904,0.32853, &
 0.32788,0.32728,0.32664,0.32595,0.32511,0.32431,0.32335,0.32242, &
 0.32143,0.32026,0.31910,0.31777,0.31640,0.31491,0.31329, &
 0.31155,0.30978,0.30780,0.30568,0.30342,0.30101,0.29844,0.29561, &
 0.29270,0.28950,0.28621,0.28262,0.27881,0.27479,0.27053, &
 0.26603,0.26120,0.25612,0.25078,0.24517,0.23930,0.23308,0.22659, &
 0.21976,0.21275,0.20541,0.19777,0.18989,0.18181,0.17347, &
 0.16491,0.15622,0.14738,0.13847,0.12949,0.12048,0.11151,0.10261, &
 0.09384,0.08527,0.07694/
!  DECRL=ALOG(DECR)
DATA X2/1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10, &
2.20,2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40, &
3.50,3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,4.50,4.60,4.70, &
4.80,4.90,5.00,5.10,5.20,5.30,5.40,5.50,5.60,5.70,5.80,5.90,6.00, &
6.10,6.20,6.30,6.40,6.50,6.60,6.70,6.80,6.90,7.00,7.10,7.20,7.30, &
7.40,7.50,7.60,7.70,7.80,7.90,8.00,8.10,8.20,8.30,8.40,8.50,8.60, &
8.70,8.80,8.90,9.00,9.10,9.20,9.30,9.40,9.50,9.60,9.70,9.80,9.90,10.0/
DATA Y2E2/0.14850E+00,0.12828E+00,0.11110E+00,0.96446E-01, &
0.83890E-01,0.73101E-01,0.63803E-01,0.55771E-01,0.48815E-01, &
0.42780E-01,0.37534E-01,0.32954E-01,0.28972E-01,0.25494E-01, &
0.22452E-01,0.19789E-01,0.17456E-01,0.15408E-01,0.13609E-01, &
0.12028E-01,0.10637E-01,0.94120E-02,0.83326E-02,0.73806E-02, &
0.65406E-02,0.57989E-02,0.51435E-02,0.45641E-02,0.40516E-02, &
0.35980E-02,0.31964E-02,0.28406E-02,0.25253E-02,0.22457E-02, &
0.19977E-02,0.17776E-02,0.15822E-02,0.14087E-02,0.12546E-02, &
0.11176E-02,0.99580E-03,0.88751E-03,0.79118E-03,0.70547E-03, &
0.62919E-03,0.56127E-03,0.50079E-03,0.44692E-03,0.39922E-03, &
0.35641E-03,0.31826E-03,0.28424E-03,0.25390E-03,0.22683E-03, &
0.20269E-03,0.18115E-03,0.16192E-03,0.14475E-03,0.12942E-03, &
0.11573E-03,0.10351E-03,0.92589E-04,0.82831E-04,0.74112E-04, &
0.66319E-04,0.59353E-04,0.53125E-04,0.47556E-04,0.42576E-04, &
0.38122E-04,0.34138E-04,0.30573E-04,0.27384E-04,0.24530E-04, &
0.21975E-04,0.19689E-04,0.17642E-04,0.15810E-04,0.14169E-04, &
0.12700E-04,0.11384E-04,0.10205E-04,0.91492E-05,0.82033E-05, &
0.73558E-05,0.65965E-05,0.59159E-05,0.53061E-05,0.47595E-05, &
0.42695E-05,0.38302E-05/
DATA Y2E4/0.86062E-01,0.75801E-01,0.66824E-01,0.58961E-01, &
0.52064E-01,0.46007E-01,0.40682E-01,0.35997E-01,0.31870E-01, &
0.28232E-01,0.25023E-01,0.22177E-01,0.19675E-01,0.17463E-01, &
0.15506E-01,0.13774E-01,0.12240E-01,0.10881E-01,0.96765E-02, &
0.86081E-02,0.76601E-02,0.68186E-02,0.60713E-02,0.54075E-02, &
0.48176E-02,0.42932E-02,0.38268E-02,0.34119E-02,0.30428E-02, &
0.27142E-02,0.24216E-02,0.21610E-02,0.19288E-02,0.17220E-02, &
0.15376E-02,0.13732E-02,0.12267E-02,0.10959E-02,0.97927E-03, &
0.87520E-03,0.78231E-03,0.69939E-03,0.62535E-03,0.55924E-03, &
0.50019E-03,0.44744E-03,0.40030E-03,0.35818E-03,0.32084E-03, &
0.28716E-03,0.25704E-03,0.23012E-03,0.20603E-03,0.18449E-03, &
0.16522E-03,0.14798E-03,0.13256E-03,0.11875E-03,0.10639E-03, &
0.95332E-04,0.85429E-04,0.76562E-04,0.68622E-04,0.61511E-04, &
0.55142E-04,0.49437E-04,0.44326E-04,0.39747E-04,0.35644E-04, &
0.31967E-04,0.28672E-04,0.25719E-04,0.23071E-04,0.20698E-04, &
0.18570E-04,0.16662E-04,0.14952E-04,0.13418E-04,0.12042E-04, &
0.10808E-04,0.97007E-05,0.87077E-05,0.78169E-05,0.70177E-05, &
0.63006E-05,0.56571E-05,0.50796E-05,0.45614E-05, &
0.40963E-05,0.36788E-05,0.33041E-05/
E2=0. ; E4=0. ; IF(X>9.49)RETURN
IF(X>=1.0)THEN
N=(X-1.0)/STEP2+2
C1=(X-X2(N))*(X-X2(N+1))
C2=(X-X2(N-1))*(X-X2(N+1))
C3=(X-X2(N))*(X-X2(N-1))
E2=0.5/STEP2**2*(C1*Y2E2(N-1)-2.*C2*Y2E2(N)+C3*Y2E2(N+1))
E4=0.5/STEP2**2*(C1*Y2E4(N-1)-2.*C2*Y2E4(N)+C3*Y2E4(N+1))
ELSE
IF(X<=XG(3)) THEN 
N=2
ELSE
N=1.00001+ALOG(1.+X*DECR1/STEP0)/DECRL 
END IF
C1=(X-XG(N))*(X-XG(N+1))/((XG(N-1)-XG(N))*(XG(N-1)-XG(N+1)))
C2=(X-XG(N-1))*(X-XG(N+1))/((XG(N)-XG(N-1))*(XG(N)-XG(N+1)))
C3=(X-XG(N))*(X-XG(N-1))/((XG(N+1)-XG(N))*(XG(N+1)-XG(N-1)))
E2=C1*Y_E2(N-1)+C2*Y_E2(N)+C3*Y_E2(N+1)
E4=C1*Y_E4(N-1)+C2*Y_E4(N)+C3*Y_E4(N+1)
END IF
END
SUBROUTINE Bn_creating(cosSZA,Bn,X,Y,NM)
DIMENSION Bn(NM),X(NM),Y(NM)
A0(U2,UU2)=0.75*(1.+U2*UU2+0.5*(1.-U2)*(1.-UU2))
UU2=cosSZA**2 
NMD2=NM/2
DO I=1,NMD2
U2=X(I)**2 
Bn(I)=A0(U2,UU2)
Bn(NM-I+1)=Bn(I)
END DO
END





 


