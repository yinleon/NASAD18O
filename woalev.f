      PROGRAM WOALEV
************************************************************************
* This script reads sequential big_endian binary files downloaded off of
* Columbia/LDEO's Ingrid database. The binary files are from the World  
* Ocean's Atlas 2009 and Levitus 1984 and contain objectively analyzed
* fields which can be manipulated and compared between databses.
*  
* The fields are gridded in a 1 x 1 degree grid by Longitude(360),     
* Latitude(180), and Depth(33). The entire dataset is on one record
* for WOA09 and divded by depth(z) per record in Levitus.
*
* This script also calculates potential density and fixed phosphate,
* and compares calculated values with values read from Levitus.      
*
* This code was written by Leon Yin under the Supervision of
* Dr. Allegra LeGrande at NASA GISS June 2015.
*
* Version 1.1
************************************************************************
      IMPLICIT NONE
      
      INTEGER,PARAMETER ::UNITC1=1000,P=1,im=360,jm=180,lm=33
      REAL*4, PARAMETER ::skip=1e-30,junk=9.96920997E+36,seaD=1025
      CHARACTER*100 tempLevF,tempWoaF,po4WoaF,po4LevF,o2LevF,
     *     o2WoaF,salLevF,salWoaF,po4starF,potDenF,denF,potTempF
      REAL*4 tempLev(im,jm,lm),tempWoa(im,jm,lm),
     *     po4Woa(im,jm,lm),po4Lev(im,jm,lm),
     *     o2Lev(im,jm,lm),o2Woa(im,jm,lm),
     *     salLev(im,jm,lm),salWoa(im,jm,lm),
     *     po4starCal(im,jm,lm),po4starLev(im,jm,lm),
     *     po4starD(im,jm,lm),potTemp(im,jm,lm),
     *     potdenLev(im,jm,lm),denLev(im,jm,lm),
     *     temporaryA(im,jm),RHO(im,jm,lm),RHOD(im,jm,lm),
     *     tempD(im,jm,lm),o2D(im,jm,lm),po4D(im,jm,lm),
     *     salD(im,jm,lm)
      INTEGER i,j,k
      REAL*4 T,T2,T3,T4,T5,S,S2,S32,A,A2,B,C,RHO1,RHOP,K0,KA,KB,K1,AW,A1
     *     ,BW,B1,P1,P2
********************************INPUT***********************************
*     Assign start up files...
************************************************************************
      tempLevF="annualtempLevatus94v2.s4"        !celsius
      tempWoaF="annualtempWOA09.s4"              !celsius
      po4LevF="annualpo4Levitus.s4"              !uM
      po4WoaF="annualpo4WO09.s4"                 !uM
      o2LevF="annualo2Levitus.s4"                !mL/L ~ uM
      o2WoaF="annualo2woa09.s4"                  !ml/L ~ uM
      salLevF="annualsalinityLevitus94.s4"       !psu
      salWoaF="annualsalinityWOA09.s4"           !psu
      po4starF="fixedphosphateLevitus.s4"        !uM     
      potDenF="potentialdensityLevtus94.s4"      !g/cm3  
      denF="annualdensityLevitus94.s4"           !g/cm3  
      potTempF="pottempLevitus94.s4"             !c      
************************************************************************
*     Populate arrays from the input files level by level.
************************************************************************
      PRINT*,"OPENING WOA TEMP..."
      OPEN(100,FILE=trim(tempWoaF),FORM="UNFORMATTED",access=
     * "SEQUENTIAL",convert="big_endian")
      READ(100,end=9999) tempWoa
      
      PRINT*,"OPENING LEVITUS TEMP..."
      OPEN(101,FILE=trim(tempLevF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")

      PRINT*,"OPENING WOA PO4..."
      OPEN(200,FILE=trim(po4WoaF),FORM="UNFORMATTED",access=
     * "SEQUENTIAL",convert="big_endian")
      READ(200,end=9999) po4Woa
      PRINT*,"OPENING LEVITUS PO4..."
      OPEN(201,FILE=trim(po4LevF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")

      PRINT*,"OPENING WOA O2..."
      OPEN(300,FILE=trim(o2WoaF),FORM="UNFORMATTED",access=
     * "SEQUENTIAL",convert="big_endian")
      READ(300,end=9999) o2Woa
      PRINT*,"OPENING LEVITUS o2..."
      OPEN(301,FILE=trim(o2LevF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      
      PRINT*,"OPENING WOA SAL..."
      OPEN(400,FILE=trim(salWoaF),FORM="UNFORMATTED",access=
     * "SEQUENTIAL",convert="big_endian")
      READ(400,end=9999) salWoa
      PRINT*,"OPENING LEVITUS SAL..."
      OPEN(401,FILE=trim(salLevF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      ! Read refs for calcs
      OPEN(500,FILE=trim(po4starF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      OPEN(600,FILE=trim(potDenF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      OPEN(601,FILE=trim(DenF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      OPEN(602,FILE=trim(potTempF),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      
      ! Go through each level of Levitus and transfer data
      do k=1,lm
         READ(101,end=9999) tempLev(:,:,k)
         READ(201,end=9999) po4Lev(:,:,k)
         READ(301,end=9999) o2Lev(:,:,k)
         READ(401,end=9999) salLev(:,:,k)
         READ(500,end=9999) po4starLev(:,:,k)
         READ(600,end=9999) potdenLev(:,:,k)
         READ(601,end=9999) denLev(:,:,k)
         READ(602,end=9999) potTemp(:,:,k)
      enddo

      PRINT*,"DATA TRANSFER SUCCESSFUL!"
      ! Close the files
      CLOSE(100)
      CLOSE(101)
      CLOSE(200)
      CLOSE(201)
      CLOSE(300)
      CLOSE(301)
      CLOSE(401)
      CLOSE(500)
      CLOSE(600)
      CLOSE(601)
      CLOSE(602)
      
********************************SCRUB***********************************
      do j=1,jm
         do i=1,im
            do k=1,lm
               ! For Levitus replace NaN
               if(isnan(tempLev(i,j,k)))THEN
                  tempLev(i,j,k)=skip
               endif
               if(isnan(po4Lev(i,j,k)))THEN
                  po4Lev(i,j,k)=skip
               endif
               if(isnan(o2Lev(i,j,k)))THEN
                  o2Lev(i,j,k)=skip
               endif
               if(isnan(salLev(i,j,k)))THEN
                  salLev(i,j,k)=skip
               endif
               if(isnan(po4starLev(i,j,k)))THEN
                  po4starLev(i,j,k)=skip
               endif
               if(isnan(potTemp(i,j,k)))THEN
                  potTemp(i,j,k)=skip
               endif
               if(isnan(potdenLev(i,j,k)))THEN
                  potdenLev(i,j,k)=skip
               else
                  potdenLev(i,j,k)=seaD-potdenLev(i,j,k)*10 !convert units
               endif
               if(isnan(denLev(i,j,k)))THEN
                  denLev(i,j,k)=skip
               else
                  denLev(i,j,k)=seaD-denLev(i,j,k)*10       !convert units
               endif
               ! For WOA09 replace 9.96920997E+36
               if((tempWoa(i,j,k).EQ.junk))THEN
                  tempWoa(i,j,k)=skip
               endif               
               if((po4Woa(i,j,k).EQ.junk))THEN
                  po4Woa(i,j,k)=skip
               endif
               if((o2Woa(i,j,k).EQ.junk))THEN
                  o2Woa(i,j,k)=skip
               endif
               if((salWoa(i,j,k).EQ.junk))THEN
                  salWoa(i,j,k)=skip
               endif
            enddo
         enddo
      enddo
      PRINT*,"DATA SCRUB SUCCESSFUL!"      
************************************************************************
* Calculate the corrected Phosphate using Broecker(1998)
* Note the equation for corrected Phosphate is as follows:
*                 PO4*=PO4+O2/175-1.93umol/kg
* where:
*     175 is the avg molar Redfield Ratio of O2 consumption to phosphate
*     remineralization and 1.93 is arbitrary constant Broecker (1985).
************************************************************************
      do j=1,jm
         do i=1,im
            do k=1,lm
               !Check for reasonable data
               if(po4Lev(i,j,k).EQ.skip)THEN
                  po4starD(i,j,k)=skip
                  po4starCal(i,j,k)=skip
               else if(o2Lev(i,j,k).EQ.skip)THEN
                  po4starD(i,j,k)=skip
                  po4starCal(i,j,k)=skip
               else
                  !Calculate PO4*
                  po4starCal(i,j,k)=po4Lev(i,j,k)+o2Lev(i,j,k)/175-1.95     
               endif  
            enddo
         enddo
      enddo
      PRINT*,"PO4* CALCULATED AND COMPARED!"
************************************************************************
* Calculate the (potential) density using Unesco equations cited in
* Millero and Poisson(1981) and Gill (1982).
* Adopted from Matlab code written by Raymond G. Najjar 2004.                                     
* Note: T(C),S(PSU),P(bar),RHO(kg/m3)
*       Recall 1 kg/m3 = 1000 g/cm3 (check this)
*       STANDARD ERROR ~ 3.6 X 10-3 KG/M-3
************************************************************************
      P1=P/10 ! convert the bars-> dbars
      P2=P1**2
      C=4.8314e-4
      A2 = 1.91075e-4
      do j=1,jm
         do i=1,im
            do k=1,lm
               if(potdenLev(i,j,k).EQ.skip)THEN
                  RHO(i,j,k)=skip
               ELSE   
                  if(potTemp(i,j,k).EQ.skip)THEN
                     RHO(i,j,k)=skip
                  ELSE IF(salLev(i,j,k).EQ.SKIP)THEN
                     RHO(i,j,k)=skip
                  ELSE
                     T=potTemp(i,j,k)
                     S=salLev(i,j,k)
                     T2=T**2
                     T3=T**3
                     T4=T**4
                     T5=T**5
                     S2=S**2
                     S32=S**(3/2)
                     ! Density of fresh water at atm pressure
                     RHO1 =999.842594+6.793952e-2*T-9.095290e-3*T2 +
     *                    1.001685e-4*T3-1.120083e-6*T4+6.536332e-9*T5
                     ! Efect of Salinity
                     A = 8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T2 - 
     *                    8.2467e-7*T3 + 5.3875e-9*T4
                     B = -5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T2
                     RHOP = RHO1 + A*S + B*S32 + C*S2
                     ! Effect of Pressure
                     K0 = 19652.21 + 148.4206*T - 2.327105*T2 + 1.360447
     *                    e-2*T3- 5.155288e-5*T4
                     KA = 54.6746 - 0.603459*T + 1.09987e-2*T2 - 6.1670
     *                    e-5*T3 
                     KB = 7.944e-2 + 1.6483e-2*T - 5.3009e-4*T2
                     AW = 3.239908 + 1.43713e-3*T + 1.16092D-4*T2 -
     *                    5.77905D-7*T3
                     A1 = 2.2838e-3 - 1.0981e-5*T - 1.6078e-6*T2
                     BW = 8.50935e-5 - 6.12293e-6*T + 5.2787e-8*T2
                     B1 = -9.9348e-7 + 2.0816e-8*T + 9.1697e-10*T2
                     K1 = K0 + S*KA + S32*KB + P1*(AW + S*A1 + S32*A2) +
     *                    P2*(BW + S*B1)                                                        
                     RHO(i,j,k)=RHOP/(1-P1/K1);
                  endif
               endif
            enddo
         enddo
      enddo
      PRINT*,"Potential Density Calculated!"
****************************File Prints*********************************
c      PRINT*,"Temp Lev..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,tempLev(:,:,k)
c      enddo
c      PRINT*,"Temp WOA..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,tempWoa(:,:,k)
c      enddo      
c      PRINT*,"o2 Lev..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,o2Lev(:,:,k)
c      enddo
c      PRINT*,"o2 WOA..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,o2Woa(:,:,k)
c      enddo
c      PRINT*,"sal Lev..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,salLev(:,:,k)
c      enddo
c      PRINT*,"sal WOA..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,salWoa(:,:,k)
c      enddo
c      PRINT*,"po4 Lev..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,po4Lev(:,:,k)
c      enddo
c      PRINT*,"po4 WOA..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,po4Woa(:,:,k)
c      enddo
c      PRINT*,"po4* Read..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,po4starLev(:,:,k)
c      enddo   
c      PRINT*,"po4* Calculated..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,po4starCal(:,:,k)
c      enddo
c      PRINT*,"Density Read..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,denLev(60,60,k)
c      enddo           
c      PRINT*,"Potential Desnity Read..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,potdenLev(60,60,k)
c      enddo   
c      PRINT*,"Potential Density Calculated..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,RHO(60,60,k)
c      enddo   
      
*************************Difference Calcs*******************************
* Calculate the difference between read values in WOA09 and Levitus94.
* Also print the caluclate difference between calculated and read values
* of PO4* and Potential Density.
************************************************************************     
      do j=1,jm
         do i=1,im
            do k=1,lm
               ! Only calc diff for valid values.
               if(tempLev(i,j,k).EQ.skip)then
                  tempD(i,j,k)=skip
               else if(tempWoa(i,j,k).EQ.skip)then
                  tempD(i,j,k)=skip
               else
                  tempD(i,j,k)=tempLev(i,j,k)-tempWoa(i,j,k)
               endif
               
               if(po4Lev(i,j,k).EQ.skip)then
                  po4D(i,j,k)=skip
               else if(po4Woa(i,j,k).EQ.skip)then
                  po4D(i,j,k)=skip
               else
                  po4D(i,j,k)=po4Lev(i,j,k)-po4Woa(i,j,k)
               endif

               if(o2Lev(i,j,k).EQ.skip)then
                  o2D(i,j,k)=skip
               else if(o2Woa(i,j,k).EQ.skip)then
                  o2D(i,j,k)=skip
               else
                  o2D(i,j,k)=o2Lev(i,j,k)-o2Woa(i,j,k)
               endif               

               if(salLev(i,j,k).EQ.skip)then
                  salD(i,j,k)=skip
               else if(salWoa(i,j,k).EQ.skip)then
                  salD(i,j,k)=skip
               else
                  salD(i,j,k)=salLev(i,j,k)-salWoa(i,j,k)
               endif

               if(potdenLev(i,j,k).EQ.skip)then
                  RHOD(i,j,k)=skip
               else if(RHO(i,j,k).EQ.skip)then
                  RHOD(i,j,k)=skip
               else
                 RHOD(i,j,k)=RHO(i,j,k)-potdenLev(i,j,k)
               endif

               if(po4starLev(i,j,k).EQ.skip)then
                  po4starD(i,j,k)=skip
               else if(po4starCal(i,j,k).EQ.skip)then
                  po4starD(i,j,k)=skip
               else
                  po4starD(i,j,k)=(po4starCal(i,j,k)-po4starLev(i,j,k))
               endif               
            enddo
         enddo
      enddo
      ! Print the differneces line by line per quality
c      PRINT*,"Difference in Temp..."
c      do k=1,lm
c         PRINT*,"Level ",k
c         PRINT*,tempD(:,:,k)
c      enddo
c      PRINT*,"Difference in Salinity..."
c      do k=1,lm  
c         PRINT*,"Level ",k
c         Print*,tempD(:,:,k)
c      enddo
c      PRINT*,"Difference in Dissovled Oxygen..."
c      do k=1,lm  
c         PRINT*,"Level ",K
c         PRINT*,o2D(:,:,k)
c     enddo
c      PRINT*,"Difference in PO4..."
c      do k=1,lm  
c         PRINT*,"Level ",k
c         Print*,po4D(:,:,k)
c      enddo      
c      PRINT*,"Difference in Po4*..."
c      do k=1,lm  
c         PRINT*,"Level ",k
c         Print*,po4starD(:,:,k)
c      enddo
      PRINT*,"Difference in Potential Density..."
      do k=1,lm  
         PRINT*,"Level ",k
         Print*,RHOd(:,:,k)
      enddo
*********************************END************************************
      PRINT*,"TASK COMPLETE: THANK YOU!"
      
 9999 PRINT*,"Runtime error, file has reached an abrupt end!"

      END
