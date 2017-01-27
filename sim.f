c
c This program performs an MD simulation of neat SPC/F water.
c
      PROGRAM spcf
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*4 RAND
      PARAMETER nmax=1000,nam=6
      PARAMETER boltz=1.38065e-23
      PARAMETER dmh=1.673533d-27,dmo=2.658019d-26
      PARAMETER dmch2=2.324756d-26,dmch3=2.49081d-26
      PARAMETER dmc=1.9950d-26,dmcl=5.8887d-26
      DIMENSION rat(3,nam,nmax),vat(3,nam,nmax),fat(3,nam,nmax)
      DIMENSION fattor(3,nam,nmax)
      DIMENSION nat(nmax),nmt(nmax),natm(nmax)
      DIMENSION sig(nam,nmax),eps(nam,nmax),dm(nam,nmax)
      DIMENSION qat(nam,nmax),req(nam,nam,nmax),dk(nam,nam,nmax)
      DIMENSION rij(3),fij(3),boxl(3),Vtor1(nam),Vtor2(nam),Vtor3(nam)
      DIMENSION torpot(nam,nam,nam,nam,3,nmax)
      DIMENSION w(3),vvib(3,3,nmax)
      CHARACTER*40 posfile
      CHARACTER*3 int2char
      CHARACTER*1 atype(nam,nmax),atypem
      LOGICAL equil,intlj(nam,nam,nmax),scatsim
c
c Zeroing out arrays.
c
      DO iat=1,nam
         DO jat=1,nam
            DO i=1,nmax
               intlj(iat,jat,i)=0
            END DO
         END DO
      END DO
      DO iat=1,nam
         DO jat=1,nam
            DO kat=1,nam
               DO lat=1,nam
                  DO idim=1,3
                     DO i=1,nmax
                        torpot(iat,jat,kat,lat,idim,i)=0
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO     
c
c Read in simulation input.
c
      OPEN (10,file='sim.inp',status='old')
      READ (10,*) nmol          !total number of molecules
      READ (10,*) ndm           !Number of different molecules
      DO im=1,ndm
         READ (10,*) nmt(im)    !Number of molecules of each type
      END DO
      DO im=1,ndm
         READ (10,*) natm(im)   !Number of atoms per molecule
         DO iat=1,natm(im)      !for ea. atom of one molecule
            READ (10,*) atypem
            READ (10,*) dmat
            READ (10,*) epsat
            READ (10,*) sigat
            READ (10,*) chat
            nmi=1               !number of molecules initial
            DO i=1,im-1
               nmi=nmi+nmt(i)
            END DO
            nmf=nmi+nmt(i)-1
            DO i=nmi,nmf
               atype(iat,i)=atypem
               dm(iat,i)=dmat   !*1.66053886d-27 converts mass to Kg
               eps(iat,i)=epsat*boltz
               sig(iat,i)=sigat
               qat(iat,i)=chat*1.6022177d-19*DSQRT(8.987551812d9)
               nat(i)=natm(im)
            END DO
         END DO
            
         READ (10,*) nb
         DO ib=1,nb
            READ (10,*) iat,jat,dkpair,reqpair
            DO i=nmi,nmf
               dk(iat,jat,i)=dkpair
               req(iat,jat,i)=reqpair
            END DO
         END DO
         READ (10,*) nilji
         DO ilj=1,nilji
            READ(10,*)iat,jat
            DO i=nmi,nmf
               intlj(iat,jat,i)=1
            END DO
         END DO
         READ(10,*)ntors
         DO itor=1,ntors
            READ(10,*)iat,jat,kat,lat,v1,v2,v3
            DO i=nmi,nmf
               torpot(iat,jat,kat,lat,1,i)=v1
               torpot(iat,jat,kat,lat,2,i)=v2
               torpot(iat,jat,kat,lat,3,i)=v3
            END DO
         END DO
      END DO

      READ (10,*) posfile
      READ (10,*) equil
      READ (10,*) temp
      READ (10,*) boxl(1)
      READ (10,*) boxl(2)
      READ (10,*) boxl(3)
      READ (10,*) nstep
      READ (10,*) dt
      READ (10,*) scatsim
      IF (scatsim) READ (10,*) nscat
      IF (scatsim) iscat=0
      IF (scatsim) READ (10,*) iscatmol
      IF (scatsim) READ (10,*) iseed
      IF (scatsim) READ (10,*) dmtot
      IF (scatsim) READ (10,*) dix
      IF (scatsim) READ (10,*) diy
      IF (scatsim) READ (10,*) diz
      CLOSE (10)
c
c For scatter simulation, loop back here to restart.
 5    iscat=iscat+1
c
      OPEN (10,file=posfile,status='old')
      IF (scatsim) THEN
         DO i=1,nmol
            DO iat=1,nat(i)
               IF (i.NE.iscatmol) THEN
               READ (10,*) rat(1,iat,i),rat(2,iat,i),rat(3,iat,i)
               READ (10,*) vat(1,iat,i),vat(2,iat,i),vat(3,iat,i)
               END IF
            END DO
         END DO
      ELSE
         DO i=1,nmol
            DO iat=1,nat(i)
               READ (10,*) rat(1,iat,i),rat(2,iat,i),rat(3,iat,i)
               READ (10,*) vat(1,iat,i),vat(2,iat,i),vat(3,iat,i)
            END DO
         END DO
      END IF
      CLOSE(10)
c For a scattering simulation, randomly position the scattered molecule.
      IF (scatsim) THEN
      vscale=DSQRT(boltz*temp/dmtot)
      DO idim=1,3
         vat(idim,1,iscatmol)=gauss(iseed)*vscale
         vat(idim,2,iscatmol)=vat(idim,1,iscatmol)
         vat(idim,3,iscatmol)=vat(idim,1,iscatmol)
      END DO
      DO idim=1,2
         rat(idim,2,iscatmol)=RAND(,)*boxl(idim)
      END DO
      rat(1,1,iscatmol)=rat(1,2,iscatmol)+0.8166d-10
      rat(1,3,iscatmol)=rat(1,2,iscatmol)-0.8166d-10
      rat(2,1,iscatmol)=rat(2,2,iscatmol)-0.5771d-10
      rat(2,3,iscatmol)=rat(2,2,iscatmol)-0.5771d-10
      IF (vat(3,1,iscatmol).LT.0) THEN
         rat(3,1,iscatmol)=7d-9
         rat(3,2,iscatmol)=7d-9
         rat(3,3,iscatmol)=7d-9
      ELSE
         rat(3,1,iscatmol)=3d-9
         rat(3,2,iscatmol)=3d-9
         rat(3,3,iscatmol)=3d-9
      END IF
c Rotations
      w(1)=DSQRT(boltz*temp/dix)*gauss(iseed)
      w(2)=DSQRT(boltz*temp/diy)*gauss(iseed)
      w(3)=DSQRT(boltz*temp/diz)*gauss(iseed)
      rhy=-0.5132d-10
      rhx=0.8165d-10
      roy=0.0642d-10
      vat(1,1,iscatmol)=vat(1,1,iscatmol)-w(3)*rhy
      vat(1,2,iscatmol)=vat(1,2,iscatmol)-w(3)*roy
      vat(1,3,iscatmol)=vat(1,3,iscatmol)-w(3)*rhy
      vat(2,1,iscatmol)=vat(2,1,iscatmol)+w(3)*rhx
      vat(2,2,iscatmol)=vat(2,2,iscatmol)
      vat(2,3,iscatmol)=vat(2,3,iscatmol)-w(3)*rhx
      vat(3,1,iscatmol)=vat(3,1,iscatmol)+w(1)*rhy-w(2)*rhx
      vat(3,2,iscatmol)=vat(3,2,iscatmol)+w(1)*roy
      vat(3,3,iscatmol)=vat(3,3,iscatmol)+w(1)*rhy+w(2)*rhx
c Vibrations
      uoh=(1.673533D-27*2.658019D-26)/(1.673533D-27+2.658019D-26)
      uhh=(1.673533D-27**2)/(1.673533D-27+1.673533D-27)

      vib1=DSQRT(boltz*temp/uoh)*gauss(iseed)
      vib2=DSQRT(boltz*temp/uoh)*gauss(iseed)
      vib3=DSQRT(boltz*temp/uhh)*gauss(iseed) 

      xvec1=rat(1,1,iscatmol)-rat(1,2,iscatmol)
      yvec1=rat(2,1,iscatmol)-rat(2,2,iscatmol)
      zvec1=rat(3,1,iscatmol)-rat(3,2,iscatmol)
      
      xvec2=rat(1,3,iscatmol)-rat(1,2,iscatmol)
      yvec2=rat(2,3,iscatmol)-rat(2,2,iscatmol)
      zvec2=rat(3,3,iscatmol)-rat(3,2,iscatmol)
      
      xvec3=rat(1,1,iscatmol)-rat(1,3,iscatmol)
      yvec3=rat(2,1,iscatmol)-rat(2,3,iscatmol)
      zvec3=rat(3,1,iscatmol)-rat(3,3,iscatmol)
      
      rvec1=DSQRT(xvec1**2+yvec1**2+zvec1**2)
      rvec2=DSQRT(xvec2**2+yvec2**2+zvec2**2)
      rvec3=DSQRT(xvec3**2+yvec3**2+zvec3**2)
c how far is each bond stretched
      droh1=(((boltz*temp)/933.1)**.5)*gauss(iseed)
      droh2=(((boltz*temp)/933.1)**.5)*gauss(iseed)
      drhh=(((boltz*temp)/228.3)**.5)*gauss(iseed)
c dividing the directional arrow by the length of the arrow, 
c making all the arrows value = 1
      xvec1=xvec1/rvec1
      xvec2=xvec2/rvec2
      xvec3=xvec3/rvec3
      
      yvec1=yvec1/rvec1
      yvec2=yvec2/rvec2
      yvec3=yvec3/rvec3
c totaling the vibrational energy for each atom based on their bonds with 
c the other two corresponding atoms. oxygen receives the least amount of 
c pull due to its larger mass. the pull on the two hydrogens are equal 
c since they are equal to eachother in the 1-3 bond.
      vvib(1,1,iscatmol)=xvec1*((16./17.)*vib1)+xvec3*((1./2.)*vib3)
      vvib(2,1,iscatmol)=yvec1*((16./17.)*vib1)+yvec3*((1./2.)*vib3)
      vvib(3,1,iscatmol)=zvec1*((16./17.)*vib1)+zvec3*((1./2.)*vib3)
      
      vvib(1,2,iscatmol)=xvec1*((-1./17.)*vib1)+xvec2*((-1./17.)*vib2)
      vvib(2,2,iscatmol)=yvec1*((-1./17.)*vib1)+yvec2*((-1./17.)*vib2)
      vvib(3,2,iscatmol)=zvec1*((-1./17.)*vib1)+zvec2*((-1./17.)*vib2)
      
      vvib(1,3,iscatmol)=xvec2*((16./17.)*vib2)+xvec3*((-1./2.)*vib3)
      vvib(2,3,iscatmol)=yvec2*((16./17.)*vib2)+yvec3*((-1./2.)*vib3)
      vvib(3,3,iscatmol)=zvec2*((16./17.)*vib2)+zvec3*((-1./2.)*vib3)

      rat(1,1,iscatmol)=rat(1,1,iscatmol)+
     &     (xvec1*((16./17.)*droh1)+xvec3*((1./2.)*drhh))
      rat(2,1,iscatmol)=rat(2,1,iscatmol)+
     &     (yvec1*((16./17.)*droh1)+yvec3*((1./2.)*drhh))
      rat(3,1,iscatmol)=rat(3,1,iscatmol)+
     &     (zvec1*((16./17.)*droh1)+zvec3*((1./2.)*drhh))
      
      rat(1,2,iscatmol)=rat(1,2,iscatmol)+
     &     (xvec1*((-1./17.)*droh1)+xvec2*((-1./17.)*droh2))
      rat(2,2,iscatmol)=rat(2,2,iscatmol)+
     &     (yvec1*((-1./17.)*droh1)+yvec2*((-1./17.)*droh2))
      rat(3,2,iscatmol)=rat(3,2,iscatmol)+
     &     (zvec1*((-1./17.)*droh1)+zvec2*((-1./17.)*droh2))
 
      rat(1,3,iscatmol)=rat(1,3,iscatmol)+
     &     (xvec2*((16./17.)*droh2)+xvec3*((-1./2.)*drhh))
      rat(2,3,iscatmol)=rat(2,3,iscatmol)+
     &     (yvec2*((16./17.)*droh2)+yvec3*((-1./2.)*drhh))
      rat(3,3,iscatmol)=rat(3,3,iscatmol)+
     &     (zvec2*((16./17.)*droh2)+zvec3*((-1./2.)*drhh))
      
      DO idim=1,3
         DO iat=1,3
            vat(idim,iat,iscatmol)=vat(idim,iat,iscatmol)
     &           +vvib(idim,iat,iscatmol)
         END DO
      END DO
      END IF
c
      IF (scatsim) THEN
      OPEN (20,file='sim.out-energy'//int2char(iscat),status='unknown')
      OPEN (100,file='zposition'//int2char(iscat),status='unknown')
      ELSE
      OPEN (20,file='sim.out-energy',status='unknown')
      END IF
c
c Now start the simulation.
c
      IF (iscat.EQ.24)THEN
      DO istep=1,nstep
	IF (MOD(istep,100).EQ.0) PRINT *,istep,maxstep
	CALL FLUSH(6)
c Zero the forces.
         DO i=1,nmol
            DO iat=1,nam
               DO idim=1,3
                  fat(idim,iat,i)=0.
                  fattor(idim,iat,i)=0.
               END DO
            END DO
         END DO
         vlj=0.
         vc=0.
         vb=0.
         vtor=0.
         ekin=0.
c Calculate the Lennard-Jones contribution.
         DO i=1,nmol-1
         DO j=i+1,nmol
            DO iat=1,nat(i)
            DO jat=1,nat(j)
               eprod=eps(iat,i)*eps(jat,j)
               IF (eprod.NE.0.) THEN
                  DO idim=1,3
                     rij(idim)=rat(idim,jat,j)-rat(idim,iat,i)
                     rij(idim)=rij(idim)
     *                    -IDNINT(rij(idim)/boxl(idim))*boxl(idim)
                  END DO
c                  PRINT*,rij(1),rij(2),rij(3)
c                  PAUSE
                  rsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
                  PRINT*,rsq
                  PAUSE
                  sigma=(sig(iat,i)+sig(jat,j))/2.
                  epsil=4.*DSQRT(eprod)
                  sr2=sigma*sigma/rsq
                  sr6=sr2*sr2*sr2
                  sr12=sr6*sr6
                  DO idim=1,3
                     fij(idim)=epsil*(-12*sr12+6*sr6)*rij(idim)/rsq
                     fat(idim,iat,i)=fat(idim,iat,i)+fij(idim)
                     fat(idim,jat,j)=fat(idim,jat,j)-fij(idim)
                  END DO
                  vij=epsil*(sr12-sr6)
                  vlj=vlj+vij
                  PRINT*,'vlj',vlj
               END IF
            END DO
            END DO
         END DO
         END DO
c Internal Lennard-Jones interactions.
         DO i=1,nmol
            DO iat=1,nat(i)
            DO jat=1,nat(i)
               IF (intlj(iat,jat,i)) THEN
               eprod=eps(iat,i)*eps(jat,i)
               IF (eprod.NE.0.) THEN
                  DO idim=1,3
                     rij(idim)=rat(idim,jat,i)-rat(idim,iat,i)
                     rij(idim)=rij(idim)
     *                    -IDNINT(rij(idim)/boxl(idim))*boxl(idim)
                  END DO
                  rsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
                  sigma=(sig(iat,i)+sig(jat,i))/2.
                  epsil=4.*DSQRT(eprod)
                  sr2=sigma*sigma/rsq
                  sr6=sr2*sr2*sr2
                  sr12=sr6*sr6
                  DO idim=1,3
                     fij(idim)=epsil*(-12*sr12+6*sr6)*rij(idim)/rsq
                     fat(idim,iat,i)=fat(idim,iat,i)+fij(idim)
                     fat(idim,jat,i)=fat(idim,jat,i)-fij(idim)
                  END DO
                  vij=epsil*(sr12-sr6)
                  vlj=vlj+vij
               END IF
               END IF
            END DO
            END DO
         END DO   
c Coulomb interactions.
         pi=DACOS(-1d0)
         rcut=boxl(1)/2d0
         alpha=0.2d10
         arcut=alpha*rcut
         twoalphasqpi=2.*alpha/DSQRT(pi)
         vc1=DERFC(arcut)/rcut
         vc2=(vc1+twoalphasqpi*DEXP(-arcut*arcut))/rcut
         DO i=1,nmol-1
         DO j=i+1,nmol
            DO iat=1,nat(i)
            DO jat=1,nat(j)
               qq=qat(iat,i)*qat(jat,j)
               IF (qq.NE.0.) THEN
                  DO idim=1,3
                     rij(idim)=rat(idim,jat,j)-rat(idim,iat,i)
                     rij(idim)=rij(idim)
     *                    -IDNINT(rij(idim)/boxl(idim))*boxl(idim)
                  END DO
                  rsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
                  dij=DSQRT(rsq)
                  IF (dij.LT.rcut) THEN
                     ad=alpha*dij
                     derfcad=DERFC(ad)/dij
                     psi=derfcad-vc1+vc2*(dij-rcut)
                     DO idim=1,3
                       fij(idim)=(derfcad+twoalphasqpi*DEXP(-ad*ad))/dij
     *                       -vc2
                       fij(idim)=-qq*fij(idim)*rij(idim)/dij
                       fat(idim,iat,i)=fat(idim,iat,i)+fij(idim)
                       fat(idim,jat,j)=fat(idim,jat,j)-fij(idim)
                    END DO
                    vc=vc+psi*qq
                  END IF
               END IF
               DO idim=1,3
c                  fij(idim)=-qq/rsq*rij(idim)/dij
c                  fat(idim,iat,i)=fat(idim,iat,i)+fij(idim)
c                  fat(idim,jat,j)=fat(idim,jat,j)-fij(idim)
               END DO
c               vij=qq/dij
c               vc=vc+vij
            END DO
            END DO
         END DO
         END DO
c Bond interactions.
         DO i=1,nmol
            DO iat=1,nat(i)-1
            DO jat=iat+1,nat(i)
               DO idim=1,3
                  rij(idim)=rat(idim,jat,i)-rat(idim,iat,i)
                  rij(idim)=rij(idim)
     *                 -IDNINT(rij(idim)/boxl(idim))*boxl(idim)
               END DO
               rsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
               dij=DSQRT(rsq)
               DO idim=1,3
                 fij(idim)=dk(iat,jat,i)*(dij-req(iat,jat,i))
     *                 *rij(idim)/dij
                 fat(idim,iat,i)=fat(idim,iat,i)+fij(idim)
                 fat(idim,jat,i)=fat(idim,jat,i)-fij(idim)
               END DO
               vb=vb+0.5*dk(iat,jat,i)*(dij-req(iat,jat,i))**2
            END DO
            END DO
         END DO
c Torsion potential.
         DO i=1,nmol
            DO iat=1,nat(i)-3
               DO jat=iat+1,nat(i)-2
                  DO kat=jat+1,nat(i)-1
                     DO lat=kat+1,nat(i)
                        V_1=torpot(iat,jat,kat,lat,1,i)
                        V_2=torpot(iat,jat,kat,lat,2,i)
                        V_3=torpot(iat,jat,kat,lat,3,i)
                        IF(V_1.NE.0.OR.V_2.NE.0.OR.V_3.NE.0)THEN
                           xia=rat(1,iat,i)
                           yia=rat(2,iat,i)
                           zia=rat(3,iat,i)
                           xib=rat(1,jat,i)
                           yib=rat(2,jat,i)
                           zib=rat(3,jat,i)
                           xic=rat(1,kat,i)
                           yic=rat(2,kat,i)
                           zic=rat(3,kat,i)
                           xid=rat(1,lat,i)
                           yid=rat(2,lat,i)
                           zid=rat(3,lat,i)
            CALL Torsion(xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid
     &  ,V_1,V_2,V_3,e,dedxia,dedyia,dedzia,dedxib,dedyib,dedzib,dedxic,
     &  dedyic,dedzic,dedxid,dedyid,dedzid)
c
                           vtor=vtor+e
                           fat(1,iat,i)=fat(1,iat,i)-dedxia
                           fat(2,iat,i)=fat(2,iat,i)-dedyia
                           fat(3,iat,i)=fat(3,iat,i)-dedzia
                           fat(1,jat,i)=fat(1,jat,i)-dedxib
                           fat(2,jat,i)=fat(2,jat,i)-dedyib
                           fat(3,jat,i)=fat(3,jat,i)-dedzib
                           fat(1,kat,i)=fat(1,kat,i)-dedxic
                           fat(2,kat,i)=fat(2,kat,i)-dedyic
                           fat(3,kat,i)=fat(3,kat,i)-dedzic
                           fat(1,lat,i)=fat(1,lat,i)-dedxid
                           fat(2,lat,i)=fat(2,lat,i)-dedyid
                           fat(3,lat,i)=fat(3,lat,i)-dedzid
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
c Leapfrog algorithm.
         DO i=1,nmol
            DO iat=1,nat(i)
               vsq=0.
               DO idim=1,3
                  vkin=vat(idim,iat,i)+fat(idim,iat,i)/dm(iat,i)*dt/2
                  vat(idim,iat,i)=vat(idim,iat,i)
     *                 +fat(idim,iat,i)/dm(iat,i)*dt
                  rat(idim,iat,i)=rat(idim,iat,i)+vat(idim,iat,i)*dt
                  vsq=vsq+vkin*vkin
               END DO
               ekin=ekin+0.5*vsq*dm(iat,i)
            END DO
         END DO
c
	 nattot=0
	 DO i=1,nmol
	    nattot=nattot+nat(i)
	 END DO
         tc=ekin*2./3./nattot/boltz
         vscale=1.
         IF (equil) vscale=DSQRT(temp/tc)
         DO i=1,nmol
            DO iat=1,nat(i)
               DO idim=1,3
                  vat(idim,iat,i)=vat(idim,iat,i)*vscale
               END DO
            END DO
         END DO
c
c Move back into the center of the box.
c
         DO i=1,nmol
            DO idim=1,3
               IF (rat(idim,2,i).LT.0.) THEN
                  DO iat=1,nat(i)
                     rat(idim,iat,i)=rat(idim,iat,i)+boxl(idim)
                  END DO
               END IF
               IF (rat(idim,2,i).GT.boxl(idim)) THEN
                  DO iat=1,nat(i)
                     rat(idim,iat,i)=rat(idim,iat,i)-boxl(idim)
                  END DO
               END IF
            END DO
         END DO
c
c For scattering simulation, track scattered molecule.
c
         IF (scatsim) THEN
c            OPEN (100+iscat,
c     *           file='zposition'//int2char(iscat),status='unknown')
c           OPEN (100,file='zposition',status='unknown')
           WRITE (100,*) istep*dt,rat(3,2,iscatmol)
         END IF
c
c For scattering simulation, set timer.
c
         IF (scatsim) THEN
            IF (istep.EQ.1) THEN
               maxstep=1000000
               iflag=0
            END IF
            DO i=1,nmol
               IF (i.NE.iscatmol) THEN
                  rij(1)=rat(1,2,iscatmol)-rat(1,2,i)
                  rij(2)=rat(2,2,iscatmol)-rat(2,2,i)
                  rij(3)=rat(3,2,iscatmol)-rat(3,2,i)
                  r=DSQRT(rij(1)**2+rij(2)**2+rij(3)**2)
                  IF (r.LE.4/1D10.AND.iflag.EQ.0) THEN
                     maxstep=istep+200000
                     iflag=1
                  END IF
               END IF
            END DO
            IF (istep.EQ.maxstep) GOTO 15
         END IF

                  
c
      vtot=vc+vlj+vb+vtor
      WRITE (20,*) istep*dt,vtot,ekin+vtot

c
      END DO
      END IF
 15   CLOSE(20)
      CLOSE(100)
      IF (iscat.LT.nscat) GOTO 5
c
      STOP
      END
      
      SUBROUTINE Torsion(xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid
     &  ,V_1,V_2,V_3,e,dedxia,dedyia,dedzia,dedxib,dedyib,dedzib,dedxic,
     &  dedyic,dedzic,dedxid,dedyid,dedzid)
      IMPLICIT REAL*8 (a-h,o-z)

       xba = xib - xia
       yba = yib - yia
       zba = zib - zia
       xcb = xic - xib
       ycb = yic - yib
       zcb = zic - zib
       xdc = xid - xic
       ydc = yid - yic
       zdc = zid - zic
c
       dot1=xba*xcb+yba*ycb+zba*zcb
       dot2=xcb*xdc+ycb*ydc+zcb*zdc
c
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
             IF (cosine.LE.-1d0) THEN
                ang=acos(-1d0)
             ELSE IF (cosine.GE.1d0) THEN
                ang=acos(1d0)
             ELSE
c                ang=dacos(cosine)
                ang=asin(sine)
             END IF
             END IF
c     set the torsional parameters for this angle
c
               v1 = V_1
               v2 = V_2
               v3 = V_3

c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
c
               phi1 = 1.0d0 + cosine
               phi2 = 1.0d0 - cosine2
               phi3 = 1.0d0 + cosine3
c
               dphi1 = -sine
               dphi2 = 2.0d0 *sine2
               dphi3 = -3.0d0 *sine3

c     calculate torsional energy and master chain rule term
c
               e =(v1*phi1 + v2*phi2 + v3*phi3)

               dedphi =(v1*dphi1 + v2*dphi2 + v3*dphi3)

c     chain rule terms for first derivative components
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib

               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute first derivative components for this angle
c
               dedxia = zcb*dedyt - ycb*dedzt
               dedyia = xcb*dedzt - zcb*dedxt
               dedzia = ycb*dedxt - xcb*dedyt
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     increment the total torsional angle energy and gradient
c

               RETURN
               END


      FUNCTION rnd(seed)
      IMPLICIT REAL*4(A-H,O-Z)
      REAL*8 rnd
      INTEGER*4 seed,nseed,tick
      SAVE tick,nseed
      DATA tick/0/
c     
      IF (tick.EQ.0) THEN
         tick=1
         nseed=seed
      ENDIF
      ix=nseed
      k1=ix/127773
      ix=16807*(ix-k1*127773)-k1*2836
c
      IF (ix .LT. 0) ix=ix+2147483647
c
      nseed=ix
      xx=ix*4.6566128e-10
      rnd=DBLE(xx)
c
      RETURN
      END
c
c     This function is copied from W H Press, et al., "Numerical
c     Recipes", Cambidge. 1986. pp. 202 - 203. See this reference for
c     deatils on this function. It generates a pseudo-random number in a
c     Gaussian distribution. To improve the randomness of the generator,
c     use a function other than the simple RND function for the
c     initial source of the random number.
c
      FUNCTION gauss(seed)
      IMPLICIT NONE
c
      REAL*8 gauss,num1,num2,v1,v2,tmp,r2,rnd
      INTEGER*4 seed,called
      SAVE num2,called
      DATA called,num2/0,0.e0/
c
c     If the function was called and num2 unused, then use num2.
c
      IF (called.EQ.1) THEN
         gauss=num2
         called=0
         GOTO 100
      ENDIF
c
c     Otherwise, calculate two new random deviates.
c
 10   v1=2.d0*DBLE(rnd(seed))-1.d0
      v2=2.d0*DBLE(rnd(seed))-1.d0
      r2=v1*v1 + v2*v2
      IF(r2.GE.1.d0) GOTO 10
c
      tmp =DSQRT( -2.d0 * LOG(r2)/r2 )
      num1=tmp*v1
      num2=tmp*v2
      gauss=num1
      called=1
 100  RETURN
      END
c
c  Converts an integer to a character... (written by JR Schmidt)
c
      character*3 function int2char(n)
      int2char = ""
      do npow = 0, 9
         if (10**npow .le. n) then
            int2char = char(48 + mod(n/10**npow, 10))//int2char
         endif
      enddo
      end
