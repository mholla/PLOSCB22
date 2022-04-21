************************************************************************
! 
! User element for cell migration, and large deformation in 2D or 3D.
!  This is for plane strain, axisymetric, and 3D.
! 
! Solution variables (or nodal variables) are the displacements and the
! cell density 
! 
! This subroutine is for the following element types
!  > two-dimensional 4 node isoparametric element as shown below
!       with 1pt (reduced) or 4pt (full) gauss integration.
!  > three-dimensional 8 node isoparametric element as shown below
!       with 1pt (reduced) or 8pt (full) gauss integration.
! 
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
! 
!  Mechanical, traction- and pressure-type boundary conditions 
!   may be applied to the dummy mesh using the Abaqus built-in 
!   commands *Dload or *Dsload.
! 
! 
! 
! 
!     
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
! 
! 
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
! 
! Shuolun Wang, May 2021 @ND 
! Version 4 has 3 species and in full 3D 
! 
! Work along with: 
! 
! 2D input:
! 2DBarMultiCell_Deform_v3.inp
! 2DCircleMultiCell_Deform_v3.inp
! 
! 3D input:
! 3DBarMultiCell_Deform_v4.inp
! 3DBallMultiCell_Deform_v4.inp
***********************************************************************
! 
! User element statement in the input file (set ? values as needed):
! 
!  2D elements
!  *User Element,Nodes=4,Type=U?,Iproperties=3,Properties=18,Coordinates=2,Variables=<nVars>,Unsymm
!  1,2,11
! 
!  3D elements
!  *User Element,Nodes=8,Type=U3,Iproperties=3,Properties=18,Coordinates=3,Variables=<nVars>,Unsymm
!  1,2,3,11
! 
! 
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!       1) alpha
! 
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nIntPt
!          svars(1+j) = alpha_tau ---- extent of reaction at integ pt k
!          j = j + nlSdv
!       end loop over k
! 
!     In the input file, set 'User output variables'= number of global SDV's
! 
!     In the input file, set 'ngSdv'= number of global SDV's
! 
!     In the input file, set 'nlSdv'= number of local SDV's
! 
!     In the input file, set 'nInt'= number of volume integration points
! 
! 
!     Material Properties Vector
!     --------------------------------------------------------------
! 
! 
!***********************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=6000) !2500=2D job, 10000=3D job
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=10000) 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      end module global

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.
     
      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
      uvar(1) = globalSdv(noel-ElemOffset,npt,1) ! source
      uvar(2) = globalSdv(noel-ElemOffset,npt,2) ! velocity
      uvar(3) = globalSdv(noel-ElemOffset,npt,3) ! stiffness
      uvar(4) = globalSdv(noel-ElemOffset,npt,4) ! growth parameter
      uvar(5) = globalSdv(noel-ElemOffset,npt,5) ! growth parameter
      uvar(6) = globalSdv(noel-ElemOffset,npt,6) ! elastic Jacobian
      uvar(7) = globalSdv(noel-ElemOffset,npt,7) ! growth Jacobian     
      uvar(8) = globalSdv(noel-ElemOffset,npt,8) ! Effective total stretch    
      uvar(9) = globalSdv(noel-ElemOffset,npt,9) ! Effective elastic stretch                        
      uvar(10) = globalSdv(noel-ElemOffset,npt,10) ! Trace of Cauchy stress            
      uvar(11) = globalSdv(noel-ElemOffset,npt,11) ! total density  
      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim,nInt,nIntS,i,j,k,l
      character*256 jobName,outDir,fileName


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter(nIntS=1) ! number of surface integration points
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !----------------------------------------------------------------
      ! 
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')


      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !  which are any of the following (64,65,72,73)
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.62).or.(lflags(1).eq.63).or.
     +     (lflags(1).eq.71).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and chekc the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         write(80,*) 'Abaqus does not have the right procedure'
         write(80,*) 'go back and chekc the procedure type'
         write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear purturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
      !
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! 
      ! Call the paricular element to perform the analysis
      !
      ! Obtain the number of volume integration points
      nInt = jprops(3)
      !
      if(jtype.eq.1) then
         !
         ! This is a plane strain analysis
         !
         nDim = 2
         call UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)
         !
         !
      elseif(jtype.eq.2) then
         !
         ! This is an axisymmetric analysis
         !
         nDim = 2
         call UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)
         !
         !
      elseif(jtype.eq.3) then
         !
         ! This is a 3D analysis
         !
         nDim = 3
         call U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)
         !
         !
      else
         !
         ! We have a problem...
         !
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
         !
      endif
      !
      ! Done with this element, RHS and AMATRX already returned
      !  as output from the specific element routine called
      !
      !----------------------------------------------------------------


      return
      end subroutine uel

************************************************************************

      subroutine UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),INew(nNode)
      real*8 IOld(nNode),dI(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nIntS

      logical crosslink,cleave

      real*8 Iden(3,3),Le,theta0,Ru(2*nNode,1),RI(nNode,1),body(3),C0
      real*8 Kuu(2*nNode,2*nNode),KII(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),ds,flux
      real*8 sh(nNode),detMapJ,dsh(nNode,2),detMapJC,umeror,epsilon
      real*8 dshC(nNode,2),I_tau,I_t,dIdX(nDim,1),F_tau(3,3),normal(2,1)
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3),F_starInv(3,3),sigma,direction(3,1)
      real*8 Smat(3,1),Bmat(3,2*nNode),BodyForceRes(2*nNode,1),Qmat(4,4)
      real*8 Gmat(4,2*nNode),G0mat(4,2*nNode),Amat(4,4),wS(nIntS)
      real*8 xLocal(nIntS),yLocal(nIntS),KIu(nNode,2*nNode),detF_t
      real*8 KuI(2*nNode,nNode),Nvec(1,nNode),ResFac,AmatUI(3,1),TanFac
      real*8 SpUIMod(3,3),SpIUMod2(3,3,3),SpIUMod1(3,3),AmatIU2(2,4)
      real*8 Lmt,AmatIU1(1,4),SpIUModFac2(nDim,nDim),SpIUMod(3,3)
      real*8 DsigmaDI,alpha0,alpha_t,alpha_tau,test,SUPG
      real*8 directionVec(nDim,1),DiffMat(nDim,nDim)
      real*8 AmatIU(1,4),dII(nNode),IInew(nNode)
      real*8 Rc(nNode,1),Kcc(nNode,nNode), KIc(nNode,nNode),KcI(nNode,nNode)
      real*8 cNew(nNode),dc(nNode),cOld(nNode)
      real*8 c_tau,c_t,dcdX(nDim,1),dcdt,diff,eps_r,sigma_c,DcDI,DIDc
      real*8 Kcu(nNode,2*nNode),Kuc(2*nNode,nNode)
      real*8 mono_t,mono_tau,DcdotDI,mono_t0,Gshear,I_tau_here
      real*8 DegreeX
      real*8 qflux(nDim,1),Source,dqdc(nDim,1),Diffusivity
      real*8 coordintC(nDim,1),coordint(nDim,1),velocity,velocity2
      real*8 dJdt,AmatUC(3,1),SpUCMod(3,3),AmatCU(2,4),SpCUMod(3,3,3)
      real*8 mu,k_normal,k_para,Je_tau,detFg,EffStretch,EffStretche,TrT
      real*8 Kcc2(nNode,nNode),Kc2u(nNode,2*nNode),Kc2c(nNode,nNode)
      real*8 Kc2c2(nNode,nNode),c2New(nNode),dc2(nNode),c2Old(nNode)
      real*8 c2_tau,c2_t,dc2dX(nDim,1),dc2dt,Rc2(nNode,1)
      real*8 qflux2(nDim,1),source2,AmatUC2(3,1),SpUCMod2(3,3)
      real*8 kuc2(2*nNode,nNode),AmatCU2(2,4),SpCUMod2(3,3,3)
      real*8 Diffusivity2,dqdc2(nDim,1)
      real*8 Rc3(nNode,1),Kuc3(2*nNode,nNode),Kcc3(nNode,nNode)
      real*8 Kc2c3(nNode,nNode),Kc3u(nNode,2*nNode)
      real*8 Kc3c(nNode,nNode),Kc3c2(nNode,nNode),Kc3c3(nNode,nNode)
      real*8 c3New(nNode),dc3(nNode),c3Old(nNode)
      real*8 c3_tau,c3_t,dc3dX(nDim,1),dc3dt,qflux3(nDim,1),Source3
      real*8 AmatCU3(2,4),SpCUMod3(3,3,3),Diffusivity3,dqdc3(nDim,1)
      real*8 velocity3,err_arry12(4,1),err_12,err_arry13(4,1),err_13
      real*8 test1(4,1),test2(4,1),test3(4,1),c_total






      
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      ! 
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point
      

      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- UPE4 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '---------- nlSdv=',nlSdv
         write(*,*) '-------------------------------------------------'
      endif

      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain initial conditions
      !
c      theta0 = props(4)


      

      ! Initialize the residual and tangent matrices to zero.
      ! 
      Ru  = zero
      Rc = zero
      Rc2 = zero ! update
      Rc3 = zero ! update

      test1 = zero
      test2 = zero
      test3 = zero 

      Kuu = zero
      Kuc = zero
      Kuc2 = zero ! update 
      Kuc3 = zero ! update 



      Kcu = zero
      Kcc = zero
      Kcc2 = zero ! update 
      Kcc3 = zero

      Kc2u = zero ! update 
      Kc2c = zero ! update
      Kc2c2 = zero ! update 
      Kc2c3 = zero ! update 

      Kc3u = zero ! update 
      Kc3c = zero ! update 
      Kc3c2 = zero ! update 
      Kc3c3 = zero ! update 


      Energy = zero


      ! Body forces
      !
      body(1:3) = zero


      ! Obtain nodal displacements and cell density 
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         ! the neuron cell density goes in here 
         k = k + 1
         cNew(i) = Uall(k)
         dc(i) = DUall(k,1)
         cOld(i) = cNew(i) - dc(i) 
         ! the neuron cell2 density goes in here
         k = k + 1
         c2New(i) = Uall(k)
         dc2(i) = DUall(k,1)
         c2Old(i) = c2New(i) - dc2(i)
         ! the neuron cell2 density goes in here
         k = k + 1
         c3New(i) = Uall(k)
         dc3(i) = DUall(k,1)
         c3Old(i) = c3New(i) - dc3(i)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


c      !
c      ! displacement increment, based on element diagonal
c      !
c      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
c     +     ((coordsC(2,1)-coordsC(2,3))**two))
c      !
c      do i=1,nNode
c         do j=1,nDim
c            if(dabs(du(i,j)).gt.10.0*Le) then
c               pnewdt = 0.5
c               return
c            endif
c         enddo
c      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the
      !  `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif



      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method. `Tau' represents the end of the increment
      !  and `t' the previous increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for plane-strain
      !
      Fc_tau(3,3) = one
      Fc_t(3,3) = one
      !
      ! 2D plane-strain implementation detF
      !
      detFc_t = Fc_t(1,1)*Fc_t(2,2) - Fc_t(1,2)*Fc_t(2,1)
      detFc_tau = Fc_tau(1,1)*Fc_tau(2,2) - Fc_tau(1,2)*Fc_tau(2,1)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration
         elseif(nInt.eq.9) then
            call xint2D9pt(xi,w,nIntPt) ! 9-pt integration
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt
         
         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif





         ! Obtain the neuron cell density and its derivative's at 
         !  this intPt at the begining and end of the incrment
         !
         c_tau = zero
         c_t = zero
         dcdX = zero
         dcdt = zero 
         do k=1,nNode
            c_tau = c_tau + cNew(k)*sh(k)
            c_t   = c_t + cOld(k)*sh(k)
            do i=1,nDim
               dcdX(i,1) = dcdX(i,1) + cNew(k)*dshC(k,i)
            enddo
         enddo
         dcdt = (c_tau - c_t)/dtime 


         !cell2         
         c2_tau = zero
         c2_t = zero
         dc2dX = zero
         dc2dt = zero 
         do k=1,nNode
            c2_tau = c2_tau + c2New(k)*sh(k)
            c2_t   = c2_t + c2Old(k)*sh(k)
            do i=1,nDim
               dc2dX(i,1) = dc2dX(i,1) + c2New(k)*dshC(k,i)
            enddo
         enddo
         dc2dt = (c2_tau - c2_t)/dtime

         !cell3       
         c3_tau = zero
         c3_t = zero
         dc3dX = zero
         dc3dt = zero 
         do k=1,nNode
            c3_tau = c3_tau + c3New(k)*sh(k)
            c3_t   = c3_t + c3Old(k)*sh(k)
            do i=1,nDim
               dc3dX(i,1) = dc3dX(i,1) + c3New(k)*dshC(k,i)
            enddo
         enddo
         dc3dt = (c3_tau - c3_t)/dtime


         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for plane-strain 
         !
         F_tau(3,3) = one
         F_t(3,3) = one
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            !  2D plane-strain implementation
            !
            detF_t = F_t(1,1)*F_t(2,2) - F_t(1,2)*F_t(2,1)
            detF_tau = F_tau(1,1)*F_tau(2,2) - F_tau(1,2)*F_tau(2,1)
            do i=1,nDim
               do j=1,nDim
                  F_tau(i,j) =((detFc_tau/detF_tau)**half)*F_tau(i,j)
                  F_t(i,j) = ((detFc_t/detF_t)**half)*F_t(i,j)
               enddo
            enddo
         endif
         call mdet(F_tau,detF)






c         ! Obtain state variables from previous increment
c         !
c         !
c         if((kstep.eq.1).and.(kinc.le.1)) then
c            !
c            ! this is the first increment, of the first step
c            !  give initial conditions
c            !
c            mono_t   = mono_t0 
c            !
c            !
c         else
c            !
c            ! read old values
c            !
c            mono_t  = svars(1+jj)
c         endif


         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo     

         ! obtain coordinates at this integration point
         ! 
         coordintC = matmul(coordsC,Nvec) ! current
         coordint  = matmul(coords,Nvec) ! old 

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         ! 
         ! Perform the constitutive time integration at this integ. point 
         ! 
   
         call integ(props,nprops,dtime,TIME,F_tau,F_t,
     +               c_tau,c2_tau,c3_tau,c_total,
     +               dcdX,dc2dX,dc3dX,coordintC,coordint,
     +               qflux,qflux2,qflux3,
     +               source,source2,source3,
     +               velocity,velocity2,velocity3,
     +               Diffusivity,Diffusivity2,Diffusivity3,
     +               dqdc,dqdc2,dqdc3,
     +               T_tau,dJdt,
     +               SpTanMod,SpUCMod,
     +               SpCUMod,SpCUMod2,SpCUMod3,
     +               mu,k_normal,k_para,Je_tau,detFg,EffStretch,
     +               EffStretche,TrT)                   

         ! I need  qflux, source,Diffusivity,dqdc
         ! I need  T_tau, dJdt, SpTanMod, SpUCMod, SpCUMod
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




c         ! Save the state variables at this integ point
c         !  at the end of the increment
c         !
c         svars(1+jj) =  mono_tau
         jj = jj + nlSdv ! setup for the next intPt
         

         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         ! 
         globalSdv(jelem,intPt,1) = source ! source term profile
         globalSdv(jelem,intPt,2) = velocity ! velocity profile 
         globalSdv(jelem,intPt,3) = mu ! stiffness profile          
         globalSdv(jelem,intPt,4) = k_normal ! growth parameter          
         globalSdv(jelem,intPt,5) = k_para ! growth parameter             
         globalSdv(jelem,intPt,6) = Je_tau ! elastic Jacobian             
         globalSdv(jelem,intPt,7) = detFg ! growth Jacobian             
         globalSdv(jelem,intPt,8) = EffStretch !  Effective total Stretch   
         globalSdv(jelem,intPt,9) = EffStretche !  Effective elastic Stretch            
         globalSdv(jelem,intPt,10) = TrT !  trace of Cauchy stress                           
         globalSdv(jelem,intPt,11) = c_total !  total cell density                       







         ! Compute/update the displacement residual vector
         ! I need: T_tau
         ! 
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(1,2)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,2+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        +matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )
 

c         write(*,*),'Ru=',Ru





         ! Compute/update the neuron cell density residual
         ! I need: qflux, source
         ! I need: dJdt




c         Rc = Rc + detmapJC*w(intpt)*
c     +        (
c     +       +transpose(Nvec)*(dcdt + c_tau*dJdt/detF)
c     +       +matmul(dshC,qflux)
c     +       -transpose(Nvec)*Source
c     +        )  

         test1 = test1 + detmapJC*w(intpt)*
     +        (
     +       +transpose(Nvec)*(dcdt + c_tau*dJdt/detF)
     +       +matmul(dshC,qflux)
     +       -transpose(Nvec)*Source
     +        )       

c         Rc2 = Rc2 + detmapJC*w(intpt)*
c     +        (
c     +       +transpose(Nvec)*(dc2dt + c2_tau*dJdt/detF)
c     +       +matmul(dshC,qflux2)
c     +       -transpose(Nvec)*Source2
c     +        ) 

         test2 = test2 + detmapJC*w(intpt)*
     +        (
     +       +transpose(Nvec)*(dc2dt + c2_tau*dJdt/detF)
     +       +matmul(dshC,qflux2)
     +       -transpose(Nvec)*Source2
     +        ) 

c         Rc3 = Rc3 + detmapJC*w(intpt)*
c     +        (
c     +       +transpose(Nvec)*(dc2dt + c2_tau*dJdt/detF)
c     +       +matmul(dshC,qflux2)
c     +       -transpose(Nvec)*Source2
c     +        ) 

         test3 = test3 + detmapJC*w(intpt)*
     +        (
     +       +transpose(Nvec)*(dc3dt + c3_tau*dJdt/detF)
     +       +matmul(dshC,qflux3)
     +       -transpose(Nvec)*Source3
     +        )


c         err_arry23 = test1 - test2
c         err_23 = sqrt(err_arry23(1,1)**2.0+
c     +                 err_arry23(2,1)**2.0+
c     +                 err_arry23(3,1)**2.0+
c     +                 err_arry23(4,1)**2.0)


c         write(*,*),'err_23=',err_23













         ! Compute/update the displacement tangent matrix
         ! I need:SpTanMod
         ! 
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(4,2+nDim*(kk-1)) = dshC(kk,2)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(4,2+nDim*(kk-1)) = dshC0(kk,2)
         enddo

         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,1,2)
         Amat(1,4) = SpTanMod(1,1,2,2)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,1,2)
         Amat(2,4) = SpTanMod(2,1,2,2)
         Amat(3,1) = SpTanMod(1,2,1,1)
         Amat(3,2) = SpTanMod(1,2,2,1)
         Amat(3,3) = SpTanMod(1,2,1,2)
         Amat(3,4) = SpTanMod(1,2,2,2)
         Amat(4,1) = SpTanMod(2,2,1,1)
         Amat(4,2) = SpTanMod(2,2,2,1)
         Amat(4,3) = SpTanMod(2,2,1,2)
         Amat(4,4) = SpTanMod(2,2,2,2)


         Qmat = zero
         Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
         Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            

         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! This is the tangent using the F-bar method with the
            !  4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           - matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           - matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +           )
         else 
            ! 
            ! This is the tangent not using the F-bar method with all
            !  other elements
            ! 
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           -matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
         endif



         ! Compute/update the displacement - cell density tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         ! I need: SpUCMod
         ! 
         AmatUC = zero
         AmatUC(1,1) = SpUCMod(1,1)
         AmatUC(2,1) = SpUCMod(2,2)
         AmatUC(3,1) = SpUCMod(1,2)
         ! 
         kuc = Kuc + detMapJC*w(intpt)*
     +        (
     +        -matmul(matmul(transpose(Bmat),AmatUC),Nvec)
     +        )


c         write(*,*),'kuc=', kuc



c         AmatUC2 = zero
c         AmatUC2(1,1) = SpUCMod2(1,1)
c         AmatUC2(2,1) = SpUCMod2(2,2)
c         AmatUC2(3,1) = SpUCMod2(1,2)
c         !
c         kuc2 = Kuc2 + detMapJC*w(intpt)*
c     +        (
c     +        -matmul(matmul(transpose(Bmat),AmatUC2),Nvec)
c     +        )

c         kuc2 = kuc
c         kuc2 = zero
c         write(*,*),'kuc2=', kuc2



         ! Compute/update the cell density  - displacement tangent matrix
         ! The F-bar method will have some effect, however we neglect that here.
         ! I need SpCUMod
         ! 
         AmatCU = zero
         AmatCU(1,1) = SpCUMod(1,1,1)
         AmatCU(1,2) = SpCUMod(1,2,1)
         AmatCU(1,3) = SpCUMod(1,1,2)
         AmatCU(1,4) = SpCUMod(1,2,2)
         AmatCU(2,1) = SpCUMod(2,1,1)
         AmatCU(2,2) = SpCUMod(2,2,1)
         AmatCU(2,3) = SpCUMod(2,1,2)
         AmatCU(2,4) = SpCUMod(2,2,2)
         !
         Kcu = Kcu + detMapJC*w(intpt)*
     +        (
     +        +matmul(matmul(dshC,AmatCU),Gmat)
     +        )


c         write(*,*),'kcu=', kc
         AmatCU2 = zero
         AmatCU2(1,1) = SpCUMod2(1,1,1)
         AmatCU2(1,2) = SpCUMod2(1,2,1)
         AmatCU2(1,3) = SpCUMod2(1,1,2)
         AmatCU2(1,4) = SpCUMod2(1,2,2)
         AmatCU2(2,1) = SpCUMod2(2,1,1)
         AmatCU2(2,2) = SpCUMod2(2,2,1)
         AmatCU2(2,3) = SpCUMod2(2,1,2)
         AmatCU2(2,4) = SpCUMod2(2,2,2)
         !
         Kc2u = Kc2u + detMapJC*w(intpt)*
     +        (
     +        +matmul(matmul(dshC,AmatCU2),Gmat)
     +        )


c         write(*,*),'kc2u=', kc2u

         AmatCU3 = zero
         AmatCU3(1,1) = SpCUMod3(1,1,1)
         AmatCU3(1,2) = SpCUMod3(1,2,1)
         AmatCU3(1,3) = SpCUMod3(1,1,2)
         AmatCU3(1,4) = SpCUMod3(1,2,2)
         AmatCU3(2,1) = SpCUMod3(2,1,1)
         AmatCU3(2,2) = SpCUMod3(2,2,1)
         AmatCU3(2,3) = SpCUMod3(2,1,2)
         AmatCU3(2,4) = SpCUMod3(2,2,2)
         !
         Kc3u = Kc3u + detMapJC*w(intpt)*
     +        (
     +        +matmul(matmul(dshC,AmatCU3),Gmat)
     +        )





         ! Compute/update the neuron cell density tangent matrix
         ! I need Diffusivity,dqdc
         ! I need dJdt

         Kcc = Kcc + detmapJC*w(intPt)*
     +        (
     +        -(1.0/dtime + dJdt/detF)*matmul(transpose(Nvec),Nvec)
     +        -Diffusivity*matmul(dshC,transpose(dshC))
     +        -matmul(matmul(dshC,dqdc),Nvec)
     +        )

         Kc2c2 = Kc2c2 + detmapJC*w(intPt)*
     +        (
     +        -(1.0/dtime + dJdt/detF)*matmul(transpose(Nvec),Nvec)
     +        -Diffusivity2*matmul(dshC,transpose(dshC))
     +        -matmul(matmul(dshC,dqdc2),Nvec)
     +        )

         Kc3c3 = Kc3c3 + detmapJC*w(intPt)*
     +        (
     +        -(1.0/dtime + dJdt/detF)*matmul(transpose(Nvec),Nvec)
     +        -Diffusivity3*matmul(dshC,transpose(dshC))
     +        -matmul(matmul(dshC,dqdc3),Nvec)
     +        )






      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface terms here
      !
      !
      ! This is stil experimental and if improved can be used as an 
      ! adaptive essential boundary condition
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the flux/traction
            !  acts on is the flux/traction ``label''
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! flux magnitude
            
            if((face.ge.1).and.(face.le.4)) then
               !
               ! fluid flux applied
               !
               if(face.eq.1) then
                  faceFlag = 1
               elseif(face.eq.2) then
                  faceFlag = 2
               elseif(face.eq.3) then
                  faceFlag = 3
               elseif(face.eq.4) then
                  faceFlag = 4
               endif
               !
               if(nIntS.eq.1) then
                  call xintSurf2D1pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.2) then
                  call xintSurf2D2pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.3) then
                  call xintSurf2D3pt(faceFlag,xLocal,yLocal,wS)
               else
                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
                  write(80,*) 'Invalid nIntS points, nIntS=',nIntS
                  call xit
               endif
               !
               ! loop over integ points
               !
               do ii=1,nIntS
                  !
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),faceFlag,
     +                 coordsC,sh,ds,normal)




                  ! Obtain the light intensity at this integ point
                  !
                  I_tau = zero
                  do n=1,nNode
                     I_tau = I_tau + INew(n)*sh(n)
                  enddo


                  ! Check for light hitting this surface, if test is
                  !  less than zero, then light is hitting this surface
                  !
                  test = direction(1,1)*normal(1,1) 
     +                 + direction(2,1)*normal(2,1)
                  if(test.lt.zero) then
c                     test = (1 - tanh(test/.2)) * 1.5
                      test = 1.0
                  else
                     test = zero
                  endif



                  !
                  ! Modify the residual, loop over nodes, recall
                  !  sh(n)=0 when n is not on this face
                  !
                  do n=1,nNode
                     RI(n,1) = RI(n,1) 
     +                    + wS(ii)*ds*sh(n)*(I_tau - flux)*test
                  enddo


                  ! Modify the tangent
                  !
                  do n=1,nNode
                     do m=1,nNode
                        KII(n,m) = KII(n,m) - w(ii)*ds*sh(n)*sh(m)*test
                     enddo
                  enddo

               enddo ! loop over nIntS
               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over flux and traction terms
      !----------------------------------------------------------------


      ! I turn off the deformation related residual and tangents 
      ! 
c      Ru  = zero
c      Rc = zero
c      Rc2 = zero ! update
c      Rc3 = zero ! update
      

c      Kuu = zero
c      Kuc = zero
c      Kuc2 = zero ! update 
c      Kuc3 = zero ! update 



c      Kcu = zero
c      Kcc = zero
c      Kcc2 = zero ! update 
c      Kcc3 = zero

c      Kc2u = zero ! update 
c      Kc2c = zero ! update
c      Kc2c2 = zero ! update 
c      Kc2c3 = zero ! update 

c      Kc3u = zero ! update 
c      Kc3c = zero ! update 
c      Kc3c2 = zero ! update 
c      Kc3c3 = zero ! update 








      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      ! 
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,test1,test2,test3,
     +     Kuu,Kuc,Kuc,Kuc,
     +     Kcu,Kcc,Kcc2,kcc3,
     +     Kc2u,Kc2c,Kc2c2,kc2c3,
     +     Kc3u,Kc3c,Kc3c2,Kc3c3,
     +     rhs,amatrx)     
      ! 
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------



      return
      end subroutine UPE4

************************************************************************

      subroutine UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)
     
      use global

      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),INew(nNode)
      real*8 IOld(nNode),dI(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,faceFlag,nIntS

      real*8 Iden(3,3),Le,theta0,Ru(2*nNode,1),RI(nNode,1),body(3)
      real*8 Kuu(2*nNode,2*nNode),KII(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,AR0
      real*8 ARc,AR_t,Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),AR
      real*8 sh(nNode),detMapJ,dsh(nNode,2),detMapJC,Lmt,umeror
      real*8 dshC(nNode,2),I_tau,I_t,dIdX(2,1),F_tau(3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3)
      real*8 Smat(3,1),BodyForceRes(2*nNode,1),SpIUMod(3,3,3)
      real*8 SmatAx(4,1),BodyForceResAx(2*nNode,1),dTRdF(3,3,3,3)
      real*8 BmatAx(4,2*nNode),Gmat(4,2*nNode),G0mat(4,2*nNode),flux,ds
      real*8 Amat(4,4),Qmat(4,4),AmatAx(5,5),QmatAx(5,5),yLocal(nIntS)
      real*8 G0matAx(5,2*nNode),GmatAx(5,2*nNode),xLocal(nIntS),detF_t
      real*8 wS(nIntS),KuI(2*nNode,nNode),KIu(nNode,2*nNode),ResFac
      real*8 TanFac,Nvec(1,nNode),AmatUI(4,1),SpIUModFac(3,3)
      real*8 AmatIU(2,5),SpUIMod(3,3)

      real*8 sigma,direction(2,1),normal(2,1),test,SUPG,hXi(2,1)
      real*8 hEta(2,1),eXi(2,1),eEta(2,1)
      
      real*8 Fc1_int(3,3),Fc1_int_inv(3,3),F_c1(3,3),detF_c1,alpha0
      real*8 alpha_t,alpha_tau
      
      real*8 CUB_t,CUB_tau,CB_t,CB_tau,C0
      real*8 epsilon

      real*8 zero,one,two,half,Pi,three,third
      
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- UAX4 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt= ',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '---------- nlSdv=',nlSdv
         write(*,*) '-------------------------------------------------'
      endif


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain initial conditions
      !
      theta0 = props(5)
      C0 = props(11)
      epsilon = props(17)

      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      RI = zero
      Kuu = zero
      KII = zero
      KuI = zero
      KIu = zero
      Energy = zero


      ! Body forces
      !
      body(1:3) = zero


      ! Obtain nodal displacements and light intensities
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         INew(i) = Uall(k)
         dI(i) = DUall(k,1)
         IOld(i) = INew(i) - dI(i)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose any time-stepping changes on the increments of
      !  light intensity or displacement if you wish
      !
      ! light intensity increment
      !
      do i=1,nNode
         if(dabs(dI(i)).gt.1.d6) then
            pnewdt = 0.5
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo



      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the
      !  `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! For an axisymmetric problem, find the ``r'' that
      !  shows up in the integrals for axisymmetric
      !  and in the F(3,3), the factors of 2Pi are for integrals
      !  i.e., dV = 2 pi r dr dz
      !
      AR0  = zero
      ARc  = zero
      AR_t = zero
      do i=1,nNode
         ! radial coord in ref config at centroid
         AR0  = AR0 + sh0(i)*coords(1,i)
         ! radial coord in current config at centroid
         ARc  = ARc + sh0(i)*(coords(1,i) + u(i,1))
         ! radial coord in current config at centroid in previous step
         AR_t = AR_t + sh0(i)*(coords(1,i) + uOld(i,1))
      enddo



      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method. `Tau' represents the end of the increment
      !  and `t' the previous increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for axisymmetric
      !
      Fc_tau(3,3) = ARc/AR0
      Fc_t(3,3) = AR_t/AR0
      !
      ! axisymmetric implementation detF
      !
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif



      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt

         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! For an axisymmetric problem, find the ``r'' that
         !  shows up in the integrals for axisymmetric
         !  and in the F(3,3), the factors of 2Pi are for integrals
         !  i.e., dV = 2 pi r dr dz
         !
         !
         AR0  = zero
         AR   = zero
         AR_t = zero
         do i=1,nNode
            AR0 = AR0 + sh(i)*coords(1,i)
            AR  = AR  + sh(i)*(coords(1,i) + u(i,1))
            AR_t = AR_t + sh(i)*(coords(1,i) + uOld(i,1))
         enddo
         AR0  = two*Pi*AR0
         AR   = two*Pi*AR
         AR_t = two*Pi*AR_t


         ! Obtain the light intensity and its derivative's at 
         !  this intPt at the begining and end of the incrment
         !
         I_tau = zero
         I_t = zero
         dIdX = zero
         do k=1,nNode
            I_tau = I_tau + INew(k)*sh(k)
            I_t   = I_t + IOld(k)*sh(k)
            do i=1,nDim
               dIdX(i,1) = dIdX(i,1) + INew(k)*dshC(k,i)
            enddo
         enddo



         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for axisymetric, give R/R0
         !
         F_tau(3,3) = AR/AR0
         F_t(3,3) = AR_t/AR0
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)


         ! Obtain state variables from previous increment
         !
         if((kstep.eq.1).and.(kinc.le.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions
            !
            alpha_t  = alpha0
            Fc1_int_inv = Iden
            CUB_t = C0
            CB_t = zero
            !
         elseif((kstep.eq.2).and.(kinc.le.1)) then
            alpha_t = svars(1+jj)
            CUB_t = svars(11+jj)
            CB_t = svars(12+jj)
            !
            ! save the deformation gradient now
            !
            call matInv3D(F_tau,Fc1_int_inv,detF,stat)
            if(stat.eq.0) then
               write(*,*) 'Problem: detF.lt.zero'
               call xit
            endif
            !
         else
            !
            ! read old values
            !
            alpha_t  = svars(1+jj)
            Fc1_int_inv(1,1) = svars(2+jj)
            Fc1_int_inv(1,2) = svars(3+jj)
            Fc1_int_inv(1,3) = svars(4+jj)
            Fc1_int_inv(2,1) = svars(5+jj)
            Fc1_int_inv(2,2) = svars(6+jj)
            Fc1_int_inv(2,3) = svars(7+jj)
            Fc1_int_inv(3,1) = svars(8+jj)
            Fc1_int_inv(3,2) = svars(9+jj)
            Fc1_int_inv(3,3) = svars(10+jj)
            CUB_t = svars(11+jj)
            CB_t = svars(12+jj)
            
         endif

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
         call integ(props,nprops,dtime,
     +        F_tau,I_tau,theta0,
     +        T_tau,SpTanMod,
     +        SpUIMod,SpIUModFac,
     +        sigma,direction,
     +        Fc1_int_inv,
     +        alpha_t,alpha_tau,
     +        time,kstep,CUB_t,CUB_tau,CB_t,CB_tau)
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = alpha_tau
         svars(2+jj) =  Fc1_int_inv(1,1)
         svars(3+jj) =  Fc1_int_inv(1,2)
         svars(4+jj) =  Fc1_int_inv(1,3)
         svars(5+jj) =  Fc1_int_inv(2,1)
         svars(6+jj) =  Fc1_int_inv(2,2)
         svars(7+jj) =  Fc1_int_inv(2,3)
         svars(8+jj) =  Fc1_int_inv(3,1)
         svars(9+jj) =  Fc1_int_inv(3,2)
         svars(10+jj) =  Fc1_int_inv(3,3)
         svars(11+jj) = CUB_tau
         svars(12+jj) = CB_tau
         jj = jj + nlSdv ! setup for the next intPt


         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = alpha_tau   ! extent of reaction
         globalSdv(jelem,intPt,2) = sigma       ! extinction field

         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the chemical reaction.
         !
         Lmt = 0.05d0
         umeror = dabs((alpha_tau - alpha_t)/Lmt)
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0
         else
            pnewdt = half
         endif


         ! Compute/update the displacement residual vector
         !
         SmatAx(1,1) = T_tau(1,1)
         SmatAx(2,1) = T_tau(2,2)
         SmatAx(3,1) = T_tau(1,2)
         SmatAx(4,1) = T_tau(3,3)
         !
         BmatAx = zero
         do kk=1,nNode
            BmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(2,2+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,2+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(4,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo
         !
         BodyForceResAX = zero
         do kk=1,nNode
            BodyForceResAx(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceResAx(2+nDim*(kk-1),1) = sh(kk)*body(2)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*AR*
     +        (
     +        -matmul(transpose(BmatAx),SmatAx)
     +        + BodyForceResAx
     +        )


         ! Compute the stabilization parameter for this element
         !
         ! I did a lot of work by hand, this is Figure 3.4, on page
         !  215 of that article I sent you by Brooks and Hughes
         !
         hXi(1,1) = coordsC(1,3)+coordsC(1,4)-coordsC(1,1)-coordsC(1,2)
         hXi(2,1) = coordsC(2,3)+coordsC(2,4)-coordsC(2,1)-coordsC(2,2)
         hEta(1,1) = coordsC(1,2)+coordsC(1,3)-coordsC(1,1)-coordsC(1,4)
         hEta(2,1) = coordsC(2,2)+coordsC(2,3)-coordsC(2,1)-coordsC(2,4)
         !
         eXi = hXi/dsqrt(hXi(1,1)**2.d0 + hXi(2,1)**2.d0)
         eEta = hEta/dsqrt(hEta(1,1)**2.d0 + hEta(2,1)**2.d0)
         !
         SUPG = (
     +        dabs(direction(1,1)*hXi(1,1) + direction(2,1)*hXi(2,1))
     +        *dabs(direction(1,1)*eXi(1,1) + direction(2,1)*eXi(2,1))
     +        +dabs(direction(1,1)*hEta(1,1) + direction(2,1)*hEta(2,1))
     +        *dabs(direction(1,1)*eEta(1,1) + direction(2,1)*eEta(2,1))
     +        )/dsqrt(15.d0)
     
     

         ! Compute/update the light intensity residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         RI = RI + detmapJC*w(intpt)*
     +        (
         !
         ! standard terms
         !
     +        (
     +        transpose(Nvec)*(matmul(transpose(direction),dIdX)
     +        + sigma*I_tau)
     +        )
         !
         ! numerical diffusion and stabilization terms
         !
     +        + (
     +        epsilon*matmul(dshC,dIdX)
     +        + SUPG*matmul(dshC,direction)*(
     +        matmul(transpose(direction),dIdX) + sigma*I_tau)
     +        )
     +        )


         ! Compute/update the displacement tangent matrix
         !
         GmatAx = zero
         do kk=1,nNode
            GmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(2,2+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(4,2+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(5,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo

         G0matAx = zero
         do kk=1,nNode
            G0matAx(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0matAx(4,2+nDim*(kk-1)) = dshC0(kk,2)
            G0matAX(5,1+nDim*(kk-1)) = sh0(kk)/ARc
         enddo

         AmatAx = zero
         AmatAx(1,1) = SpTanMod(1,1,1,1)
         AmatAx(1,2) = SpTanMod(1,1,2,1)
         AmatAx(1,3) = SpTanMod(1,1,1,2)
         AmatAx(1,4) = SpTanMod(1,1,2,2)
         AmatAx(1,5) = SpTanMod(1,1,3,3)
         AmatAx(2,1) = SpTanMod(2,1,1,1)
         AmatAx(2,2) = SpTanMod(2,1,2,1)
         AmatAx(2,3) = SpTanMod(2,1,1,2)
         AmatAx(2,4) = SpTanMod(2,1,2,2)
         AmatAx(2,5) = SpTanMod(2,1,3,3)
         AmatAx(3,1) = SpTanMod(1,2,1,1)
         AmatAx(3,2) = SpTanMod(1,2,2,1)
         AmatAx(3,3) = SpTanMod(1,2,1,2)
         AmatAx(3,4) = SpTanMod(1,2,2,2)
         AmatAx(3,5) = SpTanMod(1,2,3,3)
         AmatAx(4,1) = SpTanMod(2,2,1,1)
         AmatAx(4,2) = SpTanMod(2,2,2,1)
         AmatAx(4,3) = SpTanMod(2,2,1,2)
         AmatAx(4,4) = SpTanMod(2,2,2,2)
         AmatAx(4,5) = SpTanMod(2,2,3,3)
         AmatAx(5,1) = SpTanMod(3,3,1,1)
         AmatAx(5,2) = SpTanMod(3,3,2,1)
         AmatAx(5,3) = SpTanMod(3,3,1,2)
         AmatAx(5,4) = SpTanMod(3,3,2,2)
         AmatAx(5,5) = SpTanMod(3,3,3,3)

         QmatAx = zero
         QmatAx(1,1) = third*(AmatAx(1,1)+AmatAx(1,4)+AmatAx(1,5)) 
     +        - (two/three)*T_tau(1,1)
         QmatAx(2,1) = third*(AmatAx(2,1)+AmatAx(2,4)+AmatAx(2,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(3,1) = third*(AmatAx(3,1)+AmatAx(3,4)+AmatAx(3,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(4,1) = third*(AmatAx(4,1)+AmatAx(4,4)+AmatAx(4,5))
     +        - (two/three)*T_tau(2,2)
         QmatAx(5,1) = third*(AmatAx(5,1)+AmatAx(5,4)+AmatAx(5,5))
     +        - (two/three)*T_tau(3,3)
         QmatAx(1,4) = QmatAx(1,1)
         QmatAx(2,4) = QmatAx(2,1)
         QmatAx(3,4) = QmatAx(3,1)
         QmatAx(4,4) = QmatAx(4,1)
         QmatAx(5,4) = QmatAx(5,1)
         QmatAx(1,5) = QmatAx(1,1)
         QmatAx(2,5) = QmatAx(2,1)
         QmatAx(3,5) = QmatAx(3,1)
         QmatAx(4,5) = QmatAx(4,1)
         QmatAx(5,5) = QmatAx(5,1)
            

         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! This is the tangent using the F-bar method with the
            !  4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           + matmul(transpose(GmatAx),matmul(QmatAx,
     +           (G0matAx-GmatAx)))
     +           )
         else
            !
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           )
         endif



         ! Compute/update the light intensity tangent matrix
         !
         KII = KII - detmapJC*w(intPt)*
     +        (
         !
         ! standard terms
         !
     +        (
     +        matmul(matmul(dshC,direction),Nvec)
     +        + sigma*matmul(transpose(Nvec),Nvec)
     +        )
         !
         ! numerical diffusion and stabilization terms
         !
     +        + (
     +        epsilon*matmul(dshC,transpose(dshC))
     +        + SUPG*matmul(matmul(dshC,direction),
     +        (transpose(matmul(dshC,direction)) + sigma*Nvec))
     +        )
     +        )



      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface fluid flux terms here
      !
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the flux/traction
            !  acts on is the flux/traction ``label''
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! flux magnitude
            
            if((face.ge.1).and.(face.le.4)) then
               !
               ! fluid flux applied
               !
               if(face.eq.1) then
                  faceFlag = 1
               elseif(face.eq.2) then
                  faceFlag = 2
               elseif(face.eq.3) then
                  faceFlag = 3
               elseif(face.eq.4) then
                  faceFlag = 4
               endif
               !
               if(nIntS.eq.1) then
                  call xintSurf2D1pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.2) then
                  call xintSurf2D2pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.3) then
                  call xintSurf2D3pt(faceFlag,xLocal,yLocal,wS)
               else
                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
                  write(80,*) 'Invalid nIntS points, nIntS=',nIntS
                  call xit
               endif
               !
               ! loop over integ points
               !
               do ii=1,nIntS
                  !
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),faceFlag,
     +                 coordsC,sh,ds,normal)




                  ! Obtain the light intensity at this integ point
                  !
                  I_tau = zero
                  do n=1,nNode
                     I_tau = I_tau + INew(n)*sh(n)
                  enddo


                  ! Check for light hitting this surface, if test is
                  !  less than zero, then light is hitting this surface
                  !
                  test = direction(1,1)*normal(1,1) 
     +                 + direction(2,1)*normal(2,1)
                  if(test.lt.zero) then
c                     test = (1 - tanh(test/.2)) * 1.5
                      test = 1.0
                  else
                     test = zero
                  endif



                  !
                  ! Modify the residual, loop over nodes, recall
                  !  sh(n)=0 when n is not on this face
                  !
                  do n=1,nNode
                     RI(n,1) = RI(n,1) 
     +                    + wS(ii)*ds*sh(n)*(I_tau - flux)*test
                  enddo


                  ! Modify the tangent
                  !
                  do n=1,nNode
                     do m=1,nNode
                        KII(n,m) = KII(n,m) - w(ii)*ds*sh(n)*sh(m)*test
                     enddo
                  enddo




               enddo ! loop over nIntS
               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over flux and traction terms
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,RI,Kuu,KuI,KIu,KII,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


     
      return
      end subroutine UAX4

************************************************************************

      subroutine U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      use global

      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


      real*8 u(nNode,3),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),INew(nNode)
      real*8 IOld(nNode),dI(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,3)
      real*8 coordsC(mcrd,nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,nIntS,faceFlag

      real*8 Iden(3,3),Le,theta0,Ru(3*nNode,1),RI(nNode,1),body(3)
      real*8 Kuu(3*nNode,3*nNode),KII(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,3),dsh0(nNode,3),dshC0(nNode,3),detMapJ0C
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),SUPG
      real*8 sh(nNode),detMapJ,dsh(nNode,3),detMapJC,Lmt,umeror,test
      real*8 dshC(nNode,3),I_tau,I_t,dIdX(3,1),F_tau(3,3),DsigmaDI
      real*8 F_t(3,3),detF_tau,xi(nInt,3),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3),ds,alpha0,alpha_t,alpha_tau,C0,epsilon
      real*8 Smat(6,1),Bmat(6,3*nNode),BodyForceRes(3*nNode,1),flux
      real*8 Gmat(9,3*nNode),G0mat(9,3*nNode),Amat(9,9),Qmat(9,9),dA
      real*8 xLocal(nIntS),yLocal(nIntS),zLocal(nIntS),wS(nIntS),detF_t
      real*8 KuI(3*nNode,nNode),KIu(nNode,3*nNode),Nvec(1,nNode),ResFac
      real*8 AmatUI(6,1),TanFac,AmatIU(1,9),SpUIMod(3,3),SpIUMod(3,3)
      real*8 SpIUModFac(3,3),sigma,direction(3,1),normal(3,1)
      real*8 F_starInv(3,3),DiffMat(nDim,nDim)

      real*8 Rc(nNode,1),Kcc(nNode,nNode),KIc(nNode,nNode)
      real*8 KcI(nNode,nNode),Kcu(nNode,3*nNode),Kuc(3*nNode,nNode)
      real*8 cNew(nNode),dc(nNode),cOld(nNode)
      real*8 c_tau,c_t,dcdX(nDim,1),dcdt,diff,eps_r,sigma_c,DcDI,DIDc
      real*8 mono_t0,mono_t,mono_tau
      real*8 DcdotDI,Gshear,I_tau_here
      real*8 directionVec(nDim,1)
      real*8 DegreeX,test1(nNode,1),test2(nNode,1),test3(nNode,1)
      real*8 Rc2(nNode,1),Rc3(nNode,1)
      real*8 Kuc2(3*nNode,nNode),Kuc3(3*nNode,nNode)
      real*8 Kcc2(nNode,nNode),Kcc3(nNode,nNode)
      real*8 Kc2u(nNode,3*nNode),Kc2c(nNode,nNode)
      real*8 Kc2c2(nNode,nNode),Kc2c3(nNode,nNode)
      real*8 Kc3u(nNode,3*nNode),Kc3c(nNode,nNode)
      real*8 Kc3c2(nNode,nNode),Kc3c3(nNode,nNode) 
      real*8 c2New(nNode),dc2(nNode),c2Old(nNode)
      real*8 c3New(nNode),dc3(nNode),c3Old(nNode)
      real*8 c2_tau,c2_t,dc2dX(nDim,1),dc2dt
      real*8 c3_tau,c3_t,dc3dX(nDim,1),dc3dt
      real*8 coordintC(nDim,1),coordint(nDim,1)
      real*8 c_total,qflux(nDim,1),qflux2(nDim,1)
      real*8 qflux3(nDim,1),source,source2,source3
      real*8 velocity,velocity2,velocity3
      real*8 Diffusivity,Diffusivity2,Diffusivity3
      real*8 dqdc(nDim,1),dqdc2(nDim,1),dqdc3(nDim,1)
      real*8 dJdt,SpUCMod(3,3)
      real*8 SpCUMod(3,3,3),SpCUMod2(3,3,3),SpCUMod3(3,3,3)
      real*8 mu,k_normal,k_para,Je_tau,detFg,EffStretch
      real*8 EffStretche,TrT,AmatUC(6,1),AmatCU(3,9)
      real*8 AmatCU2(3,9),AmatCU3(3,9)







      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)

      logical crosslink,cleave

      character*256 jobName,outDir,fileName

c      write(*,*),'line2280'

      ! Get element parameters
      !
      nlSdv  = jprops(1) ! number of local sdv's per integ point
      ngSdv  = jprops(2) ! number of global sdv's per integ point


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- U3D8 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '---------- nlSdv=',nlSdv
         write(*,*) '-------------------------------------------------'
      endif


c      write(*,*),'line2345'


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      ! 
      Ru  = zero
      Rc = zero
      Rc2 = zero ! update
      Rc3 = zero ! update

      test1 = zero
      test2 = zero
      test3 = zero 

      Kuu = zero
      Kuc = zero
      Kuc2 = zero ! update 
      Kuc3 = zero ! update 



      Kcu = zero
      Kcc = zero
      Kcc2 = zero ! update 
      Kcc3 = zero

      Kc2u = zero ! update 
      Kc2c = zero ! update
      Kc2c2 = zero ! update 
      Kc2c3 = zero ! update 

      Kc3u = zero ! update 
      Kc3c = zero ! update 
      Kc3c2 = zero ! update 
      Kc3c3 = zero ! update 


      Energy = zero


      ! Body forces
      !
      body(1:3) = zero


      ! Obtain nodal displacements and cell density 
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         ! the neuron cell density goes in here 
         k = k + 1
         cNew(i) = Uall(k)
         dc(i) = DUall(k,1)
         cOld(i) = cNew(i) - dc(i) 
         ! the neuron cell2 density goes in here
         k = k + 1
         c2New(i) = Uall(k)
         dc2(i) = DUall(k,1)
         c2Old(i) = c2New(i) - dc2(i)
         ! the neuron cell2 density goes in here
         k = k + 1
         c3New(i) = Uall(k)
         dc3(i) = DUall(k,1)
         c3Old(i) = c3New(i) - dc3(i)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


c      ! Impose any time-stepping changes on the increments of
c      !  light intensity or displacement if you want
c      !
c      ! light intensity increment
c      !
c      do i=1,nNode
c         if(dabs(dI(i)).gt.1.d6) then
c            pnewdt = 0.5
c            return
c         endif
c      enddo
c      !
c      ! displacement increment, based on element diagonal
c      !
c      Le = dsqrt(((coordsC(1,1)-coordsC(1,7))**two) + 
c     +     ((coordsC(2,1)-coordsC(2,7))**two) +
c     +     ((coordsC(3,1)-coordsC(3,7))**two))
c      !
c      do i=1,nNode
c         do j=1,nDim
c            if(dabs(du(i,j)).gt.10.d0*Le) then
c               pnewdt = 0.5
c               return
c            endif
c         enddo
c      enddo



      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Here, check for hourglass stabilization and get
      !  the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         write(80,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif



      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nInt.eq.1) then
         call xint3D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
      elseif(nInt.eq.8) then
         call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
      else
         write(*,*) 'Invalid number of int points, nInt=',nInt
         write(80,*) 'Invalid number of int points, nInt=',nInt
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt

         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            write(80,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Obtain the neuron cell density and its derivative's at 
         !  this intPt at the begining and end of the incrment
         !
         c_tau = zero
         c_t = zero
         dcdX = zero
         dcdt = zero 
         do k=1,nNode
            c_tau = c_tau + cNew(k)*sh(k)
            c_t   = c_t + cOld(k)*sh(k)
            do i=1,nDim
               dcdX(i,1) = dcdX(i,1) + cNew(k)*dshC(k,i)
            enddo
         enddo
         dcdt = (c_tau - c_t)/dtime 

         !cell2         
         c2_tau = zero
         c2_t = zero
         dc2dX = zero
         dc2dt = zero 
         do k=1,nNode
            c2_tau = c2_tau + c2New(k)*sh(k)
            c2_t   = c2_t + c2Old(k)*sh(k)
            do i=1,nDim
               dc2dX(i,1) = dc2dX(i,1) + c2New(k)*dshC(k,i)
            enddo
         enddo
         dc2dt = (c2_tau - c2_t)/dtime

         !cell3       
         c3_tau = zero
         c3_t = zero
         dc3dX = zero
         dc3dt = zero 
         do k=1,nNode
            c3_tau = c3_tau + c3New(k)*sh(k)
            c3_t   = c3_t + c3Old(k)*sh(k)
            do i=1,nDim
               dc3dX(i,1) = dc3dX(i,1) + c3New(k)*dshC(k,i)
            enddo
         enddo
         dc3dt = (c3_tau - c3_t)/dtime






         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 8 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)



         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo  

         ! obtain coordinates at this integration point
         ! 
         coordintC = matmul(coordsC,Nvec) ! current
         coordint  = matmul(coords,Nvec) ! old 


c         write(*,*),'check before integ'

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         ! 
         ! Perform the constitutive time integration at this integ. point 
         ! 
         call integ(props,nprops,dtime,TIME,F_tau,F_t,
     +               c_tau,c2_tau,c3_tau,c_total,
     +               dcdX,dc2dX,dc3dX,coordintC,coordint,
     +               qflux,qflux2,qflux3,
     +               source,source2,source3,
     +               velocity,velocity2,velocity3,
     +               Diffusivity,Diffusivity2,Diffusivity3,
     +               dqdc,dqdc2,dqdc3,
     +               T_tau,dJdt,
     +               SpTanMod,SpUCMod,
     +               SpCUMod,SpCUMod2,SpCUMod3,
     +               mu,k_normal,k_para,Je_tau,detFg,EffStretch,
     +               EffStretche,TrT) 
         ! 
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



c         ! Save the state variables at this integ point
c         !  at the end of the increment
c         !
c         svars(1+jj) =  mono_tau
         jj = jj + nlSdv ! setup for the next intPt
         

         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         ! 
         globalSdv(jelem,intPt,1) = source ! source term profile
         globalSdv(jelem,intPt,2) = velocity ! velocity profile 
         globalSdv(jelem,intPt,3) = mu ! stiffness profile          
         globalSdv(jelem,intPt,4) = k_normal ! growth parameter          
         globalSdv(jelem,intPt,5) = k_para ! growth parameter             
         globalSdv(jelem,intPt,6) = Je_tau ! elastic Jacobian             
         globalSdv(jelem,intPt,7) = detFg ! growth Jacobian             
         globalSdv(jelem,intPt,8) = EffStretch !  Effective total Stretch   
         globalSdv(jelem,intPt,9) = EffStretche !  Effective elastic Stretch            
         globalSdv(jelem,intPt,10) = TrT !  trace of Cauchy stress                           
         globalSdv(jelem,intPt,11) = c_total !  total cell density   


         ! Compute/update the displacement residual vector
         ! 
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         ! 
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,3+nDim*(kk-1)) = dshC(kk,3)
            Bmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(4,2+nDim*(kk-1)) = dshC(kk,1)
            Bmat(5,2+nDim*(kk-1)) = dshC(kk,3)
            Bmat(5,3+nDim*(kk-1)) = dshC(kk,2)
            Bmat(6,1+nDim*(kk-1)) = dshC(kk,3)
            Bmat(6,3+nDim*(kk-1)) = dshC(kk,1)
         enddo
         ! 
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
            BodyForceRes(3+nDim*(kk-1),1) = sh(kk)*body(3)
         enddo
         ! 
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )


         test1 = test1 + detmapJC*w(intpt)*
     +        (
     +       +transpose(Nvec)*(dcdt + c_tau*dJdt/detF)
     +       +matmul(dshC,qflux)
     +       -transpose(Nvec)*Source
     +        )  


         test2 = test2 + detmapJC*w(intpt)*
     +        (
     +       +transpose(Nvec)*(dc2dt + c2_tau*dJdt/detF)
     +       +matmul(dshC,qflux2)
     +       -transpose(Nvec)*Source2
     +        )


         test3 = test3 + detmapJC*w(intpt)*
     +        (
     +       +transpose(Nvec)*(dc3dt + c3_tau*dJdt/detF)
     +       +matmul(dshC,qflux3)
     +       -transpose(Nvec)*Source3
     +        )


         ! Compute/update the displacement tangent matrix
         ! 
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,3+nDim*(kk-1)) = dshC(kk,1)
            Gmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(5,2+nDim*(kk-1)) = dshC(kk,2)
            Gmat(6,3+nDim*(kk-1)) = dshC(kk,2)
            Gmat(7,1+nDim*(kk-1)) = dshC(kk,3)
            Gmat(8,2+nDim*(kk-1)) = dshC(kk,3)
            Gmat(9,3+nDim*(kk-1)) = dshC(kk,3)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,3+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(4,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(5,2+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(6,3+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(7,1+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(8,2+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(9,3+nDim*(kk-1)) = dshC0(kk,3)
         enddo

         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,3,1)
         Amat(1,4) = SpTanMod(1,1,1,2)
         Amat(1,5) = SpTanMod(1,1,2,2)
         Amat(1,6) = SpTanMod(1,1,3,2)
         Amat(1,7) = SpTanMod(1,1,1,3)
         Amat(1,8) = SpTanMod(1,1,2,3)
         Amat(1,9) = SpTanMod(1,1,3,3)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,3,1)
         Amat(2,4) = SpTanMod(2,1,1,2)
         Amat(2,5) = SpTanMod(2,1,2,2)
         Amat(2,6) = SpTanMod(2,1,3,2)
         Amat(2,7) = SpTanMod(2,1,1,3)
         Amat(2,8) = SpTanMod(2,1,2,3)
         Amat(2,9) = SpTanMod(2,1,3,3)
         Amat(3,1) = SpTanMod(3,1,1,1)
         Amat(3,2) = SpTanMod(3,1,2,1)
         Amat(3,3) = SpTanMod(3,1,3,1)
         Amat(3,4) = SpTanMod(3,1,1,2)
         Amat(3,5) = SpTanMod(3,1,2,2)
         Amat(3,6) = SpTanMod(3,1,3,2)
         Amat(3,7) = SpTanMod(3,1,1,3)
         Amat(3,8) = SpTanMod(3,1,2,3)
         Amat(3,9) = SpTanMod(3,1,3,3)
         Amat(4,1) = SpTanMod(1,2,1,1)
         Amat(4,2) = SpTanMod(1,2,2,1)
         Amat(4,3) = SpTanMod(1,2,3,1)
         Amat(4,4) = SpTanMod(1,2,1,2)
         Amat(4,5) = SpTanMod(1,2,2,2)
         Amat(4,6) = SpTanMod(1,2,3,2)
         Amat(4,7) = SpTanMod(1,2,1,3)
         Amat(4,8) = SpTanMod(1,2,2,3)
         Amat(4,9) = SpTanMod(1,2,3,3)
         Amat(5,1) = SpTanMod(2,2,1,1)
         Amat(5,2) = SpTanMod(2,2,2,1)
         Amat(5,3) = SpTanMod(2,2,3,1)
         Amat(5,4) = SpTanMod(2,2,1,2)
         Amat(5,5) = SpTanMod(2,2,2,2)
         Amat(5,6) = SpTanMod(2,2,3,2)
         Amat(5,7) = SpTanMod(2,2,1,3)
         Amat(5,8) = SpTanMod(2,2,2,3)
         Amat(5,9) = SpTanMod(2,2,3,3)
         Amat(6,1) = SpTanMod(3,2,1,1)
         Amat(6,2) = SpTanMod(3,2,2,1)
         Amat(6,3) = SpTanMod(3,2,3,1)
         Amat(6,4) = SpTanMod(3,2,1,2)
         Amat(6,5) = SpTanMod(3,2,2,2)
         Amat(6,6) = SpTanMod(3,2,3,2)
         Amat(6,7) = SpTanMod(3,2,1,3)
         Amat(6,8) = SpTanMod(3,2,2,3)
         Amat(6,9) = SpTanMod(3,2,3,3)
         Amat(7,1) = SpTanMod(1,3,1,1)
         Amat(7,2) = SpTanMod(1,3,2,1)
         Amat(7,3) = SpTanMod(1,3,3,1)
         Amat(7,4) = SpTanMod(1,3,1,2)
         Amat(7,5) = SpTanMod(1,3,2,2)
         Amat(7,6) = SpTanMod(1,3,3,2)
         Amat(7,7) = SpTanMod(1,3,1,3)
         Amat(7,8) = SpTanMod(1,3,2,3)
         Amat(7,9) = SpTanMod(1,3,3,3)
         Amat(8,1) = SpTanMod(2,3,1,1)
         Amat(8,2) = SpTanMod(2,3,2,1)
         Amat(8,3) = SpTanMod(2,3,3,1)
         Amat(8,4) = SpTanMod(2,3,1,2)
         Amat(8,5) = SpTanMod(2,3,2,2)
         Amat(8,6) = SpTanMod(2,3,3,2)
         Amat(8,7) = SpTanMod(2,3,1,3)
         Amat(8,8) = SpTanMod(2,3,2,3)
         Amat(8,9) = SpTanMod(2,3,3,3)
         Amat(9,1) = SpTanMod(3,3,1,1)
         Amat(9,2) = SpTanMod(3,3,2,1)
         Amat(9,3) = SpTanMod(3,3,3,1)
         Amat(9,4) = SpTanMod(3,3,1,2)
         Amat(9,5) = SpTanMod(3,3,2,2)
         Amat(9,6) = SpTanMod(3,3,3,2)
         Amat(9,7) = SpTanMod(3,3,1,3)
         Amat(9,8) = SpTanMod(3,3,2,3)
         Amat(9,9) = SpTanMod(3,3,3,3)


         Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         

         if((nNode.eq.8).and.(nInt.eq.8)) then
            ! 
            ! This is the tangent using the F-bar method with the
            !  8 node fully integrated linear element
            ! 
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +           )
         else
            ! 
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
            ! 
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
         endif



         ! Compute/update the displacement - cell density tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         ! I need: SpUCMod
         ! 
         AmatUC = zero
         AmatUC(1,1) = SpUCMod(1,1)
         AmatUC(2,1) = SpUCMod(2,2)
         AmatUC(3,1) = SpUCMod(3,3)
         AmatUC(4,1) = SpUCMod(1,2)
         AmatUC(5,1) = SpUCMod(2,3)
         AmatUC(6,1) = SpUCMod(1,3)
         ! 
         Kuc = Kuc + detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(transpose(Bmat),AmatUC),Nvec)
     +        )



         ! Compute/update the cell density  - displacement tangent matrix
         ! The F-bar method will have some effect, however we neglect that here.
         ! I need SpCUMod
         ! 
         AmatCU = zero
         AmatCU(1,1) = SpCUMod(1,1,1)
         AmatCU(1,2) = SpCUMod(1,2,1)
         AmatCU(1,3) = SpCUMod(1,3,1)
         AmatCU(1,4) = SpCUMod(1,1,2)
         AmatCU(1,5) = SpCUMod(1,2,2)
         AmatCU(1,6) = SpCUMod(1,3,2)
         AmatCU(1,7) = SpCUMod(1,1,3)
         AmatCU(1,8) = SpCUMod(1,2,3)
         AmatCU(1,9) = SpCUMod(1,3,3)
         AmatCU(2,1) = SpCUMod(2,1,1)
         AmatCU(2,2) = SpCUMod(2,2,1)
         AmatCU(2,3) = SpCUMod(2,3,1)
         AmatCU(2,4) = SpCUMod(2,1,2)
         AmatCU(2,5) = SpCUMod(2,2,2)
         AmatCU(2,6) = SpCUMod(2,3,2)
         AmatCU(2,7) = SpCUMod(2,1,3)
         AmatCU(2,8) = SpCUMod(2,2,3)
         AmatCU(2,9) = SpCUMod(2,3,3)
         AmatCU(3,1) = SpCUMod(3,1,1)
         AmatCU(3,2) = SpCUMod(3,2,1)
         AmatCU(3,3) = SpCUMod(3,3,1)
         AmatCU(3,4) = SpCUMod(3,1,2)
         AmatCU(3,5) = SpCUMod(3,2,2)
         AmatCU(3,6) = SpCUMod(3,3,2)
         AmatCU(3,7) = SpCUMod(3,1,3)
         AmatCU(3,8) = SpCUMod(3,2,3)
         AmatCU(3,9) = SpCUMod(3,3,3)
         ! 
         Kcu = Kcu - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatCU),Gmat)
     +        )

         AmatCU2 = zero
         AmatCU2(1,1) = SpCUMod2(1,1,1)
         AmatCU2(1,2) = SpCUMod2(1,2,1)
         AmatCU2(1,3) = SpCUMod2(1,3,1)
         AmatCU2(1,4) = SpCUMod2(1,1,2)
         AmatCU2(1,5) = SpCUMod2(1,2,2)
         AmatCU2(1,6) = SpCUMod2(1,3,2)
         AmatCU2(1,7) = SpCUMod2(1,1,3)
         AmatCU2(1,8) = SpCUMod2(1,2,3)
         AmatCU2(1,9) = SpCUMod2(1,3,3)
         AmatCU2(2,1) = SpCUMod2(2,1,1)
         AmatCU2(2,2) = SpCUMod2(2,2,1)
         AmatCU2(2,3) = SpCUMod2(2,3,1)
         AmatCU2(2,4) = SpCUMod2(2,1,2)
         AmatCU2(2,5) = SpCUMod2(2,2,2)
         AmatCU2(2,6) = SpCUMod2(2,3,2)
         AmatCU2(2,7) = SpCUMod2(2,1,3)
         AmatCU2(2,8) = SpCUMod2(2,2,3)
         AmatCU2(2,9) = SpCUMod2(2,3,3)
         AmatCU2(3,1) = SpCUMod2(3,1,1)
         AmatCU2(3,2) = SpCUMod2(3,2,1)
         AmatCU2(3,3) = SpCUMod2(3,3,1)
         AmatCU2(3,4) = SpCUMod2(3,1,2)
         AmatCU2(3,5) = SpCUMod2(3,2,2)
         AmatCU2(3,6) = SpCUMod2(3,3,2)
         AmatCU2(3,7) = SpCUMod2(3,1,3)
         AmatCU2(3,8) = SpCUMod2(3,2,3)
         AmatCU2(3,9) = SpCUMod2(3,3,3)
         ! 
         Kc2u = Kc2u - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatCU2),Gmat)
     +        )


         AmatCU3 = zero
         AmatCU3(1,1) = SpCUMod3(1,1,1)
         AmatCU3(1,2) = SpCUMod3(1,2,1)
         AmatCU3(1,3) = SpCUMod3(1,3,1)
         AmatCU3(1,4) = SpCUMod3(1,1,2)
         AmatCU3(1,5) = SpCUMod3(1,2,2)
         AmatCU3(1,6) = SpCUMod3(1,3,2)
         AmatCU3(1,7) = SpCUMod3(1,1,3)
         AmatCU3(1,8) = SpCUMod3(1,2,3)
         AmatCU3(1,9) = SpCUMod3(1,3,3)
         AmatCU3(2,1) = SpCUMod3(2,1,1)
         AmatCU3(2,2) = SpCUMod3(2,2,1)
         AmatCU3(2,3) = SpCUMod3(2,3,1)
         AmatCU3(2,4) = SpCUMod3(2,1,2)
         AmatCU3(2,5) = SpCUMod3(2,2,2)
         AmatCU3(2,6) = SpCUMod3(2,3,2)
         AmatCU3(2,7) = SpCUMod3(2,1,3)
         AmatCU3(2,8) = SpCUMod3(2,2,3)
         AmatCU3(2,9) = SpCUMod3(2,3,3)
         AmatCU3(3,1) = SpCUMod3(3,1,1)
         AmatCU3(3,2) = SpCUMod3(3,2,1)
         AmatCU3(3,3) = SpCUMod3(3,3,1)
         AmatCU3(3,4) = SpCUMod3(3,1,2)
         AmatCU3(3,5) = SpCUMod3(3,2,2)
         AmatCU3(3,6) = SpCUMod3(3,3,2)
         AmatCU3(3,7) = SpCUMod3(3,1,3)
         AmatCU3(3,8) = SpCUMod3(3,2,3)
         AmatCU3(3,9) = SpCUMod3(3,3,3)
         ! 
         Kc3u = Kc3u - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatCU3),Gmat)
     +        )



         ! Compute/update the neuron cell density tangent matrix
         ! I need Diffusivity,dqdc
         ! I need dJdt

         Kcc = Kcc + detmapJC*w(intPt)*
     +        (
     +        -(1.0/dtime + dJdt/detF)*matmul(transpose(Nvec),Nvec)
     +        -Diffusivity*matmul(dshC,transpose(dshC))
     +        -matmul(matmul(dshC,dqdc),Nvec)
     +        )

         Kc2c2 = Kc2c2 + detmapJC*w(intPt)*
     +        (
     +        -(1.0/dtime + dJdt/detF)*matmul(transpose(Nvec),Nvec)
     +        -Diffusivity2*matmul(dshC,transpose(dshC))
     +        -matmul(matmul(dshC,dqdc2),Nvec)
     +        )

         Kc3c3 = Kc3c3 + detmapJC*w(intPt)*
     +        (
     +        -(1.0/dtime + dJdt/detF)*matmul(transpose(Nvec),Nvec)
     +        -Diffusivity3*matmul(dshC,transpose(dshC))
     +        -matmul(matmul(dshC,dqdc3),Nvec)
     +        )




      enddo
      ! 
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface flux terms
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the flux
            !  acts on is the flux ``label''
            !
            face = jdltyp(i,1)
            flux = adlmag(i,1)

            
            if((face.ge.1).and.(face.le.6)) then
               !
               ! fluid flux applied
               !
               if(face.eq.1) then
                  faceFlag = 1
               elseif(face.eq.2) then
                  faceFlag = 2
               elseif(face.eq.3) then
                  faceFlag = 3
               elseif(face.eq.4) then
                  faceFlag = 4
               elseif(face.eq.5) then
                  faceFlag = 5
               else
                  faceFlag = 6
               endif
               !
               if(nIntS.eq.1) then
                  call xintSurf3D1pt(faceFlag,xLocal,yLocal,zLocal,wS)
               elseif(nIntS.eq.4) then
                  call xintSurf3D4pt(faceFlag,xLocal,yLocal,zLocal,wS)
               else
                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
                  write(80,*) 'Invalid nIntS points, nIntS=',nIntS
                  call xit
               endif
               !
               ! loop over integ points on this element face
               !
               do ii=1,nIntS
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (dA)
                  !
                  call computeSurf3D(xLocal(ii),yLocal(ii),zLocal(ii),
     +                 faceFlag,coordsC,sh,dA)
                  !
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  ! Obtain the light intensity at this integ point
                  !
                  I_tau = zero
                  do n=1,nNode
                     I_tau = I_tau + INew(n)*sh(n)
                  enddo


                  ! Check for light hitting this surface, if test is
                  !  less than zero, then light is hitting this surface
                  !
                  test = direction(1,1)*normal(1,1) 
     +                 + direction(2,1)*normal(2,1)
     +                 + direction(3,1)*normal(3,1)
     
                  if(test.lt.zero) then
c                     test = (1 - tanh(test/.2)) * 1.5
                      test = 1.0
                  else
                     test = zero
                  endif



                  !
                  ! Modify the residual, loop over nodes, recall
                  !  sh(n)=0 when n is not on this face
                  !
                  do n=1,nNode
                     RI(n,1) = RI(n,1) 
     +                    + wS(ii)*ds*sh(n)*(I_tau - flux)*test
                  enddo


                  ! Modify the tangent
                  !
                  do n=1,nNode
                     do m=1,nNode
                        KII(n,m) = KII(n,m) - w(ii)*ds*sh(n)*sh(m)*test
                     enddo
                  enddo




               enddo ! end loop over integ points
               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over surface flux terms
      !----------------------------------------------------------------



      ! I turn off the deformation related residual and tangents 
      ! 
c      Ru  = zero
c      Rc = zero
c      Rc2 = zero ! update
c      Rc3 = zero ! update
      

c      Kuu = zero
c      Kuc = zero
c      Kuc2 = zero ! update 
c      Kuc3 = zero ! update 



c      Kcu = zero
c      Kcc = zero
c      Kcc2 = zero ! update 
c      Kcc3 = zero

c      Kc2u = zero ! update 
c      Kc2c = zero ! update
c      Kc2c2 = zero ! update 
c      Kc2c3 = zero ! update 

c      Kc3u = zero ! update 
c      Kc3c = zero ! update 
c      Kc3c2 = zero ! update 
c      Kc3c3 = zero ! update 



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      ! 
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,test1,test2,test3,
     +     Kuu,Kuc,Kuc,Kuc,
     +     Kcu,Kcc,Kcc2,kcc3,
     +     Kc2u,Kc2c,Kc2c2,kc2c3,
     +     Kc3u,Kc3c,Kc3c2,Kc3c3,
     +     rhs,amatrx)     
      ! 
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      
      return
      end subroutine U3D8

************************************************************************

      subroutine computeSUPG(nDim,nNode,coordsC,direction,epsilon,SUPG)

      ! This subroutine computes the SUPG stabilization parameter often
      !  called ``tau'' in the literature.  Since we have coth(a)-1/a as
      !  part of this calculation, we use the taylor expansion near zero
      !  to avoid the divide by zero.

      implicit none
      
      integer nDim,nNode

      real*8 hXi,hEta,hZeta,eXi(nDim,1),eEta(nDim,1),eZeta(nDim,1),dXi
      real*8 dEta,dZeta,alphaXi,alphaEta,alphaZeta,xiBar,etaBar,zetaBar
      real*8 coordsC(nDim,nNode),direction(nDim,1),epsilon,tmpX,tmpY
      real*8 tmpZ,SUPG

      if(nDim.eq.2) then
         !
         ! This is a 2D problem
         !
         hXi = 0.5d0*dsqrt(
     +       (coordsC(1,3)+coordsC(1,2)-coordsC(1,1)-coordsC(1,4))**2.d0
     +      +(coordsC(2,3)+coordsC(2,2)-coordsC(2,1)-coordsC(2,4))**2.d0
     +        )
         hEta = 0.5d0*dsqrt(
     +       (coordsC(1,3)+coordsC(1,4)-coordsC(1,1)-coordsC(1,2))**2.d0
     +      +(coordsC(2,3)+coordsC(2,4)-coordsC(2,1)-coordsC(2,2))**2.d0
     +        )
      
         eXi(1,1) = 
     +        (coordsC(1,3)+coordsC(1,2)-coordsC(1,1)-coordsC(1,4))
     +        /(2.d0*hXi)
         eXi(2,1) = 
     +        (coordsC(2,3)+coordsC(2,2)-coordsC(2,1)-coordsC(2,4))
     +        /(2.d0*hXi)
         eEta(1,1) = 
     +        (coordsC(1,3)+coordsC(1,4)-coordsC(1,1)-coordsC(1,2))
     +        /(2.d0*hEta)
         eEta(2,1) = 
     +        (coordsC(2,3)+coordsC(2,4)-coordsC(2,1)-coordsC(2,2))
     +        /(2.d0*hEta)
         

         dXi = direction(1,1)*eXi(1,1) + direction(2,1)*eXi(2,1)
         dEta = direction(1,1)*eEta(1,1) + direction(2,1)*eEta(2,1)
      

         if(epsilon.lt.1.d-10) then
            !
            ! we have zero diffusion in the PDE, but we assume
            !  a very small value only to compute the SUPG parameter
            !
            epsilon = 1.d-9
            !
         endif


         alphaXi = (dXi*hXi)/(2.d0*epsilon)
         alphaEta = (dEta*hEta)/(2.d0*epsilon)
         !
         if(abs(alphaXi).le.0.1d0) then
            xiBar = (1.d0/3.d0)*alphaXi 
     +           - (1.d0/45.d0)*alphaXi**3.d0 
     +           + (2.d0/945.d0)*alphaXi**5.d0
     +           - (1.d0/4725.d0)*alphaXi**7.d0
         else
            xiBar = 1.d0/dtanh(alphaXi) - 1.d0/alphaXi   
         endif
         !
         if(abs(alphaEta).le.0.1d0) then
            etaBar = (1.d0/3.d0)*alphaEta 
     +           - (1.d0/45.d0)*alphaEta**3.d0 
     +           + (2.d0/945.d0)*alphaEta**5.d0
     +           - (1.d0/4725.d0)*alphaEta**7.d0
         else
            etaBar = 1.d0/dtanh(alphaEta) - 1.d0/alphaEta
         endif
         !
         SUPG = (xiBar*dXi*hXi + etaBar*dEta*hEta)/2.d0
         !
      elseif(nDim.eq.3) then
         !
         ! This is a 3D problem
         !
         tmpX = coordsC(1,2)+coordsC(1,3)+coordsC(1,6)+coordsC(1,7)
     +        - coordsC(1,1)-coordsC(1,4)-coordsC(1,5)-coordsC(1,8)
         tmpY = coordsC(2,2)+coordsC(2,3)+coordsC(2,6)+coordsC(2,7)
     +        - coordsC(2,1)-coordsC(2,4)-coordsC(2,5)-coordsC(2,8)
         tmpZ = coordsC(3,2)+coordsC(3,3)+coordsC(3,6)+coordsC(3,7)
     +        - coordsC(3,1)-coordsC(3,4)-coordsC(3,5)-coordsC(3,8)
         hXi = 0.25d0*dsqrt(tmpX**2.d0 + tmpY**2.d0 + tmpZ**2.d0)
         eXi(1,1) = 0.25d0*tmpX/hXi
         eXi(2,1) = 0.25d0*tmpY/hXi
         eXi(3,1) = 0.25d0*tmpZ/hXi


         tmpX = coordsC(1,3)+coordsC(1,4)+coordsC(1,7)+coordsC(1,8)
     +        - coordsC(1,1)-coordsC(1,2)-coordsC(1,5)-coordsC(1,6)
         tmpY = coordsC(2,3)+coordsC(2,4)+coordsC(2,7)+coordsC(2,8)
     +        - coordsC(2,1)-coordsC(2,2)-coordsC(2,5)-coordsC(2,6)
         tmpZ = coordsC(3,3)+coordsC(3,4)+coordsC(3,7)+coordsC(3,8)
     +        - coordsC(3,1)-coordsC(3,2)-coordsC(3,5)-coordsC(3,6)
         hEta = 0.25d0*dsqrt(tmpX**2.d0 + tmpY**2.d0 + tmpZ**2.d0)
         eEta(1,1) = 0.25d0*tmpX/hEta
         eEta(2,1) = 0.25d0*tmpY/hEta
         eEta(3,1) = 0.25d0*tmpZ/hEta

         tmpX = coordsC(1,5)+coordsC(1,6)+coordsC(1,7)+coordsC(1,8)
     +        - coordsC(1,1)-coordsC(1,2)-coordsC(1,3)-coordsC(1,4)
         tmpY = coordsC(2,5)+coordsC(2,6)+coordsC(2,7)+coordsC(2,8)
     +        - coordsC(2,1)-coordsC(2,2)-coordsC(2,3)-coordsC(2,4)
         tmpZ = coordsC(3,5)+coordsC(3,6)+coordsC(3,7)+coordsC(3,8)
     +        - coordsC(3,1)-coordsC(3,2)-coordsC(3,3)-coordsC(3,4)
         hZeta = 0.25d0*dsqrt(tmpX**2.d0 + tmpY**2.d0 + tmpZ**2.d0)
         eZeta(1,1) = 0.25d0*tmpX/hZeta
         eZeta(2,1) = 0.25d0*tmpY/hZeta
         eZeta(3,1) = 0.25d0*tmpZ/hZeta


         dXi = direction(1,1)*eXi(1,1) + direction(2,1)*eXi(2,1) 
     +        + direction(3,1)*eXi(3,1) 
         dEta = direction(1,1)*eEta(1,1) + direction(2,1)*eEta(2,1)
     +        + direction(3,1)*eEta(3,1)
         dZeta = direction(1,1)*eZeta(1,1) + direction(2,1)*eZeta(2,1)
     +        + direction(3,1)*eZeta(3,1)



         if(epsilon.lt.1.d-10) then
            !
            ! we have zero diffusion in the PDE, but we assume
            !  a very small value only to compute the SUPG parameter
            !
            epsilon = 1.d-9
            !
         endif

         alphaXi = (dXi*hXi)/(2.d0*epsilon)
         alphaEta = (dEta*hEta)/(2.d0*epsilon)
         alphaZeta = (dZeta*hZeta)/(2.d0*epsilon)
            
         if(abs(alphaXi).le.0.1d0) then
            xiBar = (1.d0/3.d0)*alphaXi 
     +           - (1.d0/45.d0)*alphaXi**3.d0 
     +           + (2.d0/945.d0)*alphaXi**5.d0
     +           - (1.d0/4725.d0)*alphaXi**7.d0
         else
            xiBar = 1.d0/dtanh(alphaXi) - 1.d0/alphaXi   
         endif
         
         if(abs(alphaEta).le.0.1d0) then
            etaBar = (1.d0/3.d0)*alphaEta 
     +           - (1.d0/45.d0)*alphaEta**3.d0 
     +           + (2.d0/945.d0)*alphaEta**5.d0
     +           - (1.d0/4725.d0)*alphaEta**7.d0
         else
            etaBar = 1.d0/dtanh(alphaEta) - 1.d0/alphaEta
         endif
         
         if(abs(alphaZeta).le.0.1d0) then
            zetaBar = (1.d0/3.d0)*alphaZeta 
     +           - (1.d0/45.d0)*alphaZeta**3.d0 
     +           + (2.d0/945.d0)*alphaZeta**5.d0
     +           - (1.d0/4725.d0)*alphaZeta**7.d0
         else
            zetaBar = 1.d0/dtanh(alphaZeta) - 1.d0/alphaZeta
         endif
         !
         SUPG = (xiBar*dXi*hXi + etaBar*dEta*hEta + zetaBar*dZeta*hZeta)
     +        /2.d0
         !
      else
         write(*,*) 'SUPG, bad nDim=',nDim
         write(80,*) 'SUPG: bad nDim=',nDim
         call xit
      endif

      return
      end subroutine computeSUPG
        
************************************************************************ 
      subroutine integ(props,nprops,dtime,TIME,F_tau,F_t,
     +               c_tau,c2_tau,c3_tau,c_total,
     +               dcdX,dc2dX,dc3dX,coordintC,coordint,
     +               qflux,qflux2,qflux3,
     +               source,source2,source3,
     +               velocity,velocity2,velocity3,
     +               Diffusivity,Diffusivity2,Diffusivity3,
     +               dqdc,dqdc2,dqdc3,
     +               T_tau,dJdt,
     +               SpTanMod,SpUCMod,
     +               SpCUMod,SpCUMod2,SpCUMod3,
     +               mu,k_normal,k_para,Je_tau,detFg,EffStretch,
     +               EffStretche,TrT) 
      ! This subroutine computes everything required for the time integration
      ! of the problem.


      implicit none

      integer i,j,k,l,m,n,nprops,stat,kstep,a,b,c,d

      logical crosslink,cleave

      real*8 Iden(3,3),props(nprops),F_tau(3,3),I_tau,dtime,theta,tmpA
      real*8 T_tau(3,3),dTRdF(3,3,3,3),Gshear1,Kbulk1,Gshear2,kRate,tmpB
      real*8 spTanMod(3,3,3,3),detF,FinvT(3,3),B_tau(3,3),trB_tau,deltaI
      real*8 Finv(3,3),SpUIMod(3,3),SpIUmod(3,3),Bdis(3,3)
      real*8 Bdis0(3,3),trBdis,sigma,direction(3,1),dx,dy,dz,Avogadro
      real*8 Plank,alpha0,C0,sigmaMatrix,sigmaUB,sigmaB,alpha_tau,phiB
      real*8 freqCrosslink,freqCleave,factor_UB,factor_B,alpha_t,phiUB
      real*8 F_starInv(3,3),F_tilde(3,3),detF_tilde,T1_tau(3,3),epsilon
      real*8 B_tilde(3,3),trB_tilde,Bdis_tilde(3,3),trBdis_tilde,tmpC
      real*8 Bdis0_tilde(3,3),T2_tau(3,3),F_tildeInv(3,3),epsilonUB
      real*8 F_tildeInvT(3,3),dTRdF_A(3,3,3,3),dTRdF_B(3,3,3,3),DalphaDI
      real*8 spTanMod_A(3,3,3,3),spTanMod_B(3,3,3,3),factorA,factorB
      real*8 DsigmaDI,lamBar
      real*8 diff,eps_r,sigma_c,k_t,k_p,qyield,c_tau,fac
      real*8 DcDI,DcdotDI,DIDc,mono_t,mono_tau,mono_t0
      real*8 Gshear,Gshear0,Gshearc,xi
      real*8 I_tau_here
      real*8 dcdX(3,1),coordint(3,1),coordintC(3,1)
      real*8 qflux(3,1),source,Diffusivity
      real*8 alpha,alpha_Gc,alpha_v,Gc,Gc2,vel
      real*8 tmp,nx,ny,nz,cdiff,xcoord,ycoord,zcoord,coordiff
      real*8 Hhatc,Hhats,Hhatv,dqdc(3,1),fac1,fac2,fac3,dHdc
      real*8 velocity,radius,angletheta
      real*8 F_t(3,3),dJdt,SpUCMod(3,3),SpCUMod(3,3,3)
      real*8 direction_cur(3,1),theta_normal,k_normal,alpha_normal
      real*8 theta_para,k_para,alpha_para,Fg_tau(3,3),Fginv(3,3)
      real*8 detFg,FginvT(3,3),Fe_tau(3,3),Je_tau,lambda,mu
      real*8 J_tau,J_t,Kirk_tau(3,3),Ce_mat(3,3,3,3),C_sp(3,3,3,3)
      real*8 C_mat(3,3,3,3),dqdFFT(3,3,3),qIden(3,3,3),dFgdc(3,3)
      real*8 fac4,fac5,dTdc(3,3),term1(3,3),term2(3,3),term3(3,3)
      real*8 alpha_d,Hhatd,alpha_mu,lambda_c,lambda_s,mu_c,mu_s,Hhatmu
      real*8 alpha_k,k_s,k_c_normal,k_c_para,Hhatk
      real*8 AR_mat1(3,3,3,3),AR_mat2(3,3,3,3),AR_mat3(3,3,3,3)
      real*8 AR_mat4(3,3,3,3)
      real*8 AR_mat(3,3,3,3),A_mat(3,3,3,3)
      real*8 dTdFg1(3,3,3,3),dTdFg2(3,3,3,3),dTdFg(3,3,3,3),Be_tau(3,3)
      real*8 Cten_tau(3,3),TrC,EffStretch,Ceten_tau(3,3),TrCe,EffStretche
      real*8 TrT
      real*8 c_t
      real*8 epsilon1,Gc_time1,TIME
      real*8 cdiff2,c02,Hhatc2,c2_tau,alpha2,coordiff2,alpha_v2,Hhatv2
      real*8 alpha_d2,Hhatd2
      real*8 alpha_Gc2,Hhats2,velocity2,vel2,Diffusivity2,diff2
      real*8 qflux2(3,1),dc2dX(3,1),Gc_time2,epsilon2,source2
      real*8 fac12,fac22,fac32,dHdc2,dqdc2(3,1),c_total
      real*8 dqdFFT2(3,3,3),qIden2(3,3,3)
      real*8 SpCUMod2(3,3,3)
      real*8 c3_tau,dc3dX(3,1),qflux3(3,1),source3,velocity3
      real*8 Diffusivity3,dqdc3(3,1),SpCUMod3(3,3,3)
      real*8 alpha3,alpha_Gc3,alpha_v3,alpha_d3,Gc3,vel3,diff3,c03
      real*8 epsilon3,cdiff3,Hhatc3,coordiff3,Hhatv3
      real*8 Hhats3,Gc_time3,fac13,fac23,fac33,dHdc3
      real*8 dqdFFT3(3,3,3),qIden3(3,3,3),Hhatd3
      real*8 radius_3d,theta_3d,phi
      real*8 delta_v1,delta_v2,delta_v3
 

      real*8 zero,one,two,three,third,half,nine,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0,nine=9.d0,four=4.d0)


c      write(*,*),'F_tau = ',F_tau



      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain material properties
      ! 
      alpha         = props(1) ! smooth parameter for neuron cell migration  
      alpha2        = props(2) ! 
      alpha3        = props(3)



      alpha_Gc      = props(4) ! smooth parameter for source term 
      alpha_Gc2     = props(5) ! 
      alpha_Gc3     = props(6)


      alpha_v       = props(7) ! smooth parameter for velocity
      alpha_v2      = props(8)  ! 
      alpha_v3      = props(9)  ! 

      alpha_d       = props(10) ! smooth parameter for diffusicity 
      alpha_d2      = props(11)
      alpha_d3      = props(12)


      alpha_mu      = props(13) ! smooth parameter for stiffness 



      alpha_k       = props(14) ! smooth parameter for growth parameter

      Gc            = props(15) ! base value for source term
      Gc2           = props(16) !  
      Gc3           = props(17) !  



      vel           = props(18) ! base value for velocity
      vel2          = props(19) ! 
      vel3          = props(20) ! 


      diff          = props(21) ! diffusivity
      diff2         = props(22) ! 
      diff3         = props(23) ! 


      c0            = props(24) ! cell migration threshold 
      c02           = props(25) ! 
      c03           = props(26) ! 


      k_s           = props(27) ! shape parameter for growth 
      k_c_normal    = props(28) ! shape parameter for growth
      k_c_para      = props(29) ! shape parameter for growth      
      alpha_normal  = props(30) ! shape parameter for growth
      alpha_para    = props(31) ! shape parameter for growth
      lambda_c      = props(32) ! Lame constant 
      lambda_s      = props(33) ! Lame constant 
      mu_c          = props(34) !Lame constant
      mu_s          = props(35) !Lame constant

      epsilon1      = props(36) !delta function parameter for Gc as function of time
      epsilon2      = props(37) ! 
      epsilon3      = props(38) ! 

      delta_v1      = props(39) ! final destinations of cells 
      delta_v2      = props(40)
      delta_v3      = props(41)

      ! extract coordinates at this integration point 
      ! 
      xcoord = coordint(1,1)
      ycoord = coordint(2,1)
      zcoord = coordint(3,1)


      ! calculate the radius based on the polar coordinates 
      ! 
      call CartoPolar(xcoord,ycoord,radius,angletheta)
      call CartoShere(xcoord,ycoord,zcoord,radius_3d,theta_3d,phi)


      ! Ensure the velocity vector is a unit vector
      ! 
      ! 1-D velocity 
       nx = 1.0
       ny = 0.0
       nz = 0.0
      ! 2-D half circle
c      nx = xcoord
c      ny = ycoord
c      nz = 0.0
      ! 3D ball  
c      nx = xcoord
c      ny = ycoord
c      nz = zcoord



      ! initial direction vector 
      direction(1,1) = nx
      direction(2,1) = ny
      direction(3,1) = nz
      ! make sure the vector is a unit vector
      tmp = sqrt(direction(1,1)**2.0 + direction(2,1)**2.0 + direction(3,1)**2.0)
      direction = direction/tmp
      
      ! map the vector to the current one 
      direction_cur = matmul(F_tau,direction)
      
      
      
      ! make sure the vector is a unit vector
      tmp = sqrt(direction_cur(1,1)**2.0 + direction_cur(2,1)**2.0 + direction_cur(3,1)**2.0)
      direction_cur = direction_cur/tmp

      ! construct the qflux 
      ! 
      cdiff  = c_tau - c0
      cdiff2 = c2_tau - c02 ! update
      cdiff3 = c3_tau - c03 ! update

      call Hhat(cdiff,alpha,Hhatc)
      call Hhat(cdiff2,alpha2,Hhatc2) ! update 
      call Hhat(cdiff3,alpha3,Hhatc3) ! update 


    
      ! 1-D 
      coordiff  = -(xcoord - delta_v1)
      coordiff2 = -(xcoord - delta_v2) ! update 
      coordiff3 = -(xcoord - delta_v3) ! update 
      ! 2-D
c      coordiff = -(radius - 0.95)
c      coordiff2 = -(radius - 0.95) ! update 
c      coordiff3 = -(radius - 0.95) ! update 
      ! 3-D ball
c      coordiff  = -(radius_3d - 1.0)
c      coordiff2 = -(radius_3d - 1.0)
c      coordiff3 = -(radius_3d - 1.0)

      call Hhat(coordiff,alpha_v,Hhatv) ! Heaveistep function for velocity
      call Hhat(coordiff2,alpha_v2,Hhatv2) ! update
      call Hhat(coordiff3,alpha_v3,Hhatv3) ! update

      call Hhat(coordiff,alpha_d,Hhatd) ! Heaveistep function for diffusivity
      call Hhat(coordiff2,alpha_d2,Hhatd2) ! update
      call Hhat(coordiff3,alpha_d3,Hhatd3) ! update


      ! 1-D 
      coordiff = (xcoord - 0.98*239.5)
      ! 2-D 
c      coordiff = (radius - 0.98)
      ! 3-D
c      coordiff = (radius_3d - 0.98)
c      call Hhat(coordiff,alpha_mu,Hhatmu) ! Heaveistep function for stiffness
c      call Hhat(coordiff,alpha_k,Hhatk) ! Heaveistep function for growth parameter




      ! 1-D
      coordiff = -(xcoord - 0.2*239.5)
      coordiff2 = -(xcoord - 0.2*239.5)
      coordiff3 = -(xcoord - 0.2*239.5)
      ! 2-D 
c      coordiff = -(radius - 0.2)   
c      coordiff2 = -(radius - 0.2)    
c      coordiff3 = -(radius - 0.2)    
      ! 3-D
c      coordiff  = -(radius_3d - 0.2) 
c      coordiff2 = -(radius_3d - 0.2) 
c      coordiff3 = -(radius_3d - 0.2) 
          
      call Hhat(coordiff,alpha_Gc,Hhats) ! Gc
      call Hhat(coordiff2,alpha_Gc2,Hhats2) ! update 
      call Hhat(coordiff3,alpha_Gc3,Hhats3) ! update 



      ! velocity 
      ! 
      velocity  = vel*Hhatv   - 0.5*vel
      velocity2 = vel2*Hhatv2 - 0.5*vel2 ! update 
      velocity3 = vel3*Hhatv3 - 0.5*vel3 ! update 

      
      ! diffusicity 
      ! 
c      Diffusivity = diff*(1.5 - 0.95*Hhatd)
      Diffusivity = diff 
      Diffusivity2 = diff2 ! update 
      Diffusivity3 = diff3 ! update 



      ! Neuron cell flux 
      ! 
      qflux = -c_tau*Hhatc*velocity*direction_cur + Diffusivity*dcdX
      qflux2 = -c2_tau*Hhatc2*velocity2*direction_cur + Diffusivity2*dc2dX ! update 
      qflux3 = -c3_tau*Hhatc3*velocity3*direction_cur + Diffusivity3*dc3dX ! update 
    


      ! source term 
      ! 


      Gc_time1 = (epsilon1**2.0)/((TIME - 3.0)**2.0 + epsilon1**2.0)
      Gc_time2 = (epsilon2**2.0)/((TIME - 3.0)**2.0 + epsilon2**2.0) ! update 
      Gc_time3 = (epsilon3**2.0)/((TIME - 6.0)**2.0 + epsilon3**2.0) ! update 




      source = Gc_time1*Gc*Hhats
      source2 = Gc_time2*Gc2*Hhats2 ! update 
      source3 = Gc_time3*Gc3*Hhats3 ! update 

c      source  = Gc*Hhats ! dubug 
c      source2 = Gc2*Hhats2 !debug
c      source3 = Gc3*Hhats3 !debug

c      source  = Gc ! dubug 
c      source2 = Gc2 !debug
c      source3 = Gc3 !debug


      ! dqdc term 
      ! 
      fac1 = alpha*exp(alpha*(c_tau-c0))
      fac12 = alpha2*exp(alpha2*(c2_tau-c02)) !update 
      fac13 = alpha3*exp(alpha3*(c3_tau-c03)) !update 


      fac2 = alpha*exp(2.0*alpha*(c_tau-c0))
      fac22 = alpha2*exp(2.0*alpha2*(c2_tau-c02)) ! update 
      fac23 = alpha3*exp(2.0*alpha3*(c3_tau-c03)) ! update 

      fac3 = exp(alpha*(c_tau-c0)) + 1.0
      fac32 = exp(alpha2*(c2_tau-c02)) + 1.0 ! update 
      fac33 = exp(alpha3*(c3_tau-c03)) + 1.0 ! update 

      dHdc = fac1/fac3 - fac2/(fac3**2.0)
      dHdc2 = fac12/fac32 - fac22/(fac32**2.0)! update 
      dHdc3 = fac13/fac33 - fac23/(fac33**2.0)! update 


      dqdc = -(Hhatc + c_tau*dHdc)*vel*Hhatv*direction_cur
      dqdc2 = -(Hhatc2 + c2_tau*dHdc2)*vel2*Hhatv2*direction_cur ! update       
      dqdc3 = -(Hhatc3 + c3_tau*dHdc3)*vel3*Hhatv3*direction_cur ! update       
      
      ! stiffness 
      ! 
      mu = mu_s + (mu_c - mu_s)*Hhatmu
      lambda = lambda_s + (lambda_c - lambda_s)*Hhatmu
      
      
      ! k values 
      ! 
c      k_normal = k_s + (k_c_normal - k_s)*Hhatk
      k_normal = 0.0


c      k_para = k_s + (k_c_para - k_s)*Hhatk
      k_para = k_s
      
      ! c_total 
      c_total = c_tau + c2_tau + c3_tau

      ! growth paraemter
      !  
      theta_normal  = (1.0 + k_normal*c_total)**alpha_normal
      theta_para    = (1.0 + k_para*c_total)**alpha_para
      

c      write(*,*),'theta_normal=', theta_normal
c      write(*,*),'theta_para=', theta_para


      ! growth deformation gradient 
      ! 
      Fg_tau = theta_normal*(Iden - matmul(direction,transpose(direction)))
     +        +theta_para*matmul(direction,transpose(direction))


           
     
     
      ! Compute the inverse of Fg, its determinant, and its transpose
      ! 
      call matInv3D(Fg_tau,Fginv,detFg,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         call xit
      endif
      FginvT = transpose(Fginv)     
      
      
      ! Compute the inverse of F_tau, its determinant, and its transpose
      !     
      call matInv3D(F_tau,Finv,detF,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         call xit
      endif


      ! calculate the elastic deformation gradient and elastic Jacobian 
      ! 
      Fe_tau = matmul(F_tau,Fginv)

      call mdet(Fe_tau,Je_tau)
      
      Be_tau = matmul(Fe_tau,transpose(Fe_tau))


      ! calculate the total effective stretch 
      ! 
      Cten_tau = matmul(transpose(F_tau),F_tau)
      TrC = Cten_tau(1,1)+Cten_tau(2,2)+Cten_tau(3,3)
      EffStretch = sqrt((1.0/3.0)*TrC)
      
      
      ! calculate the elastic effetive stretch 
      ! 
      Ceten_tau = matmul(transpose(Fe_tau),Fe_tau)      
      TrCe = Ceten_tau(1,1)+Ceten_tau(2,2)+Ceten_tau(3,3)
      EffStretche = sqrt((1.0/3.0)*TrCe)      


      ! calculate the Cauchy stress 
      ! 
      T_tau = (1.0/Je_tau)*((lambda*dlog(Je_tau)-mu)*Iden  
     +        +mu*matmul(Fe_tau,transpose(Fe_tau)))


c      write(*,*),'T_tau = ',T_tau



      ! trace T 
      ! 
      TrT = T_tau(1,1) + T_tau(2,2) + T_tau(3,3)



      ! calculate total jacobian 
      ! 
      call mdet(F_tau,J_tau)
      call mdet(F_t,J_t)

      ! calculate the Kirkoff stress 
      ! 
      Kirk_tau = J_tau * T_tau  





      ! calculate the AR_mat  
      ! 
      
      AR_mat1 = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                  AR_mat1(i,j,k,l) = AR_mat1(i,j,k,l)
     +                 + detFg*(lambda*dlog(Je_tau)-mu+lambda)
     +                 *Finv(j,i)*Finv(l,k)
     +                 + detFg*(mu - lambda*dlog(Je_tau))
     +                 *Finv(l,i)*Finv(j,k)
               enddo
            enddo
         enddo
      enddo 


      AR_mat2 = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  do a = 1,3
                     do b = 1,3          
                       AR_mat2(i,j,k,l) = AR_mat2(i,j,k,l)
     +                 + detFg*mu*F_tau(i,a)*Finv(l,k)*Fginv(a,b)*Fginv(j,b)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo 

      
      AR_mat3 = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3 
                  do b = 1,3         
                 AR_mat3(i,j,k,l) = AR_mat3(i,j,k,l)
     +                 + detFg*mu*Iden(i,k)*Fginv(l,b)*Fginv(j,b)
                  enddo 
               enddo
            enddo
         enddo
      enddo      
      
       
     
      
      
      AR_mat = AR_mat1 + AR_mat2 + AR_mat3


            
      
      
      ! map AR_mat to spatial one 
      
      A_mat = zero 

      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  do m = 1,3
                     do n = 1,3          
                       A_mat(i,j,k,l) = A_mat(i,j,k,l)
     +                 + (1.0/J_tau)*F_tau(j,m)*F_tau(l,n)
     +                 *AR_mat(i,m,k,n)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo      
          
      
      ! spatial modulus tangent
      spTanMod = A_mat


      ! calculate the dqdFFT
      ! 
      dqdFFT = zero
      do i=1,3
         do j = 1,3
            do k = 1,3          
                  dqdFFT(i,j,k) = dqdFFT(i,j,k)
     +            - c_tau*velocity/tmp
     +            *(Iden(i,j)*direction_cur(k,1)  
     +            - (1.0/(tmp**2.0))*direction_cur(i,1)
     +            * direction_cur(j,1)*direction_cur(k,1))     
            enddo
         enddo
      enddo

      dqdFFT2 = zero ! update 
      do i=1,3
         do j = 1,3
            do k = 1,3          
                  dqdFFT2(i,j,k) = dqdFFT2(i,j,k)
     +            - c2_tau*velocity2/tmp
     +            *(Iden(i,j)*direction_cur(k,1)  
     +            - (1.0/(tmp**2.0))*direction_cur(i,1)
     +            * direction_cur(j,1)*direction_cur(k,1))     
            enddo
         enddo
      enddo
      
      dqdFFT3 = zero ! update 
      do i=1,3
         do j = 1,3
            do k = 1,3          
                  dqdFFT3(i,j,k) = dqdFFT3(i,j,k)
     +            - c3_tau*velocity3/tmp
     +            *(Iden(i,j)*direction_cur(k,1)  
     +            - (1.0/(tmp**2.0))*direction_cur(i,1)
     +            * direction_cur(j,1)*direction_cur(k,1))     
            enddo
         enddo
      enddo
      
      
      


      ! calculate the qIden
      ! 
      
      
      qIden = zero
      do i=1,3
         do j = 1,3
            do k = 1,3          
                  qIden(i,j,k) = qIden(i,j,k)
     +            + qflux(i,1)*Iden(j,k)   
            enddo
         enddo
      enddo

      qIden2 = zero ! update 
      do i=1,3
         do j = 1,3
            do k = 1,3          
                  qIden2(i,j,k) = qIden2(i,j,k)
     +            + qflux2(i,1)*Iden(j,k)   
            enddo
         enddo
      enddo

      qIden3 = zero ! update 
      do i=1,3
         do j = 1,3
            do k = 1,3          
                  qIden3(i,j,k) = qIden3(i,j,k)
     +            + qflux3(i,1)*Iden(j,k)   
            enddo
         enddo
      enddo



      ! calculate SpCUMod
      ! 
      
      SpCUMod = qIden - dqdFFT
      SpCUMod2 = qIden2 - dqdFFT2 ! update 
      SpCUMod3 = qIden3 - dqdFFT3 ! update 
c      SpCUMod = - dqdFFT
      
            

      
      
      ! dTdFg
      ! 
      dTdFg1 = zero     
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  dTdFg1(i,j,k,l) = dTdFg1(i,j,k,l)
     +            +(1.0/Je_tau)
     +            *(
     +             (lambda*dlog(Je_tau)-mu-lambda)*Iden(i,j)*Fginv(l,k)
     +            + mu*Be_tau(i,j)*Fginv(k,l)
     +             ) 
               enddo             
            enddo
         enddo
      enddo   
      


      
      dTdFg2 = zero     
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  do a = 1,3
                  dTdFg2(i,j,k,l) = dTdFg2(i,j,k,l)
     +            +(1.0/Je_tau)
     +            *(
     +            -mu*Fe_tau(i,k)*Fe_tau(j,a)*Fginv(l,a)
     +            -mu*Fe_tau(i,a)*Fe_tau(j,k)*Fginv(l,a)
     +            )
                  enddo
               enddo             
            enddo
         enddo
      enddo        
      
      dTdFg = dTdFg1 + dTdFg2
      
      
      ! calculate dFgdc
      ! 
      
      fac4 = k_normal*alpha_normal*(1.0 + k_normal*c_total)**(alpha_normal-1.0)
      fac5 = k_para*alpha_para*(1.0 + k_para*c_total)**(alpha_para-1.0)
      dFgdc = fac4*(Iden - matmul(direction,transpose(direction)))
     +       +fac5*matmul(direction,transpose(direction))      
      
      
      ! dTdc 
      ! 
      
      dTdc = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3         
                  dTdc(i,j) = dTdc(i,j)
     +            + dTdFg(i,j,l,k)*dFgdc(k,l)
               enddo             
            enddo
         enddo
      enddo      
      
     
      
      SpUCMod = dTdc
      
      
      ! time derivative of Jacobian 
      ! 
      dJdt = (J_tau - J_t)/dtime
 

      return
      end subroutine integ

****************************************************************************
      subroutine Hhat(x,alpha,y)

      implicit none

      real*8 x,alpha,y
      
      y = (exp(alpha*x))/(1.0 + exp(alpha*x))

      end subroutine Hhat
****************************************************************************
      subroutine CartoPolar(x,y,r,theta)

      implicit none

      real*8 x,y,r,theta
 
      r = sqrt(x**2.0 + y**2.0)

      theta = atan(y/x)
      
      end subroutine CartoPolar
****************************************************************************
      subroutine solveAlpha(root,args,nargs,rootOld)

      implicit none

      integer maxit,j,nargs

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin

      parameter(maxit=50)
      parameter(xacc=1.d-6,zero=0.d0,one=1.d0)

      rootMax = 1.d0
      rootMin = 0.d0

      x1 = rootMin
      x2 = rootMax
      call phiFunc(x1,FL,DF,args,nargs)
      call phiFunc(x2,FH,DF,args,nargs)

      if(fl*fh.ge.zero) then
         root = rootOld
         write(*,*) 'FYI, root not bracketed on phi'
         write(*,*) 'fl=',fl
         write(*,*) 'fh=',fh
         write(*,*) 'rootOld=',rootOld
         write(80,*) 'FYI, the root is not bracketed on phi'
         write(80,*) 'fl=',fl
         write(80,*) 'fh=',fh
         write(80,*) 'rootOld=',rootOld

         write(*,*) 'args(1)=',args(1)
         write(*,*) 'args(2)=',args(2)
         write(*,*) 'args(3)=',args(3)

         call xit
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      ROOT = rootOld !0.5D0 *( X1 + X2)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD
      
      call phiFunc(root,F,DF,args,nargs)

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
         IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call phiFunc(root,F,DF,args,nargs)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'solveAlpha EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solveAlpha EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solveAlpha

****************************************************************************

      subroutine phiFunc(alpha,f,df,args,nargs)

      implicit none

      integer nargs,NeoHookean,Langevin,material
      parameter(NeoHookean=1,Langevin=2)

      real*8 args(nargs),f,df,alpha,alphaOld,dtime,k

      real*8 zero,one,two,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)

      
      ! Obtain relevant quantities
      !
      alphaOld = args(1)
      dtime    = args(2)
      k        = args(3)


      ! Compute the residual
      !
      f = alpha - alphaOld - dtime*k*((one - alpha)**two)


      ! Compute the tangent
      !
      df = one + two*dtime*k*(one - alpha)

      return
      end subroutine phiFunc

     
************************************************************************
************************************************************************
************************************************************************
************************************************************************   

      subroutine AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rc,Rc2,Rc3,
     +     Kuu,Kuc,Kuc2,Kuc3,
     +     Kcu,Kcc,Kcc2,kcc3,
     +     Kc2u,Kc2c,Kc2c2,kc2c3,
     +     Kc3u,Kc3c,Kc3c2,Kc3c3,
     +     rhs,amatrx)  


      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),RI(nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 KII(nNode,nNode),KuI(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 KIu(nNode,nDim*nNode),amatrx(ndofel,ndofel)
      real*8 Kuc(nDim*nNode,nNode),KIc(nNode,nNode)
      real*8 Kcu(nNode,nDim*nNode),KcI(nNode,nNode)
      real*8 Kcc(nNode,nNode)
      real*8 Rc(nNode,1),Rc2(nNode,1),Kuc2(nDim*nNode,nNode)
      real*8 Kcc2(nNode,nNode),Kc2u(nNode,nDim*nNode),Kc2c(nNode,nNode)
      real*8 Kc2c2(nNode,nNode),kcc3(nNode,nNode),kc2c3(nNode,nNode)
      real*8 Kc3u(nNode,nDim*nNode),Kc3c(nNode,nNode)
      real*8 Kc3c2(nNode,nNode),Kc3c3(nNode,nNode)
      real*8 Kuc3(nDim*nNode,nNode),Rc3(nNode,1)


      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode


      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            ! 
            ! displacement
            ! 
            rhs(A11,1) = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            ! 
            ! neuron cell density 
            ! 
            rhs(A11+2,1) = Rc(i,1)
            rhs(A11+3,1) = Rc2(i,1)     
            rhs(A11+4,1) = Rc3(i,1)     
            !  
         enddo
         ! 
         ! Assemble the element level tangent matrix
         ! 
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               ! 
               ! c1 - c1 , c2 - c2, c3-c3
               ! 
               amatrx(A11+2,B11+2) = Kcc(i,j)
               amatrx(A11+3,B11+3) = Kc2c2(i,j)
               amatrx(A11+4,B11+4) = Kc3c3(i,j)
               ! 
               ! 
               ! displacement - cell1
               ! 
               ! 
               amatrx(A11,B11+2) = Kuc(A12,j)
               amatrx(A11+1,B11+2) = Kuc(A12+1,j)
               ! 
               ! displacement - cell 2 
               !  
               amatrx(A11,B11+3) = Kuc2(A12,j)
               amatrx(A11+1,B11+3) = Kuc2(A12+1,j)
               ! 
               ! displacement - cell 3 
               !  
               amatrx(A11,B11+4) = Kuc3(A12,j)
               amatrx(A11+1,B11+4) = Kuc3(A12+1,j)
               ! 
               ! cell1 - displacement
               ! 
               amatrx(A11+2,B11) = Kcu(i,B12)
               amatrx(A11+2,B11+1) = Kcu(i,B12+1)
               ! 
               ! cell2 - displacement
               ! 
               amatrx(A11+3,B11) = Kc2u(i,B12)
               amatrx(A11+3,B11+1) = Kc2u(i,B12+1)
               ! 
               ! cell3 - displacement
               ! 
               amatrx(A11+4,B11) = Kc3u(i,B12)
               amatrx(A11+4,B11+1) = Kc3u(i,B12+1)
               !  
               ! cell1 - cell2 
               ! 
               amatrx(A11+2,B11+3) = Kcc2(i,j)
               !  
               ! cell1 - cell3 
               ! 
               amatrx(A11+2,B11+4) = Kcc3(i,j)
               !  
               ! cell2 - cell1 
               ! 
               amatrx(A11+3,B11+2) = Kc2c(i,j) 
               !  
               ! cell2 - cell3
               ! 
               amatrx(A11+3,B11+4) = Kc2c3(i,j) 
               !  
               ! cell3 - cell1 
               ! 
               amatrx(A11+4,B11+2) = Kc3c(i,j)
               !  
               ! cell3 - cell2 
               ! 
               amatrx(A11+4,B11+3) = Kc3c2(i,j)
               ! 
            enddo
         enddo
         !
      elseif(nDim.eq.3) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            ! 
            ! displacement
            ! 
            rhs(A11,1)   = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            rhs(A11+2,1) = Ru(A12+2,1)
            ! 
            ! neuron cell density 
            ! 
            rhs(A11+3,1) = Rc(i,1)
            rhs(A11+4,1) = Rc2(i,1)     
            rhs(A11+5,1) = Rc3(i,1)
            ! 
         enddo
         !
         ! Assembly the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               ! 
               ! displacement
               ! 
               amatrx(A11,B11)     = Kuu(A12,B12)
               amatrx(A11,B11+1)   = Kuu(A12,B12+1)
               amatrx(A11,B11+2)   = Kuu(A12,B12+2)
               amatrx(A11+1,B11)   = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
               amatrx(A11+2,B11)   = Kuu(A12+2,B12)
               amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
               amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
               ! 
               ! c1 - c1 , c2 - c2, c3-c3
               ! 
               amatrx(A11+3,B11+3) = Kcc(i,j)
               amatrx(A11+4,B11+4) = Kc2c2(i,j)
               amatrx(A11+5,B11+5) = Kc3c3(i,j) 
               ! 
               ! displacement - cell1
               !                
               amatrx(A11,B11+3) = Kuc(A12,j)
               amatrx(A11+1,B11+3) = Kuc(A12+1,j)
               amatrx(A11+2,B11+3) = Kuc(A12+2,j)
               ! 
               ! displacement - cell 2 
               ! 
               amatrx(A11,B11+4) = Kuc2(A12,j)
               amatrx(A11+1,B11+4) = Kuc2(A12+1,j)
               amatrx(A11+2,B11+4) = Kuc2(A12+2,j)
               ! 
               ! displacement - cell 3
               ! 
               amatrx(A11,B11+5) = Kuc3(A12,j)
               amatrx(A11+1,B11+5) = Kuc3(A12+1,j)
               amatrx(A11+2,B11+5) = Kuc3(A12+2,j)
               ! 
               ! cell1 - displacement
               ! 
               amatrx(A11+3,B11) = Kcu(i,B12)
               amatrx(A11+3,B11+1) = Kcu(i,B12+1)
               amatrx(A11+3,B11+2) = Kcu(i,B12+2)
               ! 
               ! cell2 - displacement
               ! 
               amatrx(A11+4,B11) = Kc2u(i,B12)
               amatrx(A11+4,B11+1) = Kc2u(i,B12+1)
               amatrx(A11+4,B11+2) = Kc2u(i,B12+2)
               ! 
               ! cell3 - displacement
               ! 
               amatrx(A11+5,B11) = Kc3u(i,B12)
               amatrx(A11+5,B11+1) = Kc3u(i,B12+1)
               amatrx(A11+5,B11+2) = Kc3u(i,B12+2)
               !  
               ! cell1 - cell2 
               ! 
               amatrx(A11+3,B11+4) = Kcc2(i,j)
               !  
               ! cell1 - cell3 
               ! 
               amatrx(A11+3,B11+5) = Kcc3(i,j)
               !  
               ! cell2 - cell1 
               ! 
               amatrx(A11+4,B11+3) = Kc2c(i,j) 
               !  
               ! cell2 - cell3
               ! 
               amatrx(A11+4,B11+5) = Kc2c3(i,j)               
               !  
               ! cell3 - cell1 
               ! 
               amatrx(A11+5,B11+3) = Kc3c(i,j)
               !  
               ! cell3 - cell2 
               ! 
               amatrx(A11+5,B11+4) = Kc3c2(i,j)
               !  
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt


************************************************************************

      subroutine xint2D9pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 9 gauss ponits for integration as shown
      !
      !              A eta (=xi_2)
      !              |
      !              |
      !        o-----o-----o
      !        |  7  8  9  |
      !        |     |     |
      !        o  4  5--6--o---> xi (=xi_1)
      !        |           |
      !        |  1  2  3  |
      !        o-----o-----o
      !
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(9,2), w(9)


      ! Number of Gauss points
      !
      nIntPt = 9


      ! Gauss weights
      !
      w(1) = (5.d0/9.d0)*(5.d0/9.d0)
      w(2) = (5.d0/9.d0)*(8.d0/9.d0)
      w(3) = (5.d0/9.d0)*(5.d0/9.d0)
      w(4) = (5.d0/9.d0)*(8.d0/9.d0)
      w(5) = (8.d0/9.d0)*(8.d0/9.d0)
      w(6) = (5.d0/9.d0)*(8.d0/9.d0)
      w(7) = (5.d0/9.d0)*(5.d0/9.d0)
      w(8) = (5.d0/9.d0)*(8.d0/9.d0)
      w(9) = (5.d0/9.d0)*(5.d0/9.d0)
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(3.d0/5.d0)
      xi(1,2) = -dsqrt(3.d0/5.d0)
      xi(2,1) = 0.d0
      xi(2,2) = -dsqrt(3.d0/5.d0)
      xi(3,1) = dsqrt(3.d0/5.d0)
      xi(3,2) = -dsqrt(3.d0/5.d0)
      xi(4,1) = -dsqrt(3.d0/5.d0)
      xi(4,2) = 0.d0
      xi(5,1) = 0.d0
      xi(5,2) = 0.d0
      xi(6,1) = dsqrt(3.d0/5.d0)
      xi(6,2) = 0.d0
      xi(7,1) = -dsqrt(3.d0/5.d0)
      xi(7,2) = dsqrt(3.d0/5.d0)
      xi(8,1) = 0.d0
      xi(8,2) = dsqrt(3.d0/5.d0)
      xi(9,1) = dsqrt(3.d0/5.d0)
      xi(9,2) = dsqrt(3.d0/5.d0)

      return
      end subroutine xint2D9pt
      
************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      end subroutine xint3D1pt
     
************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

************************************************************************

      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),w(1),zero,one,two
      parameter(zero=0.d0,one=1.d0,two=2.d0)


      ! Gauss weights
      !
      w(1) = two
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = zero
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D1pt

************************************************************************

      subroutine xintSurf2D2pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(2),yLocal(2),w(2),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = -dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D2pt

************************************************************************

      subroutine xintSurf2D3pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(3),yLocal(3),w(3),zero,one,two,three,five,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,five=5.d0,
     +     eight=8.d0,nine=9.d0)


      ! Gauss weights
      !
      w(1) = five/nine
      w(2) = eight/nine
      w(3) = five/nine
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = -one
         xLocal(2) = zero
         yLocal(2) = -one
         xLocal(2) = dsqrt(three/five)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(three/five)
         xLocal(2) = one
         yLocal(2) = zero
         xLocal(3) = one
         yLocal(3) = dsqrt(three/five)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = one
         xLocal(2) = zero
         yLocal(2) = one
         xLocal(3) = dsqrt(three/five)
         yLocal(3) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(three/five)
         xLocal(2) = -one
         yLocal(2) = zero
         xLocal(3) = -one
         yLocal(3) = -dsqrt(three/five)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D3pt

************************************************************************

      subroutine xintSurf3D1pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 1 gauss point for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  zLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),zLocal(1),w(1),zero,one,four
      parameter(zero=0.d0,one=1.d0,four=4.d0)


      ! Gauss weights
      !
      w(1) = four
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = one
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = -one
         zLocal(1) = zero
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = zero
         zLocal(1) = zero
      elseif(face.eq.5) then
         xLocal(1) = zero
         yLocal(1) = one
         zLocal(1) = zero
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = zero
         zLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D1pt

************************************************************************

      subroutine xintSurf3D4pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 4 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  yLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(4),yLocal(4),zLocal(4),w(4),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = -one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = -one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = -one
      elseif(face.eq.2) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = one
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = -one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = -one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.5) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = -one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D4pt
     
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i,j

      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3)
      real*8 d2shxi(8,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)
      
      return
      end subroutine calcShape3DLinear

!************************************************************************


      subroutine computeSurf(xLocal,yLocal,face,coords,sh,ds,normal)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements

      implicit none

      integer face

      real*8 xLocal,yLocal,ds,dshxi(4,2),sh(4),dXdXi,dXdEta,dYdXi
      real*8 dYdEta,one,coords(2,4),fourth,shape,normal(2,1)
      parameter(one=1.d0,fourth=1.d0/4.d0)

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)
      
      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


      ! Jacobian of the mapping
      !
      if((face.eq.2).or.(face.eq.4)) then
         ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif((face.eq.1).or.(face.eq.3)) then
         ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      else
         write(*,*) 'never should get here'
         call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.2).or.(face.eq.4)) then
         normal(1,1) = dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         normal(2,1) = -dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         if(face.eq.4) normal = -normal
      elseif((face.eq.1).or.(face.eq.3)) then
         normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         if(face.eq.3) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif

      return
      end subroutine computeSurf

************************************************************************

      subroutine computeSurf3D(xLocal,yLocal,zLocal,face,coords,sh,dA)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 8-node brick elements

      implicit none

      integer face,stat,i,j,k

      real*8 xLocal,yLocal,zLocal,dA,dshxi(8,3),sh(8),zero,dsh(8,3),one
      real*8 coords(3,8),two,eighth,mapJ(3,3),mag,normal(3,1)

      real*8 dXdXi,dXdEta,dXdZeta,dYdXi,dYdEta,dYdZeta,dZdXi,dZdEta
      real*8 dZdZeta

      parameter(one=1.d0,two=2.d0,eighth=1.d0/8.d0,zero=0.d0)

      ! The shape functions
      !
      sh(1) = eighth*(one - xLocal)*(one - yLocal)*(one - zLocal)
      sh(2) = eighth*(one + xLocal)*(one - yLocal)*(one - zLocal)
      sh(3) = eighth*(one + xLocal)*(one + yLocal)*(one - zLocal)
      sh(4) = eighth*(one - xLocal)*(one + yLocal)*(one - zLocal)
      sh(5) = eighth*(one - xLocal)*(one - yLocal)*(one + zLocal)
      sh(6) = eighth*(one + xLocal)*(one - yLocal)*(one + zLocal)
      sh(7) = eighth*(one + xLocal)*(one + yLocal)*(one + zLocal)
      sh(8) = eighth*(one - xLocal)*(one + yLocal)*(one + zLocal)


      ! Shape function derivatives
      !
      dshxi(1,1) = -eighth*(one - yLocal)*(one - zLocal)
      dshxi(1,2) = -eighth*(one - xLocal)*(one - zLocal)
      dshxi(1,3) = -eighth*(one - xLocal)*(one - yLocal)
      dshxi(2,1) = eighth*(one - yLocal)*(one - zLocal)
      dshxi(2,2) = -eighth*(one + xLocal)*(one - zLocal)
      dshxi(2,3) = -eighth*(one + xLocal)*(one - yLocal)
      dshxi(3,1) = eighth*(one + yLocal)*(one - zLocal)
      dshxi(3,2) = eighth*(one + xLocal)*(one - zLocal)
      dshxi(3,3) = -eighth*(one + xLocal)*(one + yLocal)
      dshxi(4,1) = -eighth*(one + yLocal)*(one - zLocal)
      dshxi(4,2) = eighth*(one - xLocal)*(one - zLocal)
      dshxi(4,3) = -eighth*(one - xLocal)*(one + yLocal)
      dshxi(5,1) = -eighth*(one - yLocal)*(one + zLocal)
      dshxi(5,2) = -eighth*(one - xLocal)*(one + zLocal)
      dshxi(5,3) = eighth*(one - xLocal)*(one - yLocal)
      dshxi(6,1) = eighth*(one - yLocal)*(one + zLocal)
      dshxi(6,2) = -eighth*(one + xLocal)*(one + zLocal)
      dshxi(6,3) = eighth*(one + xLocal)*(one - yLocal)
      dshxi(7,1) = eighth*(one + yLocal)*(one + zLocal)
      dshxi(7,2) = eighth*(one + xLocal)*(one + zLocal)
      dshxi(7,3) = eighth*(one + xLocal)*(one + yLocal)
      dshxi(8,1) = -eighth*(one + yLocal)*(one + zLocal)
      dshxi(8,2) = eighth*(one - xLocal)*(one + zLocal)
      dshxi(8,3) = eighth*(one - xLocal)*(one + yLocal)


      dXdXi = zero
      dXdEta = zero
      dXdZeta = zero
      dYdXi = zero
      dYdEta = zero
      dYdZeta = zero
      dZdXi = zero
      dZdEta = zero
      dZdZeta = zero
      do k=1,8
         dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
         dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
         dXdZeta = dXdZeta + dshxi(k,3)*coords(1,k)
         dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
         dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
         dYdZeta = dYdZeta + dshxi(k,3)*coords(2,k)
         dZdXi = dZdXi + dshxi(k,1)*coords(3,k)
         dZdEta = dZdEta + dshxi(k,2)*coords(3,k)
         dZdZeta = dZdZeta + dshxi(k,3)*coords(3,k)
      enddo


      ! Jacobian of the mapping
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdEta - dYdEta*dZdXi)**two
     +        + (dXdXi*dZdEta - dXdEta*dZdXi)**two
     +        + (dXdXi*dYdEta - dXdEta*dYdXi)**two
     +        )
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdZeta - dYdZeta*dZdXi)**two
     +        + (dXdXi*dZdZeta - dXdZeta*dZdXi)**two
     +        + (dXdXi*dYdZeta - dXdZeta*dYdXi)**two
     +        )
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         dA = dsqrt(
     +          (dYdEta*dZdZeta - dYdZeta*dZdEta)**two
     +        + (dXdEta*dZdZeta - dXdZeta*dZdEta)**two
     +        + (dXdEta*dYdZeta - dXdZeta*dYdEta)**two
     +        )
         else
            write(*,*) 'never should get here'
            call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         normal(1,1) = dYdXi*dZdEta - dYdEta*dZdXi
         normal(2,1) = dXdXi*dZdEta - dXdEta*dZdXi
         normal(3,1) = dXdXi*dYdEta - dXdEta*dYdXi
         if(face.eq.1) normal = -normal
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         normal(1,1) = dYdXi*dZdZeta - dYdZeta*dZdXi
         normal(2,1) = dXdXi*dZdZeta - dXdZeta*dZdXi
         normal(3,1) = dXdXi*dYdZeta - dXdZeta*dYdXi
         if(face.eq.5) normal = -normal
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         normal(1,1) = dYdEta*dZdZeta - dYdZeta*dZdEta
         normal(2,1) = dXdEta*dZdZeta - dXdZeta*dZdEta
         normal(3,1) = dXdEta*dYdZeta - dXdZeta*dYdEta
         if(face.eq.6) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif
      mag = dsqrt(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
      normal(1,1) = normal(1,1)/mag
      normal(2,1) = normal(2,1)/mag
      normal(3,1) = normal(3,1)/mag

      end subroutine computeSurf3D

************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv

      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

****************************************************************************
      subroutine CartoShere(x,y,z,r,theta,phi)

      implicit none

      real*8 x,y,z,r,theta,phi

      r = sqrt(x**2.0 + y**2.0 + z**2.0)

      theta = acos(z/r)

      phi = atan(y/x)


      end subroutine CartoShere
****************************************************************************
