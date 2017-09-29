module mod_edge
implicit none
SAVE
! constantes usadas no programa
  integer,    parameter :: dp     = kind(1.d0)            ! precisão dupla
  real(dp),   parameter :: pi     = 3.141592653589793_dp     ! constante pi
  complex(dp),parameter :: ci     = (0.d0,1.d0)           ! i complexo
  real(dp),   parameter :: mu0    = 4.d0*pi*1e-7     
  real(dp),   parameter :: e0     = 8.8541878176*1e-12

! VARIÁVEIS DO PROBLEMA
!=====================================================================
! vetor fonte
  complex(dp) :: kP, zethaP, zetha, dzetha, k0, IMP1, IMP2, RTM
  real(dp)    :: ethaP, etha, detha, resP   
  complex(dp) :: Ex(3), Ez(3), Et(3)
  real(dp),allocatable, dimension(:) :: mont_i, mont_f, sol_i, sol_f
  real(dp), allocatable:: freq(:)
  integer:: grau, el_=0, ll 

! VARIÁVEIS DE ENTRADA 
!=====================================================================
    
! variaveis relacionadas à malha

  integer:: nnos,    &  
            nelem,   &  
            nedges 
                                    
  integer, allocatable:: fronteira(:),         &   
                         indice_local(:,:),    &
                         edges(:,:),           &             
                         indice_local_e(:,:),  &
                         fronteira_e(:)   
  
  real(dp), allocatable:: coordenadas(:,:),    &
                          resistividade(:)  


! VARIAVEIS RELACIONADAS À SOLUÇÃO NUMÉRICA
!====================================================================                                                    
  integer:: numb,    &   
            nbanda
  
  integer, allocatable::  ind_b(:)
  real(dp), allocatable:: ordem(:,:)         
  
  real(dp):: a(3),   & 
             b(3),   &
             c(3),   &
             Area,   &
             Amn(3), & 
             Bmn(3), &
             Cmn(3), &
             Dmn(3), &
             lmn(3), &
             cos1(3), &
             cos2(3) 

  real(dp):: intx,   & 
             intz,   &
             intx2,  &
             intz2,  &
             intxz,  &
             inte,   &
             intez
  
  real(dp):: xe(3),  &
             ze(3) 
                                                                           
  complex(dp), allocatable:: M(:,:), &
                             F(:)
                            
!=======================================================================
 

  Contains


  subroutine campoE(e)
  implicit none
    integer, intent(in):: e 
    integer:: v(4), i
    real(dp):: prof(3)

    ethaP=1.d0/resP
    kP = (1.d0-ci)*dsqrt(pi*freq(ll)*mu0*ethaP)
    zethaP = ci*2.d0*pi*freq(ll)*mu0
   
    etha=1.d0/resistividade(e)
    zetha = zethaP
    dzetha = 0.d0   !zetha-zethaP

    k0  = 2.d0*pi*freq(ll)*sqrt(e0*mu0)
                            
    ! Devido às condições de fronteira ser sobre o campo elétrico: 
    IMP1 = ci*k0/zethaP    !   Corresponde à adimitancia do MEF nodal
    IMP2 = (ci*kP/zethaP)  !                  II


    RTM = (IMP1-IMP2)/(IMP1+IMP2)
    
    v = (/1, 2, 3, 1/)

    do i=1,3
       prof(i)=(ze(v(i+1))+ze(v(i)))/2.d0    !centro da aresta
    end do

    if(media(ze,3) < 0.d0) then   
           detha = 0.d0
           if(grau==1) Ex = exp(ci*k0*prof) + RTM*exp(ci*k0*prof)
           if(grau==2) Ex = exp(ci*k0*ze) + RTM*exp(ci*k0*ze)
    else
           if(grau==1) Ex = exp(-ci*kP*prof)
           if(grau==2) Ex = exp(-ci*kP*ze) 
           detha = etha-ethaP
           el_=el_+1
    end if
           
     
    Ez = 0.d0

    Et(:) = (Ex*cos1(:) + Ez*cos2(:))
    

  end subroutine campoE 



  subroutine calc_banda()
  implicit none
  
    integer:: i, aux

        
    nbanda=0
    do i=1, nelem
       aux = max( abs( indice_local_e(i,1)-indice_local_e(i,2) ), &
                  abs( indice_local_e(i,2)-indice_local_e(i,3) ), &
                  abs( indice_local_e(i,1)-indice_local_e(i,3) ) )
       if ( aux > nbanda ) then
             nbanda = aux
       end if
    end do

  end subroutine calc_banda

  subroutine solver_edge()
  implicit none

    integer     :: e                ! elemento 
    integer     :: i, j, v(4)             ! contadores
    complex(dp) :: Me(3,3), Ke(3,3), Qe(3,3)
    complex(dp) :: Fe(3)
    complex(dp):: alfax(3), alfaz(3), betax(3), betaz(3), gamax(3), gamaz(3)
!    complex(dp):: alfax, alfaz, betax, betaz, gamax, gamaz

    call cpu_time(mont_i(ll)) 

    v=(/1,2,3,1/)

!    write(*,*) 'Digite o grau das bases de interpolação do campo primário (1 ou 2)'
!    read(*,*) grau

    grau=1

    M = (0.d0, 0.d0)         ! Matriz global
    F = (0.d0, 0.d0)         ! Vetor fonte & solução nas arestas

    ! laco de elementos       
    do e = 1, nelem          
 
       call coef_bases(e)

       call int_quad()        
 
       call campoE(e)
          

       ! Elemento Local
       !===================================================================================
       
       Ke = (0.d0,0.d0); Qe = (0.d0,0.d0); Me = (0.d0,0.d0); Fe = (0.d0,0.d0);

       do i=1,3
            
          do j=1,3

             Ke(i,j) = ( lmn(i)*lmn(j) ) *Bmn(i)*Bmn(j)    ! int [dot( curl(N_m), curl(N_n) )]dxdy
  
             Qe(i,j) = ( lmn(i)*lmn(j) )                                 &   
                              * ( Amn(i)*Amn(j)*Area                     &
                              + ( Amn(i)*Bmn(j) + Amn(j)*Bmn(i))*intz    &    
                              +   Bmn(i)*Bmn(j)*intz2                    &  
                                                +                        &  ! int [dot( N_m, N_n )]dxdy
                                  Cmn(i)*Cmn(j)*Area                     &
                              + ( Cmn(i)*Dmn(j) + Cmn(j)*Dmn(i) )*intx   &   
                              +   Dmn(i)*Dmn(j)*intx2                  )
         
 
                 Fe(i) = Fe(i) + Et(j)*lmn(i)*lmn(j)                   &   
                              * ((Amn(i)*Amn(j)*Area                      &
                              + ( Amn(i)*Bmn(j) + Amn(j)*Bmn(i))*intz     &    
                              +   Bmn(i)*Bmn(j)*intz2)                    &       
                                                +                         &           ! int [dot( E, N_m )]dxdy
                                ( Cmn(i)*Cmn(j)*Area                      &
                              + ( Cmn(i)*Dmn(j) + Cmn(j)*Dmn(i) )*intx    &   
                              +   Dmn(i)*Dmn(j)*intx2)                  )
                             
          end do 
         
      
       
       end do
 
       Ke =  Ke / (4.d0*Area**3) 
       Qe =  Qe / (16.d0*Area**4)  
       Fe =  Fe / (16.d0*Area**4)
      
       Me = Ke/zetha  + etha*Qe;       
       Fe = -(dzetha/zetha *ethaP + detha)*Fe 


       call CASSEMBG(M,nbanda+1,nbanda,indice_local_e(e,:),3,Me,3) 
  
       call CASSEMBB(indice_local_e(e,:),3,Fe,F)
       !====================================================================================

    end do
    

    call cpu_time(mont_f(ll))     

    write(*,*) 'SISTEMA GLOBAL MONTADO' 
    
    call cpu_time(sol_i(ll))

    ! condicoes de contorno

    call CDIRICHB(F,ind_b,numb,dcmplx(0.0d0+0.00*ci)*ind_b)

    call CDIRICHG(M,nbanda+1,nbanda,ind_b,numb)
    
    call CREDUCEG(M,nbanda+1,nedges,nbanda)
  
    call CSOLVEG(M,nbanda+1,nedges,nbanda,F)

  call cpu_time(sol_f(ll))

  write(*,*) 'SISTEMA LINEAR RESOLVIDO'

  end subroutine solver_edge


  subroutine coef_bases(e)
  implicit none

    integer, intent(in):: e   
    integer:: j(3), k(3), w(3), t3(3), t4(3)
    integer:: i
                  
    t3=(/1,2,1/)
    t4=(/2,3,3/) 

    xe(:) = coordenadas(indice_local(e,:),1)
    ze(:) = coordenadas(indice_local(e,:),2)
    
    Area = ((xe(2)-xe(1))*(ze(3)-ze(1)) - ((xe(3)-xe(1))*(ze(2)-ze(1))))/2.d0
    
    ! BASE NODAL  
    ! phi(i) = (a(i) + b(i)x + c(i)y) / (2 Area)
    ! indice "i" referente ao nó associado à base nodal
    j = (/2, 3, 1/); k = (/3, 1, 2/)    
    do i=1,3
       a(i) = xe(j(i))*ze(k(i)) - xe(k(i))*ze(j(i)) 
       b(i) = ze(j(i)) - ze(k(i))                      ! Volakis(1998), pg. 102    
       c(i) = xe(k(i)) - xe(j(i)) 
    end do
    
   ! BASE DE WHITNEY
   ! N(m) = l(m) ( phi(i) grad( phi(j) ) - phi(j)( grad(phi(i)) ) ) 
    
   ! Nx(m) = l(m) ( A(m) + B(m) ze ) / (4  Area**2)
   ! Ny(m) = l(m) ( C(m) + D(m) xe ) / (4  Area**2)                                                       
   ! indice m referente à aresta associada à base de Whitney
    do i=1, 3 
     
       lmn(i) = dsqrt( ( xe(t4(i))  - xe( t3(i) ) )**2      &
                     + ( ze( t4(i) ) - ze( t3(i) ) )**2 )
       Amn(i) = a(t3(i))*b(t4(i)) - a(t4(i))*b(t3(i))           ! Volakis(1998), pg. 143    
       Bmn(i) = c(t3(i))*b(t4(i)) - c(t4(i))*b(t3(i))
       Cmn(i) = a(t3(i))*c(t4(i)) - a(t4(i))*c(t3(i))

       cos1(i) = (xe(t4(i)) - xe(t3(i))) / lmn(i)              
       cos2(i) = (ze(t4(i)) - ze(t3(i))) / lmn(i)
                   
    end do
       
    do i=1,3   
       if ( indice_local( e, t3(i) ) > indice_local( e, t4(i) ) ) then

                  Amn(i)=-Amn(i)
                  Bmn(i)=-Bmn(i)          ! Inverte o sentido se necessário
                  Cmn(i)=-Cmn(i)

                  cos1(i) = - cos1(i) 
                  cos2(i) = - cos2(i) 
                                   
       end if

    end do

    Dmn = -Bmn
                
  end subroutine coef_bases


  real(dp) function media(x, n)
  implicit none

    integer, intent(in):: n 
    real(dp), intent(in):: x(n)

    media=sum(x)/n

  end function media
  

  subroutine int_quad()
  implicit none

    real(dp), dimension(3):: xi, zi, wi, fxz
    real(dp):: IQ(5)
    integer:: p(5), q(5)
    integer:: i, j, prec
    real(dp)::aux
    complex:: IQQ(2)

    prec=3
 
    wi(1:3)=(/ 0.333333333333333333_dp, 0.33333333333333333333_dp, 0.333333333333333333_dp/)
    xi(1:3)=(/ 0.5d0, 0.5d0, 0.0d0/)
    zi(1:3)=(/ 0.0d0, 0.5d0, 0.5d0/) 
    
!     wi(1:7)=(/ 0.225000d0, 0.125939d0, 0.125939d0, &
!                0.125939d0, 0.132394d0, 0.132394d0, &
!                0.132394d0 /)

!     xi(1:7)=(/ 0.333333333333333333d0, 0.797427d0, 0.101287d0, 0.101287d0, 0.059716d0, 0.470142d0, 0.470142d0/)
!     zi(1:7)=(/ 0.333333333333333333d0, 0.101287d0, 0.797427d0, 0.101287d0, 0.470142d0, 0.059716d0, 0.470142d0/) 
    

    call Trans(wi,xi,zi, prec)

    p = (/1, 0, 2, 0, 1/)
    q = (/0, 1, 0, 2, 1/)

    do j=1, 5
       fxz   = ( xi**p(j) ) * ( zi**q(j) )
       IQ(j) = 0.d0
       do i=1, prec
          IQ(j) = IQ(j) + wi(i)*fxz(i) / 2.d0
       end do
    end do

    intx=IQ(1)
    intz=IQ(2)
    intx2=IQ(3)
    intz2=IQ(4)
    intxz=IQ(5)

  end subroutine int_quad

  subroutine Trans(wi,xi,zi, n)
  implicit none
 
    integer, intent(in):: n 
    real(dp), intent(inout):: wi(n), xi(n), zi(n)
    real(dp):: x(n), z(n)

    x = xi; z = zi

    xi = ( 1.d0-x-z )*xe(1) + x*xe(2) + z*xe(3)
    zi = ( 1.d0-x-z )*ze(1) + x*ze(2) + z*ze(3)
    wi = 2.d0*wi*Area

  end subroutine Trans


  subroutine buble_sort(x,y,n)

    implicit none
    
    integer, intent(in):: n 
    real(dp), intent(inout):: x(n)
    complex(dp), intent(inout):: y(n)
    real(dp):: aux
    complex(dp):: aux2
    integer:: i, j, troca
   
    j=0.
    do
         troca=0
         j=j+1    
         do i=1, n-j
        
            if( x(i) > x(i+1) ) then
                 aux = x(i+1);  aux2 = y(i+1)
                 x(i+1) = x(i); y(i+1) = y(i)
                 x(i) = aux;    y(i) = aux2
                 troca=1
            end if
                   
         end do
         if(troca == 0) exit    

     end do


  end subroutine buble_sort 


!================================================================
! SUBROTINAS PARA 
!                   SOLUÇÃO DO SISTEMA LINEAR
!================================================================

		
		SUBROUTINE CASSEMBB(node,nnode,velem,b)
		Implicit none
      INTEGER(4),Intent(In)                    :: nnode
		Integer(4)                               :: i,j
		Integer(4),Dimension(:),Intent(In)       :: node
      COMPLEX(8),Dimension(:),Intent(Inout)    :: b,velem

!     CASSEMBB/EGS 2000 Versão 1.0 05/89, Versão 2.0 10/2000.
!     Autor: Luiz Rijo
!
!     Constrói o vetor fonte b
!
!     ENTRADA
!
!        node   INTEGER (*)
!               Vetor dos n¢s do elemento
!
!        nnode  INTEGER
!               N£mero de n¢s do elemento
!
!        velem  COMPLEX*8 (*)
!               Vetor fonte do elemento               
!
!     SAIDA
!
!        b      COMPLEX*8 (*)
!               Vetor fonte 

      DO i = 1,nnode
         j = node(i)
         b(j) = b(j)+velem(i)       
      ENDDO
      RETURN
      END SUBROUTINE CASSEMBB
!==========================================================
      
		SUBROUTINE CASSEMBG(g,ldg,m,node,nnode,melem,ldm)
		Implicit none
      INTEGER(4),Intent(In)                           :: m,ldg,ldm,nnode
      INTEGER(4)                                      :: i,j,ij,ig,jg,m1,m2
		Integer(4),Dimension(:),Intent(in)              :: node
      COMPLEX(8),Dimension(:,:),Intent(in)            :: melem
		Complex(8),Dimension(:,:),Intent(inout)         :: g

!     CASSEMBB/EGS 2000 Versão 1.0 05/89, Versão 2.0 10/2000.
!     Autor: Luiz Rijo
!     
!     Constr¢i a matriz global g na forma semi-bandeada horizontal
!
!             *  * 13 24 35 46 57 68 79
!             * 12 23 34 45 56 67 78 89
!            11 22 33 44 55 66 77 88 99
!
!
!     ENTRADA
!
!        g      COMPLEX*8 (ldg,n)
!               Matriz global na forma semi-badeada horizontal
!
!        ldg    INTEGER
!               DimensÆo principal de g
!               O valor de ldg deve ser maior ou igual a m+1
!
!        m      INTEGER
!               N£mero de diagonais acima da diagonal principal
!
!        node   INTEGER (ldn,*)
!               Vetor dos n¢s do elemento
!
!        nnode  INTEGER
!               N£mero de n¢s do elemento
!
!        melem  COMPLEX*8 (ldn,*)
!               Matriz do elemento
!
!     SAIDA
!
!         g     Matriz global na forma semi-bandeada horizontal

      m1 = m+1
      m2 = m+2
      DO i = 1,nnode
         ig = node(i)
         DO j = 1,nnode
            jg = m1-node(j)+ig
            IF(jg.LT.m2) THEN
              ij = ig+m1-jg
              g(jg,ij) = g(jg,ij)+melem(j,i)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CASSEMBG
!==========================================================
      
		SUBROUTINE CDIRICHB(b,noded,nnoded,bd)
		Implicit none
      INTEGER(4),Intent(In)                 :: nnoded
		Integer(4),Dimension(:),Intent(in)    :: noded
		COMPLEX(8),Dimension(:),Intent(in)    :: bd
		COMPLEX(8),Dimension(:),Intent(inout) :: b
      COMPLEX(8)                            :: big
		Integer(4)                            :: i,j

!     CASSEMBB/EGS 2000 Versão 1.0 05/89, Versão 2.0 10/2000.
!      Autor Luiz Rijo
!      
!      Impäe as condiäes de fronteira Dirichlet no vetor b
! 
!      ENTRADA
! 
!         b      COMPLEX*8 (*)
!                Vetor fonte
! 
!         noded  INTEGER (*)
!                Vetor com os n¢s da fronteira de Dirichlet
! 
!         nnoded INTEGER
!                N£mero de n¢s na fronteira Dirichlet
! 
!         bd     COMPLEX*8 (*)
!                Vetor dos valores de b na fronteira de Dirichlet
! 
!      SAIDA
! 
!         b      COMPLEX*8 (*)
!                Vetor fonte modificado
      big = (1.0E+30,0.0)
      DO i = 1,nnoded
         j = noded(i)
         b(j) = big*bd(i)
      ENDDO
      RETURN
      END SUBROUTINE CDIRICHB

!==========================================================

      SUBROUTINE CDIRICHG(g,ldg,m,noded,nnoded)
		Implicit none
      INTEGER(4)                              :: i,j,m1
		Integer(4),Intent(in)                   :: m,ldg,nnoded
		Integer(4),Dimension(:), Intent(In)     :: noded
      COMPLEX(8)                              :: big
		Complex(8),Dimension(:,:),Intent(inout) :: g
!
!     CASSEMBB/EGS 2000 Versão 1.0 05/89, Versão 2.0 10/2000.
!     Autor: Luiz Rijo
!     
!     Impäe as condiäes de fronteira Dirichlet na matriz g armazenada
!     na forma semi-bandeada horizontal
!
!             *  * 13 24 35 46 57 68 79
!             * 12 23 34 45 56 67 78 89
!            11 22 33 44 55 66 77 88 99
!
!
!     ENTRADA
!
!        g      COMPLEX*8 (ldg,n)
!               Matriz global na forma semi-badeada horizontal
!
!        ldg    INTEGER
!               DimensÆo principal da array g
!               O valor de ldg deve ser maior ou igual a m+1.
!
!        m      INTEGER
!               N£mero de diagonais acima da diagonal principal.
!
!        noded  INTEGER
!               Vetor com os n¢s da fronteira do tipo Dirichlet
!
!        nnoded INTEGER
!               Comprimento do vetor noded
!
!     SAIDA
!
!        g      Matriz global modificada
!   
      m1 = m+1
      big = (1.0E+30,0.0)
      DO i = 1,nnoded
         j = noded(i)
         G(m1,j) = big
!   G(m1,j) = big * G(j)
      ENDDO
      RETURN
      END SUBROUTINE CDIRICHG
!**********************************************************

      SUBROUTINE CREDUCEG(g,ldg,n,m)
		Implicit none
      INTEGER(4)                                :: i,j,k,l,ij,kl,ml,m1
      INTEGER(4),Intent(in)                     :: m,n,ldg
      COMPLEX(8),Dimension(:,:),Intent(inOut)   :: g
		Complex(8)                                :: c

!     CASSEMBB/EGS 2000 Versão 1.0 05/89, Versão 2.0 10/2000.
!     Autor: Luiz Rijo
!     
!     Triangulariza a matriz global g na forma semi-badeada horizontal
!
!             *  * 13 24 35 46 57 68 79
!             * 12 23 34 45 56 67 78 89
!            11 22 33 44 55 66 77 88 99
!
!
!     ENTRADA
!
!        g     COMPLEX*8 (ldg,n)
!              Matriz global na forma semi-badeada horizontal
!
!        ldg   INTEGER
!              DimensÆo principal da array g
!              O valor de ldg deve ser maior ou igual a m+1
!
!        n     INTEGER
!              Ordem da matriz global g
!
!        m     INTEGER
!              N£mero de diagonais acima da diagonal principal
!
!     SAIDA
!
!        g     Matriz global triangularizada 
   
      m1 = m+1
      DO j = 1,n
         ml = MAX0(1,j-n+m1)
         DO i = m,ml,-1
            k = j-i+m1
            IF(CDABS(g(i,k)).NE.0.0) THEN
              c = g(i,k)/g(m1,j)
              ij = m1+1
              kl = k-1
              DO l = i,ml,-1
                 ij = ij-1
                 kl = kl+1
                 g(ij,kl) = g(ij,kl)-c*g(l,kl)
              ENDDO
              g(i,k) = c
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CREDUCEG
!***********************************************************

      SUBROUTINE CSOLVEG(g,ldg,n,m,b)
		Implicit none
      INTEGER(4),Intent(in)                   :: m,n,ldg
      INTEGER(4)                              :: i,j,k,l,ml,m1
		COMPLEX(8),Dimension(:,:),Intent(in)    :: g 
		COMPLEX(8),Dimension(:),  Intent(inOut) :: b

!      CSOLVEG/EGS 2000 Versão 1.0 05/89, Versão 2.0 10/2000.
!      Autor: Luiz Rijo
!      
!      Reduz o vetor fonte b e soluciona o sistema de equaäes na forma
!      semi-bandeada horizontal.
! 
!      ENTRADA
! 
!         g     COMPLEX*8 (ldg,*)
!               Matriz global triangularizada pela rotina CREDUCEG
! 
!         ldg   INTEGER
!               DimensÆo principal da array g
!               O valor de ldg deve ser maior ou igual a m+1
! 
!         n     INTEGER
!               Ordem da matriz global g
! 
!         m     INTEGER
!               N£umero de diagonais acima da diagonal principal
! 
!         b     COMPLEX*8 (*)
!              Vetor fonte
! 
!      SAIDA
! 
!         b     Vetor com a soluÆo do sistema
  
      m1 = m+1
      DO j = 1,n
         ml = MAX0(1,j-n+m1)
         DO i = m,ml,-1
              k = j-i+m1
              b(k) = b(k)-g(i,k)*b(j)
         ENDDO
         b(j) = b(j)/g(m1,j)
      ENDDO

      DO j = 2,n
         k = n+1-j
         ml = MAX0(1,k-n+m1)
         DO i = m,ml,-1
            l = k-i+m1
            b(k) = b(k)-g(i,l)*b(l)
         ENDDO
      ENDDO
      RETURN
      END Subroutine CSOLVEG

end module mod_edge
