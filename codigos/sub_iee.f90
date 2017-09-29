module mod_edge_iee
implicit none
SAVE
! constantes usadas no programa
  integer,    parameter :: dp     = kind(1.d0)            ! precisão dupla
    
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

                                
!=======================================================================
 

  Contains


  subroutine iee_edge
  implicit none

    integer:: arestas(3*nelem,4)
    integer:: i, j, k, aux, l, m
   
    arestas(:,4) = 1; k = 1
    do i=1,nelem
       arestas(k,2:3)=indice_local(i,1:2); arestas(k,1)=k
       arestas(k+1,2:3)=indice_local(i,2:3); arestas(k+1,1)=k+1
       arestas(k+2,2:3)=indice_local(i,3:1:-2);arestas(k+2,1)=k+2
       k=k+3;
    end do

    do i=1,3*nelem
       if ( arestas(i,2) > arestas(i,3) ) then
                 aux = arestas(i,2)
                 arestas(i,2) = arestas(i,3)
                 arestas(i,3) = aux
                 arestas(i,4) = -1
       end if
    end do

    do i=4, 3*nelem
       do j=1,i-1
          if( arestas(i,2) == arestas(j,2) .and. arestas(i,3) == arestas(j,3) ) then
               arestas(i,1) = arestas(j,1);
               arestas(i+1:3*nelem,1) = arestas(i+1:3*nelem,1)-1;
          end if
       end do
    end do
  
    k = 1
    do i=1, nelem
       indice_local_e(i,1) = arestas(k,1)
       indice_local_e(i,2) = arestas(k+1,1)
       indice_local_e(i,3) = arestas(k+2,1)
       k=k+3
    end do
  
    fronteira_e=0
   
    j=1
    do i=1,3*nelem
       if(arestas(i,1)==j)then
          edges(j,1:2)=arestas(i,2:3)
          if(fronteira(edges(j,1))==1 .and. fronteira(edges(j,2))==1)then
               fronteira_e(j)=1;

       !     9__________10
       !     |    15    | 
       !     |          | 
       !   16|          |18
       !     |          |  
       !     |__________|
       !    19    17    20
       ! 
          else if(fronteira(edges(j,1))==9 .and. fronteira(edges(j,2))==10)then
               fronteira_e(j)=15
          !else if(fronteira(edges(j,1))==10 .and. fronteira(edges(j,2))==9)then
          !     fronteira_e(j)=15 
          else if(fronteira(edges(j,1))==19 .and. fronteira(edges(j,2))==20)then
               fronteira_e(j)=17 
          !else if(fronteira(edges(j,1))==20 .and. fronteira(edges(j,2))==19)then
          !     fronteira_e(j)=17 
          else if(fronteira(edges(j,1))==9 .and. fronteira(edges(j,2))==19)then
              fronteira_e(j)=16 
          !else if(fronteira(edges(j,1))==19 .and. fronteira(edges(j,2))==9)then
          !    fronteira_e(j)=16 
          else if(fronteira(edges(j,1))==10 .and. fronteira(edges(j,2))==20)then
              fronteira_e(j)=18 
          !else if(fronteira(edges(j,1))==20 .and. fronteira(edges(j,2))==10)then
          !    fronteira_e(j)=18 
          end if

          j=j+1;
       end if

    end do
  

  end subroutine iee_edge

end module mod_edge_iee

