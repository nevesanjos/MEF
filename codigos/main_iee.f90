program edge_iee
use mod_edge_iee
implicit none  

  integer:: i, j
  character (len=20):: nome

  call getarg(1,nome)  

  open(unit=101, file= '.././malhas/' // trim(nome) // '.node', status='old', action='read')
  open(unit=202, file= '.././malhas/' // trim(nome) // '.ele' , status='old', action='read')
  open(unit=303, file= '.././malhas/' // trim(nome) // '.edge' , status='old', action='read')
 
  !========================================================================================

  ! LEITURA DOS ARQUIVOS
  !========================================================================================
  ! dimensões   
  read(101,*) nnos     
  read(202,*) nelem
  read(303,*) nedges

  ! Alocacao das variaveis de entrada
  allocate(coordenadas(nnos,2), fronteira(nnos))
  allocate(indice_local(nelem,3), resistividade(nelem))
  allocate(edges(nedges,2), fronteira_e(nedges))
  
  ! Leitura 
  do i=1, nnos
     read(101,*) j, coordenadas(i,1:2), fronteira(i)
  end do

  do i=1,nelem
     read(202,*) j, indice_local(i,1:3), resistividade(i)
  end do


  ! Apesar de ler, não utiliza estas informações
  ! A matriz edges e o vetor fronteira_e serão construidos na subrotina iee_edge
  ! Este procedimento pode ser otmizado! 
  do i=1,nedges
     read(303,*) j, edges(i,1:2), fronteira_e(i)
  end do

  close(101); close(202); close(303);
  !========================================================================================

  ! + Alocação 
  allocate(indice_local_e(nelem,3))
  
  
  ! Define o indice local das arestas
  !==========================================================================
  call iee_edge()
  !==========================================================================

  ! Imprime
  !==========================================================================
  ! Arquivos de saída
  open(unit=1201, file='.././malhas/' // trim(nome) // '_iee.ele', status='replace', action='write') 
  open(unit=1202, file='.././malhas/' // trim(nome) // '_iee.edge', status='replace', action='write') 

  ! Artifício encontrado para utilizar a rotina RCM do John Burkardt:
  ! Esta rotina tem como entrada a matriz local dos nós (.ele) e coordenadas dos nós (.node),
  ! porém para obter uma matriz com banda reduzida associada às arestas,
  ! forneceremos a matriz local e as "coordenadas" das arestas (nós pertencentes às mesmas)     

  write(1201,*) nelem 
  do i=1,nelem
     write(1201,*) i, indice_local_e(i,:), resistividade(i)
  end do

  write(1202,*) nedges 
  do i=1,nedges
     write(1202,*) i, edges(i,:), fronteira_e(i) 
  end do
  
end program edge_iee
