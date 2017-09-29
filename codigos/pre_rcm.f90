PROGRAM pre_rcm
!
!**************************************************************************************
!  Autor: modificado de Walleson Santos
!**************************************************************************************
IMPLICIT NONE
INTEGER(4):: i, j, n, n_elem
INTEGER(4), ALLOCATABLE:: nosElem(:,:), nos_1e0(:)
REAL(8), ALLOCATABLE:: coordXZ(:,:), resistividade(:)
REAL(8), ALLOCATABLE:: nos_2e0(:)
CHARACTER (LEN=20):: NOME
!
call getarg(1,nome)

!
!----------------Dados da malha

!write(*,*) 'Por favor, digite o prefixo'
!read(*,*) nome
open(10, file='.././malhas/' // trim(nome) // '_iee.ele', status='old', action='read')
open(20, file='.././malhas/' // trim(nome) // '_iee.edge', status='old', action='read')
open(30, file='.././malhas/' // trim(nome) //'_elements.txt', status='replace', action='write')
open(40, file='.././malhas/' // trim(nome) //'_edges.txt', status='replace', action='write')
!
read(10,*) n_elem !numero de elementos
read(20,*) n   !numero de nós

allocate(nosElem(n_elem,3), coordXZ(n,2), nos_1e0(n), nos_2e0(n), resistividade(n_elem))
i=0
DO WHILE (i .LT. n_elem)
   read(10,*) i, (nosElem(i,j),j=1,3), resistividade(i)
   write(30,*) (nosElem(i,j),j=1,3)!, resistividade(i)
ENDDO
i=0
DO WHILE (i .LT. n)
   read(20,*) i, (coordXZ(i,j),j=1,2), nos_1e0(i)
   write(40,*) (coordXZ(i,j),j=1,2), nos_1e0(i)
ENDDO
!
close(10)
close(20)
close(30)
close(40)
!
END PROGRAM pre_rcm
