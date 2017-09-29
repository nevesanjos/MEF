PROGRAM pos_rcm
!
!**************************************************************************************
!  Autor: Walleson Santos
!**************************************************************************************
IMPLICIT NONE
INTEGER(4):: i, j, n, n_elem, label1_ele, label2_ele, label1_node, label2_node, label3_node
INTEGER(4), ALLOCATABLE:: nosElem(:,:), nos_1e0(:)
REAL(8), ALLOCATABLE:: coordXZ(:,:), resistividade(:)
REAL(8), ALLOCATABLE:: nos_2e0(:)
INTEGER(4):: banda, aux
CHARACTER (LEN=20):: NOME
!
!----------------Dados da malha

!write(*,*) 'Por favor, digite o prefixo'
!read(*,*) nome
call getarg(1,nome)
open(10, file='.././malhas/' // trim(nome) // '_iee.ele', status='old', action='read')
open(20, file='.././malhas/' // trim(nome) // '_iee.edge', status='old', action='read')
!
!----------------Dados da malha
read(10,*) n_elem!, label1_ele, label2_ele !numero de elementos
read(20,*) n!, label1_node, label2_node, label3_node   !numero de nós

allocate(nosElem(n_elem,3), coordXZ(n,2), nos_1e0(n), nos_2e0(n), resistividade(n_elem))
!
open(50, file='.././malhas/' // trim(nome) // '_iee_rcm_elements.txt', status='old', action='read')
open(60, file='.././malhas/' // trim(nome) // '_iee_rcm_edges.txt', status='old', action='read')
open(70, file='.././malhas/' // trim(nome) // '_iee_rcm.ele', status='replace', action='write')
open(80, file='.././malhas/' // trim(nome) // '_iee_rcm.edge', status='replace', action='write')
!

DO i=1,n_elem
   read(10,*) j, (nosElem(i,j),j=1,3), Resistividade(i)
ENDDO

DO i=1,n_elem
   read(50,*) (nosElem(i,j),j=1,3)
ENDDO

banda=0
DO i=1,n_elem
      aux=max(abs(nosElem(i,1)-nosElem(i,2)),abs(nosElem(i,1)-nosElem(i,3)),abs(nosElem(i,2)-nosElem(i,3)))
      if (aux>banda)then
         banda=aux
      end if
ENDDO

!
DO i=1,n
   read(60,*)(coordXZ(i,j),j=1,2), nos_2e0(i)
ENDDO
!
write(70,*) n_elem, banda !numero de elementos e a banda
write(80,*) n             !numero de nós
DO i=1,n_elem
   write(70,'(4i8,f20.3)') i, (nosElem(i,j),j=1,3), 1*resistividade(i)
ENDDO

nos_1e0=1*nos_2e0
DO i=1,n
   write(80,'(i5,2i10,i5)') i, (int(coordXZ(i,j)),j=1,2), (nos_1e0(i))
ENDDO

close(50)
close(60)
close(70)
close(80)

END PROGRAM pos_rcm
