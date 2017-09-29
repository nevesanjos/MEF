program edge
use mod_edge
implicit none  

  integer:: i, j, k, ii, hh

  complex(dp), allocatable:: ex_plot(:), hy_plot(:), campos_tangentes(:,:), est_etha(:), est_zetha(:)
  real(dp), allocatable:: xx(:,:), rho(:,:), fase(:,:), zz(:,:), rho2(:), fase2(:)
  integer:: nfreq

  real(dp):: inicio, fim

  character (len=20):: nome

  call getarg(1,nome)  

  open(unit=101, file= '.././malhas/' // trim(nome) // '.node', status='old', action='read')
  open(unit=202, file= '.././malhas/' // trim(nome) // '.ele' , status='old', action='read')!
  open(unit=303, file= '.././malhas/' // trim(nome) // '_iee_rcm.edge', status='old', action='read')
  open(unit=404, file= '.././malhas/' // trim(nome) // '_iee_rcm.ele', status='old', action='read')


  ! LEITURA DO MODELO PRIMÁRIO
  !========================================================================================

  write(*,*) 'Digite a resistividade do meio encaixante: '
  read(*,*) resP

  write(*,*) 'Digite a quantidade de frequencias: '
  read(*,*) nfreq

  allocate(freq(nfreq+1))

  write(*,*) 'Digite a primeira frequência: '
  read(*,*) freq(1)

 call cpu_time(inicio) 

 allocate(mont_i(nfreq), mont_f(nfreq), sol_i(nfreq), sol_f(nfreq)) 

  ! LEITURA DOS ARQUIVOS
  !========================================================================================
  ! dimensões   
  read(101,*) nnos     
  read(202,*) nelem
  read(303,*) nedges
 
  write(*,*) ' '
  write(*,*) 'INFORMAÇÕES DA MALHA UTILIZADA: '
  write(*,*) ' '
  write(*,*) 'Número de elementos: ', nelem
  write(*,*) 'Número de nós: ', nnos
  write(*,*) 'Número de arestas: ', nedges

  read(404,*) nelem  ! releitura devido ao formato do arquivo

  ! Alocacao das variaveis de entrada
  allocate(coordenadas(nnos,2), fronteira(nnos))
  allocate(indice_local(nelem,3), indice_local_e(nelem,3), resistividade(nelem))
  allocate(edges(nedges,2), fronteira_e(nedges))
  
  ! Leitura do arquivos de entrada
  do i=1, nnos
     read(101,*) j, coordenadas(i,1:2), fronteira(i)
  end do

  do i=1,nelem
     read(202,*) j, indice_local(i,1:3), resistividade(i)
     read(404,*) j, indice_local_e(i,1:3)
  end do

  do i=1,nedges
     read(303,*) j, edges(i,1:2), fronteira_e(i)
  end do

  close(101); close(202); close(303); close(404)
  !========================================================================================
 
  j=0
  do i=1,nnos
     if(fronteira(i)==9)then
          j=j+1
     end if
  end do

  write(*,*) 'Quantidade de estações MT: ', j

  allocate(xx(j,2), zz(j,2), campos_tangentes(j,4), rho(nfreq,j), fase(nfreq,j), ordem(j,4))

  allocate(ex_plot(j), hy_plot(j), est_etha(j), est_zetha(j))


  j=1; k=1
  do i=1,nnos
     if(fronteira(i)==9)then
          xx(j,1)=coordenadas(i,1)
          zz(j,1)=coordenadas(i,2)
          j=j+1 
     end if
     if(fronteira(i)==20)then
          xx(k,2)=coordenadas(i,1)
          zz(k,2)=coordenadas(i,2)
          k=k+1 
     end if
  end do
    
    ! Calcula a quantidade de arestas das bordas  
    !==========================================================================
    numb=0
    do i=1,nedges
       if (fronteira_e(i)==1)then
            numb=numb+1
       end if
    end do


    write(*,*) 'Número de arestas na borda: ', numb    

    allocate(ind_b(numb))

    j=1
    do i=1,nedges
       if (fronteira_e(i)==1)then
            ind_b(j)=i
            j=j+1
       end if
     end do
    !========================================================================== 


  ! Calcula a semibanda da matriz global, pelo indice local das arestas
  !========================================================================== 
  call calc_banda() 
  !==========================================================================

  write(*,*) 'Posições de memória armazenada ([nbanda+1]*nedges):', nbanda+1, 'x',nedges,'=',(nbanda+1)*nedges
  write(*,*) '-----------------------------------------------------------------------------------------'
  write(*,*) ' '

  ! + Alocação
  allocate(M(nbanda+1,nedges), F(nedges))

  do ll=1,nfreq

     write(*,*) ll, '/', nfreq 


     ! MEF DE ARESTAS
     !==========================================================================
     call solver_edge()
     !==========================================================================

     j=1; k=1; ii=1; hh=1
     do i=1,nedges
        if(fronteira_e(i)==15)then
             campos_tangentes(j,1)=F(i)
             if(edges(i,1)>edges(i,2))  campos_tangentes(j,1) = - campos_tangentes(j,1)
             ordem(j,1) = (coordenadas(edges(i,1),1) + coordenadas(edges(i,2),1))/2.d0 
             j=j+1 
        else if(fronteira_e(i)==17)then
             campos_tangentes(k,2)=F(i)
             if(edges(i,1)>edges(i,2))  campos_tangentes(k,2) = - campos_tangentes(k,2)
             ordem(k,2) = (coordenadas(edges(i,1),1) + coordenadas(edges(i,2),1))/2.d0
             k=k+1
        else if(fronteira_e(i)==16)then
             campos_tangentes(ii,3)=F(i)
             if(edges(i,1)>edges(i,2))  campos_tangentes(ii,3) = - campos_tangentes(ii,3)
             ordem(ii,3) = (coordenadas(edges(i,1),1) + coordenadas(edges(i,2),1))/2.d0
             ii=ii+1
        else if(fronteira_e(i)==18)then
             campos_tangentes(hh,4)=F(i)
             if(edges(i,1)>edges(i,2))  campos_tangentes(hh,4) = - campos_tangentes(hh,4)
             ordem(hh,4) = (coordenadas(edges(i,1),1) + coordenadas(edges(i,2),1))/2.d0
             hh=hh+1
        end if
     end do

     call buble_sort(ordem(:,1), campos_tangentes(:,1), j-1 )
     call buble_sort(ordem(:,2), campos_tangentes(:,2), k-1 )
     call buble_sort(ordem(:,3), campos_tangentes(:,3), ii-1)
     call buble_sort(ordem(:,4), campos_tangentes(:,4), hh-1)

  
     hy_plot = (1.d0/zethaP) * ( ( (campos_tangentes(:,4) - campos_tangentes(:,3))/(xx(:,2)-xx(:,1))  &
                             -     (campos_tangentes(:,2) - campos_tangentes(:,1))/(zz(:,2)-zz(:,1))) &
                             +     (ci*kP)*exp(-ci*kP*(zz(:,1)+zz(:,2))/2.d0) )      ! Primario

     ex_plot = ( campos_tangentes(:,1) + campos_tangentes(:,2) ) / 2.d0 &
                + exp(-ci*kP*(zz(:,1)+zz(:,2))/2.d0)  ! Primário

  
     rho(ll,:)  = (1.d0/(2.d0*pi*freq(ll)*mu0))*(abs(ex_plot/hy_plot))**2
     fase(ll,:) = 180.d0*atan2(aimag(ex_plot/hy_plot),real(ex_plot/hy_plot))/pi  

     freq(ll+1) = freq(ll)*10**0.1
   
  end do 

  call cpu_time(fim)

  ! Imprime 1
  !==========================================================================
  ! Arquivos de saída
  open(unit=1201, file='.././saidas/sondagem_rho.dat', status='replace', action='write')  
  open(unit=1202, file='.././saidas/sondagem_fase.dat', status='replace', action='write')
  open(unit=1301, file='.././saidas/tempo.dat', status='replace', action='write')

  open(unit=1404, file='.././saidas/Ex2_est.dat', status='replace', action='write')
  open(unit=1505, file='.././saidas/Ez2_est.dat', status='replace', action='write')
  
  ! Formatos dos arquivos de saida 
  !-----------------------------------------
  !  n_est    x_1      x_2    ...   x_n
  ! freq_1  valor_1  valor_2  ... valor_n
  ! freq_2  valor_1  valor_2  ... valor_n
  !   .       .        .            .  
  !   .       .        .            .
  !   .       .        .            .
  ! freq_m  valor_1  valor_2  ... valor_n
  !-----------------------------------------

  write(1201,*) j-1, ((xx(:,1)+xx(:,2))/2000.d0) 
  write(1202,*) j-1, ((xx(:,1)+xx(:,2))/2000.d0) 
  write(1404,*) j-1, ((xx(:,1)+xx(:,2))/2000.d0)
  write(1505,*) j-1, ((xx(:,1)+xx(:,2))/2000.d0)
  do i=1,nfreq 
      write(1201,*) freq(i), rho(i,:)
      write(1202,*) freq(i), fase(i,:)
      write(1404,*) freq(i),  real( campos_tangentes(:,1) + campos_tangentes(:,2) )/2.d0
      write(1404,*) freq(i), aimag( campos_tangentes(:,1) + campos_tangentes(:,2) )/2.d0
      write(1505,*)  freq(i),  real( campos_tangentes(:,3) + campos_tangentes(:,4) )/2.d0
      write(1505,*)  freq(i), aimag( campos_tangentes(:,3) + campos_tangentes(:,4) )/2.d0

      write(1301,*) freq(i), mont_f(i)-mont_i(i), sol_f(i)-sol_i(i) 
  end do
  
write(*,*) 'tempo de processamento total: ', fim-inicio

end program edge
