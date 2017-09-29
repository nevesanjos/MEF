#!/bin/bash

echo $1 
echo ola $USER

cd codigos

echo $1
entrada=$1
if [ $entrada=='' ];
then
echo Por favor, informe o prefixo da malha.
read entrada
fi

echo $1
echo -n "Deseja reenumerar a malha? (1=sim, 2=nao) " 

read sn 

if [ "$sn" -eq "1" ]; 


then 
  
  gfortran sub_iee.f90 main_iee.f90 -o iee.x    
  ./iee.x $entrada                          # Constroi o indice local das arestas

  gfortran pre_rcm.f90 -o pre.x
  ./pre.x $entrada                          # Gera a entrada para o rcm 

  gfortran triangulation_rcm.f90 -o rcm.x
  ./rcm.x $entrada

  gfortran pos_rcm.f90 -o pos.x
  ./pos.x $entrada                         # Gera a entrada para o codigo principal  

fi  

  gfortran sub_edge.f90 main_edge.f90 -o mt.x
  ./mt.x $entrada                          # Porgrama principal  

  rm *.mod*
  
  cd ../malhas

  rm *.txt*

echo Obrigado!
