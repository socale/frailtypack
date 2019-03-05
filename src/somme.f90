subroutine somme(ab,s,nboot,nbr,vOut)
  !$ use OMP_LIB
    implicit none
    integer,intent(in)::nbr,nboot
    double precision, dimension(2), intent(in):: ab
    double precision, intent(out)::s
    integer:: i,num_thread,num_proc,nb_thread,d
    double precision::c,a,b
    double precision, dimension(3), intent(inout)::vOut ! pour tester les parametre de sortie, vecteur allocatable
    
    !allocate(vOut(3))
    vOut(1)= 5.d0
    vOut(2)=10.d0
    vOut(3)=15.d0
    !deallocate(vout)
    a=ab(1)
    b=ab(2)
    s=a+b
    !open(2,file="testResult.txt")
    !write(2,*)"resultat de la somme=", s
    
    num_thread=0
    !call intpr("contenu de la variable contenant le nombre de thread avant:", -1, num_thread, 1)
    !$ call omp_set_num_threads(nbr)
    !$ num_thread = OMP_GET_NUM_THREADS()
    !$ num_proc = omp_get_num_procs()
    !call intpr("Nombre de processus diponibles sur la machine: ", -1, num_proc, 1)
    !call intpr("Nombre de threads utilises: ", -1, num_thread, 1)
    c=0.d0
    !$OMP PARALLEL DO shared (a,b,nb_thread,nboot) private (d,num_thread,i) reduction(+:c)
      do i=1,nboot
        !$ num_thread = omp_get_thread_num()
        !$ nb_thread = OMP_GET_NUM_THREADS() ! doit etre place dans la section parallele 
                                             ! sinon donnd automatiquement 1
        if(num_thread==0)then
          !d = d + 1
          !call intpr("Nombre de threads section parallele: ", -1, nb_thread, 1)
          !call intpr("nombre itteration: ", -1, d, 1)
          !call dblepr("i vaut: ", -1, dble(i), 1)
        endif
        !if(num_thread==0)then
        !  call dblepr("c avant: ", -1, c, 1)
        !endif
        c = c + 1
        !if(num_thread==0)then
        !  call dblepr("c apres: ", -1, c, 1)
        !endif
      end do
    !$OMP end parallel DO
    !call dblepr("la somme totale vaut", -1, c, 1) ! pour afficher un reel double precision
    
  !close(2)
endsubroutine somme