program griglia
    implicit none

    !integer:: dim !dimensione della griglia di spaziatura non costante
    integer:: i 

    real*8, parameter:: pc=3.086d18 !parsec (cgs)
    real*8, parameter:: max=3*pc*1d6 !dimensione della griglia -> 3Mpc
    real*8, parameter:: f=1.005d0 !fattore scala nella spaziatura dei punti della griglia
    real*8:: h !spaziatura NON COSTANTE tra i punti della griglia
    !eal*8:: test !distanza di controllo per la dimensione della griglia
    real*8, dimension(0:1005):: r, rr !griglie delle coordinate radiali

    !test=0.d0 !primo punto della griglia
    !dim=1 !inizializzazione della dimensione della griglia
    !space=h !prima spaziatura della griglia
    
    !do while (test<max)
        !print*, test/pc
    !    test=test+space
    !    space=f*space !incremento della spaziatura per l'eventuale prossimo punto di griglia
    !    dim=dim+1
    !enddo

    !print*, dim-1 !1006

    open(1, file="griglia.txt")
    h=100*pc
    r(0)=0.d0
    rr(0)=0.d0+h/2
    do i=1, 1005
        r(i)=r(i-1)+h
        h=f*h !vedi "griglia.f90"
        rr(i)=r(i)+h/2
        write(1,*) i, log10(r(i)/(pc*1d3)), log10(rr(i)/(pc*1d3))
    enddo
    close(1)

end program griglia