!DIFFUSIONE DEL FERRO NEL MEZZO INTRA-AMMASSO
!--------------------------------------------

module algoritmi
    implicit none
    
    contains
    subroutine eulero(old, new, rhs, h)
        real*8, intent(in):: rhs !rhs dell'equazione da risolvere
        real*8, intent(in):: h !passo spaziale di integrazione (costante o variabile)
        real*8, intent(in):: old !soluzione al passo precendente
        real*8, intent(out):: new !soluzione al passo successivo

        new=old+h*rhs

    end subroutine eulero

    subroutine ftcs(old, new, rhs, dt)
        real*8, intent(in):: rhs !rhs dell'equazione da risolvere
        real*8, intent(in):: dt !passo temporale di integrazione (costante o variabile)
        real*8, intent(in):: old !soluzione al tempo precendente
        real*8, intent(out):: new !soluzione al tempo successivo

        new=old+dt*rhs
    
    end subroutine ftcs

end module algoritmi

module equazioni
    implicit none
    
    contains
    subroutine eq_idro1(lhs, m, r) !equazione dell'equilibrio idrostatico nel caso di un gas perfetto isotermo
        real*8:: grav !campo gravitazionale
        real*8, intent(in):: r !raggio r, griglia per la massa  
        real*8, intent(in):: m !massa entro il raggio r
        real*8, intent(out):: lhs !lhs dell'equazione: derivata rispetto al raggio del logaritmo naturale della densità
        
        real*8, parameter:: k=1.38066d-16 !costante di Boltzmann (cgs)
        real*8, parameter:: T=8.9d7 !temperatura costante (Kelvin)
        real*8, parameter:: mu=0.61d0 !massa molecolare media del gas (assunta quella solare)
        real*8, parameter:: mh=1.67265d-24 !massa di un protone (cgs)
        real*8, parameter:: g=6.6720d-8 !costante di gravità (cgs)
        
        grav=-g*m/r**2
        lhs=grav*mu*mh/(k*T)

    end subroutine eq_idro1

    subroutine eq_idro2(lhs, m, t, cost, r, h) !equazione dell'equilibrio idrostatico nel caso di un gas perfetto con un profilo di temperatura
        real*8:: grav !campo gravitazionale
        real*8, intent(in):: h !passo d'integrazione
        real*8, intent(in):: r !raggio r  
        real*8, intent(in):: m !massa entro il raggio r
        real*8, intent(in):: t !temperatura al raggio r, ottenuta come media tra rr(i) e rr(i-1) 
        real*8, intent(in):: cost !differenza dei logaritmi naturali tra rr(i) e rr(i-1)
        real*8, intent(out):: lhs !lhs dell'equazione: derivata rispetto al raggio del logaritmo naturale della densità
        
        real*8, parameter:: k=1.38066d-16 !costante di Boltzmann (cgs)
        real*8, parameter:: mu=0.61d0 !massa molecolare media del gas (assunta quella solare)
        real*8, parameter:: mh=1.67265d-24 !massa di un protone (cgs)
        real*8, parameter:: g=6.6720d-8 !costante di gravità (cgs)
        
        grav=-g*m/r**2

        lhs=grav*mu*mh/(k*t)-cost/h

    end subroutine eq_idro2

    subroutine eq_idro3(lhs, m, t, cost, r, h) !equazione dell'equilibrio idrostatico nel caso di un gas perfetto con un profilo di temperatura e considerando la turbolenza
        real*8:: grav !campo gravitazionale
        real*8:: alpha !fattore dipendente dalla temperatura nell'equazione
        real*8, intent(in):: h !passo d'integrazione
        real*8, intent(in):: r !raggio r  
        real*8, intent(in):: m !massa entro il raggio r
        real*8, intent(in):: t !temperatura al raggio r, ottenuta come media tra rr(i) e rr(i-1) 
        real*8, intent(in):: cost !differenza dei logaritmi naturali tra rr(i) e rr(i-1)
        real*8, intent(out):: lhs !lhs dell'equazione: derivata rispetto al raggio del logaritmo naturale della densità
        
        real*8, parameter:: k=1.38066d-16 !costante di Boltzmann (cgs)
        real*8, parameter:: mu=0.61d0 !massa molecolare media del gas (assunta quella solare)
        real*8, parameter:: mh=1.67265d-24 !massa di un protone (cgs)
        real*8, parameter:: g=6.6720d-8 !costante di gravità (cgs)
        real*8, parameter:: vt=3d7 !velocità di riferimento della turbolenza, (cgs) cm/s 300 km/s
        
        alpha=1d0+mu*mh*vt**2/(k*t)
        grav=-g*m/r**2

        lhs=(1d0/alpha)*(grav*mu*mh/(k*t)-cost/h)

    end subroutine eq_idro3

    subroutine diff(lhs, div, source) !equazione di diffusione per la densità dei metalli nell'ICM
        real*8, intent(in):: div !divergenza del vettore flusso di metalli
        real*8, intent(in):: source !eventuale sorgente di metalli
        real*8, intent(out):: lhs !derivata rispetto al tempo della densità di metalli

        lhs=div+source
        
    end subroutine diff

end module equazioni

module funzioni
    implicit none
    
    contains
    subroutine NFW(rho, r, n) !dicretizzazione del profilo di densità di NFW
        integer, intent(in):: n !dimensione della lista
        integer:: i !indice di scorrimento della lista

        real*8, parameter:: pc=3.086d18 !pc nelle unità cgs
        real*8, parameter:: rho0=7.35d-26 !cgs, parametro di densità
        real*8, parameter:: rs=435.7*pc*10**3 !cgs, paramentro di distanza
        real*8, dimension(n), intent(in):: r
        real*8, dimension(n), intent(out):: rho 
        real*8:: x !distanza adimensionale
    
        do i=1, n
            x=r(i)/rs
            rho(i)=rho0/(x*(1+x)**2)
        enddo

    end subroutine NFW

    subroutine massa(m, rho, r, n) !calcolo della massa di una data distribuzione radiale entro un certo raggio
        integer, intent(in):: n
        integer:: i

        real*8:: guscio
        real*8, parameter:: pi=acos(-1.d0)
        real*8, dimension(n), intent(in) :: r 
        real*8, dimension(n), intent(in):: rho 
        real*8, dimension(n):: vol
        real*8, dimension(n), intent(out):: m !massa entro il raggio r

        do i=1, n-1 !calcolo del volume di ognuno dei n-1 gusci sferici
            vol(i)=4.d0/3.d0*pi*(r(i+1)**3-r(i)**3)
        enddo

        m(1)=0.d0
        do i=2, n 
            guscio=rho(i-1)*vol(i-1) !calcolo della massa del guscio sferico i-1
            m(i)=m(i-1)+guscio
        enddo
        
    end subroutine massa

    subroutine t_r(t, r, n)
        integer, intent(in):: n !dimensione della griglia dei raggi
        integer:: i !indice di scorrimento

        real*8, parameter:: pc=3.086d18 !pc nelle unità cgs
        real*8, parameter:: tcost=8.9d7 !temperatura scala (Kelvin)
        real*8, dimension(n), intent(in) :: r
        real*8, dimension(n), intent(out) :: t !profilo di temperatura
        real*8:: x !raggio adimensionale
        real*8:: fact1, fact2 !due fattori del profilo di temperatura
        
        do i=1, n
            x=r(i)/(1.4*pc*1d6)
            fact1=((x/0.045d0)**(1.9d0)+0.45d0)/((x/0.045d0)**(1.9d0)+1d0)
            fact2=1d0/(1d0+(x/0.6d0)**2)**(0.45d0)
            t(i)=1.35d0*tcost*fact1*fact2
        enddo

    end subroutine t_r

end module funzioni




program DIFFUSIONE_METALLI
    use algoritmi
    use equazioni
    use funzioni

    implicit none
    integer, parameter:: n=1006 !punti nella griglia, vedi "griglia.f90"
    integer:: i !indice di scorrimento nelle liste

    real*8, parameter:: e=2.71828d0 !numero di Nepero, usato per la conversione da logarirmo naturale a logaritmo in base 10
    real*8, parameter:: msun=1.989d33
    real*8, parameter:: pc=3.086d18 !pc nelle unità cgs
    real*8, parameter:: rmin=0d0, rmax=3*pc*10**6 !massima distanza dal centro dell'ammasso studiata, 3Mpc
    real*8, parameter:: zsol=1.8d-3 !abbondanza del ferro nel sole

    real*8:: h !passo d'integrazione sulla coordinata radiale, assunto non costante, vedi "griglia.f90"
    real*8:: lnrho0 !condizione iniziale per la densità del gas all'equilibrio idrostatico (NFW)
    real*8:: lnrho0_bcg !condizione iniziale per la densità del gas all'equilibrio idrostatico  (NFW+BCG)
    real*8:: lnrho0_t !condizione iniziale per la densità del gas all'equilibrio idrostatico  (NFW+BCG+T non costante)
    real*8:: lnrho0_turb !condizione iniziale per la densità del gas all'equilibrio idrostatico  (NFW+BCG+T non costante+turbolenza)
    real*8:: lhs !lhs delle equazioni
    !real*8:: fb !frazione di massa barionica che deve essere 0.16
    real*8:: dcost !coefficiente di diffusione costante
    real*8:: dt !passo temporale
    real*8:: time !variabile di controllo del tempo 
    real*8:: fi, fip1 !funzione f, vedi relazione nella posizione i e nella posizione i+1 
    real*8:: div !divergenza del flusso del ferro
    real*8:: source !sorgente di ferro
    real*8:: snu !quantificatore del contributo nell'arricchimento di ferro da parte delle SNIa 
    real*8:: a !parametro a nel profilo di Hernquist

    real*8, dimension(0:n-1):: r, rr !griglie delle coordinate radiali
    real*8, dimension(0:n-1):: rho !profilo di densità di NFW
    real*8, dimension(0:n-1):: bcg !prfilo di densità di Hernquist
    real*8, dimension(0:n-1):: rhogas !profilo di densità del gas (considerando solo NFW)
    real*8, dimension(0:n-1):: rhogas_bcg !profilo di densità del gas (considerando NFW+BCG)
    real*8, dimension(0:n-1):: rhogas_t !profilo di densità del gas (considerando NFW+BCG+T non costante)
    real*8, dimension(0:n-1):: rhogas_turb !profilo di densità del gas (considerando NFW+BCG+T non costante+turbolenza)
    real*8, dimension(0:n-1):: mdm !massa della Materia Oscura entro un raggio r
    real*8, dimension(0:n-1):: mbcg !massa della Bright Central Galaxy entro un raggio r, profilo di Hernquist (1990)
    real*8, dimension(0:n-1):: mgas !massa della gas entro un raggio r (considerando solo NFW)
    real*8, dimension(0:n-1):: mgas_bcg !massa della gas entro un raggio r (considerando NFW+BCG)
    real*8, dimension(0:n-1):: mgas_t !massa della gas entro un raggio r (considerando NFW+BCG+T non costante)
    real*8, dimension(0:n-1):: mgas_turb !massa della gas entro un raggio r (considerando NFW+BCG+T non costante+turbolenza)
    real*8, dimension(0:n-1):: lnrho !logaritmo NATURALE della densità del gas all'equilibrio idrostatico (NFW)
    real*8, dimension(0:n-1):: lnrho_bcg !logaritmo NATURALE della densità del gas all'equilibrio idrostatico (NFW+BCG)
    real*8, dimension(0:n-1):: lnrho_t !logaritmo NATURALE della densità del gas all'equilibrio idrostatico (NFW+BCG+T non costante)
    real*8, dimension(0:n-1):: lnrho_turb !logaritmo NATURALE della densità del gas all'equilibrio idrostatico (NFW+BCG+T non costante+turbolenza)
    real*8, dimension(0:n-1):: t !profilo radiale di temperatura all'interno dell'ammasso
   
    real*8, dimension(0:n-1):: zfe !profilo radiale dell'abbondanza di ferro osservata nell'ammasso di Perseo
    !real*8, dimension(0:n-1):: d !eventuale coefficente di diffusione dipendente dal tempo
    real*8, dimension(0:n-1):: rhofe !profilo radiale della densità dei ferro nell'ammasso
    real*8, dimension(0:n-1):: gradfe !gradiente del profilo di abbondanza di ferro 

!EQUILIBRIO IDROSTATICO IN UN AMMASSO DI GALASSIE
    !CREAZIONE DELLA GRIGLIE DELLE DISTANZE, assunta di SPAZIATURA COSTANTE
    !h=(rmax-rmin)/n
    !do i=0, n-1
    !    r(i)=rmin+i*h
    !    rr(i)=r(i)+h/2
    !enddo

    !CREAZIONE DELLA GRIGLIE DELLE DISTANZE, assunta di SPAZIATURA NON COSTANTE
    h=100*pc !prima distanza tra i punti della griglia
    r(0)=rmin
    rr(0)=rmin+h/2
    do i=1, n-1
        r(i)=r(i-1)+h
        h=1.005d0*h !vedi "griglia.f90"
        rr(i)=r(i)+h/2
    enddo

    !open(1, file="r.txt") !file con la griglia delle distanze in Kparsec
    !do i=0, n-1
    !    write(1, *) r(i)/(pc*10**3)
    !enddo
    !close(1)

    !CREAZIONE DEL PROFILO DI DENSITA' DI NFW DISCRETIZZATO
    call NFW(rho, rr, n)
    open(2, file="NFW.txt") !file con l distanze in log10(Kparsec) e la NFW in log10(cgs)
    do i=0, n-1 
        write(2, *) log10(rr(i)/(pc*10**3)), log10(rho(i))
    enddo
    close(2)

    !CREAZIONE DEL PROFILO DI MASSA DI NFW ENTRO UNA CERTA DISTANZA DAL CENTRO
    call massa(mdm, rho, r, n)
    open(3, file="massa_NFW.txt") !file con le distanze in log10(Kparsec) e la NFW log10(masse solari)
    do i=0, n-1 
        write(3, *) log10(r(i)/(pc*10**3)), log10(mdm(i)/msun)
    enddo
    close(3)

    !RISOLUZIONE DELL'EQUAZIONE DELL'EQUILIBRIO IDROSTATICO PER UN GAS PERFETTO A TEMPERTATURA COSTANTE
    !IL CAMPO GRAVITAZIONALE E' GENERATO SOLO DAL PROFILO NFW
    lnrho0=log(3.95d-26) !cgs, CONDIZIONE INIZIALE SCELTA VERIFICANDO CHE LA FRAZIONE DI MASSA BARIONICA SIA 0.16
    lnrho(0)=lnrho0
    do i=1, n-1
        call eq_idro1(lhs, mdm(i), r(i))
        call eulero(lnrho(i-1), lnrho(i), lhs, r(i)-r(i-1)) !PASSO NON COSTANTE
    enddo
    open(4, file="gas_NFW.txt") !file con le distanze in log10(Kparsec) e in log10(cgs)
    do i=0, n-1 
        write(4, *) log10(rr(i)/(pc*10**3)), lnrho(i)/log(10d0)
    enddo
    close(4)
    do i=0, n-1 !densità del gas in cgs
        rhogas(i)=e**(lnrho(i))
    enddo
    call massa(mgas, rhogas, r, n) !massa del gas entro il raggio r 
    open(5, file="massa_gas_NFW.txt") !file con le distanze in log10(Kparsec) e la massa di gas log10(masse solari), con solo NFW
    do i=0, n-1 
        write(5, *) log10(r(i)/(pc*10**3)), log10(mgas(i)/msun)
    enddo
    close(5)

    !CREAZIONE DEL PROFILO DI MASSA DELLA BRIGHT CENTRAL GALAXY ENTRO UNA CERTA DISTANZA DAL CENTRO
    do i=0, n-1
        mbcg(i)=(msun*1d12*r(i)**2)/(r(i)+(12d3*pc/(1+sqrt(2d0))))**2
    enddo
    open(6, file="massa_BCG.txt") !file con le distanze in log10(Kparsec) e la BCG log10(masse solari)
    do i=0, n-1 
        write(6, *) log10(r(i)/(pc*10**3)), log10(mbcg(i)/msun)
    enddo
    close(6)

    !RISOLUZIONE DELL'EQUAZIONE DELL'EQUILIBRIO IDROSTATICO PER UN GAS PERFETTO A TEMPERTATURA COSTANTE
    !IL CAMPO GRAVITAZIONALE E' GENERATO DAI PROFILI NFW ED HERNQUIST
    lnrho0_bcg=log(8d-26) !cgs, CONDIZIONE INIZIALE SCELTA VERIFICANDO CHE LA FRAZIONE DI MASSA BARIONICA SIA 0.16
    lnrho_bcg(0)=lnrho0_bcg
    do i=1, n-1
        call eq_idro1(lhs, mdm(i)+mbcg(i), r(i))
        call eulero(lnrho_bcg(i-1), lnrho_bcg(i), lhs, r(i)-r(i-1)) !PASSO NON COSTANTE
    enddo
    open(7, file="gas_NFW_BCG.txt") !file con le distanze in log10(Kparsec) e in log10(cgs)
    do i=0, n-1 
        write(7, *) log10(rr(i)/(pc*10**3)), lnrho_bcg(i)/log(10d0)
    enddo
    close(7)
    do i=0, n-1 !densità del gas in cgs
        rhogas_bcg(i)=e**(lnrho_bcg(i))
    enddo
    call massa(mgas_bcg, rhogas_bcg, r, n) !massa del gas entro il raggio r 
    open(8, file="massa_gas_NFW_BCG.txt") !file con le distanze in log10(Kparsec) e la massa di gas log10(masse solari), NFW+BCG
    do i=0, n-1 
        write(8, *) log10(r(i)/(pc*10**3)), log10(mgas_bcg(i)/msun)
    enddo
    close(8)

    !RISOLUZIONE DELL'EQUAZIONE DELL'EQUILIBRIO IDROSTATICO PER UN GAS PERFETTO A TEMPERTATURA NON COSTANTE
    !IL CAMPO GRAVITAZIONALE E' GENERATO DAI PROFILI NFW ED HERNQUIST
    call t_r(t, rr, n) !creazione del profilo di temperatura
    open(9, file="temp.txt") !file con il raggio in log10(raggio/Kpc) e log10(T/Kelvin)
    do i=0, n-1
        write(9, *) log10(rr(i)/(pc*1d3)), log10(t(i))
    enddo
    close(9)
    lnrho0_t=log(1.54d-25) !cgs, CONDIZIONE INIZIALE SCELTA VERIFICANDO CHE LA FRAZIONE DI MASSA BARIONICA SIA 0.16
    lnrho_t(0)=lnrho0_t
    do i=1, n-1
        call eq_idro2(lhs, mdm(i)+mbcg(i), (0.5d0)*(t(i)+t(i-1)), log(t(i))-log(t(i-1)), r(i), r(i)-r(i-1))
        call eulero(lnrho_t(i-1), lnrho_t(i), lhs, r(i)-r(i-1)) !PASSO NON COSTANTE
    enddo
    open(10, file="gas_temp.txt") !file con le distanze in log10(Kparsec) e le densità in log10(cgs)
    do i=0, n-1 
        write(10, *) log10(rr(i)/(pc*10**3)), lnrho_t(i)/log(10d0)
    enddo
    close(10)
    do i=0, n-1 !densità del gas in cgs
        rhogas_t(i)=e**(lnrho_t(i))
    enddo
    call massa(mgas_t, rhogas_t, r, n) !massa del gas entro il raggio r 
    open(11, file="massa_gas_temp.txt") !file con le distanze in log10(Kparsec) e la massa di gas log10(masse solari), NFW+BCG+T non costante
    do i=0, n-1 
        write(11, *) log10(r(i)/(pc*10**3)), log10(mgas_t(i)/msun)
    enddo
    close(11)

    !RISOLUZIONE DELL'EQUAZIONE DELL'EQUILIBRIO IDROSTATICO PER UN GAS PERFETTO A TEMPERTATURA NON COSTANTE E CON TURBOLENZA
    !IL CAMPO GRAVITAZIONALE E' GENERATO DAI PROFILI NFW ED HERNQUIST
    lnrho0_turb=log(9.15d-26) !cgs, CONDIZIONE INIZIALE VERIFICANDO CHE LA FRAZIONE DI MASSA BARIONICA SIA 0.16
    lnrho_turb(0)=lnrho0_turb
    do i=1, n-1
        call eq_idro3(lhs, mdm(i)+mbcg(i), (0.5d0)*(t(i)+t(i-1)), log(t(i))-log(t(i-1)), r(i), r(i)-r(i-1))
        call eulero(lnrho_turb(i-1), lnrho_turb(i), lhs, r(i)-r(i-1)) !PASSO NON COSTANTE
    enddo
    open(12, file="gas_turb.txt") !file con le distanze in log10(Kparsec) e le densità in log10(cgs)
    do i=0, n-1 
        write(12, *) log10(rr(i)/(pc*10**3)), lnrho_turb(i)/log(10d0)
    enddo
    close(12)
    do i=0, n-1 !densità del gas in cgs
        rhogas_turb(i)=e**(lnrho_turb(i))
    enddo
    call massa(mgas_turb, rhogas_turb, r, n) !massa del gas entro il raggio r 
    open(13, file="massa_gas_turb.txt") !file con le distanze in log10(Kparsec) e la massa di gas log10(masse solari), NFW+BCG+T non costante
    do i=0, n-1 
        write(13, *) log10(r(i)/(pc*10**3)), log10(mgas_turb(i)/msun)
    enddo
    close(13)

    !CONTROLLO DELLA FRAZIONE DI MASSA BARIONICA PER SETTARE LE CORRETTE DENSITA' CENTRALI INIZIALI NELL'INTEGRAZIONE
    !massa totale=massa(1467-1)
    !fb=mgas(n-1)/(mgas(n-1)+mdm(n-1)) !NFW e temperatura costante
    !fb=(mgas_bcg(n-1)+mbcg(n-1))/(mgas_bcg(n-1)+mbcg(n-1)+mdm(n-1)) !NFW+BCG e temperatura costante
    !fb=(mgas_t(n-1)+mbcg(n-1))/(mgas_t(n-1)+mbcg(n-1)+mdm(n-1)) !NFW+BCG e temperatura NON costante
    !fb=(mgas_turb(n-1)+mbcg(n-1))/(mgas_turb(n-1)+mbcg(n-1)+mdm(n-1)) !NFW+BCG e temperatura NON costante+turbolenza
    !print*, "massa barionica", fb
    

!EQUAZIONE DI DIFFUSIONE DELLA DENSITA' DI FERRO
    !LA DENSITA' DEL GAS E' ASSUNTA ESSERE QUELLA CALCOLATA CONSIDERANDO NFW+BCG ED IL PROFILO DI TEMPERATURA (rhogas_t)
    !CREAZIONE DEL PROFILO DI METALLICITA' OSSERVATO NELL'AMMASSO DI PERSEO 
    open(14, file="zobs.txt")
    do i=0, n-1
        zfe(i)=0.3d0*1.4d0*zsol*(2.2d0+(rr(i)/(pc*80d3))**3)/(1d0+(rr(i)/(pc*80d3))**3)
        rhofe(i)=(1d0/1.4d0)*zfe(i)*rhogas_t(i)
        write(14,*) log10(rr(i)/(pc*1d3)), zfe(i)/zsol
    enddo 
    close(14)

    dcost=1.32d29 !coefficiente di diffusione costante (cgs, velocità per spazio tipiche della turbolenza 10^2km/s e 10kpc)
    dt=0.5d0*(100*pc)**2/(2d0*dcost) !step temporale scelto minore del rapporto tra il quadrato del minimo step spaziale (100pc) e il minimo del coefficiente di diffusione (in questo caso costante)
    !print*, dt/(3.14d7), 5d9/(dt/3.14d7)
    !ASSUMO COME PROFILO INIZIALE QUELLO OSSERVATO NELL'AMMSSO DI PERSEO E TRASCURO LA SORGENTE DI METALLI
    time=0.d0
    do while (time.le.5d9*3.14d7) !integra per 5 miliardi di anni FARE CINQUE SALVATAGGI DEL PROFILO DI DENSITA' IN CORRISPONDENZA DEI PASSI DA 1Gyr!!!
        gradfe(0)=0.d0
        gradfe(n-1)=0.d0

        do i=1, n-2
            gradfe(i)=(zfe(i)-zfe(i-1))/(rr(i)-rr(i-1)) !gradiente del profilo d'abbondanza di ferro centrato in i
        enddo

        do i=1, n-2 !ciclo sulla distanza radiale per la costruzione del profilo radiale al tempo t+dt
            !if(i.eq.1) print*, r(i-1), r(i), rr(i-1), r(i)
            !print*, gradfe(i)
            fi=dcost*r(i)**2*(rhogas_t(i)+rhogas_t(i-1))/2d0*gradfe(i)
            fip1=dcost*r(i+1)**2*(rhogas_t(i+1)+rhogas_t(i))/2d0*gradfe(i+1)
            div=(1d0/1.4d0)*(fip1-fi)/((r(i+1)**3-r(i)**3)/3d0)
            !print*, log10(rr(i)/(), gradfe(i), fi, fip1, div
            !print*, div
            source=0.d0 !ALCUNA SORGENTE DI METALLI
            !lhs=div+source
            call diff(lhs, div, source) 
            call ftcs(rhofe(i), rhofe(i), lhs, dt) !costruzione del nuovo profilo di densità al tempo t+dt
            !rhofe(i)=rhofe(i)+dt*lhs
        enddo

        do i=1, n-2 !creazione del dato iniziale per il passo temporale successivo
            zfe(i)=1.4d0*rhofe(i)/rhogas_t(i)
        enddo
        !condizioni di outflow
        rhofe(0)=rhofe(1)
        rhofe(n-1)=rhofe(n-2)
        zfe(0)=zfe(1)
        zfe(n-1)=zfe(n-2)

        time=time+dt !faccio un passo temporale
    enddo 
    
    open(15, file="diff.txt")
    do i=2, n-1
        write(15,*) log10(rr(i)/(pc*1d3)), zfe(i)/zsol
    enddo
    close(15)

    !ASSUMO L'ASSENZA DI UN PROFILO INIZIALE DI FERRO E L'ASSENZA DI DIFFUSIONE
    !CONSIDERO UNA SORGENTE DI FERRO DOVUTA ALLE SNIa E AI VENTI STELLARI
    dcost=1.32d29 !coefficiente di diffusione costante (cgs, velocità per spazio tipiche della turbolenza 10^2km/s e 10kpc)
    rhofe=0.d0
    snu=0.15d0 !

    !creazione del profilo di densità di Hernquist
    a=12d3*pc/(1d0+sqrt(2d0))   
    do i=0, n-1
        bcg(i)=1d12*msun/(2d0*3.14)*a/(rr(i)*(rr(i)+a)**3)
    enddo

    time=0.d0
    do while (time.le.1d9*3.14d7) !integra per 5 miliardi di anni FARE CINQUE SALVATAGGI DEL PROFILO DI DENSITA' IN CORRISPONDENZA DEI PASSI DA 1Gyr!!!

        do i=0, n-1 !ciclo sulla distanza radiale per la costruzione del profilo radiale al tempo t+dt
            div=0d0 !assenza di diffusione
            source=bcg(i)*(6d-23+snu*3.13d-21)
            call diff(lhs, div, source) 
            call ftcs(rhofe(i), rhofe(i), lhs, dt) !costruzione del nuovo profilo di densità al tempo t+dt
        enddo

        do i=0, n-1 !creazione del dato iniziale per il passo temporale successivo
            zfe(i)=1.4d0*rhofe(i)/rhogas_t(i)
        enddo
        
        !condizioni di outflow
        rhofe(0)=rhofe(1)
        rhofe(n-1)=rhofe(n-2)
        zfe(0)=zfe(1)
        zfe(n-1)=zfe(n-2)

        time=time+dt !faccio un passo temporale
    enddo

    open(16, file="source.txt")
    do i=2, n-1
        write(16,*) log10(rr(i)/(pc*1d3)), zfe(i)/zsol, log10(bcg(i)) !prova del profilo di Hernquist
    enddo
    close(16)

    !ASSUMO L'ASSENZA DI UN PROFILO INIZIALE DI FERRO E LA PRESENZA DELLA DIFFUSIONE
    !CONSIDERO UNA SORGENTE DI FERRO DOVUTA ALLE SNIa E AI VENTI STELLARI
    dcost=1.32d29 !coefficiente di diffusione costante (cgs, velocità per spazio tipiche della turbolenza 10^2km/s e 10kpc)
    zfe=0.d0
    rhofe=0.d0
    snu=0.15d0 !
    time=0.d0
    do while (time.le.5d9*3.14d7) !integra per 5 miliardi di anni FARE CINQUE SALVATAGGI DEL PROFILO DI DENSITA' IN CORRISPONDENZA DEI PASSI DA 1Gyr!!!
        gradfe(0)=0.d0
        gradfe(n-1)=0.d0

        do i=1, n-2
            gradfe(i)=(zfe(i)-zfe(i-1))/(rr(i)-rr(i-1)) !gradiente del profilo d'abbondanza di ferro centrato in i
        enddo

        do i=1, n-2 !ciclo sulla distanza radiale per la costruzione del profilo radiale al tempo t+dt
            fi=dcost*r(i)**2*(rhogas_t(i)+rhogas_t(i-1))/2d0*gradfe(i)
            fip1=dcost*r(i+1)**2*(rhogas_t(i+1)+rhogas_t(i))/2d0*gradfe(i+1)
            div=(1d0/1.4d0)*(fip1-fi)/((r(i+1)**3-r(i)**3)/3d0)
            source=bcg(i)*(6d-23+snu*3.13d-21)
            call diff(lhs, div, source) 
            call ftcs(rhofe(i), rhofe(i), lhs, dt) !costruzione del nuovo profilo di densità al tempo t+dt
        enddo

        do i=1, n-2 !creazione del dato iniziale per il passo temporale successivo
            zfe(i)=1.4d0*rhofe(i)/rhogas_t(i)
        enddo
        !condizioni di outflow
        rhofe(0)=rhofe(1)
        rhofe(n-1)=rhofe(n-2)
        zfe(0)=zfe(1)
        zfe(n-1)=zfe(n-2)

        time=time+dt !faccio un passo temporale
    enddo

    open(17, file="source_diff.txt")
    do i=2, n-1
        write(17,*) log10(rr(i)/(pc*1d3)), zfe(i)/zsol
    enddo
    close(17)
end program DIFFUSIONE_METALLI