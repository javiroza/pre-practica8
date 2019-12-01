! ---------------------------------- Pre-pràctica 8 ------------------------------------- !
! Autor: Javier Rozalén Sarmiento
! Grup: B1B
! Data: 03/12/2019
!
! Funcionalitat: es resol l'equació de Scdxrödinger independent del temps per trobar els
! autovalors i autovectors d'una partícula en una caixa 1D de mida L emprant el mètode 
! artiller amb Runge-Kutta 3
!
! Comentaris: els apartats 1, 2 i 3 es resumeixen en les subrutines RLSTN3,derivades i 
! artiller, respectivament.
!
! Nota: el programa principal i les subrutines tenen indentacions tals que, en alguns
! editors de text (SublimeText), es poden plegar i desplegar per facilitar la lectura.

program pre_practica8
    implicit none
    double precision V,E,E1,E2,dx
    double precision integral
    double precision, allocatable :: vectphi(:)
    integer i,j,N,nequs
    common/cts/V,E

    V=-1.2d0
    nequs=2

    ! -------------------------------- Apartat 4 --------------------------------------- !
    N=640 ! Nombre de passos per Ralston, i també dimensions de l'autovector
    dx=1.d0/dble(N)
    allocate(vectphi(N))
    open(11,file="P8-1920-res.dat")
    open(12,file="aux.dat")

    ! ---------------- Primer estat (n=1, estat fonamental) ---------------- !
    E1=0.1d0
    E2=0.09d0
    call artiller(E1,E2,nequs,N,vectphi)
    call write(11)
    call trapezoids(0.d0,1.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(12,*) i*dx,vectphi(i)
    enddo
    call write(12)

    ! ----------------------- Segon estat (n=2) ---------------------------- !
    E1=16.d0
    E2=15.d0
    call artiller(E1,E2,nequs,N,vectphi)
    call write(11)
    call trapezoids(0.d0,1.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(12,*) i*dx,vectphi(i)
    enddo
    call write(12)

    ! ----------------------- Tecer estat (n=3) ---------------------------- !
    E1=45.d0
    E2=44.d0
    call artiller(E1,E2,nequs,N,vectphi)
    call write(11)
    call trapezoids(0.d0,1.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(12,*) i*dx,vectphi(i)
    enddo
    call write(12)

    ! ----------------------- Quart estat (n=4) ---------------------------- !
    E1=65.d0
    E2=64.d0
    call artiller(E1,E2,nequs,N,vectphi)
    call write(11)
    call trapezoids(0.d0,1.d0,N,abs(vectphi)**2.d0,integral)
    vectphi=vectphi/dsqrt(integral)
    ! Escrivim el vector propi en un arxiu
    do i=1,N
        write(12,*) i*dx,vectphi(i)
    enddo
    call write(12)

    close(11)
    close(12)
end program pre_practica8

! Subrutina RLSTN3 --> Calcula un pas del mètode de Ralston de tercer ordre per un sistema
! d'n equacions de primer ordre acoblades
subroutine RLSTN3(x,dx,nequs,yyin,yyout)
    ! x --> Punt a partir del qual es fa el pas
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt anterior
    ! yyout --> Vector amb les dades del punt següent
    implicit none
    integer nequs,i
    double precision x,dx,yyin(nequs),yyout(nequs)
    double precision k1(nequs),k2(nequs),k3(nequs)
    double precision dyout(nequs) ! Vector mut necessari per cridar la subrutina derivades

    do i=1,nequs 
        ! Càlcul dels vectors k1,k2,k3 
        call derivades(nequs,x,yyin,dyout)
        k1(i)=dyout(i)
        call derivades(nequs,x,yyin+dx/2.d0*k1,dyout)
        k2(i)=dyout(i)
        call derivades(nequs,x,yyin+(3.d0*dx/4.d0)*k2,dyout)
        k3(i)=dyout(i)

        ! Càlcul del vector yyout 
        yyout(i)=yyin(i)+dx/9.d0*(2.d0*k1(i)+3.d0*k2(i)+4.d0*k3(i))
    enddo

    return 
end subroutine RLSTN3

! Subrutina derivades --> Calcula la derivada de la funció (vectorial) a trobar
subroutine derivades(nequ,x,yin,dyout)
    ! nequ --> nombre d'equacions
    ! x --> Punt on es calcula la derivada
    ! yin --> Vector amb les dades del punt on es calcula la derivada
    ! dyout --> Vector amb les derivades
    implicit none
    integer nequ
    double precision x,yin(nequ),dyout(nequ)
    double precision V,E
    common/cts/V,E

    dyout(1)=yin(2)
    dyout(2)=2.d0*(V-E)*yin(1)

    return
end subroutine derivades

! Subrutina artiller --> Resol una equació diferencial amb condicions de contorn pel mètode de tir
subroutine artiller(E1,E2,nequs,npassos,vectphi)
    ! nequs,npassos --> nombre d'equacions del sistema a resoldre i nombre de passos pel solver
    ! E1,E2 --> Valors inicials aleatoris per inicialitzar el mètode de la secant
    ! yyin,yyout --> vectors amb les conidions inicials i resultats finals
    ! vectphi --> vector amb la solució a l'equació 
    implicit none
    integer i,j
    double precision E1,E2,yyin(nequs),yyout(nequs)
    integer nequs,npassos,N
    double precision phiE1,phiE2,phiE3,vectphi(npassos),E,E3,V,dx,x
    common/cts/V,E

    N=npassos
    dx=1.d0/dble(N)

    do j=1,10000 
        ! Trobem una primera solució per E=E1
        yyin(1)=0.d0
        yyin(2)=0.25d0
        E=E1
        do i=1,N
            call RLSTN3(x,dx,nequs,yyin,yyout)
            vectphi(i)=yyout(1)
            ! Reescrivim variables
            yyin=yyout
        enddo
        phiE1=vectphi(N)

        ! Trobem una segona solució per E=E2
        yyin(1)=0.d0
        yyin(2)=0.25d0
        E=E2
        do i=1,N
            call RLSTN3(x,dx,nequs,yyin,yyout)
            vectphi(i)=yyout(1)
            ! Reescrivim variables
            yyin=yyout
        enddo
        phiE2=vectphi(N)

        ! Amb una nova E calculada amb el mètode de la secant, trobem una solució 
        ! aproximada
        yyin(1)=0.d0
        yyin(2)=0.25d0
        E3=(E1*phiE2-E2*phiE1)/(phiE2-phiE1)
        E=E3
        do i=1,N
            call RLSTN3(x,dx,nequs,yyin,yyout)
            vectphi(i)=yyout(1)
            ! Reescrivim variables
            yyin=yyout
        enddo
        phiE3=vectphi(N)
        write(11,*) E3,phiE3

        if (abs(phiE3).lt.0.00005d0) then
            print*,E3
            exit
        else
            E1=E2
            E2=E3
        endif 
    enddo
    return
end subroutine artiller

! Subrutina trapezoids --> Calcula una integral 1-D per trapezis
subroutine trapezoids(x1,x2,ndim,funci,integral)
    ! x1,x2 --> Punts inicial i Final 
    ! ndim --> nombre de dimensions del vector funció
    ! funci --> vector funció (conté totes les imatges)
    ! integral --> pues eso
    implicit none
    double precision x1,x2,integral,funci(ndim)
    integer i,ndim

    integral = 0.d0
    do i=1,ndim-1
        integral = integral + funci(i)
    enddo
    integral = (integral + funci(1)/2.d0 + funci(ndim)/2.d0)*((x2-x1)/dble(ndim))

    return
end subroutine trapezoids

! Subrutina write --> Escriu dues línies en blanc en un arxiu
subroutine write(arxiu)
    implicit none
    integer arxiu
    write(arxiu,*) ""
    write(arxiu,*) ""
    return
end subroutine