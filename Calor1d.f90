! CV_Jacobi.f90 - Solución numérica con Jacobi y comparación con solución analítica
! Conducción de calor en estado estacionario sin generación
! 13-Feb-2025

program calor1d_jacobi
    implicit none
    integer, parameter :: max_iter = 10000  ! Máximo de iteraciones
    real, parameter :: tol = 1.0e-6  ! Tolerancia de convergencia

    ! Declaración de variables
    real k, Tc, Th, dx, L, error, q
    real Se, Sw
    integer i, nx, iter
    real, allocatable :: T(:), T_old(:), x(:), aP(:), aE(:), aW(:), Sp(:), b(:), T_analytical(:)

    ! Inicialización de variables
    nx = 25
    Tc = 100.0; Th = 500.0  ! Temperaturas en °C
    L = 0.5                 ! Longitud de la barra en m
    k = 1000.0              ! Conductividad térmica en W/mK
    Se = 0.012; Sw = 0.012  ! Área transversal
    dx = L / float(nx-1)    ! Espaciado de la malla
    q = 1000.0e3            ! Generación de calor en W/m^3 (1000 kW/m^3)

    ! Asignación de memoria
    allocate(aE(nx), aW(nx), aP(nx), Sp(nx), b(nx), T(0:nx+1), T_old(0:nx+1), x(0:nx+1), T_analytical(0:nx+1))

    ! Inicialización de temperatura
    T = 0.0
    T(0) = Tc
    T(nx+1) = Th

    ! Definir nodos de la malla
    do i = 0, nx+1
        x(i) = i * dx
    end do

    ! #############################
    ! Solución analítica
    ! #############################
    do i = 0, nx+1
        T_analytical(i) = Tc + (Th - Tc) / L * x(i) - q / (2.0 * k) * x(i) * (L - x(i))
    end do

    ! #############################
    ! Cálculo de coeficientes
    ! #############################
    do i = 1, nx
        aE(i) = k * Se / dx
        aW(i) = k * Sw / dx
        aP(i) = aE(i) + aW(i)
        Sp(i) = 0.0
        b(i) = 0.0
    end do

    ! Condiciones de frontera
    aP(1) = aP(1) + aW(1)
    b(1) = b(1) + 2.0 * aW(1) * Tc
    aW(1) = 0.0

    aP(nx) = aP(nx) + aE(nx)
    b(nx) = b(nx) + 2.0 * aE(nx) * Th
    aE(nx) = 0.0

    ! #############################
    ! Método de Jacobi
    ! #############################
    iter = 0
    error = 1.0
    do while (error > tol .and. iter < max_iter)
        iter = iter + 1
        T_old = T
        error = 0.0

        do i = 1, nx
            T(i) = (aE(i) * T_old(i+1) + aW(i) * T_old(i-1) + b(i)) / aP(i)
            error = max(error, abs(T(i) - T_old(i)))
        end do
    end do

    ! Guardar resultados para Gnuplot
    open(1, file="solucion.dat", status="replace")
do i = 1, nx
    write(1,*) x(i), T_numerica(i), T_analitica(i)
end do
close(1)


end program calor1d_jacobi
