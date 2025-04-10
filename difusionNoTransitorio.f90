! CV.f90
! Conductividad variable
! 6-feb-2025
!
! Este programa resuelve la ecuacion de energía sin generación de calor, con volumenes finitos,
! en estado permanente

program calor1d
	implicit none
	real k, Tc, Th
	real Se, Sw, dx, dv
	real x0, xl,L
	integer i, j, nx
	real, allocatable :: T(:),x(:), ap(:), ae(:), aw(:), sp(:)

	!Inicialización de variables
nx = 25
Tc = 100.0; Th = 200.0		![°C]
L = 0.02			![m]
x0=0.0;xl=L
k = 0.5				![W/mK]
Se = 1.0; Sw = 1.0		!Unidad de área
dx = L/float(nx); dv = dx	!Operación de variales reales-->Resultado real!
q = 1000.0

!Alojamiento de memoria dinámica para los arreglos
allocate(aE(nx),aW(nx),aP(nx),Sp(nx),T(0:nx+1),x(0:nx+1))

T = 0.0 !Valor inicial de memoria
T(0) = Tc; T(nx+1)=Th

!#####################################################
!Solución numérica
!#####################################################

!Cálculo de coeficientes
do i=1,nx
	aE(i) = k*Se/dx
	aW(i) = k*Sw/dx
	aP(i) = aE(i)+aW(i)
	Sp(i) = q*dv
end do

!Condiciones de frontera
!*OESTE
aP(1) = aP(1) + aW(1)
Sp(1) = Sp(1) + 2.0*aW(1)*T(0)
aW(1) = 0.0

!*ESTE
aP(nx) = aP(nx) + aE(nx)
Sp(nx) = Sp(nx) + 2.0*aE(nx)*T(nx+1)
aE(nx) = 0.0

!https://matrixcalc.org/es/slu.html
!Método de Jacobi
do j=1,1000
	do i=1,nx
		T(i) = (aE(i)*T(i+1) + aW(i)*T(i-1) + Sp(i)) / aP(i)
	end do
end do

	!Imprime datos archivo .dat
	open(1,file='./CV.dat', status='replace')
	
	write(1,*) x0, T(0)
	do i=1, nx
		write(1,*)((float(i)-0.5)*dx+x0), T(i)
	enddo
	write(1,*) L, T(nx+1)
	
	!Imprime datos archivo .gp
	open(2,file='./CV.gp', status='replace')
	write(2,*) 'unset key'
	write(2,*)'set terminal png size 1024,768 enhanced font "Helvetica,20"'	
	write(2,*) 'set sty da l'
	write(2,*) 'set xlabel "x" '
	write(2,*) 'set ylabel "T" norotate'
	write(2,*) 'set style data points'
	write(2,*) 'set pointsize 1'
	write(2,*) 'Tc=', Tc
	write(2,*) 'Th=', Th
	write(2,*) 'k=', k
    write(2,*) 'L=', L
	write(2,*) 'a= (Th-Tc)/L'
	write(2,*) 'f(x)= Tc+a*x'
    write(2,*) 'set xrange [0:L]'
    write(2,*)'set output "output.png" '
	write(2,*) 'p "CV.dat" lc 3,f(x) with lines lt -1'


		
end program
