%Represento la matriz por 3 vectores. Uno para la diagonal central, otra para la primera supradiagonal, y otra para las E. Las inferiores salen por simetría, no hace falta guardarlas.

RTOL = 0.0001;
n = 9;
%Se muestran 3 decimales.
output_max_field_width(3);

function [D1,D2,D3] = construirDiagonales(n)
	if (n<6)
		error("La dimension de la matriz tiene que ser 6 o mas.");
	endif
	if (mod(n,3) != 0)
		error("La dimension de la matriz debe ser divisible por 3 para construirlas segun la consigna.");
	endif
	%Diagonal central.
	D1 = [[2,2,2],repmat([4,4,4],1,n/3-2),[3,3,3]]; 
	%Supra-Diagonal 1.
	D2 = [[-1/2,-1/2],repmat([0,-1,-1],1,n/3-2),[0,-1/2,-1/2]];
	%Supra-Diagonal 3.
	D3 =repmat([-1,-1,-1],1,n/3-1);	
endfunction

function M = construirMatriz(D1,D2,D3)
	n = length(D1);
	for i = [1:n]
		for j = [1:n]
			if (i == j)
				M(i,j) = D1(i);
			elseif ((j == i + 1) || (j == i - 1))
				M(i,j) = D2(min(i,j));
			elseif ((j == i + 3) || (j == i - 3))
				M(i,j) = D3(min(i,j));
			endif;
		endfor
	endfor			
endfunction
	
function s = suma(D2,D3, n, x)
	s(1) = D2(1)*x(2) + D3(1)*x(4);
	for i = [2,3]
		s(i) =  D2(i)*x(i+1) + D2(i-1)*x(i-1) + D3(i)*x(i+3);
	endfor
	for i = [4:n-3]
		s(i) = D2(i)*x(i+1) + D2(i-1)*x(i-1) + D3(i)*x(i+3) + D3(i-3)*x(i-3);
	endfor
	for i = [n-2,n-1]
		s(i) = D2(i)*x(i+1) + D2(i-1)*x(i-1) + D3(i-3)*x(i-3);
	endfor
	s(n) = D2(n-1)*x(n-1) + D3(n-3)*x(n-3);
endfunction

function [x, iteraciones, diferencias] = jacobi(D1,D2,D3,b,semilla,RTOL)
	n = length(D1);
	tolerable = false;
	x = semilla;
	iteraciones = 0;
	while (!tolerable)
		iteraciones++;
		%s es suma productos en cada fila.
		s  = suma(D2,D3,n,x);
		for i = [1:n]
			xNuevo(i,1) = (1/D1(i))*(b(i)-s(i));
		endfor;
		diferencia = max(abs(xNuevo-x));
		diferencias(iteraciones) = diferencia;
		tolerable =  diferencia/(max(abs(xNuevo))) < RTOL;
		x = xNuevo;
	endwhile
endfunction

function xNuevo = sumaGS(D1,D2,D3,b,x)
	n = length(D1);
	s = D2(1)*x(2) + D3(1)*x(4);
	xNuevo(1,1) = (1/D1(1))*(b(1)-s(1));
	for i = [2,3]
		s =  D2(i)*x(i+1) + D2(i-1)*xNuevo(i-1) + D3(i)*x(i+3);
		xNuevo(i,1) = (1/D1(i))*(b(i)-s);
	endfor
	for i = [4:n-3]
		s = D2(i)*x(i+1) + D2(i-1)*xNuevo(i-1) + D3(i)*x(i+3) + D3(i-3)*xNuevo(i-3);
		xNuevo(i,1) = (1/D1(i))*(b(i)-s);
	endfor
	for i = [n-2,n-1]
		s = D2(i)*x(i+1) + D2(i-1)*xNuevo(i-1) + D3(i-3)*xNuevo(i-3);
		xNuevo(i,1) = (1/D1(i))*(b(i)-s);
	endfor 
	s = D2(n-1)*xNuevo(n-1) + D3(n-3)*xNuevo(n-3);
	xNuevo(n,1) = (1/D1(n))*(b(n)-s);
endfunction

function [x, iteraciones, diferencias] = GS(D1,D2,D3,b,semilla,RTOL)
	n = length(D1);
	tolerable = false;
	x = semilla;
	iteraciones = 0;
	while (!tolerable)
		iteraciones++;
		xNuevo = sumaGS(D1,D2,D3,b,x);
		diferencia = max(abs(xNuevo-x));
		diferencias(iteraciones) = diferencia;
		tolerable =  diferencia/(max(abs(xNuevo))) < RTOL;
		x = xNuevo;
	endwhile
endfunction

function b = construirB(D1,D2,D3)
	n = length(D1);
	%Mi padrón es 94881.
	x = repmat([8;8;1],n/3,1);
	s = suma(D2,D3,n,x);
	for i = [1:n]
		b(i,1) = s(i) + D1(i)*x(i);
	endfor;
endfunction

%Construcciones previas
[D1,D2,D3]= construirDiagonales(n);
b = construirB(D1,D2,D3)
semilla = repmat(1,n,1);

%Resolución del sistema con cada método.
[xJacobi, iteracionesJacobi, diferenciasJacobi] = jacobi(D1,D2,D3,b,semilla,RTOL);
[xGS, iteracionesGS, diferenciasGS] = GS(D1,D2,D3,b,semilla,RTOL);
m = construirMatriz(D1,D2,D3);
xPC = m\b;

printf("Resultado por Jacobi:\n");
disp(xJacobi);

printf("Resultado por Gauss-Seidel:\n");
disp(xGS);

printf("Resultado de Octave:\n");
disp(xPC);

%Calculo de Radios espectrales.
puntosJacobi = [1:iteracionesJacobi];
rEspectralJacobi = e^(polyfit(puntosJacobi,log(diferenciasJacobi),1)(1))

puntosGS = [1:iteracionesGS];
rEspectralGS =  e^(polyfit(puntosGS,log(diferenciasGS),1)(1))

%Parte grafica.
figure(1)
plot(puntosJacobi, log(diferenciasJacobi),puntosGS,log(diferenciasGS));
title('diferencias en funcion de las iteraciones');
xlabel('Numero de Iteracion');
ylabel('Logaritmo de la norma infinito de la diferencia')
legend('Jacobi','Gauss-Seidel');
print -dpng 'grafico.png'
printf("FIN")