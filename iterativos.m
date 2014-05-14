clear;
RTOL = 0.0001;
%listaDimensiones = [6,18,24,30];
listaDimensiones = [6,18,24,30];
%Se muestran 3 decimales.
output_max_field_width(4);

%Represento la matriz por 3 vectores. Uno para la diagonal central, otra para la primera supradiagonal, y otra para las E. Las inferiores salen por simetría, no hace falta guardarlas.

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

function [x, xk] = jacobi(D1,D2,D3,b,semilla,RTOL)
	n = length(D1);
	tolerable = false;
	x = semilla;
	k = 0;
	while (!tolerable)
		k++;
		%s es suma productos en cada fila.
		s  = suma(D2,D3,n,x);
		for i = [1:n]
			xNuevo(i,1) = (1/D1(i))*(b(i)-s(i));
		endfor;
		xk{k} = xNuevo;
		tolerable =  max(abs(xNuevo-x))/(max(abs(xNuevo))) < RTOL;
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

function [x, xk] = GS(D1,D2,D3,b,semilla,RTOL)
	n = length(D1);
	tolerable = false;
	x = semilla;
	k = 0;
	while (!tolerable)
		k++;
		xNuevo = sumaGS(D1,D2,D3,b,x);
		xk{k} = xNuevo;
		tolerable =  max(abs(xNuevo-x))/(max(abs(xNuevo))) < RTOL;
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

exp = 0;
for n = listaDimensiones
	exp++;
	printf("**** DIMENSION %i ****\n\n",n);
	%Construcciones previas
	printf("Se obtiene el vector x a partir del padron y b a partir de este\n");
	x = repmat([8;8;1],n/3,1);
	
	[D1,D2,D3]= construirDiagonales(n);
	
	b = construirB(D1,D2,D3);
	printf ("x               b\n");
	for indice = [1:n]
		printf("%i            %.4f\n", x(indice),b(indice));
	endfor
	
	semilla = repmat(1,n,1);
	
	%Resolución del sistema con cada método.
	tic
		[xJacobi, xkJacobi] = jacobi(D1,D2,D3,b,semilla,RTOL);
	tJacobi(exp) = toc;
	tic
		[xGS, xkGS] = GS(D1,D2,D3,b,semilla,RTOL);
	tGS(exp) = toc;
	
	m = construirMatriz(D1,D2,D3);
		xPC = m\b;
	printf("\nResultados para cada metodo\n");
	printf("Jacobi            GS            Octave\n");
	for indice = [1:n]
		printf("%.4f          %.4f          %.4f\n",xJacobi(indice),xGS(indice),xPC(indice));
	endfor
	
	printf("\nIteraciones\n  %i              %i",length(xkJacobi),length(xkGS));
	printf("\nTiempo\n%.4f          %.4f\n",tJacobi(exp),tGS(exp));
	
	
	%Calculo de Radios espectrales.
	
	puntosJacobi = [1:length(xkJacobi)];
	for i = puntosJacobi
		diferenciasJacobi(i) = log(max(abs(xkJacobi{i}-x)));
	endfor
	rEspectralJacobi(exp) = e^(polyfit(puntosJacobi,diferenciasJacobi,1)(1));
	
	printf("\n\nEl radio espectral de Jacobi para la dimension %i es %.3f.\n",n,rEspectralJacobi(exp));
	puntosGS = [1:length(xkGS)];
	for i = puntosGS
		diferenciasGS(i) = log(max(abs(xkGS{i}-x)));
	endfor
	rEspectralGS(exp) =  e^(polyfit(puntosGS,diferenciasGS,1)(1));
	printf("El radio espectral de Gauss-Seidel para la dimension %i es %.3f.\n\n",n,rEspectralGS(exp));
	
	%Parte grafica.
	figure(1);
	plot(puntosJacobi, diferenciasJacobi,puntosGS,diferenciasGS);
	titulo = ["n = " num2str(n)];
	title(titulo);
	set(gca,'XTick',[[1:5:50],[length(puntosJacobi) length(puntosGS)]]);
	xlabel('Numero de Iteracion');
	xlim([0 50]);
	ylabel('Logaritmo de la norma infinito de la diferencia');
	legend('Jacobi','Gauss-Seidel');
	fName = sprintf("grafico%i.png",n);
	print ('-dpng', fName);

endfor;

%Grafico de radios espectrales.
plot(listaDimensiones,rEspectralJacobi,listaDimensiones,rEspectralGS);
titulo = ["n = " strrep(num2str(listaDimensiones),'   ',',')];
title(titulo);
xlabel('Dimension de la matriz del sistema');
ylabel('Radios espectrales para cada metodo');
legend('Jacobi','Gauss-Seidel');
print('-dpng',"rEspectral.png");

%Grafico de tiempos

plot(listaDimensiones,tJacobi,listaDimensiones,tGS);
titulo = ["n = " strrep(num2str(listaDimensiones),'   ',',')];
title(titulo);
xlabel('Dimension de la matriz del sistema');
ylabel('tiempo tardado por cada metodo.');
legend('Jacobi','Gauss-Seidel');
print('-dpng',"tiempos.png");

printf("FIN")