#La idea es pasarle una matriz ampliada (R^(n*(n+1))). Lo que devuelva cada metodo debe ser un vector de R^n, solucion del sistema. Por esto, para todos los métodos se supone un sistema representado por una matriz cuadrada.

#A los iterativos se les pasa la cantidad de iteraciones. Despues sera una condicion de corte.

#Directos. Se recibe matriz diagonal ampliada. Resolución de Lower (diagonal inferior).

function x = despejarInferior(L)
	x = 0;
	n = rows(L);
	for i = [1:n]
		s = sum( L( i,[1: (i-1)] ) .* x([1:(i-1)]) );
		#Los b están incluidos en L. b(i) = L(i,end)
		x(i) = (L(i,n + 1) - s)/L(i,i);
	endfor
endfunction

function x = despejarSuperior(U)
	x = 0;
	n = rows(U);
	for i = [n: -1: 1]
		s = sum( U( i,[n:-1:i+1] ) .* x([n:-1:i+1]) );
		x(i) = (U(i,n + 1) - s)/U(i,i);
	endfor
endfunction

#Con pivote parcial.
function LU = factorizarLU(M)
	#Eliminación de Gauss, pero aprovecha los inferiores para guardar los multiplicadores.
	n = rows(M);
	for i = [2:n]
		huboCambios = false;
		#pivoteo parcial. Buscar forma de no intercambiar filas, sino guardar los cambios.
		
		
		#eliminación
		if ( M(i, i - 1) != 0 )
			huboCambios  = true;
			multiplicador  = M(i,i-1)/M(i-1,i-1);
			#guardo multiplicador en la inferior.
			M(i,i - 1) = multiplicador;
			M(i,[i:n]) -= M(i - 1,[i:n]) * multiplicador;
		endif
	endfor	
endfunction

L = [1,0,0,4;1,1,0,7;1,1,1,8];
despejarInferior(L);

U = [1,2,3,4;0,1,2,3;0,0,1,2];
despejarSuperior(U);
