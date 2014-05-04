%La idea es pasarle una matriz ampliada (R^(n*(n+1))). Lo que devuelva cada metodo debe ser un vector de R^n, solucion del sistema. Por esto, para todos los métodos se supone un sistema representado por una matriz cuadrada.

%Directos. Se recibe matriz diagonal ampliada. Resolución de Lower (diagonal inferior).

function x = despejarInferior(L)
	x = 0;
	n = rows(L);
	for i = [1:n]
		s = sum( L( i,[1: (i-1)] ) .* x([1:(i-1)]) );
		%Los b están incluidos en L. b(i) = L(i,end)
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

%Con pivote parcial.
function LU = factorizarLU(M)
	n = rows(M);
	allCols = [1:n+1];
	i = 1;
	ld = false;
	
	while (i <= n) && (!ld)
		%Pivoteo parcial. OPT: Buscar forma de no intercambiar filas, sino guardar los cambios.
		yMax = i; vMax = M(i,i);
		
		for j = [i,n]
			if ( M(j,i) < yMax )
				yMax = j;
				vMax = M(j,i);
			endif
		endfor
		
		%OPT: Habría que revisar de que la igualdad sea estricta o no.
		ld = (vMax == 0);
		
		if (!ld)
			%Intercambio de filas i,j por pivoteo.
			aux = M(j,allCols);
			M(j,allCols) = M(i,allCols);
			M(i,allCols) = aux;
			
			%eliminación
			for j = [i + 1: n]
				if M(j,i) != 0
					multiplicador  = M(j, i)/M(i,i);
					%guardo multiplicador en la inferior.
					M(j,i) = multiplicador;
					M(j,[i + 1: n + 1]) -= M(i,[i + 1: n + 1]) * multiplicador;
				endif
			endfor
		endif
		
		i++;
	endwhile	
	if (ld) 
		printf("es linealmente dependiente")
	endif
	LU = M;
endfunction

function x = Gauss(M)
	 LU = factorizarLU(M);
	 x = despejarSuperior(LU);
endfunction

L = [1,0,0,4;1,1,0,7;1,1,1,8];
despejarInferior(L);

U = [1,2,3,4;0,1,2,3;0,0,1,2];
despejarSuperior(U);

M = [1,1,1,1;1,1,1,1;1,1,1,1];
LU = factorizarLU(M)