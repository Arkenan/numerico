%Empezamos con un polinomio de ejemplo. Luego se hará una entrada más general.
function r = polinomio (coefs,x)
	res = 0.0;
	i = 0;
	for c = coefs;
		res += c*(x.^i);
		i++;
	endfor
	r = res;
endfunction

%Segun los coeficientes crea un polinomio listo para evaluar
function p = polEval(coeficientes)
	p = @(a) polinomio(coeficientes,a);
endfunction
	
%Funcion principal de biseccion. Dada una funcion de evaluación y un rango([inicio,fin]), devuelve una raiz.
function raiz = biseccion(f,rango)
	a = rango(1);
	b = rango(2);
	%comienzo con una cantidad arbitraria de iteraciones, luego hago algun criterio mejor.
	for i = [1:5]
		m = (a + b)/2;
		if (f(m)*f(a) < 0)
			b = m;
		else
			a = m;
		endif
	endfor
	raiz = m;
endfunction

function raiz = RegulaFalsi(f,rango)
	a = rango(1);
	b = rango(2);
	%comienzo con una cantidad arbitraria de iteraciones, luego hago algun criterio mejor.
	for i = [1:5]
		%esto es lo unico que cambia en RS. El punto medio se busca de otra forma.
		m = a - f(a)*((b-a)/(f(b)-f(a)));
		if (f(m)*f(a) < 0)
			b = m;
		else
			a = m;
		endif
	endfor
	raiz = m;
endfunction

%Se obtienen los coeficientes de un
function nc = coeficientesDerivados(coefs)
	%saco el primer coeficiente.
	res = coefs(2:length(coefs));
	if (!isempty(res))
		for (i = 1:length(res))
			res(i) *= i;
		endfor
	else
		res = [0];
	endif 
	nc = res;
endfunction

%Se pasan los coeficientes de un polinomio y se obtiene el polinomio derivado.
function pd = PolinomioDerivado(coefs)
	pd = @(x) polinomio(coeficientesDerivados(coefs),x);
endfunction

%
function pf = NewtonRaphson(f,df,seed)
	x = seed;
	for i  = 1:5
		x = x - f(x)/df(x);
	endfor
	pf = x;
endfunction

pol = polEval([-2,0,1]);
dpol = PolinomioDerivado([-2,0,1])

%fplot(pol,[-5,5])
%fplot(dpol,[-5,5])

r = biseccion(pol,[0,10])
r = RegulaFalsi(pol,[0,10])
r = NewtonRaphson(pol,dpol,10)