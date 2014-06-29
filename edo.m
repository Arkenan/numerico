%Arjovsky-94881-Tarela-18/6/2014

clear();

%Espesor: 25,4 milimetros.
global rExt = 122.23750;
global rInt = 96.83750;

%padron: 94881. 
global TExt = 81.0;
global TInt = 948.0/6.0;

%Euler. f(y,t) es un function handler, h es el paso. t es un vector con la grilla de la variable independiente.
function y = Euler (f, h, y0,t)
  y(1) = y0; 
  n = columns(t);
  for i = [1:n-1]
    y(i+1) = y(i) + h* f(y(i),t(i));
  endfor
endfunction

%Se diferencia en que en lugar de recibir una funcion f a evaluar para la derivada recibe los valores de la derivada para la grilla.
function y = EulerDiscreto (f, h, y0,t)
  y(1) = y0; 
  n = columns(t);
  for i = [1:n-1]
    y(i+1) = y(i) + h* f(i);
  endfor
endfunction

%Esto se encarga de hacer los tiros y sacar una solucion aproximada.
function ret = tirarEuler (dvdr,paso,grilla)
  global TInt;
  global TExt;
  %Tiro1.
  vEuler1 = Euler(dvdr, paso, 0.0, grilla);
  TEuler1 = EulerDiscreto(vEuler1, paso, TInt, grilla);
  %Tiro 2.
  vEuler2 = Euler(dvdr, paso, 1.0, grilla);
  TEuler2 = EulerDiscreto(vEuler2, paso, 0.0, grilla);
  %Solucion aprox.
  ret = TEuler1 + (TExt-TEuler1(end))/(TEuler2(end))*TEuler2;
endfunction

%RK2
function y = Heun (f,h,y0,t)
  y(1) = y0;
  n = columns(t);
  for i = [1:n-1]
    yAux = y(i) + (2.0/3.0) * h * f(y(i),t(i));
    tAux = t(i) + h*(2.0/3.0);
    y(i+1) = y(i) + (h/4.0) * ( f(y(i),t(i)) + 3.0*f(yAux,tAux) );
  endfor
endfunction

%RK2. Como se usa para la segunda vuelta, es discreto (f es un vector).
function y = EulerMod (f,h,y0,t)
  y(1) = y0;
  n = columns(t);
  for i = [1:n-1]
    % No hace falta usar yAux = y(i) + h*f(i), ya que la f solo depende de la posicion. Ya esta calculada para todo punto en la grilla.
    y(i+1) = y(i) + (h/2)*(f(i)+f(i+1));
  endfor
endfunction

function ret = tirarRK(dvdr,paso,grilla)
  global TInt;
  global TExt;
  %Tiro1.
  vRK1 = Heun(dvdr, paso, 0.0, grilla);
  TRK1 = EulerMod(vRK1, paso, TInt, grilla);
  %Tiro2.
  vRK2 = Heun(dvdr, paso, 1.0, grilla);
  TRK2 = EulerMod(vRK2, paso, 0.0, grilla);
  %Solucion aprox.
  ret = TRK1 + (TExt-TRK1(end))/(TRK2(end))*TRK2;
endfunction;

function T = solAnalitica(r)
  global TExt;
  global TInt;
  global rExt;
  global rInt;
  T = (TInt-TExt)*(log(rExt/r))/(log(rExt/rInt))+TExt;
endfunction

%Funcion derivada de la velocidad de transmision de temperatura segun el radio.
dvdr = @(v,r) -v/r;

%Se inician diferente las iteraciones ya que las estimadas dejan de contar antes.
iteraciones = 0;
itEstimadas = 1;

%Cantidad de partes a dividir para paso inicial cercano a 1 milimetro.
partes = ceil(rExt-rInt);

pasos = [];
erroresE = [];
erroresRK = [];

do
  paso = (rExt - rInt)/partes;
  pasos = [pasos paso];
  
  %grilla discretizada de radios.
  r = [rInt:paso:rExt];
  grillas{iteraciones + 1} = r;
  
  %Euler.
  TEuler = tirarEuler(dvdr, paso, r);
  
  %Runge-Kutta. 
  TRK = tirarRK(dvdr,paso,r);
  
  solAnaliticaDiscreta = [0];
  for i = [1: columns(r)]
    solAnaliticaDiscreta(i) = solAnalitica(r(i));
  endfor
  
  errorEstimado = max(abs(TEuler - TRK));
  if (errorEstimado > 0.02)
    itEstimadas += 1;
  endif
  errorEuler = max(abs(TEuler - solAnaliticaDiscreta));
  erroresE = [erroresE errorEuler];
  errorRK = max(abs(TRK - solAnaliticaDiscreta));
  erroresRK = [erroresRK errorRK];
  errorMax = max([errorEuler, errorRK]);
  
  %Output.
  printf('\nPaso: %.5f\n',paso);
  printf ('Grilla     T(Euler)    T(RK)      T(Analitico) \n');
  for i = [1:columns(r)]
    printf('%.5f  %.5f   %.5f   %.5f\n', r(i),TEuler(i),TRK(i),solAnaliticaDiscreta(i));
  endfor
  printf('\nError(Euler) Error(RK) Error estimado\n');
  printf('%.5f       %.5f     %.5f\n\n',errorEuler,errorRK,errorEstimado);
  
  partes += 1;
  iteraciones += 1.0;
until (errorMax < 0.02)

pasos
erroresE
erroresRK

plot(pasos,erroresRK)
title('Error en Runge-Kutta segun el paso')
xlabel('Paso (h)')
ylabel('Error')
print('-dpng',"errRK.png");

plot(pasos,erroresE)
title('Error en Euler segun el paso')
xlabel('Paso (h)')
ylabel('Error')
print('-dpng',"errE.png");

%Estimativo del orden.
pasos = [];
erroresE = [];
erroresRK = [];
for j=[1:10]
  %Para el orden lo resuelvo otra vez pero reduciendo el paso a la mitad cada vez.
  paso = (rExt - rInt)/(2^j);
  
  %grilla discretizada de radios.
  r = [rInt:paso:rExt];
  
  %Euler. 
  TEuler = tirarEuler(dvdr,paso,r);  
  %Runge-Kutta. 
  TRK = tirarRK(dvdr,paso,r);
  
  solAnaliticaDiscreta = [0];
  for i = [1: columns(r)]
    solAnaliticaDiscreta(i) = solAnalitica(r(i));
  endfor
  
  pasos = [pasos paso];
  
  errorEuler = max(abs(TEuler - solAnaliticaDiscreta));
  erroresE = [erroresE errorEuler];
  errorRK = max(abs(TRK - solAnaliticaDiscreta));
  erroresRK = [erroresRK errorRK];
  errorMax = max([errorEuler, errorRK]);
endfor

printf('\nValores para estimacion del orden:\n')

erroresE
erroresRK
pasos