%Arjovsky-94881-Tarela-18/6/2014

clear();

%Radio exterior: 244,475 milimetros. Radio interior: 219,075 milimetros. Espesor: 25,4 milimetros.
global rExt = 244.47500;
global rInt = 219.07500;

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

function T = solAnalitica(r)
  global TExt;
  global TInt;
  global rExt;
  global rInt;
  T = (TInt-TExt)*(log(rExt/r))/(log(rExt/rInt))+TExt;
endfunction

function dTdr = obtenerV(paso,v)
  global rInt;
  indice = (r- rInt)/paso
  dTdr = @(T,r) v(indice);
endfunction

%funcion derivada de la velocidad de transmision de temperatura segun el radio.
dvdr = @(v,r) -v/r;

%Se inician diferente las iteraciones ya que las estimadas dejan de contar antes.
iteraciones = 0;
itEstimadas = 1;

%Paso inicial: 1 milimetro.
paso = (rExt-rInt)/3;
pasos = [];
erroresE = [];
erroresRK = [];

do
  
  %grilla discretizada de radios.
  r = [rInt:paso:rExt];
  grillas{iteraciones + 1} = r;
  %Euler. 

  %Tiro1.
  vEuler1 = Euler(dvdr, paso, 0.0, r);
  TEuler1 = EulerDiscreto(vEuler1, paso, TInt, r);
  %Tiro 2.
  vEuler2 = Euler(dvdr, paso, 1.0, r);
  TEuler2 = EulerDiscreto(vEuler2, paso, 0.0, r);
  %Solucion aprox.
  TEuler = TEuler1 + (TExt-TEuler1(end))/(TEuler2(end))*TEuler2;
  TEs{iteraciones + 1} = TEuler;
  %Runge-Kutta. 

  %Tiro1.
  vRK1 = Heun(dvdr, paso, 0.0, r);
  TRK1 = EulerMod(vRK1, paso, TInt, r);
  %Tiro2.
  vRK2 = Heun(dvdr, paso, 1.0, r);
  TRK2 = EulerMod(vRK2, paso, 0.0, r);
  %Solucion aprox.
  TRK = TRK1 + (TExt-TRK1(end))/(TRK2(end))*TRK2;
  TRKs{iteraciones + 1} = TRK;
  solAnaliticaDiscreta = [0];
  for i = [1: columns(r)]
    solAnaliticaDiscreta(i) = solAnalitica(r(i));
  endfor
  
  iteraciones += 1.0;
  pasos = [pasos paso];
  printf('\nPaso: %.5f\n',paso)
  paso = (rExt - rInt)/(iteraciones+4);
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
  
  printf ('Grilla     T(Euler)    T(RK)      T(Analitico) \n')
  for i = [1:columns(r)]
    printf('%.5f  %.5f   %.5f   %.5f\n', r(i),TEuler(i),TRK(i),solAnaliticaDiscreta(i));
  endfor
  printf('\nError(Euler) Error(RK) Error estimado\n');
  printf('%.5f       %.5f     %.5f\n\n',errorEuler,errorRK,errorEstimado);
until (errorMax < 0.02)

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

%Para el orden lo resuelvo otra vez pero reduciendo el paso a la mitad cada vez.
paso = (rExt-rInt)/2;

pasos = [];
erroresE = [];
erroresRK = [];
iteraciones = 0;
for j=[1:10]
  %grilla discretizada de radios.
  r = [rInt:paso:rExt];
  
  %Euler. 

  %Tiro1.
  vEuler1 = Euler(dvdr, paso, 0.0, r);
  TEuler1 = EulerDiscreto(vEuler1, paso, TInt, r);
  %Tiro 2.
  vEuler2 = Euler(dvdr, paso, 1.0, r);
  TEuler2 = EulerDiscreto(vEuler2, paso, 0.0, r);
  %Solucion aprox.
  TEuler = TEuler1 + (TExt-TEuler1(end))/(TEuler2(end))*TEuler2;
  %Runge-Kutta. 

  %Tiro1.
  vRK1 = Heun(dvdr, paso, 0.0, r);
  TRK1 = EulerMod(vRK1, paso, TInt, r);
  %Tiro2.
  vRK2 = Heun(dvdr, paso, 1.0, r);
  TRK2 = EulerMod(vRK2, paso, 0.0, r);
  %Solucion aprox.
  TRK = TRK1 + (TExt-TRK1(end))/(TRK2(end))*TRK2;
  solAnaliticaDiscreta = [0];
  for i = [1: columns(r)]
    solAnaliticaDiscreta(i) = solAnalitica(r(i));
  endfor
  
  iteraciones += 1.0;
  
  pasos = [pasos paso];
  paso = (rExt - rInt)/(2^(iteraciones));
  
  errorEuler = max(abs(TEuler - solAnaliticaDiscreta));
  erroresE = [erroresE errorEuler];
  errorRK = max(abs(TRK - solAnaliticaDiscreta));
  erroresRK = [erroresRK errorRK];
  errorMax = max([errorEuler, errorRK]);
endfor

erroresE
erroresRK
pasos

function str = csv(vec)
  str = num2str(vec(1));
  for i = [2: columns(vec)]
    str = [str,' ', num2str(vec(i))];
  endfor
endfunction 