#Empezamos con un polinomio de ejemplo. Luego se hará una entrada más general.

function polinomio (coefs,x)
	res = 0.0
	i = 1
	for c = coefs
		res += c*(x^i)
		i++
	end
end

coefs = [-2,1]
polinomio(a,0)

#biseccion funciona cortando intervalos
