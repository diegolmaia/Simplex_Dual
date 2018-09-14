function [m_piv] = teste_razao(A, n_piv)
[m,n]=size(A);  %retorna tamanho de linha e colunas da matriz A
teste(1)=1000; % vetor que recebe o LD para análise de linha pivô
				% Valor alto para marcar LD do Z 

for i=2:m
	if (A(i,n_piv)~=0)
		teste(i)=A(i,n)/A(i,n_piv); % vetor que recebe o LD para análise de linha pivô 
	else
		teste(i)=10000;
	end
end

[menor, m_piv] = min(teste);
	