%montar simplex
function [p_otimo] = simplex(A,flag_sinal, B, C, flag_tipo)
	% The simplex algorithm for the LP problem

%                    min(max) z = c*x
%                 Subject to: Ax <= b
%                              x >= 0
%	matriz A -> coeficientes
%	matriz B -> lado direito  
%	matriz C ->	coeficientes de Z (equacao)


%	funcao z = 3x1 + 5x2
%	Maximizar Z
%	sujeito a 	x1        <=4
%					x2	  <=2
%				3 x1 + 2 x2	  <=18
%				x1,x2	  >=0
%
% 


%A= [1 0;	0 2; 3 2];

%B= [4;	12; 18];

%C= [-3 -5]; 

%flag_tipo = 1;    % 1 - Maximizar    
				  % 2 - Minimizar


if flag_tipo==1 
    [m_C, n_C] = size(C);
	[m n] = size(A);				%identifica o tamanho da matriz coef. de A
	A=[A eye(m)];					%insere as variáveis de folga
	A=[A B];						%Insere o Lado Direito (B)
	D=[C zeros(1, m+1)];			%Modela os coeficientes de Z
	A=[D;A];
	[m2,n2]=size(A);

	%Insere 'Z' na primeira linha de 'A '
	%vb=zeros[m,1];
	%coefi=zeros[1,n2-1];
    var_aux=zeros(m,2*m);			%cria matriz de var (excedente e artificial) nessa ordem

	linha = 1;
	coluna =1;

	% define as var (excedente e artificial) nessa ordem
	for i=1:m	
        if flag_sinal(i,1)==-1	% ta pronto
            var_aux(linha,coluna)=1;
			%var_extra(i)=coluna+n_C+1;
			coluna=coluna+1;
			%nada	
        end
        linha=linha+1;
    end	
	[m_var_aux, n_var_aux] = size(var_aux); %pega num de linhas e colunas da var_aux
    linha=1;
	coluna=1;
    
	for i=1:m_var_aux
			if var_aux(linha,coluna)==1 && var_aux(linha,coluna+1)==0;
				var_extra(i,1)=n_C+coluna;
				coluna=coluna+1;%soma-se mais 1 pq pula soh um numero	
            end
            linha=linha+1;
    end
	[A, vb_sai]=pivotamento(A,var_extra);
	
	A=[vb_sai A];
	
	var_otimo=n2-m2;
	j=1;
	for i=1:m2
		if ((A(i,1)<=var_otimo) && (A(i,1)~=0))
			p_otimo(j,1)=A(i,1);
			p_otimo(j,2)=A(i,n2+1);
			j=j+1;
		end
    end
    var_otimo=A(1,n2+1)
elseif flag_tipo==2 
	[m_C, n_C] = size(C);
	[m n] = size(A);				%identifica o tamanho da matriz coef. de A
	C=-C;							%minimizar
	M_grande=9999;					%M_grande pra Z
	var_aux=zeros(m,2*m);			%cria matriz de var (excedente e artificial) nessa ordem

	linha = 1;
	coluna =1;
	
	% define as var (excedente e artificial) nessa ordem
	for i=1:m			
		if flag_sinal(i,1)==1 % insere 1 var_excedente + 1 var_artificial
			var_aux(linha,coluna)=-1;
			var_aux(linha,coluna+1)=1;		%var_excedente é negativa
			%var_extra(i)=coluna+n_C;
			coluna=coluna+2;
			
		
		elseif flag_sinal(i,1)==0 % insere 1 var_artificial
			var_aux(linha,coluna)=1;
			%var_extra(i)=coluna+n_C+1;
			coluna=coluna+1;
		
		elseif flag_sinal(i,1)==-1	% ta pronto
            var_aux(linha,coluna)=1;
			%var_extra(i)=coluna+n_C+1;
			coluna=coluna+1;
			%nada		
		end
		linha=linha+1;
	end	
	
	%Incrementando variaveis artificial a Z
 	[m_var_aux, n_var_aux] = size(var_aux); %pega num de linhas e colunas da var_aux
	[m_C, n_C] = size(C);					% num de linhas  e colunas da matriz C ( coefi de Z)
	D=[C zeros(1, n_var_aux+1)];			%Modela os coeficientes de Z
	linha=1;
	coluna=1;
	for i=1:m_var_aux
			if var_aux(linha,coluna)==1 && var_aux(linha,coluna+1)==0;
				D(1,n_C+coluna)=M_grande;
				var_extra(i,1)=n_C+coluna;
				coluna=coluna+1;%soma-se mais 1 pq pula soh um numero	
			elseif var_aux(linha,coluna)==-1 && var_aux(linha,coluna+1)==1;
				D(1,n_C+coluna+1)=M_grande;
				var_extra(i,1)=n_C+coluna+1;
                coluna=coluna+2;%soma-se mais 2 pq pula dois numeros
				
			end
		linha=linha+1;
    end
   	
	%maximizar padrão
	%A=[A eye(m)];	
	A=[A var_aux];
	A=[A B];						%Insere o Lado Direito (B)
	%D=[C zeros(1, m+1)];			%Modela os coeficientes de 
   
	A=[D;A];						% insere os coeficientes de D em A
	[m2,n2]=size(A);
	
	%Insere 'Z' na primeira linha de 'A '
	%vb=zeros[m,1];
	%coefi=zeros[1,n2-1];

	%retirando variavel artificial de Z
	A_aux=A;
	for i=2:m2
		A_aux(i,:)=A_aux(i,:).*(-M_grande);
	end
	
	for i=2:m2
		for j=1:n2
			A_aux(1,j)=A_aux(1,j) + A_aux(i,j);		
		end
	end
	A(1,:)=A_aux(1,:);
	
	[A, vb_sai]=pivotamento(A,var_extra);
	

	A=[vb_sai A];

	var_otimo=A(1,n2+1)
	j=1;
	for i=1:m2
		if ((A(i,1)<=var_otimo) && (A(i,1)~=0))
			p_otimo(j,1)=A(i,1);
			p_otimo(j,2)=A(i,n2+1);
			j=j+1;
		end
    end	
end
