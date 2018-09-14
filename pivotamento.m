function [A, vb_sai] = pivotamento(A,var_extra)
[m,n]=size(A);
j=2;
vb_sai(1,1)=0;

% Criar o vetor posicao das VB
for i=2:m
	vb_sai(j,1)=var_extra(i-1,1);
	j=j+1;
end

% Define a Coluna Pivô
var=A(1,1);

for i=2:n-1
	if (var > A(1,i))
		var=A(1,i);
		n_piv=i;
	end
end
A

%repete interacoes até coeficientes de Z >=0
while var < 0
	%var
	%n_piv
	%retorna a linha pivô
	m_piv = teste_razao(A, n_piv);
	
	coefi_entra_pos = n_piv ;		%pega posicao dos Coeficientes X que entra    
	vb_sai_pos = vb_sai(m_piv);  	%pega posicao da VB que sai
	
	%ld=A(m_piv,vb_sai_pos)
	%pivo=A(m_piv,n_piv)
	multiplicador = A(m_piv,vb_sai_pos)/A(m_piv,n_piv);  % calcula o multiplicador que transforma Coeficiente X que entra em VB que sai	
	
	A(m_piv,:)=A(m_piv,:).*multiplicador;	%atualiza linha pivô

	
%a3=A	



	for i=1:m
		linha=i;
		multiplicador_linha=-A(linha,n_piv)/A(m_piv,n_piv); %calcula multiplicador para igualar coluna do Coeficiente X q entra em Coluna do VB que sai  
		if linha~=m_piv
			for j=1:n
			
				A(linha,j)=A(linha,j) + multiplicador_linha*A(m_piv,j);
			end
		end
		
	
    end	

    A	
   
	vb_sai(m_piv)=n_piv;
	var=A(1,1);

    for i=1:n-1
        if (var >= A(1,i))
            var=A(1,i);
            n_piv=i;
        end
    end	
		
	%vb_sai
	%A
	%pause 
	
end

