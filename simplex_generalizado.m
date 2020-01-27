%Implementação do Método Simplex Generalizado
%  Recebe como parâmetros uma matriz A ∈ R^(m x n), vetores b ∈ R^m vetor
%  (coluna) e c ∈ R^n (Linha). Além de valores de limite inferior (l) e
%  limite superior (u) de cada variável - a sigla inf indica infinito.
function[X,Z,I,J1,J2,it,tipo]= simplex_generalizado(A,b,c,l,u)
%limpar as variáveis de saida
clear X;
clear Z;
clear I;
clear it;
clear tipo;
clear J1;
clear J2;
clear tamanho; 
clear m;
clear n;
clear j1_aux;
clear j2_aux;
clear delta;
clear x;

%valores iniciais vazios
flag_otmfinita1=[];
flag_candidate1=[];
flag_otmfinita2=[];
flag_candidate2=[];


X=[];
Z=[];
I=[];
it=0;

%busca o tamanho da matriz A
tamanho = size(A); % retorna o tamanho da matriz Amostras
m = tamanho(1); % n1 igual a quantidade de linhas
n = tamanho(2); % n2 igual a quantidade de colunas

%validar se e infactivel
if (size(A,2)~= size(c,2)) || (size(A,1)~= size(b,1)) || (rank([A b])> rank(A))
    tipo=-1;
    ind=-1;
    return;
end


% Multiplicamos por -1 as restrições em que b(i) < 0.
% Assim, garante-se que b >= 0.
for i = 1:m
    if b(i) < 0
        A(i, :)= A(i,:)*-1;
        b(i) = b(i)*-1;
     end
end

% Constrói o problema auxiliar
    A_aux = [A, eye(m)];
    c_aux = [zeros(n, 1); ones(m, 1)];
    x = [zeros(n, 1); b];
    i = (n + 1):(n + m);
    i_inv = eye(m);
j1_aux=(1:n);
j2_aux=[]; 
    
%atualiza os limites com os dados das variáveis artificiais
l=[l;zeros(m,1)];
u=[u;Inf(m,1)];

 % Aplica a Fase 1 no problema auxiliar - encontrar a abase inicial ótima
  [X,Z,I,J1,J2,it,tipo] =   simplex_body(A_aux, b, c_aux, m, n, x, i,j1_aux,j2_aux, it,l,u);
  % Elimina as variáveis artificiais do jogo
        X = X(1:n);
        l=l(1:n);
        u=u(1:n);
        J2(find(J2>n))=[];
        J1(find(J1>n))=[];
    if Z>0 %solução infactivel
        X=[];
        Z=[];
        I=[];
        tipo=-1;
        return;
    end
    
    %validar se o Z=0 então solucao otima do PA
    if Z==0      
        %validar situação de base degenerada
        remove=find(i==I); %se tem alguma variável artificial que está em I
        if isempty(remove)
        else %se tiver algum elemento para renover
            A(remove,:)=[];
            b(remove)=[];
            I(remove)=[];
            J2(remove)=[];
        end
        %executa o simples com os valores atualizados -2a. Fase do Simplex
       [X,Z,I,J1,J2,it,tipo]= simplex_body(A, b, c', m, n, X, I,J1,J2,it,l,u);
            return;
     end
end

function[X,Z,I,J1,J2,it,tipo]= simplex_body(A, b, c, m, n, x, i,j1_aux,j2_aux,it,l,u)
   
%inicilialização das variáveis para a primeira iteração
set=[1:length(c)]; %cria um vetor com as posições do vetor c
set(find(ismember(set,i)==1))=[]; %busca somente as variáveis basicas
N_set=set;
Ai=A(:,i); %cria a matriz de variáveis basicas
Ci=c(i);   %vetor C'das variáveis basicas

Cj1=c(j1_aux);  %vetor C' das variáveis não basicas J1 lower (l)
Cj2=c(j2_aux); %vetor C' das variáveis não basicas J2 Up (u)

I=i;
%atualiza as variáveis não basicas
J1=j1_aux;
J2=j2_aux;

xi=Ai\b; %calcula o valor de X das variávies básicas
Aj1=A(:,j1_aux); %cria a matriz de variáveis não basicas J1 (l)
Aj2=A(:,j2_aux); %cria a matriz de variáveis não báscas J2 (u)

X(i,1)=xi;%vetor X de acordo com as variávei1's básicas calculadas

xj=zeros(length(N_set),1); %preenche com zeros o vetor das variáveis não básicas
X(N_set,1)=xj;%vetor X de acordo com as variáveis não básicas calculadas

Z=c'*X;

while it>=0
    it=it+1; %atualiza o numero de interações

    pie=Ai'\Ci; %calcular pi = Ci*Ai(inversa)
    custoj1=pie'*Aj1-Cj1' ;%computa o custo das variâveis não basicas Ĉ em J1
    custoj2=pie'*Aj2-Cj2' ;%computa o custo das variâveis não basicas Ĉ em J2
    
    
    flag_otmfinita1=find(custoj1<=0);  %vetor dos custos negativos para J1
    flag_otmfinita2=find(custoj2>0);  %vetor dos custos positivos para J2
     
    flag_candidate1=find(custoj1>0); %vetor de custos positivos para J1
    flag_candidate2=find(custoj2<0); %vetor de custos negativos para J2
   
    if isempty(flag_candidate1)
       maiorindice=0; 
    else
        maiorindice=custoj1(find(custoj1==max(custoj1(flag_candidate1))));
    end
    if isempty(flag_candidate2)
        menorindice=0;
    else
        menorindice=custoj2(find(custoj2==min(custoj2(flag_candidate2))));
    end 
    
    
    poscandidate1=find(custoj1==maiorindice);    
    poscandidate2=find(custoj2==menorindice);
    
    indice1=j1_aux(poscandidate1);
    indice2=j2_aux(poscandidate2);
    
    flag_otmfinita1=j1_aux(flag_otmfinita1);
    flag_otmfinita2=j2_aux(flag_otmfinita2);
    
    
    if (isempty(flag_candidate1) & isempty(flag_candidate2) & (length(flag_otmfinita1)>0 |length(flag_otmfinita2)>0)) |Z==0
        tipo=0;% solução otima finita e unica
        return;
    end
    
    if isempty(flag_candidate1) & isempty(flag_candidate2) & max(custoj1(flag_otmfinita2)==0)
        tipo=1; % solução otima finita e multipla
        return;
    end
        
    %calcular o indice k de quem entra na base
    candidates1=[];
    candidates2=[];
    
    if isempty(maiorindice)
        maiorindice=0;
    end
    if isempty(menorindice)
        menorindice=0;
    end
    if maiorindice>=abs(menorindice)
        candidates1 =indice1;
        enter = candidates1(1); %em caso de empate o primeiro indice entra na base

    end 
    if abs(menorindice)>maiorindice
         candidates2=indice2;
         enter = candidates2(1); %em caso de empate o primeiro indice entra na base

    end
    
    %extrair o vetor Ak
    Ak=A(:,enter);
    
    if size(candidates1)>0 & isempty(candidates2)  %entra na base a partir do limite inferior l    
            %%%%%%%%calculo de quem sai da base a partir do limite inferior J1 (l)
            %definir quem sai da base r - calcular os gama1, gama2 e gama3       
            %retorna os indices positivos em Ak para gama1
            indp=find(Ak>0);
            if isempty(indp)
                leave_gama1=Inf;
                r_min1=Inf;
            else
                %calcular yk positvos
                yk1=Ai\Ak; 
                if size(yk1)>0
                    d= yk1.\(xi(indp)-l(I)); %calculo do vetor bi/yk

                   r=find(d>0); %retorna os indices positivos do vetor d
                   r_min1=min(d(r)); %busca o menor elemento em d de acordo com os indices positivos    
                    if isempty(r_min1)
                         leave_gama1=Inf;
                         r_min1=Inf;
                    else
                           leave_gama1 = i(find(d==r_min1));  %indice de quem sai da base indice candidato de acordo com gama1
                           leave_gama1(1);
                    end
                end
            end
            %calcula os indices negativos para gama2

            %retorna os indices negativos em Ak
            indn=find(Ak<0);
            if isempty(indn)
                leave_gama2=Inf;
                r_min2=Inf;
            else
                %calcular yk negativos
                yk2=Ai\Ak;
                if size(yk2)<0
                     d= -yk2.\(u(I) - xi(indn)); %calculo do vetor bi/yk
                     r=find(d>0); %retorna os indices positivos do vetor d
                     r_min2=min(d(r)); %busca o menor elemento em d de acordo com os indices positivos
                     if isempty(r_min2)
                         leave_gama2=Inf;
                         r_min2=Inf;
                     else
                        leave_gama2 = i(find(d==r_min2));  %indice de quem sai da base indice candidato de acordo com gama2
                        leave_gama2(1);
                     end
                end
            end
            if u(enter)==Inf | l(enter)==Inf
                leave_gama3=Inf;
                r_min3=Inf;
            else
                %gama 3
                r_min3= u(enter)-l(enter);
                leave_gama3=enter;
            end
            
            %valida solução ilimitada
            if leave_gama1==Inf & leave_gama2==Inf & leave_gama3==Inf
                    %solução ilimitada
                    tipo=2;
                    return;
            end
                       
            %seleciona o menor entre os r_min
            vetr_min=[r_min1;r_min2;r_min3];
            delta=min(vetr_min);
            indleave=find(vetr_min==delta); %indice para saber quem foi selecionado como menor
            indleave=indleave(1); %no caso de empate buscar o primeiro indice
            
            if indleave==1
               leave=leave_gama1;%sai da base o indice de r_min1     
            elseif indleave==2
                leave=leave_gama2; %sai da base o indice de rmin2    
            else
                leave=enter;%rmin3
            end
           indice=find(j1_aux==leave);
           j1_aux(indice)=[];
           j2_aux=[j2_aux;leave];
    end
    
    
     if size(candidates2)>0 & isempty(candidates1)  %entra na base a partir do limite superior u    
             %%%%%%%%calculo de quem sai da base a partir do limite superior J2 (u)

              %definir quem sai da base r - calcular os gama1, gama2 e gama3
           
            %retorna os indices negativos em Ak para gama1
            indp=find(Ak<0);
            if isempty(indp)
                leave_gama1=Inf;
                r_min1=Inf;
            else
                %calcular yk negativos
                yk1=Ai\Ak;
                if size(yk1)<0
                    d= -yk1.\(xi(indp)-l(I)); %calculo do vetor bi/yk

                   r=find(d<0); %retorna os indices positivos do vetor d
                   r_min1=min(d(r)); %busca o menor elemento em d de acordo com os indices positivos    
                   
                    if isempty(r_min1)
                         leave_gama1=Inf;
                         r_min1=Inf;
                    else
                   
                         leave_gama1 = i(find(d==r_min1));  %indice de quem sai da base indice candidato de acordo com gama1
                         leave_gama1(1);
                    end
                end
            end
            %calcula os indices positvos para gama2

            %retorna os indices postivos em Ak
            indn=find(Ak>0);
            if isempty(indn)
                leave_gama2=Inf;
                r_min2=Inf;
            else
                %calcular yk positvos
                yk2=Ai\Ak;
                if size(yk2)>0
                     d= yk2.\(u(I) - xi(indn)); %calculo do vetor bi/yk 
                     r=find(d>0); %retorna os indices positivos do vetor d
                     r_min2=min(d(r)); %busca o menor elemento em d de acordo com os indices positivos
                     if isempty(r_min2)
                         leave_gama2=Inf;
                         r_min2=Inf;
                     else
                         leave_gama2 = i(find(d==r_min2));  %indice de quem sai da base indice candidato de acordo com gama2
                         leave_gama2(1);
                     end
                end
            end
            if u(enter)==Inf | l(enter)==Inf
                leave_gama3=Inf;
                r_min3=Inf;
            else
                %gama 3
                r_min3= u(enter)-l(enter);
                leave_gama3=enter;
            end
            
            %valida solução ilimitada
            if leave_gama1==Inf & leave_gama2==Inf & leave_gama3==Inf
                    %solução ilimitada
                    tipo=2;
                    return;
            end
                       
            %seleciona o menor entre os r_min
            vetr_min=[r_min1;r_min2;r_min3];
            delta=min(vetr_min);
            indleave=find(vetr_min==delta); %indice para saber quem foi selecionado como menor
            indleave=indleave(1); %no caso de empate buscar o primeiro indice
            
            if indleave==1
               leave=leave_gama1;%sai da base o indice de r_min1     
            elseif indleave==2
                leave=leave_gama2; %sai da base o indice de rmin2    
            else
                leave=enter;%rmin3
            end
           indice=find(j1_aux==leave);
           j2_aux(indice)=[];
           j1_aux=[j1_aux;leave];
     end            
        
    %atualiza bases i e j para a próxima iteração
    indice_i=find(i==leave(1));
    i(indice_i)=enter;
    indice_j=find(N_set==enter);
    N_set(indice_j)=leave(1);
    I=i;
    %atualiza J1 e J2
   
    j1_aux(find(j1_aux==enter))=[];
    j2_aux(find(j2_aux==enter))=[];

    J1=j1_aux;
    J2=j2_aux;
     %atualiza as variáveis não basicas
    Cj1=c(j1_aux);  %vetor C' das variáveis não basicas J1 lower (l)
    Cj2=c(j2_aux); %vetor C' das variáveis não basicas J2 Up (u)

   
    Aj1=A(:,J1); %cria a matriz de variáveis não basicas J1 (l)
    Aj2=A(:,J2); %cria a matriz de variáveis não báscas J2 (u)

    %atualiza variáveis 
    Ai=A(:,i); %cria a matriz de variáveis basicas
    Ci=c(i);   %vetor C'das variáveis basicas
    Cj=c(N_set);  %vetor C' das variáveis não basicas
    
    %%%%%%%%%atualiza xi de acordo com o método
    yk=Ai\Ak;
    
    xi=(Ai\b)-(yk*delta); %calcula o valor de X das variávies básicas
    
    
    
    Aj=A(:,N_set); %cria a matriz de variáveis não basicas
    
    %Preencher o xj de acordo com o l
    %xj=zeros(length(N_set),1); %preenche com zeros o vetor das variáveis não básicas
    xj(J1)=l(J1);
    xj(J2)=l(J2);
    
    X=[];
    
    X(J1,1)=l(J1);
    X(J2,1)=l(J2);
    
    X(i,1)=xi;%vetor X de acordo com as variáveis básicas calculadas
    X(enter,1)=l(enter)+delta;
   % X(N_set,1)=xj;%vetor X de acordo com as variáveis não básicas calculadas
    Z=c'*X;

end   
end