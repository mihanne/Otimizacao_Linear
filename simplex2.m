%Implementação do Método Simplex 2 fases
%  Recebe como parâmetros uma matriz A ∈ R^(m x n), vetores b ∈ R^m vetor (coluna) e c ∈ R^n (Linha)
function[X,Z,I,it,tipo]= simplex2(A,b,c)
%limpar as variáveis de saida
clear X;
clear Z;
clear I;
clear it;
clear tipo;
clear tamanho; 
clear m;
clear n;
clear BaseO;


%valores iniciais vazios
flag_otmfinita=[];
flag_candidate=[];
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

    % Aplica a Fase 1 no problema auxiliar
  

 [X,Z,I,it,tipo] =   simplex_body(A_aux, b, c_aux, m, n, x, i, it);
  % Elimina as variáveis artificiais do jogo
        X = X(1:n);
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
        end
        %executa o simples com os valores atualizados
       [X,Z,I,it,tipo]= simplex_body(A, b, c', m, n, X, I, it);
            return;
     end
end

%verificar%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[X,Z,I,it,tipo]= simplex_body(A, b, c, m, n, x, i,it)

    
%inicilialização das variáveis para a primeira iteração
set=[1:length(c)]; %cria um vetor com as posições do vetor c
set(find(ismember(set,i)==1))=[]; %busca somente as variáveis basicas
N_set=set;
Ai=A(:,i); %cria a matriz de variáveis basicas
Ci=c(i);   %vetor C'das variáveis basicas
Cj=c(N_set);  %vetor C' das variáveis não basicas
I=i;
xi=Ai\b; %calcula o valor de X das variávies básicas
Aj=A(:,N_set); %cria a matriz de variáveis não basicas
xj=zeros(length(N_set),1); %preenche com zeros o vetor das variáveis não básicas
X(i,1)=xi;%vetor X de acordo com as variávei1's básicas calculadas
X(N_set,1)=xj;%vetor X de acordo com as variáveis não básicas calculadas
Z=c'*X;

while it>=0
    it=it+1; %atualiza o numero de interações

    pie=Ai'\Ci; %calcular pi = Ci*Ai(inversa)
    custoj=pie'*Aj-Cj' ;%computa o custo das variâveis não basicas Ĉ
    flag_otmfinita=find(custoj<=0);  %vetor dos custos negativos
    flag_candidate=find(custoj>0); %vetor de custos positivos ou zero
    
    if isempty(flag_candidate) & length(flag_otmfinita)>0 & max(custoj(flag_otmfinita))<0
        tipo=0;% solução otima finita e unica
        return;
    end
    
    if isempty(flag_candidate) & max(custoj(flag_otmfinita))==0
        tipo=1; % solução otima finita e multipla
        return;
    end
        
    %calcular o indice k de quem entra na base
    
    indice=find(custoj==max(custoj(flag_candidate)));
    candidates =N_set(indice);
    enter = candidates(1); %em caso de empate o primeiro indice entra na base
    
    %extrair o vetor Ak
    Ak=A(:,enter);
    %calcular yk
    yk=Ai\Ak;
    %definir quem sai da base r
    d= yk.\xi; %calculo do vetor bi/yk
    
    if all(d <= 0)
       tipo=2; %solução ilimitada
       return;
    end

   r=find(d>0); %retorna os indices positivos do vetor d
   r_min=min(d(r)); %busca o menor elemento em d de acordo com os indices positivos
   leave = i(find(d==r_min));  %indice de quem sai da base
     
    %atualiza bases i e j para a próxima iteração

    indice_i=find(i==leave(1));
    i(indice_i)=enter;
    indice_j=find(N_set==enter);
    N_set(indice_j)=leave(1);
    I=i;       
    %atualiza variáveis 
    Ai=A(:,i); %cria a matriz de variáveis basicas
    Ci=c(i);   %vetor C'das variáveis basicas
    Cj=c(N_set);  %vetor C' das variáveis não basicas
    xi=Ai\b; %calcula o valor de X das variávies básicas
    Aj=A(:,N_set); %cria a matriz de variáveis não basicas
    xj=zeros(length(N_set),1); %preenche com zeros o vetor das variáveis não básicas
    X=[];
    X(i,1)=xi;%vetor X de acordo com as variáveis básicas calculadas
    X(N_set,1)=xj;%vetor X de acordo com as variáveis não básicas calculadas
    Z=c'*X;
end   
end