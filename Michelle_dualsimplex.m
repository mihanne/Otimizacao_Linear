function[X,Z,I,it,tipo]= dualsimplex(A,b,c,i)
%Implementação dual simplex
%A igual a matriz
%b igual ao vetor coluna
%c igual ao vetor linha de custos
%i igual ao vetor linha da base inicial

%limpar as variáveis de saida
clear X;
clear Z;
clear BaseO;
clear it;
clear tipo;

%valores iniciais vazios
X=[];
Z=[];
I=[];
it=0;

%busca o tamanho da matriz A
tamanho = size(A); % retorna o tamanho da matriz Amostras
m = tamanho(1); % n1 igual a quantidade de linhas
n = tamanho(2); % n2 igual a quantidade de colunas


%para min e negativos para max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%validar se e infactivel
if (size(A,2)~= size(c,2)) || (size(A,1)~= size(b,1)) | (rank([A b])> rank(A))
    tipo=-1;
    return;
end
%Monta o problema primal
%inicilialização das variáveis para a primeira iteração
set=[1:length(c)]; %cria um vetor com as posições do vetor c
set(find(ismember(set,i)==1))=[]; %busca somente as variáveis basicas
N_set=set;
I=i;

%%%%%%%% Aplicar o problema dual na tabela primal
c_PP=c.*-1;
A_PP=A.*-1;
b_PP=b.*-1;

Ai=A_PP(:,i); %cria a matriz de variáveis basicas
Ci=c_PP(i);   %vetor C'das variáveis basicas
Cj=c_PP(N_set);  %vetor C' das variáveis não basicas


xi=Ai\b_PP; %calcula o valor de X das variávies básicas
Aj=A_PP(:,N_set); %cria a matriz de variáveis não basicas
xj=zeros(length(N_set),1); %preenche com zeros o vetor das variáveis não básicas
X(i,1)=xi;%vetor X de acordo com as variávei1's básicas calculadas
X(N_set,1)=xj;%vetor X de acordo com as variáveis não básicas calculadas

Z=c*X;

while it>=0
    it=it+1; %atualiza o numero de interações
    %Encontrar quem sai da base
    negrc = [b_PP(b_PP<0)];   %calcula custo negativo
    if isempty(negrc)
         tipo=0; % solução otima finita e unica
         return;
    end
    if min(b_PP)==0
        tipo=1; % solução otima finita e multipla
    end
    val_sair = min(negrc);   %Valor de b que sairá da base o menor custo negativo
    pos_sair=find(b_PP==val_sair); %posicao da linha que sai da tabela
    base_sair= I(pos_sair);   %indice que sai da base
    
    %encontrar quem entra na base
    calc_entrar=Cj./Aj(pos_sair,:);
    val_entrar=min(calc_entrar(find(calc_entrar>0))); %valor que entra na base
    if isempty(val_entrar)
        tipo=2; %solução ilimitada 
        return;
    end
    base_entrar=N_set(find(calc_entrar==val_entrar));  %indice que entra na base
    base_entrar=base_entrar(1);
    A_PP(pos_sair,:) = A_PP(pos_sair,:)./A_PP(pos_sair,base_entrar);
    b_PP(pos_sair)=b_PP(pos_sair)/Aj(pos_sair,base_entrar);
    %Atualizar Valores de custo
    c_PP=c_PP-c_PP(base_entrar).*A_PP(pos_sair,:);
  
    %Fazer pivoteamento
    [A_PP,b_PP] = pivot_step(A_PP,pos_sair,base_entrar,b_PP);
    
    %atualiza bases i e j para a próxima iteração

    indice_i=find(i==base_sair);
    i(indice_i)=base_entrar;
    indice_j=find(N_set==base_entrar);
    N_set(indice_j)=base_sair;
    I=i;
     
     %atualiza variáveis 
    Ai=A_PP(:,i); %cria a matriz de variáveis basicas
    Ci=c_PP(i);   %vetor C'das variáveis basicas
    Cj=c_PP(N_set);  %vetor C' das variáveis não basicas


    xi=Ai\b_PP; %calcula o valor de X das variávies básicas
    Aj=A_PP(:,N_set); %cria a matriz de variáveis não basicas
    xj=zeros(length(N_set),1); %preenche com zeros o vetor das variáveis não básicas
    X=[];
    X(i,1)=xi;%vetor X de acordo com as variávei1's básicas calculadas
    X(N_set,1)=xj;%vetor X de acordo com as variáveis não básicas calculadas
%     for k=1:n
%         pos=find(Pi==k);
%         C_calc(k)=c(pos);
%     end
 
    Z=c*X;
end
end

function [A_PP,b_PP] = pivot_step(A_PP,i,j,b_PP)
 [n,m] = size(A_PP);

    for k=1:n
       if k~= i
        b_PP(k)=b_PP(k)-A_PP(k,j).*b_PP(i);  
        A_PP(k,:)=A_PP(k,:)-A_PP(k,j).*A_PP(i,:);
        
       end
    end
end

