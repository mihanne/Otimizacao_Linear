function[X,Z,I,it,tipo]= primal_dual(A,b,c,pi)
%Implementação dual simplex
%A igual a matriz
%b igual ao vetor coluna
%c igual ao vetor linha de custos
%pi igual ao vetor linha do problema dual

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
q=[];
c_aux=[];
c_dual=[];
%busca o tamanho da matriz A
tamanho = size(A); % retorna o tamanho da matriz Amostras
m = tamanho(1); % n1 igual a quantidade de linhas
n = tamanho(2); % n2 igual a quantidade de colunas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%validar se e infactivel
if (size(A,2)~= size(c,2)) || (size(A,1)~= size(b,1)) | (rank([A b])> rank(A))
    tipo=-1;
    return;
end
%Monta o problema dual e testar os valores iniciais de pi
dual=pi*b;
if (dual==0)
    b_aux=b;
    Z=0;
    Z0=0;
else %calcula os valores de lado direito da tabela
    b_aux=b.*pi
    Z=dual;
    Z0=dual;
end    
A_aux = [A, eye(m)];
c_aux = [c*-1, zeros(1, m)];
c_dual= [zeros(1, n), ones(1, m)*-1];
x = [zeros(n, 1); b_aux];
i = (n + 1):(n + m);

%inicilialização das variáveis para a primeira iteração
set=[1:length(c)]; %cria um vetor com as posições do vetor c
set(find(ismember(set,i)==1))=[]; %busca somente as variáveis basicas

N_set=set;
I=i;

%Atualiza o valor de C_dual
for j=1:m
    for k=1:n+m
    c_dual(k)=A_aux(j,k)+c_dual(k);
    end
end
%Atualiza Z0
for k=1:m
    Z0=b_aux(k)+Z0;
end

Ai=A_aux(:,i); %cria a matriz de variáveis basicas

xi=Ai\b_aux; %calcula o valor de X das variávies básicas
Aj=A_aux(:,N_set); %cria a matriz de variáveis não basicas

while it>=0
    it=it+1; %atualiza o numero de interações
    
    %Calcula quem está em q
    q=find(c_aux=0); %quem é zero no vetor dual Z0
    
    %Calcula teta =min{-z/z0}, onde i não pertence a q
    calc_teta=(c_aux./c_dual)*-1;
    if (min(calc_teta<0))
        tipo=1; % solução otima finita e multipla não tem ninguem para entrar na base
    end
    
    teta=min(calc_teta(find(calc_teta>0))); %valor que entra na base
 
    pos_teta=find(calc_teta==teta); %refere-se a posição em J(N_set) de quem vai entrar na base
   
    %Encontrar quem sai da base
    y= b_aux./A_aux(:,pos_teta);   %calcula os valores de y
    if (min(y<0))
        tipo=2; %solução ilimitada - não tem ninguem para sair da base
        return;
    end
    
    sair=min(y(find(y>0))); %valor que entra na base
    pos_sair=find(y==sair); %posição de sair para o pivoteamento
    pos_sairBase=I(pos_sair); %variável que sai da base
    
    %Atualiza o vetor c_aux e Z, sendo Z=Z+teta*Z0
    c_aux=c_aux+teta*c_dual;
    Z=Z+teta*Z0;
    
    %Fazer pivoteamento
    [A_aux,b_aux,c_dual,Z0] = pivot_step(A_aux,pos_sair,pos_teta,b_aux,c_dual,Z0);

    %Atualiza variáveis I e J
    indice_i=find(i==pos_sairBase);
    i(indice_i)=pos_teta;
    indice_j=find(N_set==pos_teta);
    N_set(indice_j)=pos_sairBase;
    I=i;
    
    %atualiza variáveis Ai 
    Ai=A_aux(:,i); %cria a matriz de variáveis basicas
       
    xi=Ai\b_aux; %calcula o valor de X das variávies básicas
    Aj=A_aux(:,N_set); %cria a matriz de variáveis não basicas
    xj=zeros(length(N_set),1); %preenche com zeros o vetor das variáveis não básicas
    X=[];
    X(i,1)=xi;%vetor X de acordo com as variávei1's básicas calculadas
    X(N_set,1)=xj;%vetor X de acordo com as variáveis não básicas calculada
   
    if (Z0==0)
         tipo=0; % solução otima finita e unica
         return;
    end
    
 end
end

function [A_aux,b_aux,c_dual,Z0] = pivot_step(A_aux,i,j,b_aux,c_dual,Z0)
 [n,m] = size(A_aux);
%atualiza o vetor c_dual
ind_c=c_dual(j);
for col=1:m
    c_dual(col)=c_dual(col)-A_aux(i,col)*ind_c;
end
%Atualiza o valor de Z0
Z0=Z0-b_aux(i)*ind_c;

    for k=1:n
       if k~= i
        b_aux(k)=b_aux(k)-A_aux(k,j).*b_aux(i);  
        A_aux(k,:)=A_aux(k,:)-A_aux(k,j).*A_aux(i,:);
       end
    end
end

