function[X,Z,I,it,tipo]= simplex1(A,b,c,i)
%limpar as variáveis de saida
clear X;
clear Z;
clear BaseO;
clear it;
clear tipo;

%valores iniciais vazios
flag_otmfinita=[];
flag_candidate=[];
X=[];
Z=[];
I=[];
it=0;

%validar se e infactivel
if (size(A,2)~= size(c,2)) || (size(A,1)~= size(b,1)) | (rank([A b])> rank(A))
    tipo=-1;
    return;
end

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
Z=c*X;

while it>=0
    it=it+1; %atualiza o numero de interações

    pie=Ai'\Ci'; %calcular pi = Ci*Ai(inversa)
    custoj=pie'*Aj-Cj ;%computa o custo das variâveis não basicas Ĉ
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
   r_min=min(d(r)); %busca o meno elemento em d de acordo com os indices positivos
   leave = i(find(d==r_min));  %indice de quem sai da base
     
    %atualiza bases i e j para a próxima iteração

    indice_i=find(i==leave);
    i(indice_i)=enter;
    indice_j=find(N_set==enter);
    N_set(indice_j)=leave;
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
    Z=c*X;
end    