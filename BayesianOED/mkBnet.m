function [bnet] = mkBnet(structure)
if(isequal(structure, 'cancer'))
    nNodes = 5;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C', 'D', 'E'}
    A = 1; B = 2; C = 3; D = 4; E =5;
    dag(A,[B C]) = 1;
    dag(B,D) = 1;
    dag(C,[D E]) = 1;
    ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
    p = .9;
%cptA = [.4 .6];
cptA = [0.1 0.9]; % must be different from [0.5 0.5] to detect intervention
cptB = [1-p p ; p 1-p];
cptC = [1-p p ; p 1-p];
cptD = zeros(2, 2, 2);
cptD(:,:,1) = [ .15 .3 ; .7 .9];
cptD(:,:,2) = 1-cptD(:,:,1);
cptE = [p 1-p ; 1-p p];

bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', cptA);
bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', cptB);
bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', cptC);
bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', cptD);
bnet.CPD{5} = tabular_CPD(bnet, 5, 'CPT', cptE);

elseif(isequal(structure, 'tree'))
    nNodes = 5;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C', 'D', 'E'}
    A = 1; B = 2; C = 3; D = 4; E =5;
    dag(A, [B C D E]) = 1;
    ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
elseif(isequal(structure, 'line'))
    nNodes = 8;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C', 'D', 'E','F','G','H'}
    A = 1; B = 2; C = 3; D = 4; E =5; F=6; G=7;H=8;
    dag(A, B) = 1;
    dag(B, C) = 1;
    dag(C, D) = 1;
    dag(D, E) = 1;
    dag(E, F) = 1;
    dag(F, G) = 1;
    dag(G, H) = 1;
    ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
elseif(isequal(structure, 'tree3'))
    nNodes = 3;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C'}
    A = 1; B = 2; C = 3;
    dag(A, [B C]) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
elseif(isequal(structure, 'line3'))
    nNodes = 3;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C'}
    A = 1; B = 2; C = 3;
    dag(A, B) = 1;
    dag(B, C) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
elseif(isequal(structure, 'tree4'))
    nNodes = 4;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C','D'}
    A = 1; B = 2; C = 3; D = 4;
    dag(A, [B C D]) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
elseif(isequal(structure, 'line4'))
    nNodes = 4;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C','D'}
    A = 1; B = 2; C = 3; D = 4;
    dag(A, B) = 1;
    dag(B, C) = 1;
    dag(C, D) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
    elseif(isequal(structure, 'tree5'))
    nNodes = 5;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C','D','E'}
    A = 1; B = 2; C = 3; D = 4; E = 5;
    dag(A, [B C D E]) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
 elseif(isequal(structure, 'tree8'))
    nNodes = 8;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C','D','E','F','G','H'}
    A = 1; B = 2; C = 3; D = 4; E = 5; F=6; G=7; H = 8;
    dag(A, [B C D E F G H]) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;  
% cptA = [0.3 0.7];
% cptB = [0.6 0.4; 0.4 0.6];
% cptC = [0.6 0.4; 0.4 0.6]; 
% cptD = [0.6 0.4; 0.4 0.6];
% cptE = [0.6 0.4; 0.4 0.6];
% cptF = [0.4 0.6; 0.6 0.4];
% cptG = [0.4 0.6; 0.6 0.4];
% cptH = [0.4 0.6; 0.6 0.4];
% bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', cptA);
%  bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', cptB);
%  bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', cptC);
%  bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', cptD);
%  bnet.CPD{5} = tabular_CPD(bnet, 5, 'CPT', cptE);
%  bnet.CPD{6} = tabular_CPD(bnet, 6, 'CPT', cptF);
%  bnet.CPD{7} = tabular_CPD(bnet, 7, 'CPT', cptG);
%  bnet.CPD{8} = tabular_CPD(bnet, 8, 'CPT', cptH);

elseif(isequal(structure, 'line5'))
    nNodes = 5;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C','D','E'};
    A = 1; B = 2; C = 3; D = 4; E = 5;
    dag(A, B) = 1;
    dag(B, C) = 1;
    dag(C, D) = 1;
    dag(D, E) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
%     cptA = [0.3 0.7];
%     cptB = [0.6 0.4; 0.4 0.6];
%     cptC = [0.6 0.4; 0.4 0.6]; 
%     cptD = [0.6 0.4; 0.4 0.6];
%     cptE = [0.6 0.4; 0.4 0.6];
%  bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', cptA);
%  bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', cptB);
%  bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', cptC);
%  bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', cptD);
%  bnet.CPD{5} = tabular_CPD(bnet, 5, 'CPT', cptE);
elseif(isequal(structure, 'line8'))
    nNodes = 8;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C','D','E','F','G','H'};
    A = 1; B = 2; C = 3; D = 4; E = 5; F = 6; G=7; H=8;
    dag(A, B) = 1;
    dag(B, C) = 1;
    dag(C, D) = 1;
    dag(D, E) = 1;
    dag(E, F) = 1;
    dag(F, G) = 1;
    dag(G, H) = 1;
     ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;
elseif(isequal(structure, 'parent3'))
    nNodes = 3;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C'}
    A = 1; B = 2; C = 3; 
    dag(A, [B C]) = 1;
    dag(B, C) = 1;
    ns = 2*ones(nNodes,1);
    bnet =  mk_bnet(dag, ns);
    bnet.nNodes = nNodes;
    bnet.dag = dag;
    bnet.names = labels;
    bnet.eligibleNodes = 1:nNodes;
    bnet.ns = ns;   
    
elseif(isequal(structure, 'toySachs11'))
    % % 2/16/19 
nNodes = 11;
dag = zeros(nNodes);
labels = {'pip3','plcy','pip2','pkc','pka','raf','mek','erk','jnk','p38','akt'};
pip3 = 1; plcy = 2; pip2 = 3; pkc = 4; pka = 5; raf = 6;
mek = 7; erk = 8; jnk = 9; p38 = 10; akt = 11;
dag(pip3, [plcy pip2 akt]) = 1;
dag(plcy, [pip2 pkc]) = 1;
dag(pip2, pkc) = 1;
dag(pkc, [jnk p38 pka raf mek]) = 1;
dag(pka, [jnk p38 akt erk mek raf]) = 1;
dag(raf, mek) = 1;
dag(mek, erk) = 1;
dag(erk, akt) = 1;
ns = 3*ones(nNodes,1);
bnet =  mk_bnet(dag, ns);
bnet.nNodes = nNodes;
bnet.dag = dag;
bnet.names = labels;
%bnet.eligibleNodes = [2,4,7,8,9];
bnet.eligibleNodes = 1:nNodes;
bnet.ns = ns;  
elseif(isequal(structure, 'asia'))
    % % 8/30/19
nNodes = 8;
dag = zeros(nNodes);
labels = {'asia','tub','smoke','lung','bronc','either','xray','dysp'};
asia = 1; tub = 2; smoke = 3; lung = 4; bronc = 5; either = 6;
xray = 7; dysp = 8;
dag(asia, tub) = 1;
dag(tub, either) = 1;
dag(smoke, [lung bronc]) = 1;
dag(lung, either) = 1;
dag(either, [dysp xray]) = 1;
dag(bronc, dysp) = 1;
ns = 2*ones(nNodes,1);
bnet =  mk_bnet(dag, ns);
bnet.nNodes = nNodes;
bnet.dag = dag;
bnet.names = labels;
bnet.eligibleNodes = 1:nNodes;
bnet.ns = ns;
elseif(isequal(structure, 'lung8'))
    % % 4/15/19
nNodes = 8;
dag = zeros(nNodes);
labels = {'asia','tub','smoke','lung','either','bronc','dysp','xray'};
asia = 1; tub = 2; smoke = 3; lung = 4; either = 5; bronc = 6;
dysp = 7; xray = 8;
dag(asia, tub) = 1;
dag(tub, either) = 1;
dag(smoke, [lung bronc]) = 1;
dag(lung, either) = 1;
dag(either, [dysp xray]) = 1;
dag(bronc, dysp) = 1;
ns = 2*ones(nNodes,1);
bnet =  mk_bnet(dag, ns);
bnet.nNodes = nNodes;
bnet.dag = dag;
bnet.names = labels;
bnet.eligibleNodes = 1:nNodes;
bnet.ns = ns;

cptA = [0.3 0.7];
 cptT = [0.8 0.2; 0.3 0.7];
 cptS = [0.4 0.6]; 
 cptL = [0.9 0.1; 0.1 0.9];
  cptE = zeros(2, 2, 2);
 cptE(:,:,1) = [ .15 .3 ; .7 .9];
 cptE(:,:,2) = 1-cptE(:,:,1);
 cptB = [0.9 0.1; 0.1 0.9];
 cptD = zeros(2, 2, 2);
 cptD(:,:,1) = [ .15 .3 ; .7 .9];
 cptD(:,:,2) = 1-cptD(:,:,1); 
 cptX = [0.9 0.1; 0.1 0.9];
 bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', cptA);
 bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', cptT);
 bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', cptS);
 bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', cptL);
 bnet.CPD{5} = tabular_CPD(bnet, 5, 'CPT', cptE);
 bnet.CPD{6} = tabular_CPD(bnet, 6, 'CPT', cptB);
 bnet.CPD{7} = tabular_CPD(bnet, 7, 'CPT', cptD);
 bnet.CPD{8} = tabular_CPD(bnet, 8, 'CPT', cptX);

elseif(isequal(structure, 'survey6'))
    % % 4/19/19
nNodes = 6;
dag = zeros(nNodes);
labels = {'a','s','e','o','r','t'};
a = 1; s = 2; e = 3; o = 4; r = 5; t = 6;
dag(a, e) = 1;
dag(s, e) = 1;
dag(e, [o r]) = 1;
dag(o, t) = 1;
dag(r, t) = 1;
ns = 2*ones(nNodes,1);
bnet =  mk_bnet(dag, ns);
bnet.nNodes = nNodes;
bnet.dag = dag;
bnet.names = labels;
bnet.eligibleNodes = 1:nNodes;
bnet.ns = ns;
elseif(isequal(structure, 'disjoint9'))
    % % 4/21/19
nNodes = 9;
dag = zeros(nNodes);
labels = {'A','B','C','D','E','F','G','H','I'};
A=1; B=2; C=3; D=4; E=5; F=6; G=7; H=8; I=9;
dag(A, [B C]) = 1;
dag(B, C) = 1;
dag(D, F) = 1;
dag(E, F) = 1;
dag(F, G) = 1;
dag(G, [H I]) = 1;
ns = 2*ones(nNodes,1);
bnet =  mk_bnet(dag, ns);
bnet.nNodes = nNodes;
bnet.dag = dag;
bnet.names = labels;
bnet.eligibleNodes = 1:nNodes;
bnet.ns = ns;
elseif(isequal(structure, 'ten2'))
    % % 7/23/19
nNodes = 10;
dag = zeros(nNodes);
labels = {'A','B','C','D','E','F','G','H','I','J'};
A=1; B=2; C=3; D=4; E=5; F=6; G=7; H=8; I=9; J=10;
dag(A, H) = 1;
dag(B, G) = 1;
dag(C, [B J]) = 1;
dag(D, [A C F]) = 1;
dag(F, E) = 1;
dag(G, E ) = 1;
dag(H,[I J]) = 1;
dag(I,[B C] ) = 1;
dag(J,F )= 1;
ns = 2*ones(nNodes,1);
bnet =  mk_bnet(dag, ns);
bnet.nNodes = nNodes;
bnet.dag = dag;
bnet.names = labels;
bnet.eligibleNodes = 1:nNodes;
bnet.ns = ns;
elseif(isequal(structure, 'ten4'))
    % % 7/23/19
nNodes = 10;
dag = zeros(nNodes);
labels = {'A','B','C','D','E','F','G','H','I','J'};
A=1; B=2; C=3; D=4; E=5; F=6; G=7; H=8; I=9; J=10;
dag(A, [C D G I]) = 1;
dag(B, [A F G H]) = 1;
dag(C, [J]) = 1;
dag(D, [C I]) = 1;
dag(E, [A B C G H]) = 1;
dag(F, [G I]) = 1;
dag(H, [A D F J]) = 1;
dag(I, [C J]) = 1;
ns = 2*ones(nNodes,1);
bnet =  mk_bnet(dag, ns);
bnet.nNodes = nNodes;
bnet.dag = dag;
bnet.names = labels;
bnet.eligibleNodes = 1:nNodes;
bnet.ns = ns;
end


% %% MZ 2/16/19 make DAG that Sachs reported in Figure 3a
% % but include the dashed adges and reverse the direction the purple arrow
% % between plcy and pip3 to better match biological consensus network
% nNodes = 11;
% dag = zeros(nNodes, nNodes);
% labels = {'pip3','plcy','pip2','pkc','pka','raf','mek','erk','jnk','p38','akt'};
% pip3 = 1; plcy = 2; pip2 = 3; pkc = 4; pka = 5; raf = 6;
% mek = 7; erk = 8; jnk = 9; p38 = 10; akt = 11;
% dag(pip3, [plcy pip2 akt]) = 1;
% dag(plcy, [pip2 pkc]) = 1;
% dag(pip2, pkc) = 1;
% dag(pkc, [jnk p38 pka raf mek]) = 1;
% dag(pka, [jnk p38 akt erk mek raf]) = 1;
% dag(raf, mek) = 1;
% dag(mek, erk) = 1;
% dag(erk, akt) = 1;
% %%
% ns = 3*ones(nNodes,1);
% bnet =  mk_bnet(dag, ns);
% bnet.nNodes = nNodes;
% bnet.dag = dag;
% bnet.names = labels;
% %bnet.eligibleNodes = 1:11;
% bnet.eligibleNodes = [3,4,5,7,11]; %Sachs intervened on pip2, pkc, pka, mek, akt
% bnet.ns = ns;
% end

% d = 11;
% G = zeros(d,d);
% labels = {'raf','mek12','plcy','pip2','pip3','erk','akt','pka','pkc','p38','jnk'};
% raf= 1; mek12= 2; plcy= 3; pip2= 4; pip3= 5; erk= 6;
% akt= 7; pka= 8; pkc= 9; p38= 10; jnk=11;
% %% MZ:DAG edges for what I am calling ground truth (Figure 6c minus one edge)
% G(pip3, pip2) = 1;
% G(plcy, [pip3 pip2]) = 1;
% % pip2 no children in Fig 6c
% G(pkc, [mek12 raf p38 pka jnk plcy pip3 erk]) = 1;
% G(pka, [raf p38 mek12 erk akt jnk ]) = 1;
% G(raf, [mek12 akt]) = 1;
% G(mek12, [erk jnk akt]) = 1;
% G(jnk, [p38 plcy]) = 1;
% %G(akt, plcy) = 1; removing this edge 2/5/19 
% G(p38, plcy) = 1;
% G(erk, akt) = 1;
%% MZ:DAG edges that Eaton and Murphy found in Figure 6c - trying again
% G(pip3, pip2) = 1;
% G(plcy, [pip3 pip2]) = 1;
% % pip2 no children in Fig 6c
% G(pkc, [mek12 raf p38 pka jnk plcy pip3 erk]) = 1;
% G(pka, [raf p38 mek12 erk akt jnk ]) = 1;
% G(raf, [mek12 akt]) = 1;
% G(mek12, [erk jnk akt]) = 1;
% G(jnk, [p38 plcy]) = 1;
% G(akt, plcy) = 1;
% G(p38, plcy) = 1;
% G(erk, akt) = 1;
%%
% %% MZ:DAG edges that Eaton and Murphy found in Figure 6c
% G(pip3, pip2) = 1;
% G(plcy, [pip3 pip2]) = 1;
% % pip2 no children in Fig 6c
% G(pkc, [mek12 raf p38 pka jnk plcy pip3 ]) = 1;
% G(pka, [raf p38 mek12 erk akt jnk ]) = 1;
% G(raf, [mek12 akt]) = 1;
% G(mek12, [erk jnk akt]) = 1;
% G(jnk, p38) = 1;
% G(akt, plcy) = 1;

%% Eaton & Murphy coded: Biologists consensus DAG (from Figure 6a)
% G(pip3,[akt plcy pip2])=1;
% G(plcy,[pkc pip2])=1;
% G(pip2,pkc)=1;
% G(pkc,[mek12 raf p38 pka jnk])=1;
% G(pka,[raf p38 mek12 erk akt jnk])=1;
% G(raf,mek12)=1;
% G(mek12,erk)=1;
% G(erk,akt)=1;



%% things the original demos2/sachsTrueDAG.m function had
% g06976 = 1; aktinh = 2;  psitect = 3; u0126 = 4; b2camp = 5; pma = 6;
% e = 6;
% F = zeros(e,d);
% % activators
% F(pma,pkc)=1;
% F(b2camp,pka)=1;
% % inhibitors
% F(g06976,pkc)=1;
% F(aktinh,akt)=1;
% F(psitect,pip2)=1;
% F(u0126,mek12)=1;
% 
% %H = [F zeros(e,d);
% %     G zeros(d,d)];
% 
% H = [zeros(e,e) F;
%      zeros(d,e) G];


