function [bnet] = mkBnet(structure)
if(isequal(structure, 'cancer'))
    nNodes = 5;
    dag = zeros(nNodes, nNodes);
    labels = {'A', 'B', 'C', 'D', 'E'};
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
    cptA = [0.1 0.9]; 
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
 elseif(isequal(structure, 'tree'))
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
elseif(isequal(structure, 'ten'))
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
end




