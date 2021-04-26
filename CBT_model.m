%% Supplementary MATLAB code for:

%Simulating the computational mechanisms of Cognitive and Behavioral Psychotherapeutic Interventions: 
%Insights from active inference

%Ryan Smith, Michael Moutoussis, Edda Bilek

%Published in Scientific Reports, 2021, DOI: https://doi.org/10.1038/s41598-021-89047-0



%% Active inference model of cognitive-behavioral therapy

%Ryan Smith

clear, close all
rng('shuffle')

%% Simulation Setting

CABi = .9;% 0-1, higher = stronger cognitive-affective-behavioral interactions
Psafe = 0.1;% 0-1, higher = stronger explicit belief that spider is now safe

Exposure = 0;% 0 = single trial; 1 = simulate exposure therapy

N = 200; % length of exposure

%% Model Specification

%--------------------------------------------------------------------------
%Initial State Priors (D Vectors)
%--------------------------------------------------------------------------

D{1} = [1 0 0 0 0 0]'; % start/stim/approach/interact/avoid/safety+cost 
D{2} = [0 1]'; %no spider/spider
D{3} = [0 1]'; % dangerous/safe


%--------------------------------------------------------------------------
%Patient's Explicit Beliefs ('d' Vectors)
%--------------------------------------------------------------------------

d = D;

d{2} = [1 1]';
d{3} = [1-Psafe Psafe]';

d{1} = d{1}*128;
d{2} = d{2}*128;
d{3} = d{3}*50;

%--------------------------------------------------------------------------
%State-Observation Mappings (A Matrices)
%--------------------------------------------------------------------------

Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
No    = [2 2 3 6]; % possible outcomes per modality: spider, arousal, affective consequences, behavior


Ng    = numel(No);
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
end

%-----------------------------------------------------------------------
% Spider observations

%state: no spider

A{1}(:,:,1,1) =   [1 1 1 1 1 1;% no spider
                   0 0 0 0 0 0];% spider
             
A{1}(:,:,1,2) =   [1 1 1 1 1 1;% no spider
                   0 0 0 0 0 0];% spider
               
%state: spider

A{1}(:,:,2,1) =   [1 0 0 0 0 0;% no spider
                   0 1 1 1 1 1];% spider
             
A{1}(:,:,2,2) =   [1 0 0 0 0 0;% no spider
                   0 1 1 1 1 1];% spider
               

%-----------------------------------------------------------------------
% Arousal observations
             
%state: no spider

A{2}(:,:,1,1) =   [1 1 1 1 0 1;% low arousal
                   0 0 0 0 1 0];% high arousal

A{2}(:,:,1,2) =   [1 1 1 1 0 1;% low arousal
                   0 0 0 0 1 0];% high arousal
               
%state: spider

A{2}(:,:,2,1) =   [1 0 1 1 0 1;% low arousal
                   0 1 0 0 1 0];% high arousal
             
A{2}(:,:,2,2) =   [1 0 1 1 0 1;% low arousal
                   0 1 0 0 1 0];% high arousal      
%-----------------------------------------------------------------------
% Affective consequences

%state: no spider

A{3}(:,:,1,1) =   [1 1 1 1 0 0;% positive affect
                   0 0 0 0 1 1;% negative affect
                   0 0 0 0 0 0];% serious harm

A{3}(:,:,1,2) =   [1 1 1 1 0 0;% positive affect
                   0 0 0 0 1 1;% negative affect
                   0 0 0 0 0 0];% serious harm

%state: spider + dangerous

A{3}(:,:,2,1) =   [1 1 0 0 0 0;% positive affect
                   0 0 1 0 1 1;% negative affect
                   0 0 0 1 0 0];% serious harm
             
%state: spider + safe

A{3}(:,:,2,2) =   [1 1 1 1 0 0;% positive affect
                   0 0 0 0 1 1;% negative affect
                   0 0 0 0 0 0];% serious harm
               
%-----------------------------------------------------------------------
% Behavioral observations

for i = 1:2
    for j = 1:2

A{4}(:,:,i,j) = eye(6);
    
    end 
end

%--------------------------------------------------------------------------
%Patient's Implicit Beliefs ('a' Matrices)
%--------------------------------------------------------------------------
         
a = A;

a{1} = a{1}*128;
a{2} = a{2}*128;
a{4} = a{4}*128;

a{3}(:,:,2,1) =   [1 1 .1 .1 0 0;% positive affect
                   0 0 .9  0 1 1;% negative affect
                   0 0  0 .9 0 0];% serious harm
               
a{3}(:,:,2,2) =   [1 1 CABi   CABi   0 0;% positive affect
                   0 0 1-CABi 0      1 1;% negative affect
                   0 0 0      1-CABi 0 0];% serious harm
               
a{3} = a{3}*5;


%--------------------------------------------------------------------------
%Action-Dependent State Transitions (B Matrices)
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end
 

    B{1}(:,:,1) = [0 0 0 0 0 0;
                   1 0 0 0 0 0;
                   0 1 0 0 0 0;
                   0 0 1 1 1 1;
                   0 0 0 0 0 0;
                   0 0 0 0 0 0]; %approach
   
    
    B{1}(:,:,2) = [0 0 0 0 0 0;
                   1 0 0 0 0 0;
                   0 0 0 0 0 0;
                   0 0 0 0 0 0;
                   0 1 0 0 0 0;
                   0 0 1 1 1 1]; %avoid



%--------------------------------------------------------------------------
%Allowable Policies (Pi)
%--------------------------------------------------------------------------
T = 4; % number of timesteps

Np        = 2; % number of policies

V         = ones(T-1,Np,Nf);
V(:,:,1) = [1 1 1; % Approach Trajectory
            2 2 2]'; % Avoid Trajectory
        

%Prior over policies (E)
    
E = [1 1]';

%--------------------------------------------------------------------------
%Preferences (C vectors)
%--------------------------------------------------------------------------
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{3}      = zeros(No(3),T);
C{4}      = zeros(No(4),T);

C{2}(2,:) = -1; % dysprefer high arousal

C{3}(2,:) = -1; % dysprefer negative affect
C{3}(3,:) = -12; % strongly dysprefer serious harm

%% Define MDP Structure

mdp.T = T;                      % number of time points
mdp.V = V;                      % allowable policies
mdp.A = A;                      % state-outcome mappings
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.E = E;                      % prior over policies

mdp.d = d;                      % explicit beliefs
mdp.a = a;                      % implicit beliefs

 
mdp.beta    = 1;                % policy precision
mdp.alpha   = 4;                % action precision (inverse temperature)
mdp.eta     = 1;                % learning rate

label.factor{1}   = 'Appoach/Avoid';   label.name{1}    = {'start','stimulus','approach','interact','avoid','safety+cost'};
label.factor{2}   = 'Clinically Relevant Stimulus';     label.name{2}    = {'no spider','spider'};
label.factor{3}   = 'Explicit Belief about Stimulus';     label.name{3}    = {'dangerous','safe'};

label.modality{1} = 'Visual';    label.outcome{1} = {'no spider','spider'};
label.modality{2} = 'Interoceptive Arousal';  label.outcome{2} = {'low arousal','high arousal'};
label.modality{3} = 'Consequences';    label.outcome{3} = {'positive affect','negative affect','serious harm'};
label.modality{4} = 'Actions';    label.outcome{4} = {'start','stimulus','approach','interact','avoid','safety+cost'};
label.action{1} = {'approach','avoid'};

mdp.label = label;

MDP         = spm_MDP_check(mdp);

MDP = spm_MDP_VB_X(MDP);

% Illustrate single trial before exposure

spm_figure('GetWin','Pre-exposure behavior'); clf
spm_MDP_VB_trial(MDP,2:3,2:3);


%% Simulate Exposure Therapy

if Exposure == 1

mdp.E = [1 0]';

MDP2(1:N) = deal(mdp);

A_startD = MDP2(1).a{3}(:,:,2,1);

MDP2 = spm_MDP_VB_X(MDP2);

A_endD = MDP2(N).a{3}(:,:,2,1);

figure
imagesc(A_startD), colormap gray % before learning
title('Implicit beliefs before therapy')

    set(gca,'XTickLabel',{'Start','Observe','Approach', 'Interact','Avoid','Safety+Cost'});
    set(gca,'YTickLabel',{'','Positive Affect', '','Negative Affect','','Serious Harm'});

figure
imagesc(A_endD), colormap gray % before learning
title('Implicit beliefs after therapy')

    set(gca,'XTickLabel',{'Start','Observe','Approach', 'Interact','Avoid','Safety+Cost'});
    set(gca,'YTickLabel',{'','Positive Affect', '','Negative Affect','','Serious Harm'});


for i = 1:N
     d_evo(1:2,i) = MDP2(i).d{1,3}(:,1);
end


figure
imagesc(spm_softmax(d_evo,.1)), colormap gray % before learning
title('Explicit beliefs over time in exposure')

    set(gca,'YTickLabel',{'','Dangerous', '','Safe',''});
    
end
    
return
%% Sample approach-avoidance behavior under beliefs at a given time point during exposure

Test_N = 200; % time point to sample (0-200)

mdp.a  = MDP2(Test_N).a;
mdp.d = MDP2(Test_N).d;
mdp.E = [1 1]';


MDP3         = spm_MDP_check(mdp);

MDP3 = spm_MDP_VB_X(MDP3);

spm_figure('GetWin','Post-Exposure Behavior'); clf
spm_MDP_VB_trial(MDP3,2:3,1:3);




