% 
% Sana Waheed
% Please contact at sana.waheed@smme.nust.edu.pk or sana.waheed@cantab.net
% 19th February, 2019
% 
% The compiled DDP code builds on top of Zebang Zheng's 'SourceCode V3.1' with the following changes:
% 1. 3-noded linear triangular elements are implemented. The error in Jacobian calculation from 6-noded elements has been fixed.
% 2. Reverse Cuthill-Mckee reordering of nodes to reduce stiffness matrix bandwidth.
% 3. Analytical solution for dislocation forces on nodes is implemeneted. 
%     See Waheed, S. (2018). 'Discrete dislocation plasticity analyses of
%     cyclic loading phenomena', (PhD thesis), Imperial College London, UK.
% 4. Dislocation velocity correction using stress gradient methodology introduced in Chakravarthy and Curtin (2010) MSMSE implemented.
% 5. Non-convex domains, including grains with holes can be partly be
% simulated, although more work needs to be done on this.
% 
% 
%  In addition the following features are also carried forward:
% . Dislocation penetration from Li et al (2009) Comp. Mat. Sci.
% . The version also carries forward the additional superposition of O'Day
% and Curtin (2004) J Appl Mech, for strain controlled and for load
% controlled loading.
% . Material heterogoneity
% . Thermal activation