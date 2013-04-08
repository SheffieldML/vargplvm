if 0 
delete TEMPbetaGradTrC.mat
 delete TEMPbetaGradNksigm.mat
 delete TEMPbetaGradPsi0.mat
 delete TEMPbetaGradTrPP.mat
 delete TEMPbetaGradLat.mat
 delete TEMPbetaGradPlm.mat
 TEMPbetaGradTrC=[]; save 'TEMPbetaGradTrC.mat' 'TEMPbetaGradTrC';
 TEMPbetaGradNksigm=[];save 'TEMPbetaGradNksigm.mat' 'TEMPbetaGradNksigm';
 TEMPbetaGradPsi0=[];save 'TEMPbetaGradPsi0.mat' 'TEMPbetaGradPsi0';
 TEMPbetaGradTrPP=[];save 'TEMPbetaGradTrPP.mat' 'TEMPbetaGradTrPP';
 TEMPbetaGradLat=[];save 'TEMPbetaGradLat.mat' 'TEMPbetaGradLat';
 TEMPbetaGradPlm=[];save 'TEMPbetaGradPlm.mat' 'TEMPbetaGradPlm';

 clear
 fixedBetaIters=0;
 experimentNo=31;
 itNo=20;
 
 demHighDimVargplvm1
 
 load TEMPbetaGradTrC.mat
 load TEMPbetaGradNksigm.mat
 load TEMPbetaGradPsi0.mat
 load TEMPbetaGradTrPP.mat
 load TEMPbetaGradLat.mat
 load TEMPbetaGradPlm.mat
 
 plot(1:length(TEMPbetaGradTrC),TEMPbetaGradTrC); figure
 plot(1:length(TEMPbetaGradNksigm),TEMPbetaGradNksigm); figure  %const
 plot(1:length(TEMPbetaGradPsi0),TEMPbetaGradPsi0); figure
 plot(1:length(TEMPbetaGradTrPP),TEMPbetaGradTrPP); figure
 plot(1:length(TEMPbetaGradLat),TEMPbetaGradLat); figure
 plot(1:length(TEMPbetaGradPlm),TEMPbetaGradPlm); 
 end
 
 
 %-----------------------------------------
 
 delete TEMPbetaLikTrC.mat
 delete TEMPbetaLikNksigm.mat
 delete TEMPbetaLikPsi0.mat
 delete TEMPbetaLikTrPP.mat
 delete TEMPbetaLikLat.mat
 TEMPbetaLikTrC=[]; save 'TEMPbetaLikTrC.mat' 'TEMPbetaLikTrC';
 TEMPbetaLikNksigm=[];save 'TEMPbetaLikNksigm.mat' 'TEMPbetaLikNksigm';
 TEMPbetaLikPsi0=[];save 'TEMPbetaLikPsi0.mat' 'TEMPbetaLikPsi0';
 TEMPbetaLikTrPP=[];save 'TEMPbetaLikTrPP.mat' 'TEMPbetaLikTrPP';
 TEMPbetaLikLat=[];save 'TEMPbetaLikLat.mat' 'TEMPbetaLikLat';

 clear
 fixedBetaIters=0;
 experimentNo=31;
 itNo=21;
 
 demHighDimVargplvm1
 
 load TEMPbetaLikTrC.mat
 load TEMPbetaLikNksigm.mat
 load TEMPbetaLikPsi0.mat
 load TEMPbetaLikTrPP.mat
 load TEMPbetaLikLat.mat
 
 plot(1:length(TEMPbetaLikTrC),TEMPbetaLikTrC); figure % closer and closer to 0
 plot(1:length(TEMPbetaLikNksigm),TEMPbetaLikNksigm); figure  %const
 plot(1:length(TEMPbetaLikPsi0),TEMPbetaLikPsi0); figure % closer and closer to 0
 plot(1:length(TEMPbetaLikTrPP),TEMPbetaLikTrPP); figure % This one changes a log depending on whetehr b is fixed or not
 % almost converges when b is fixed, decreases when not. Same for the one below
 plot(1:length(TEMPbetaLikLat),TEMPbetaLikLat); figure % Ths one changes a lot depending on whether b  is fixed or not
