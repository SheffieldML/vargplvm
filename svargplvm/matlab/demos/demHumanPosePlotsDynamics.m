% DEMHUMANPOSEPLOTSDYNAMICS
%
% SHEFFIELDML

ha = tight_subplot(3,4,[.01 .01],[.01 .01],[.01 .01]);
zview = [-90 0];
c = 1;
axes(ha(c)); c = c + 1;
xyzankurVisualise2(Z_test(78,:));
axes(ha(c)); c = c + 1;
xyzankurVisualise2(ZpredAll(78,:));
axes(ha(c)); c = c + 1;
xyzankurVisualise2(Ytr{2}(mini(78),:)); %NN_Y_sil
axes(ha(c)); c = c + 1;
xyzankurVisualise2(Ytr{2}(miniAll2(78),:)); %NN_Y_pos


axes(ha(c)); c = c + 1;
imagesc(reshape(Yim_test(78,:),height,width)), colormap('gray'), axis equal,axis tight, axis off
axes(ha(c)); c = c + 1;
imagesc(reshape(Yim(indsAllOrig(78),:), height, width)), colormap('gray'), axis equal, axis tight, axis off %from NN_X
axes(ha(c)); c = c + 1;
imagesc(reshape(Yim(mini(78),:), height, width)), colormap('gray'), axis equal, axis tight, axis off %NN_Y_sil
axes(ha(c)); c = c + 1;
imagesc(reshape(Yim(miniAll2(78),:), height, width)), colormap('gray'), axis equal, axis tight, axis off %NN_Y_pos
axes(ha(c)); c = c + 1;

xyzankurVisualise2(Z_test(78,:),zview);
axes(ha(c)); c = c + 1;
xyzankurVisualise2(ZpredAll(78,:),zview);
axes(ha(c)); c = c + 1;
xyzankurVisualise2(Ytr{2}(mini(78),:),zview);
axes(ha(c)); c = c + 1;
xyzankurVisualise2(Ytr{2}(miniAll2(78),:),zview);

%%
