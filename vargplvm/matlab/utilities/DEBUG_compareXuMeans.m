% compare X_u and means:

load demMissaVargplvm41
model = prunedModelTr;
means = sort(model.vardist.means);
xu=sort(model.X_u);

for i=1:model.q
parentHandle = bar(means(:,i),'r');
childHandle = get(parentHandle,'Children');
set(childHandle,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
hold on
parentHandle2=bar(xu(:,i),'b');
childHandle2=get(parentHandle2,'Children');
set(childHandle2,'FaceAlpha',0.5);
legend('means','X_u')
pause(0.5)
hold off
end
disp('press any key')
pause
disp('variances')
bar(var(model.vardist.means));
figure
bar(var(model.X_u),'r');