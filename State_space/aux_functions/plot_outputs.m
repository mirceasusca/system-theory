function plot_outputs(t,YY)
% Given the output matrix of a discrete LTI system, plot the individual
% outputs within a figure.
%
[n,~]=size(YY);
figure;
for i=1:n
   subplot(n,1,i);
   stairs(t, YY(i,:)); hold on;
   title('Output plot');
   ylabel(['y_', num2str(i), '[k]']);
end
xlabel('Time [sec]');

end