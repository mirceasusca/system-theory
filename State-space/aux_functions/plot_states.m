function plot_states(t,XX)
% Given the state matrix from a discrete LTI system, plot the individual
% states within a figure.
%
[n,~]=size(XX);
figure;
for i=1:n
   subplot(n,1,i);
   stairs(t, XX(i,:)); hold on;
   title('State plot');
   ylabel(['x_', num2str(i), '[k]']);
end
xlabel('Time [sec]');

end