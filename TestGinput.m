close all
clear
clc
% Sample data
x = 1:10;
y = rand(1,10);

% Plot the data
plot(x, y, '-o');
xlabel('X-axis');
ylabel('Y-axis');
title('Select a point from the plot');

% Select a point using ginput
[x_selected, y_selected] = ginput(1);  % '1' means select one point

% Display the selected point
hold on;
plot(x_selected, y_selected, 'r*', 'MarkerSize', 10);
hold off;

% Output the selected coordinates
fprintf('Selected point: (%.2f, %.2f)\n', x_selected, y_selected);
close(gcf);
