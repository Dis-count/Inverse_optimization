x = 0:80; % range for graph
y1 = max(75 - x,0); % x + y <= 75 farm land
y2 = max((4000-110*x)/30,0); % 110x + 30y <= 4000 storage
y3 = max((15000 - 120*x)/210,0); % 120x + 210y <= 15000 expenses
ytop = min([y1; y2; y3]); % array of minima
area(x,ytop); % filled area plot