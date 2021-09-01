% This m file creates the space factor and saves it in sf.m
% Once the function has been modified, RUN sfinput.m to generate
% the new data for the new function.
% Example 1:
% Rectangular SF
% Width: 45 to 135 degrees.
x=0:180;
y=zeros(1,length(x));
for i=1:length(x)
    if (x(i)>=60)&&(x(i)<=120)
        y(i)=1;
    end
end
save sf1.m x y -ascii;

% Once the function has been modified, RUN sfinput.m to generate
% the new data for the new function.
% Example 2:
% Triangular SF
% Width: 45 to 135 degrees.
x=0:180;
y=zeros(1,length(x));
for i=1:length(x)
    if (i>=60)&&(i<=120)
        y(i)=1/45*x(i)-1;
    elseif (i>=120)&&(i<=120)
        y(i)=-1/45*x(i)+3;
    end
end
save sf2.m x y -ascii;