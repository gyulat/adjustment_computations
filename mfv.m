x=[-12.5 -6.7 -2 -1.5 0.1 2.4 6.8 9.8 15 23.5 30];
M1=mean(x);
epsilon1=0.5*sqrt(3)*(max(x)-min(x));
itermax=30;
for j=1:itermax
    counter=0; counter2=0;
    name=0; name2=0;
    seg=0; seg2=0;
if j==1
        epsilon(j)=epsilon1;
        M(j)=M1;
    else
        for i=1:length(x)
            seg=(x(i)-M(j-1))^2;
            counter=counter+3*((seg)/(((epsilon(j-1)^2)+seg)^2));
            name=name+(1/((epsilon(j-1)^2)+seg)^2);
        end
        epsilon(j)=sqrt(counter/name);
        for i=1:length(x)
            seg2=(epsilon(j-1)^2)/((epsilon(j-1)^2)+((x(i)-M(j-1))^2));
            counter2=counter2+(seg2*x(i));
            name2=name2+seg2;
        end
        M(j)=counter2/name2;
    end
end

M
epsilon
