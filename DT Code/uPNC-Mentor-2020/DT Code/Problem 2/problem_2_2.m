pred = [];
modk =[];
base = [];
for i =1:232
    N = problem_2(Data, i, false);
    base = [base; N(1)];
    modk = [modk; N(2)];
    pred = [pred;N(3)];
end

    figure;
for i =1:232
    if i ==1
        scalar = compass(max(modk),0);
        scalar.Color = 'none';
    end
    x = cos(pred(i))*modk(i);
    y = sin(pred(i))*modk(i);
    hold on
    compass(x,y)
end
title('Rose Compass Plot scaled by Modulation Depth')
%%
figure;
for i =1:232
  x1 = cos(pred(i));
  y1 = sin(pred(i));
    if i ~= 1
        hold on
        compass(x1,y1)
    else 
        compass(x1,y1)
        hold on
    end
end
title('Rose Compass Plot Unit Length')
hold off