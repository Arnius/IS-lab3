%% 
clc
clear
close all


% Iejimo vektorius
x= 0.1:1/22:1;

d = ((1 + 0.6*sin(2*pi*x/0.7)) + (0.3*sin(2*pi*x)))/2;

% Parametrai
n = 0.25;
% Pirmas sluoksnis
w0 = rand(1);
w11 = rand(1);
w12 = rand(1);

% % Spindulio tipo funkciju parametrai
c1 = x(3);
r1 = x(6) - x(3);

c2 = x(19);
r2 = abs(x(19) - x(16));

% Mokymosi epochos
for index=1:2000

    % Neuronų atsako ir aktyvavimo funkcijų skaičiavimai, bei svorių
    % atnaujinimas su kiekvienų x(i) pavyzdžiu.
    for i = 1:length(x)
       F1 = gaussFunction(x(i),c1,r1);
       F2 = gaussFunction(x(i),c2,r2);
       Y = F1*w11 + F2*w12 + w0;

       e = d(i) - Y;

       w11 = w11 + n*e*F1;
       w12 = w12 + n*e*F2;
       w0 = w0 + n*e;
    end
end

% Funkcijos aproksimavimas su naujomis x vertėmis
apx_x = linspace(0.1,1,200);
apx_y = zeros(1,length(apx_x));
for i=1:length(apx_x)
       
       F1 = gaussFunction(apx_x(i),c1,r1);
       F2 = gaussFunction(apx_x(i),c2,r2);
       apx_y(i) = F1*w11 + F2*w12 + w0;
end

plot(x,d,'bo', apx_x,apx_y,'r--')
title("Funkcijos aproksimacija")
ylabel('y')
xlabel('x')

%%
clc
clear
close all

% Iėjimo vektorius
x= 0.1:1/22:1;
d = ((1 + 0.6*sin(2*pi*x/0.7)) + (0.3*sin(2*pi*x)))/2;

abs_en = 0.01; % Nominali paklaida
abs_e = 1;    % Gauta pakalaida

% Parametrai
n = 0.25;
% Pirmas sluoksnis
w0 = rand(1);
w11 = rand(1);
w12 = rand(1);

% Spindulio tipo bazinių funkcijų parametrai
% Išrenkame visus galimus lokalius maksimumus ir sudedame jų x vertes
% į centrų parametrų matricą
c = local_maximums(x,d);

% Parenkame spindulio didžius pagal savo nuožiūrą
postumis = 1;   % Postumis lokalaus maksimumo x vertės atžvilgiu
r = radius_update(x, c, postumis);

% Mokymosi epochos su automatiškai prisiderinančiu spindulio ilgiu
while 1
    for index=1:2000
        % Neuronų atsako ir aktyvavimo funkcijų skaičiavimai, bei svorių
        % atnaujinimas su kiekvienų x(i) pavyzdžiu.
        for i = 1:length(x)
           F1 = gaussFunction(x(i),c(1,1), r(1,1));
           F2 = gaussFunction(x(i),c(1,2), r(1,2));
           Y = F1*w11 + F2*w12 + w0;
    
           e = d(i) - Y;
    
           w11 = w11 + n*e*F1;
           w12 = w12 + n*e*F2;
           w0 = w0 + n*e;
        end
    end
    
    % Funkcijos aproksimavimas su naujomis x vertėmis
    apx_x = linspace(0.1,1,200);
    apx_y = zeros(1,length(apx_x));
    for i=1:length(apx_x)
    
           F1 = gaussFunction(apx_x(i),c(1,1),r(1,2));
           F2 = gaussFunction(apx_x(i),c(1,2),r(1,2));
           apx_y(i) = F1*w11 + F2*w12 + w0;
    end
    
    taskai = plot(x,d,'bo');
    hold on
    kreive = plot(apx_x, apx_y, 'r--');
    title("Funkcijos aproksimacija")
    ylabel('y')
    xlabel('x')
    
    % While ciklo nutraukimo sąlyga
    if (abs_e < abs_en) && (abs_e > -abs_en), break, end
    
    pause(0.5)
    delete(kreive)
    pause(0.5)

    % Randam lokalius maksimumus iš aproksimuotos kreivės
    apx_max = local_maximums(apx_x, apx_y);
    
    % if isempty(apx_max)
    %     apx_max(1,1) = apx_x(1,21);
    %     apx_max(1,2) = apx_x(1,182);
    % end

    % Paskaičiuojam bendrą aproksimuotos ir realios funkcijos paklaidą
    % pagal atimdami centrų reikšmes
    abs_e = 0;  % Ištrinama paskutinė paklaida
    % Ieškomi indeksai su kuriais paėmus vertes iš x matricos gaunami y
    % reikšmių maksimumai
    for indx = 1:length(apx_max)
         id1 = find(apx_max(1,indx) == apx_x);
         id2 = find(c(1,indx) == x);
         abs_e = abs_e + (apx_y(1,id1) - d(1,id2))  % Paklaidos skaičiavimas
    end

    % Parenkam naujus spindulius
    % Jei kreivė aukščiau aproksimuojamų taškų
    if abs_e > 0.5
        postumis = postumis + 1;
        r = radius_update(x, c, postumis)
    elseif (abs_e < 0.5) && (abs_e > 0)
        for indx=1:length(r)
            r(1,indx)= r(1,indx) + 0.01
        end
    % Jei kreive žemiau aproksimuojamų taškų
    elseif abs_e < -0.5
        postumis = postumis - 1;
        r = radius_update(x, c, postumis)
    elseif (abs_e > -0.5) && (abs_e < 0)
        for indx=1:length(r)
            r(1,indx)= r(1,indx) - 0.01
        end
    end

end  % While pabaiga

% Spindulio atnaujinimo funkcija naudojant poslinkį
function r = radius_update(x, c, offset)
    temp = [];
    mx_it = 1;  % Nepriklausomas c masyvo iteratorius
    for indx=1:length(x)
        if x(1,indx) == c(1,mx_it)
            X = 0;
            if indx + offset < length(x)
                X = abs(x(indx + offset) - x(indx));
            else
                X = abs(x(indx - offset) - x(indx));
            end
            temp =[temp X];
            mx_it = mx_it + 1;
        end
        if mx_it > length(c), break, end
    end
    r = temp;
end

% Lokalių maksimumų atradimo funkcija
function max = local_maximums(x,y)
    c = [];
    for indx=1:length(x)-2
        Dl = y(indx+1) - y(indx);    % Taško vertė iš kairės pusės 
        Dr = y(indx+2) - y(indx+1);  % Taško vertė iš dešinės pusės
        if Dl > 0 && Dr < 0
            c = [c x(indx+1)];
        end 
    end
    max = c;
end

% Gauso spindulio funkcija
function [F] = gaussFunction(x,c,r)
F = exp(-((x-c)^2)/(2*r^2));
end