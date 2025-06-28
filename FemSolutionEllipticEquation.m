I = input('Aralık Gir : ');
n = input('Aralığı kaç eşit parçaya böleceğiz:');
f = @(x,y) -1;
a = @(x,y) 1%exp(-1/(x^2 + y^2)); 
x = linspace(I(1),I(2),n+1);
k = @(x,y) 0;  
gD = @(x,y) 0;                % Dirichlet sınır koşulu
gN = @(x,y) 0;                      % Neumann Sınır koşulu

P = [];
for i = 1 : n+1
    for j = 1 : n+1
        P = [P,[x(j);x(i)]];
    end
end
T = [];
for i = 1:n
    for j = 1:n
        n1 = (n+1)*(i-1)+j;
        n2 = n1+1;
        n3 = n2+n;
        n4 = n3+1;
        T = [T,[n1;n2;n4],[n4;n3;n1]];
    end
end
np = length(P); % düğüm sayısı
nt = length(T); % üçgen sayısı
% Sınır kenarlarını oluştur
e = [];

% Alt kenar 
for i = 1:n
    n1 = i;
    n2 = i + 1;
    e = [e, [n1; n2]];
end

% Sağ kenar 
for i = 1:n
    n1 = (n+1)*i;
    n2 = n1 + (n+1);
    e = [e, [n1; n2]];
end

% Üst kenar 
for j = n:-1:1
    n1 = (n+1)*n + (j+1);
    n2 = n1 - 1;
    e = [e, [n1; n2]];
end

% Sol kenar 
for i = n:-1:1
    n1 = (n+1)*(i) + 1;
    n2 = n1 - (n+1);
    e = [e, [n1; n2]];
end
ne = length(e); % kenar sayısı
A = zeros(np,np);
for K = 1:nt
    loc2glb = T(1:3,K); % üçgenlerin düğüm numaraları
    x = P(1,loc2glb); 
    y = P(2,loc2glb); % üçgenin x ve y kordinatlarını tut
    alan = polyarea(x,y);
    b = [y(2)-y(3);y(3)-y(1);y(1)-y(2)]/(2*alan); % gradiyentler
    c = [x(3)-x(2);x(1)-x(3);x(2)-x(1)]/(2*alan);
    abar = a((x(1) + x(2) + x(3))/3,(y(1)+y(2)+y(3))/3);
    local_matrix = abar * alan * [b(1)^2 + c(1)^2, b(1)*b(2)+c(1)*c(2),b(1)*b(3)+c(1)*c(3);
                                  b(2)*b(1)+c(2)*c(1),b(2)^2+c(2)^2,b(2)*b(3)+c(2)*c(3);
                                  b(3)*b(1)+c(3)*c(1),b(3)*b(2)+c(3)*c(2),b(3)^2 + c(3)^2];
    A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + local_matrix;
end

b = zeros(np,1);
for K = 1:nt
    loc2glb = T(1:3,K);
    x = P(1,loc2glb);
    y = P(2,loc2glb);
    alan = polyarea(x,y);
    localb = 1/3 * alan * [f(x(1),y(1));f(x(2),y(2));f(x(3),y(3))];
    b(loc2glb) = b(loc2glb) + localb;
end

R = zeros(np,np); % sınır katkıları
for E = 1:ne
    loc2glb = e(1:2,E);
    x = P(1,loc2glb);
    y = P(2,loc2glb);
    edge_length = sqrt((y(2)-y(1))^2 + (x(2)-x(1))^2);
    R_local = (1/6) * edge_length * [2 1 ; 1 2] * k((x(1)+x(2))/2,(y(1)+y(2))/2);
    R(loc2glb,loc2glb) = R(loc2glb,loc2glb) + R_local;
end
r = zeros(np,1); % sınır katkıları
for E = 1:ne
    loc2glb = e(1:2,E);
    x = P(1,loc2glb);
    y = P(2,loc2glb);
    edge_length = sqrt((y(2)-y(1))^2 + (x(2)-x(1))^2);
    r_local = (1/2) * (k((x(1)+x(2))/2,(y(1)+y(2))/2)*gD((x(1)+x(2))/2,(y(1)+y(2))/2) + gN((x(1)+x(2))/2,(y(1)+y(2))/2))*edge_length;
    r(loc2glb) = r(loc2glb) + r_local;
end
stifness = A+R;
load = b+r;


np = size(P,2); % Toplam düğüm sayısını belirle
fixed = unique([e(1,:); e(2,:)]); % Sınır düğümlerinin indislerini al
free = setdiff([1:np],fixed); % İç düğümleri belirle

g = zeros(ne,1);
for i = 1:ne
    g(i)=g(i)+gD(P(1,fixed(i)),P(2,fixed(i)));
end

load = load(free)-stifness(free,fixed)*g; % b vektörünü güncelle

stifness = stifness(free,free); % stifness matrisini güncelle
xi = zeros(np,1); % Çözüm vektörünü tahsis et
xi(fixed) = g; % Sınır düğümlerini yerleştir
xi(free) = stifness\load; % iç düğüm değerlerini çöz

figure;
trisurf(T', P(1,:), P(2,:), xi);
colorbar;
view(3); % 3D görünüm
xlabel('x');
ylabel('y');
zlabel('Çözüm xi');
title('Sonlu Elemanlar Yöntemi ile 3D Çözüm Yüzeyi');

e = sort(eigs(stifness, 10, 'sm'));
disp(e);
condest(stifness)


