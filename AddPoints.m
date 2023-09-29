function newpts = AddPoints(samples,nadd,fun,jac,lb,ub,rad)

pl = 1;

ncells = size(samples,2);
dim    = size(samples,1);


newpts = samples;

thetas = linspace(0,2*pi,nadd+1);
thetas = thetas(1:end-1); 

coords = zeros(ncells,2);
if pl==1
figure(1)
clf
set(gcf,'DefaultAxesFontSize',18,'DefaultAxesTickLabelInterpreter','latex')
hold on
end
for i=1:ncells
        coords(i,:) = fun(samples(:,i));
        actJac = jac(samples(:,i));

        [u, s, vh] = svd(actJac);
        %actJac
        %u*s*vh.'
        vh = vh';
        
        s = diag(s);
        
        if min(s)>1e-3
        w = vh(1:length(s),:).'*diag(1./s)*u;
        %w(:,1) = w(:,1)/norm(w(:,1));
        %w(:,2) = w(:,2)/norm(w(:,2));

        
        for j = 1:length(thetas)
        	theta = thetas(j);
        	newp = samples(:,i)+rad(i)*(cos(theta)*w(:,1)+sin(theta)*w(:,2));
        	tcor = fun(newp);

        	if and(min(newp)>=lb,max(newp)<=ub)
        		newpts = [newpts newp(:)];
        		if(pl==1)
        		plot(tcor(1),tcor(2),'.r');
        		end
        	else
        		%plot(tcor(1),tcor(2),'.g');
        	end
        	%pause
        end
        
       end
end
if(pl==1)
figure(1)


plot(coords(:,1),coords(:,2),'.b');

hold off
end

axis equal

size(newpts)

