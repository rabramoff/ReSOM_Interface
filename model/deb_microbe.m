function g=deb_microbe(m0,g0,cnp,yld,je,ev)
%
%Description
%do DEB calculation

%g=deb_microbe(m0,g0,cnp,yld,je,ev)
%first the metabolic reserve turnover is used to support maintenance
%respiration. If there is excessive carbon flux then it is used to support
%growth and other activity. Were there any stresses, that should be
%effectively change the maintenance requirement
%On input it is assued the maintenance requirement has included cost to
%deal with stress
%input variables:
%m0: maintenance demand
%g0: a list of (nx1) maximum production rates
%cnp: a matrix of (nx2) elemental stoichiometry for the process indicated
%in g0
%yld: a list of (nx1) yielding rate for given processes, in the unit of carbon.
%je: the export metabolic flux to cell structure
%the returning variables should in the order of [netgrowth, activity
%investment]
%determine the number of elements
global deb;

deb.m0=m0;
deb.g0=g0;
deb.cnp=cnp;
deb.yld=yld;
deb.je=je;
deb.nelm=length(deb.cnp);
deb.ev=ev;
deb.g=zeros(size(g0));
deb.iter=0;
deb.maxiter=100;

Extra=[]; %empty structure

gb=brent_zero ( 0d0, 1d0, eps, eps, @deb_grow ,Extra);
    
    
if(deb.iter==deb.maxiter)
    gb=deb.gpm(1);
end

residual=deb_grow(gb,Extra);
%fprintf('iter=%d,residual=%f,gb0=%f,gb=%f,ge=%f\n',deb.iter,residual,gb,deb.g);
g=deb.g;

end


function residual=deb_grow(gb,Extra)
%Extra is an structure to pass parameters, but it is empty for the moment
global deb;

jc0=deb.je-gb*deb.ev;

 
%compute the excessive carbon
dc=jc0(1)-deb.m0;

if(dc>0)
    %there is carbon to support growth activity
    %compute the actual carbon flux to support growth, the yield rate here
    %is less than 1, indicating the fraction of carbon being truned into
    %the required structure after taking off the overhead
    jc=dc.*deb.yld; 

    if(deb.nelm==1)
        %carbon only
        %compute the potential growth rate for different processes
        gp=1./(1./deb.g0+1./jc);
        %compute the down-regulation factor
        scal_c=min(dc/sum(gp./deb.yld),1d0);
        %compute the actual growth rate
        deb.g=gp.*scal_c;        
    elseif(deb.nelm==2)
        %carbon and nitrogen
        %compute the potential growth rate for different processes
        gp=1./(1./deb.g0+1./jc+1./jc0(1)-1./(jc+jc0(1)));
        %compute c-based down-regulation factor
        scal_c=dc/sum(gp./deb.yld);
        %compute nitrogen based down-regulation factor
        scal_n=jc0(1)/sum(gp./deb.cnp(:,1));
        %compute the actual growth rate
        deb.g=gp.*min(scal_c,scal_n);        
    elseif(deb.nelm==3)
        %carbon and nitrogen and phosphorus
        %compute the potential growth rate for different processes
        gp=1./(1./deb.g0+1./jc+1./jc0(2)+1./jc0(3)...
            -1./(jc+jc0(2))-1./(jc+jc0(3))-1./(jc0(2)+jc0(3))+...
            1./(jc+jc0(2)+jc0(3)));
        %compute c-based down-regulation factor
        scal_c=dc/sum(gp./deb.yld);
        %compute nitrogen based down-regulation factor
        scal_n=jc0(2)/sum(gp./deb.cnp(:,1));
        %compute phosphorus based down-regulation factor
        scal_p=jc0(3)/sum(gp./deb.cnp(:,2));
        %compute the actual growth rate
        deb.g=gp.*min([scal_c,scal_n,scal_p]);           
    end
    residual=deb.je(1)-deb.g(1)*deb.ev(1)-deb.m0-sum(deb.g./deb.yld);    
else
    gp=0;
    deb.g=zeros(size(deb.yld));
    residual=deb.je(1)-deb.g(1)*deb.ev(1)-deb.m0;
    if(isnan(residual))
        fprintf('je=%f,g=%f,ev=%f,m0=%f\n',deb.je(1),deb.g(1),deb.ev(1),deb.m0);
    end
end

if(deb.iter==0)

    deb.residual=residual;
    deb.gpm=gp;

else
    if(abs(residual)<abs(deb.residual))
        deb.residual=residual;
        deb.gpm=gp;
    end
end
deb.iter=deb.iter+1;
if(deb.iter==deb.maxiter)
    residual=0d0;
end

end