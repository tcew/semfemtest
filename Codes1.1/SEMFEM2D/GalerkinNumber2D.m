
function [gnum, Nsolve, galGather, galScatter, galWinv] = GalerkinNumber2D(vmapDIR)
  
  Globals2D;

  % number local nodes consecutively
  gnum = reshape(1:K*Np, Np, K);

  % boost numbers of boundary face nodes
  %gnum(vmapDIR) = (K*Np+1):(K*Np+length(vmapDIR));

  % now use neighbor info to uniquely number nodes
  change=1;
  while(change)
    change=0;
    
    newgnum = gnum;

    for k=1:K
      for f=1:Nfaces
        ids = (k-1)*Nfaces*Nfp + (f-1)*Nfp + (1:Nfp);
        idsL = newgnum(vmapM(ids));
        idsR = newgnum(vmapP(ids));
        newgnum(vmapM(ids)) = max(idsL, idsR);
        newgnum(vmapP(ids)) = max(idsL, idsR);
      end
    end
    
    change = length(find(newgnum~=gnum));
    
    gnum = newgnum;
  end

  % create consecutive numbering
  gids = zeros(max(max(gnum)),1);
  gids(gnum) = 1;
  ids = find(gids);
  gids(ids) = 1:length(ids);
  gnum = gids(gnum);

  % find lowest boundary numbered node
  Nsolve = max(max(gnum));

  % build Galerkin gather operaot
  galNtotal = max(max(gnum));
  
  galScatter = spconvert([[1:Np*K]', gnum(:), ones(Np*K,1)]);
  if(size(galScatter,1)<Np*K | size(galScatter,2)< galNtotal)
    galScatter(Np*K,galNtotal) = 0;
  end

  galGather = spconvert([gnum(:),[1:Np*K]', ones(Np*K,1)]);
  if(size(galGather,2)<Np*K | size(galGather,1)< galNtotal)
    galGather(Np*K,galNtotal) = 0;
  end

  %wJ = diag(sum(inv(V')/V))*J;
  %galWinv = 1./(galGather*wJ(:));
    


if(0)
  PlotGrid2D();
  for k=1:K
    for n=1:Np
      text(x(n,k),y(n,k), sprintf('%d', gnum(n,k)))
    end
  end
  end
  
  
