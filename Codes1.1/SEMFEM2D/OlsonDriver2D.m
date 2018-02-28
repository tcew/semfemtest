clear all
close all

%% Driver for building FEM-SEM precon for triangles 
Globals2D;

%% 1 for WarpBlend nodes, 0 for Electrostatic nodes
useWarpBlend = 0;

%% Polynomial order used for approximation 
for N = 1:10
  
  VX = [-1  1 -1];
  VY = [-1 -1  1];

  EToV = [1 2 3];
  K = 1;
  Np = (N+1)*(N+2)/2;
  
  if(useWarpBlend==1)
    %% Warp & Blend nodes
    [r,s] = Nodes2D(N);
    [r,s] = xytors(r,s);
  else
    nodefile = sprintf('ElectrostaticNodes/um_n%d.mat', N+1); 
    foo = load(nodefile);
    r = foo.umR;
    s = foo.umS;
  end

  %% build matrices for degree N element		       
  V  = Vandermonde2D(N, r, s);
  MM = inv(V')/V;
  [Dr,Ds] = Dmatrices2D(N, r, s, V);
  
  %% build coordinates of all the nodes
  va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
  x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
  y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));
		       
  [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
  
  vmapDIR = [];
  [tnum, NsolveTri, triGather, triScatter, triWinv] = GalerkinNumber2D(vmapDIR);

  va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
  xN = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
  yN = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

  rQ = r;
  sQ = s;
  NpM = Np;
  NfpM = Nfp;
		       
  va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
  xQ = 0.5*(-(rQ+sQ)*VX(va)+(1+rQ)*VX(vb)+(1+sQ)*VX(vc));
  yQ = 0.5*(-(rQ+sQ)*VY(va)+(1+rQ)*VY(vb)+(1+sQ)*VY(vc));

  Np_bak = Np;
  Nfp_bak = Nfp;
  Np = NpM;
  Nfp = NfpM;
  vmapDIR = [];
  [tnum, Nsolve, galGather, galScatter, galWinv] = GalerkinNumber2D(vmapDIR);
  Np = Np_bak;
  Nfp = Nfp_bak;

  VXG(tnum) = xQ;
  VYG(tnum) = yQ;

  VQN = Vandermonde2D(N, rQ(:), sQ(:));
  IQN = VQN/V;

  EToVG = delaunay(VXG,VYG);

  %% visualize FEM mesh
  figure(1)
  subplot(3,4, N);
  if(N==1) trimesh(EToVG,VXG',VYG')
  else trimesh(EToVG,VXG,VYG)
  end
  axis equal
  axis off		        
  title(sprintf('N=%d', N))		       

  %% build FEM infrastructure
  KG = size(EToVG,1);
  NG = 1;
  [rG,sG] = Nodes2D(NG);
  [rG,sG] = xytors(rG,sG);
  VG = Vandermonde2D(NG,rG,sG);
  [DrG,DsG] = Dmatrices2D(NG, rG, sG, VG);
  MMG = inv(VG')/VG;

  SrrG = DrG'*MMG*DrG;
  SrsG = DrG'*MMG*DsG;
  SsrG = DsG'*MMG*DrG;
  SssG = DsG'*MMG*DsG;

  va = EToVG(:,1)'; vb = EToVG(:,2)'; vc = EToVG(:,3)';
  xG = 0.5*(-(rG+sG)*VXG(va)+(1+rG)*VXG(vb)+(1+sG)*VXG(vc));
  yG = 0.5*(-(rG+sG)*VYG(va)+(1+rG)*VYG(vb)+(1+sG)*VYG(vc));
    
  %% build C0 stiffness
  [rxG,sxG,ryG,syG,JG] = GeometricFactors2D(xG,yG,DrG,DsG);
    
  %% remove slithers and reorient negative elements

  cnt = 1;
  for k=1:KG
    if(min(abs(JG(:,k)))>1e-10)
     if(min(JG(:,k))<0)
	  EToVG(cnt,[1 3 2]) = EToVG(k,:);
     else
	  EToVG(cnt,[1 2 3]) = EToVG(k,:);
     end
     cnt = cnt + 1;
    end
  end
  KG = cnt-1;
  EToVG = EToVG(1:KG,:);
	    
  %% reset geometric factors after stripping slithers
  va = EToVG(:,1)'; vb = EToVG(:,2)'; vc = EToVG(:,3)';
  xG = 0.5*(-(rG+sG)*VXG(va)+(1+rG)*VXG(vb)+(1+sG)*VXG(vc));
  yG = 0.5*(-(rG+sG)*VYG(va)+(1+rG)*VYG(vb)+(1+sG)*VYG(vc));

  [rxG,sxG,ryG,syG,JG] = GeometricFactors2D(xG,yG,DrG,DsG);
  
  lambda = 0;

  %% loop over elements in minigrid
  AFEM = zeros(Nsolve);
  for k=1:KG

    DxG = diag(rxG(:,k))*DrG + diag(sxG(:,k))*DsG;
    DyG = diag(ryG(:,k))*DrG + diag(syG(:,k))*DsG;
    MMJG = diag(JG(:,k))*(inv(VG')/VG);

    Ak = DxG'*MMJG*DxG + DyG'*MMJG*DyG + lambda*MMJG;

    gids = EToVG(k,:);
    AFEM(gids,gids) = AFEM(gids,gids) + Ak;
    
  end

  % single element
  k = 1;
  Dx = diag(rx(:,k))*Dr + diag(sx(:,k))*Ds;
  Dy = diag(ry(:,k))*Dr + diag(sy(:,k))*Ds;
  MM = diag(J(:,k))*(inv(V')/V);

  ASEM = Dx'*MM*Dx + Dy'*MM*Dy + lambda*MM;

  %% just boost null space for all Neumann problem
  if(lambda == 0)
    ASEM = ASEM + ones(Np,1)*ones(1,Np)/Np;		       
    AFEM = AFEM + ones(Np,1)*ones(1,Np)/Np;		       
  end

  %% evaluate condition number of AFEM\ASEM
  condPreconA(N) = cond(AFEM\ASEM);
  eigPreconA{N} = sort(eig(AFEM\ASEM), 'ascend');

  outName = sprintf('olson_N%d.mat', N);

  save(outName, "r", "s", "AFEM", "ASEM")

  figure(2)
  hold on
  plot(eigPreconA{N}, N*ones(size(eigPreconA{N})), '-*')
  title('Eigenvalues of AFEM\ASEM')
  xlabel('Eigenvalues')
  ylabel('Polynomial degree')
end
K
resultCondPreconA = condPreconA'
