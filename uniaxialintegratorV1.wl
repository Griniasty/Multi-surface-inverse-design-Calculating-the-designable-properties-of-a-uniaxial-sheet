(* ::Package:: *)

(* Inverse design of uniaxial sheets package *)
(*
so the data structure is 
{uv,r,n,\[Alpha],\[Beta],\[Epsilon],b,s,p,q,dp,dq}
r[[i]] =  the 2d vector representing the position on the i^th surface
n[[i]] =  the 2d vector representing the image of the director on the i^th surface
\[Epsilon][[i]] = the i'th material response
p[[i]] = the corresponding derivative 1/\[Alpha]\!\(
\*SubscriptBox[\(\[PartialD]\), \(u\)]\(\[Epsilon][\([i]\)]\)\)
q[[i]] = the corresponding derivative 1/q\!\(
\*SubscriptBox[\(\[PartialD]\), \(v\)]\(\[Epsilon][\([i]\)]\)\) 

The integration is supplemented by
SurfaceNumber=Number of surfaces in including original
\[CapitalLambda][[i]]=List of stimuli
g[[i]][r] = the metric of the original i'th target surface. g[[i]][r]:=({g1,g2,..,gn}[[i]])&r[[i]]
*)


(* ::Section::Closed:: *)
(*Warp and weft integration*)


WWintegrateUniaxial[Lu_,Lv_,r0_,n0_,eps0_,pq0_,lam1_,lam2_,Lambda_,g_,propagationEq_,i12_,du_,umax_,signu_,dv_,vmax_,signv_]:=
Module[{umaxtab,vmaxloc,data,tu,tv,uu,vv},
vmaxloc=vmax;
umaxtab=Table[umax,{vv,vmax/dv}];
data=datatable[eps0,umax,vmax,du,dv];
(*Origin point*)
data[[1,1]]=origin[r0,n0,eps0,pq0,Lu,Lv,g[[1]],i12];
(*Integrate initial ucurve vv=1*)
tu=0;
For[uu=2,uu<=umaxtab[[1]]/du,uu++,{tu,data[[uu,1]]}=N[integrateLUstep[Lu,eps0,pq0,tu,g,data[[uu-1,1]],propagationEq,i12,du*signu]];If[ustop[data[[uu,1,5]],data[[uu,1,6]],lam1,lam2,Lambda],umaxtab[[;;]]=du(uu-1)];];
tv=0;
For[vv=2,vv<=vmaxloc/dv,vv++,
(*integrate initial condition *)
{tv,data[[1,vv]]}=N[integrateLvstep[Lv,eps0,pq0,tv,g,data[[1,vv-1]],propagationEq,i12,dv*signv]];If[vstop[data[[1,vv,4]],data[[1,vv,6]],lam1,lam2,Lambda],vmaxloc=dv(vv-1)];
(*integrate r,n,\[Alpha],\[Beta],\[Epsilon],q along v*)
For[uu=2,uu<=umaxtab[[vv]]/du,uu++,
data[[uu,vv]]=N[integrateV[g,data[[uu,vv-1]],propagationEq,i12,signv*dv]];
If[vstop[data[[uu,vv,4]],data[[uu,vv,6]],lam1,lam2,Lambda],umaxtab[[vv;;]]=du*(uu-1);
If[uu<=2,Abort]];];
(*Derive p along u*)
data[[;;,vv,9]]=pdp[data[[;;,vv,6]],data[[;;,vv,4]],signu*du][[1]];
umaxtab[[vv+1;;]]=(umaxtab[[vv]]-du);
(*Intergate {\[Beta],s} along u*)
For[uu=2,uu<=umaxtab[[vv]]/du,uu++,
data[[uu,vv,{5,8}]]+=integrateU[g,data[[uu-1,vv]],propagationEq,i12,signu*du][[{5,8}]];
If[ustop[data[[uu,vv,5]],data[[uu,vv,6]],lam1,lam2,Lambda],umaxtab[[vv;;]]=(uu-1)*du;If[uu<=2,Abort]];];
If[umaxtab[[vv]]<=du,vmaxloc=(vv-1)*dv];
];
Return[{umaxtab,vmaxloc,data}];
];



(* ::Section::Closed:: *)
(*Diagonal integration*)


diagonalintegration[Lu_,Lv_,r0_,n0_,eps0_,pq0_,lam1_,lam2_,Lambda_,g_,propagationEq_,i12_,maxdiagonal_,du_,signu_,dv_,signv_]:=
Module[{vmaxtab,umaxtab,data,tu,tv,uu,vv,diagonal,maxdiagonalloc,uvratio,tempu},
data=datatable[eps0,maxdiagonal,maxdiagonal,1,1];
(*uvratio is the ratio of how much do we account for u vs v propagaion in r,n,\[Epsilon]*)
uvratio=1/2;
(*Origin point*)
data[[1,1]]=origin[r0,n0,eps0,pq0,Lu,Lv,g[[1]],i12];
tu=0;
tv=0;
maxdiagonalloc=maxdiagonal;
vmaxtab=Table[Infinity,{uu,maxdiagonal}];
umaxtab=Table[Infinity,{vv,maxdiagonal}];
diagonal=2;
While[diagonal<=maxdiagonalloc,
vv=1;
While[vv<=Min[vmaxtab[[1+diagonal-vv]],diagonal],
If[1+diagonal-vv>umaxtab[[vv]],
vv=1+diagonal-umaxtab[[vv]];];
uu=1+diagonal-vv;
If[vv==diagonal,
(*integrate data along v*)
{tv,data[[uu,vv]]}=N[integrateLvstep[Lv,eps0,pq0,tv,g,data[[1,vv-1]],propagationEq,i12,dv*signv]];If[vstop[data[[1,vv,4]],data[[1,vv,6]],lam1,lam2,Lambda],
(*Do[umaxtab[[v]]=Max[1,diagonal-v],{v,vv,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,uu,maxdiagonal}];*)
Do[umaxtab[[v]]=Max[1,diagonal-v],{v,1,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,1,maxdiagonal}];
maxdiagonalloc=diagonal-1;
Print["vstop at "<>ToString[{uu,vv}]];];
,
If[vv==1,
(*integrate data along u*)
{tu,data[[uu,vv]]}=N[integrateLUstep[Lu,eps0,pq0,tu,g,data[[uu-1,1]],propagationEq,i12,du*signu]];If[ustop[data[[uu,1,5]],data[[uu,1,6]],lam1,lam2,Lambda],
(*Do[umaxtab[[v]]=Max[1,diagonal-v],{v,vv,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,uu,maxdiagonal}];*)
Do[umaxtab[[v]]=Max[1,diagonal-v],{v,1,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,1,maxdiagonal}];
maxdiagonalloc=diagonal
Print["ustop at "<>ToString[{uu,vv}]];
If [uu<=2,Abort[]]];
,
(*integrate r,n,\[Alpha],\[Beta],\[Epsilon],q along v*)
data[[uu,vv]]=data[[uu,vv]]=N[integrateV[g,data[[uu,vv-1]],propagationEq,i12,signv*dv]];
(*integrate \[Epsilon] another step along v to derive the relevant (i12-1)*q *)
data[[uu,vv+1]][[6]]=(i12-1)*N[integrateV[g,data[[uu,vv]],propagationEq,i12,signv*dv]][[6]];
(*derive the relevant (i12-1)*q *)
data[[uu,vv]][[10]]=(i12-1)*qdq[data[[uu,vv;;vv+1,6]],data[[uu,vv;;vv+1,5]],signv*dv][[1,1]]+(2-i12)*data[[uu,vv,10]];
If[vstop[data[[uu,vv,4]],data[[uu,vv,6]],lam1,lam2,Lambda],
(*Do[umaxtab[[v]]=Max[1,diagonal-v],{v,vv,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,uu,maxdiagonal}];*)
Do[umaxtab[[v]]=Max[1,diagonal-v],{v,1,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,1,maxdiagonal}];
maxdiagonalloc=diagonal-1;
Print["vstop at "<>ToString[{uu,vv}]];];
(*Intergate {\[Beta],s} along u*)
tempu=integrateU[g,data[[uu-1,vv]],propagationEq,i12,signu*du];
data[[uu,vv,{5,8}]]+=tempu[[{5,8}]];
data[[uu,vv,{2,3,6}]]=(1-uvratio)data[[uu,vv,{2,3,6}]]+uvratio*tempu[[{2,3,6}]];
(*integrate \[Epsilon] another step along u to derive the relevant (2-i12)*p *)
data[[uu+1,vv]][[6]]=(2-i12)*integrateU[g,data[[uu,vv]],propagationEq,i12,signu*du][[6]]+
(i12-1)*data[[uu+1,vv,6]];
(*derive the relevant (2-i12)*p *)
data[[uu,vv]][[9]]=(2-i12)*( pdp[data[[uu;;uu+1,vv,6]],data[[uu;;uu+1,vv,4]],signu*du][[1,1]])+(i12-1)*data[[uu,vv]][[9]];If[ustop[data[[uu,vv,5]],data[[uu,vv,6]],lam1,lam2,Lambda],
(*Do[umaxtab[[v]]=Max[1,diagonal-v],{v,vv,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,uu,maxdiagonal}];*)
Do[umaxtab[[v]]=Max[1,diagonal-v],{v,1,maxdiagonal}];
Do[vmaxtab[[u]]=Max[1,diagonal-u],{u,1,maxdiagonal}];
maxdiagonalloc=diagonal-1;
(*umaxtab[[vv;;]]=(uu-1);
vmaxtab[[uu;;]]=(vv-1);*)
Print["ustop at "<>ToString[{uu,vv}]];
If[uu<=2,Abort]];If[umaxtab[[vv]]<=uu,vmaxtab[[uu]]=(vv-1)];
If[vmaxtab[[uu]]<=vv,umaxtab[[vv]]=(uu-1)];
];
];
vv++;
];
diagonal++;
];
umaxtab=Table[Min[(1+maxdiagonal-vv),umaxtab[[vv]]],{vv,maxdiagonal}];
vmaxtab=Table[Min[(1+maxdiagonal-uu),vmaxtab[[uu]]],{uu,maxdiagonal}];
Return[{umaxtab,vmaxtab,data}]
];


(* ::Section::Closed:: *)
(*Initial condition integration*)


(* Given initial uv curves Lu and Lv of the form {f1[t],f2[t]}
and: initial postions, initial directors, intitila \[Epsilon], initial p and q where appropariate*)


(* ::Subsection:: *)
(*Integrate the rest*)


integrateRest[g_,lam1_,lam2_,Lambda_,dat_,propagationEq_,i12_,du_,umax_,signu_,dv_,vmax_,signv_]:=Module[{umaxtab,vmaxloc,temp,data},
data=dat;
vmaxloc=vmax;
umaxtab=Table[umax,{vv,vmax/dv}];
(*integrate r,n,\[Alpha],\[Beta],\[Epsilon],q along u*)
For[vv=2,vv<=vmaxloc/dv,vv++,
For[uu=1,uu<=umaxtab[[vv]]/du,uu++,
temp=N[integrateV[g,data[[uu,vv-1]],propagationEq,i12,signv*dv]];
If[uu==1,
data[[uu,vv,2,2;;]]=temp[[2,2;;]];
data[[uu,vv,3,2;;]]=temp[[3,2;;]];
data[[uu,vv,4;;]]+=temp[[4;;]];,
data[[uu,vv]]+=temp;];
If[vstop[data[[uu,vv,4]],data[[uu,vv,6]],lam1,lam2,Lambda],umaxtab[[vv;;]]=du*(uu-1);
If[uu<=2,Abort]];];
(*Derive p along u*)
data[[;;,vv,9]]=pdp[data[[;;,vv,6]],data[[;;,vv,4]],signu*du][[1]];
umaxtab[[vv+1;;]]=(umaxtab[[vv]]-du);
(*Intergate {\[Beta],s} along v*)
For[uu=2,uu<=umaxtab[[vv]]/du,uu++,
data[[uu,vv,{5,8}]]+=integrateU[g,data[[uu-1,vv]],propagationEq,i12,signu*du][[{5,8}]];
If[ustop[data[[uu,vv,5]],data[[uu,vv,6]],lam1,lam2,Lambda],umaxtab[[vv;;]]=(uu-1)*du;If[uu<=2,Abort]];];
If[umaxtab[[vv]]<=du,vmaxloc=(vv-1)*dv];
];
Return[{umaxtab,vmaxloc,data}];
];


(* ::Subsection::Closed:: *)
(*Integrate initial curves *)


integrateInitialCurvesWarpWarf[Lu_,Lv_,r0_,n0_,eps0_,pq0_,lam1_,lam2_,Lambda_,g_,propagationEq_,i12_,du_,umax_,signu_,dv_,vmax_,signv_]:=
Module[{data0,dat,vv,tu,tv,vmaxlocal,umaxlocal,haswarpdv},
dat=datatable[eps0,umax,vmax,du,dv];
dat[[1,1]]=origin[r0,n0,eps0,pq0,Lu,Lv,g[[1]],i12];

tv=0;
vmaxlocal=vmax;
For[vv=2,vv<=vmaxlocal/dv,vv++,
{tv,dat[[1,vv]]}=N[integrateLvstep[Lv,eps0,pq0,tv,g,dat[[1,vv-1]],propagationEq,i12,dv*signv]];
If[stops[dat[[1,vv,4]],dat[[1,vv,5]],dat[[1,vv,6]],lam1,lam2,Lambda],vmaxlocal=dv(vv-1)];
];

haswarpdv={1,{{1,1}}~Join~(0*r0[[2;;]]),1,0,1,0,0,0,(i12-1),(i12-1)};
dat[[1,2;;]]=Table[haswarpdv,{vv,(vmax/dv)-1}]*dat[[1,2;;]];

tu=0;
umaxlocal=umax;
For[uu=2,uu<=umaxlocal/du,uu++,
{tu,dat[[uu,1]]}=N[integrateLUstep[Lu,eps0,pq0,tu,g,dat[[uu-1,1]],propagationEq,i12,du*signu]];
If[stops[dat[[uu,1,4]],dat[[uu,1,5]],dat[[uu,1,6]],lam1,lam2,Lambda],umaxlocal=du(uu-1)];
];

Return[{dat,umaxlocal,vmaxlocal}];
];



(* ::Subsection::Closed:: *)
(*Origin*)


origin[r0_,n0_,\[Epsilon]0_,pq0_,Lu_,Lv_,g_,i12_]:=Module[{x,y,dat},
dat={{0,0},r0,n0,1,1,\[Epsilon]0,\[Kappa]g[Lu,0,g,{x,y}],\[Kappa]g[Lv,0,g,{x,y}]}~Join~Transpose[MapThread[If[#1==1,{#2,#3},{#3,#2}]&,{i12,D[\[Epsilon]0,t],pq0}]]/.t->0;
Return[dat];];


(* ::Subsection::Closed:: *)
(*integrate initial curves step*)


integrateLUstep[Lu_,eps0_,pq0_,t0_,g_,data_,propagationEq_,i12_,du_]:=
Module[{\[Alpha]r,datar,tr,br,epsr,pr,qr,velocity},
velocity=Sqrt[D[Lu,t].g[[1]][Lu].D[Lu,t]]/.t->t0;
tr=t0+du/velocity/.t->t0;
datar=integrateU[g,data,propagationEq,i12,du];
\[Alpha]r=1;
br=bLu[Lu,g[[1]],tr];
epsr=(2-i12)*eps0+(i12-1)*datar[[6]]/.t->tr;
pr=(2-i12)*(D[eps0,t]/velocity)+(i12-1)*datar[[9]]/.t->tr;
qr=(2-i12)pq0+(i12-1)*datar[[10]]/.t->tr;
datar[[2]][[1]]=Lu/.t->tr;
datar[[3]][[1]]=nLu[Lu,g[[1]],tr];
datar[[4]]=\[Alpha]r;
datar[[6]]=epsr;
datar[[7]]=br;
datar[[9]]=pr;
datar[[10]]=qr;
Return[{tr,datar}];]

integrateLvstep[Lv_,eps0_,pq0_,t0_,g_,data_,propagationEq_,i12_,dv_]:=
Module[{\[Beta]r,datar,tr,sr,epsr,pr,qr,velocity},
velocity=Sqrt[D[Lv,t].g[[1]][Lv].D[Lv,t]];
tr=t0+dv/velocity/.t->t0;
datar=integrateV[g,data,propagationEq,i12,dv];
\[Beta]r=1;
sr=sLv[Lv,g[[1]],tr];
epsr=(i12-1)*eps0+(2-i12)*datar[[6]]/.t->tr;
pr=(i12-1)*pq0+(2-i12)*datar[[9]]/.t->tr;
qr=(i12-1)*(D[eps0,t]/velocity)+(2-i12)*datar[[10]]/.t->tr;
datar[[2]][[1]]=Lv/.t->tr;
datar[[3]][[1]]=nLv[Lv,g[[1]],tr];
datar[[5]]=\[Beta]r;
datar[[6]]=epsr;
datar[[8]]=sr;
datar[[9]]=pr;
datar[[10]]=qr;
Return[{tr,datar}];]


(* ::Subsection::Closed:: *)
(*Stop functions*)


(*stops*)
stops[alpha_,beta_,eps_,lam1_,lam2_,Lambda_]:=
Module[{lam1list,lam2list,lamneg,isotropic,stopconditions},
lam1list=(lam1[eps,#]&/@Lambda);
lam2list=(lam2[eps,#]&/@Lambda);
lamneg=Not[AllTrue[lam1list~Join~lam2list,#>0&]];
isotropic=AllTrue[lam1list-lam2list,#==0&];
stopconditions=Not[NumberQ[alpha]]||Not[NumberQ[beta]]||alpha<=0||beta<=0||lamneg||isotropic;
Return[stopconditions];];

ustop[beta_,eps_,lam1_,lam2_,Lambda_]:=
Module[{lam1list,lam2list,lamneg,isotropic,stopconditions},
lam1list=(lam1[eps,#]&/@Lambda);
lam2list=(lam2[eps,#]&/@Lambda);
lamneg=Not[AllTrue[lam1list~Join~lam2list,#>0&]];
isotropic=AllTrue[lam1list-lam2list,#==0&];
stopconditions=Not[NumberQ[beta]]||beta<=0||lamneg||isotropic;
Return[stopconditions];];

vstop[alpha_,eps_,lam1_,lam2_,Lambda_]:=
Module[{lam1list,lam2list,lamneg,isotropic,stopconditions},
lam1list=(lam1[eps,#]&/@Lambda);
lam2list=(lam2[eps,#]&/@Lambda);
lamneg=Not[AllTrue[lam1list~Join~lam2list,#>0&]];
isotropic=AllTrue[lam1list-lam2list,#==0&];
stopconditions=Not[NumberQ[alpha]]||alpha<=0||lamneg||isotropic;
Return[stopconditions];];



(* ::Subsection::Closed:: *)
(*Initial condition set up functions*)


(*Initial condition set up functions *)
npLv[Lv_,g_,t0_]:=ReplaceAll[D[Lv,t]/Sqrt[D[Lv,t].g[Lv].D[Lv,t]],t->t0];
(* nLv = find n=-Subscript[(Subscript[n, \[Perpendicular]]), \[Perpendicular]] along Lv(t) *)
nLv[Lv_,g_,t0_]:=-fnp[npLv[Lv,g,t0],g[ReplaceAll[Lv,t->t0]]];
(* sLv = find s= s0, s0=\[Kappa]g along Lv(t) *)
sLv[Lv_,g_,t0_]:=Module[{x,y},Return[ \[Kappa]g[Lv,t0,g,{x,y}]];];

nLu[Lu_,g_,t0_]:=ReplaceAll[D[Lu,t]/Sqrt[D[Lu,t].g[Lu].D[Lu,t]],t->t0];
(* bLu = find b=\[Kappa]g along Lu(t) *)
bLu[Lu_,g_,t0_]:=Module[{x,y},Return[\[Kappa]g[Lu,t0,g,{x,y}]];];



(* ::Section::Closed:: *)
(*Integration steps*)


(* ::Subsection::Closed:: *)
(*Integration steps along u and v*)


(* Integration along u and v *)

integrateU[g_,data_,propagationEq_,i12_,du_]:=
Module[{hasdu,datar,stepu,target,variables,dp,dq},
target=Length[i12]+2;
dp=stringlist["dp",Length[i12]];
dq=stringlist["dq",Length[i12]];
variables=Flatten[{{u,v},xy[target],nxy[target],\[Alpha],\[Beta],stringlist["\[Epsilon]",target-2],b,s,stringlist["p",target-2],stringlist["q",target-2]}];
hasdu={1,1,1,0,1,1,0,1,i12-1,0};
stepu=du*{{1,0},dur,dun,du\[Alpha],du\[Beta],du\[Epsilon],db,ds,dp,dq}/.propagationEq/.MapThread[#1->#2&,{variables,Flatten[data]}];
datar=hasdu*(data+stepu);
datar[[3]]=MapThread[#3/Sqrt[(#3.#1[#2].#3)]&,{g,datar[[2]],datar[[3]]}];
Return[datar];]

integrateV[g_,data_,propagationEq_,i12_,dv_]:=
Module[{hasdv,datar,stepv,target,variables,dp,dq},
target=Length[i12]+2;
dp=stringlist["dp",Length[i12]];
dq=stringlist["dq",Length[i12]];
variables=Flatten[{u,v,xy[target],nxy[target],\[Alpha],\[Beta],stringlist["\[Epsilon]",target-2],b,s,stringlist["p",target-2],stringlist["q",target-2]}];
hasdv={1,1,1,1,0,1,1,0,0,2-i12};
stepv=dv*{{0,1},dvr,dvn,dv\[Alpha],dv\[Beta],dv\[Epsilon],db,ds,dp,dq}/.propagationEq/.MapThread[#1->#2&,{variables,Flatten[data]}];
datar=hasdv*(data+stepv);
datar[[3]]=MapThread[#3/Sqrt[(#3.#1[#2].#3)]&,{g,datar[[2]],datar[[3]]}];
Return[datar];]


(* ::Subsection::Closed:: *)
(*Propagation equations along u and v*)


(* Propagation equations *)
propagationInverseEqns[surfacemetrics_,lambda1_,lambda2_,eps_,Lambda_]:=
Module[
{\[Lambda]1times,\[Lambda]2times,n\[CapitalGamma]n,np\[CapitalGamma]n,metrics,\[CapitalGamma]s,surface,r,position,targets,p,q,surfacecurvatures,directors,directorperpendicular,\[Kappa]gu,\[Kappa]gv,result},
(*Assume surfacemetrics={g1,g2,g3,..} s.t. gi={{gi11,gi12},{gi21,gi22}}[#]&*) 

targets=Length[Lambda];
p=stringlist["p",Length[eps]];
q=stringlist["q",Length[eps]];
position=xy[targets];
directors=nxy[targets];

metrics=Table[surfacemetrics[[surface]]@position[[surface]],{surface,targets}];

directorperpendicular=MapThread[fnp,{directors,metrics}];
\[CapitalGamma]s=Table[\[CapitalGamma][metrics[[surface]],position[[surface]],i,j,k],{surface,targets},{i,2},{j,2},{k,2}];

\[Kappa]gu=Table[1/lambda2[eps,Lambda[[surface]]] (b-Sum[D[Log[lambda1[eps,Lambda[[surface]]]],eps[[ep]]]q[[ep]],{ep,targets-2}]),{surface,targets}];
\[Kappa]gv=Table[1/lambda1[eps,Lambda[[surface]]] (s+Sum[D[Log[lambda2[eps,Lambda[[surface]]]],eps[[ep]]]p[[ep]],{ep,targets-2}]),{surface,targets}];

surfacecurvatures=Table[1/2 R[metrics[[surface]],position[[surface]]],{surface,targets}];
n\[CapitalGamma]n=Table[directors[[surface]].\[CapitalGamma]s[[surface,i]] . directors[[surface]],{surface,targets},{i,2}]//Simplify;
np\[CapitalGamma]n=Table[directorperpendicular[[surface]].\[CapitalGamma]s[[surface,i]] . directors[[surface]],{surface,targets},{i,2}]//Simplify;

\[Lambda]1times[f_]:=\[Lambda]times[f,lambda1,eps,Lambda];
\[Lambda]2times[f_]:=\[Lambda]times[f,lambda2,eps,Lambda];

result=
({
dur-> \[Alpha] \[Lambda]1times[directors],
dvr-> \[Beta] \[Lambda]2times[directorperpendicular],
dun->(\[Lambda]1times[MapThread[#1 *#2&,{\[Kappa]gu,directorperpendicular}]-n\[CapitalGamma]n])\[Alpha]  ,
dvn->(\[Lambda]2times[MapThread[#1 *#2&,{\[Kappa]gv,directorperpendicular}]-np\[CapitalGamma]n])\[Beta] ,
dv\[Alpha] -> -b \[Alpha] \[Beta] ,
du\[Beta] ->s \[Alpha] \[Beta] ,
du\[Epsilon] -> \[Alpha]  p,
dv\[Epsilon] -> \[Beta]  q}
~Join~
MapThread[#1->#2&,{dipq2[lambda1,lambda2,eps,Lambda],inverseEquations[surfacecurvatures,Lambda,eps,lambda1,lambda2]}]
);
Return[result];]


(* ::Subsection:: *)
(*Given Geometry Propagation equations along u and v*)


(* Propagation equations *)
propagationInverseEqnsGivenGeometry[surfacemetrics_,surfaceK_,surface\[CapitalGamma]_,lambda1_,lambda2_,eps_,Lambda_]:=
Module[
{\[Lambda]1times,\[Lambda]2times,n\[CapitalGamma]n,np\[CapitalGamma]n,metrics,\[CapitalGamma]s,surface,r,position,targets,p,q,directors,directorperpendicular,\[Kappa]gu,\[Kappa]gv,result,surfacecurvatures},
(*Assume surfacemetrics={g1,g2,g3,..} s.t. gi={{gi11,gi12},{gi21,gi22}}[#]&*) 
(*Assume surfaceK={K1,K2,K3,..} s.t. Ki=Ki[#]&*) 
(*Assume surface\[CapitalGamma]={\[CapitalGamma]1,\[CapitalGamma]2,\[CapitalGamma]3,..} s.t. \[CapitalGamma]i1[[i,j,k]]&*) 


targets=Length[Lambda];
p=stringlist["p",Length[eps]];
q=stringlist["q",Length[eps]];
position=xy[targets];
directors=nxy[targets];

metrics=Table[surfacemetrics[[surface]]@position[[surface]],{surface,targets}];

directorperpendicular=MapThread[fnp,{directors,metrics}];
(*\[CapitalGamma]s=Table[\[CapitalGamma][metrics[[surface]],position[[surface]],i,j,k],{surface,targets},{i,2},{j,2},{k,2}];*)

\[CapitalGamma]s=Table[surface\[CapitalGamma][[surface,i,j,k]]@@position[[surface]],{surface,targets},{i,2},{j,2},{k,2}];

\[Kappa]gu=Table[1/lambda2[eps,Lambda[[surface]]] (b-Sum[D[Log[lambda1[eps,Lambda[[surface]]]],eps[[ep]]]q[[ep]],{ep,targets-2}]),{surface,targets}];
\[Kappa]gv=Table[1/lambda1[eps,Lambda[[surface]]] (s+Sum[D[Log[lambda2[eps,Lambda[[surface]]]],eps[[ep]]]p[[ep]],{ep,targets-2}]),{surface,targets}];

surfacecurvatures=Table[ surfaceK[[surface]]@position[[surface]],{surface,targets}];

n\[CapitalGamma]n=Table[directors[[surface]].\[CapitalGamma]s[[surface,i]] . directors[[surface]],{surface,targets},{i,2}]//Simplify;
np\[CapitalGamma]n=Table[directorperpendicular[[surface]].\[CapitalGamma]s[[surface,i]] . directors[[surface]],{surface,targets},{i,2}]//Simplify;

\[Lambda]1times[f_]:=\[Lambda]times[f,lambda1,eps,Lambda];
\[Lambda]2times[f_]:=\[Lambda]times[f,lambda2,eps,Lambda];

result=
({
dur-> \[Alpha] \[Lambda]1times[directors],
dvr-> \[Beta] \[Lambda]2times[directorperpendicular],
dun->(\[Lambda]1times[MapThread[#1 *#2&,{\[Kappa]gu,directorperpendicular}]-n\[CapitalGamma]n])\[Alpha]  ,
dvn->(\[Lambda]2times[MapThread[#1 *#2&,{\[Kappa]gv,directorperpendicular}]-np\[CapitalGamma]n])\[Beta] ,
dv\[Alpha] -> -b \[Alpha] \[Beta] ,
du\[Beta] ->s \[Alpha] \[Beta] ,
du\[Epsilon] -> \[Alpha]  p,
dv\[Epsilon] -> \[Beta]  q}
~Join~
MapThread[#1->#2&,{dipq2[lambda1,lambda2,eps,Lambda],inverseEquations[surfacecurvatures,Lambda,eps,lambda1,lambda2]}]
);
Return[result];]


(* ::Subsection::Closed:: *)
(*Deriving residual {pi,dpi,...} {qi,dqi,...}*)


(* Deriving residual {pi,dpi,...} {qi,dqi,...} *)
pdp[\[Epsilon]list_,\[Alpha]_,du_]:=
Module[{p,dp},
p=stringlist["p",Length[\[Epsilon]list[[1]] ]];
dp=stringlist["p",Length[\[Epsilon]list[[1]] ]];
p=MapThread[If[#3==0,#1*0,(#1-#2)/(du #3)]&,{\[Epsilon]list[[2;;]],\[Epsilon]list[[;;-2]],\[Alpha][[;;-2]]}]~Join~{(0*\[Epsilon]list[[1]])};
dp=MapThread[If[#3==0,#1*0,(#1-#2)/(du )]&,{p[[2;;]],p[[;;-2]],\[Alpha][[;;-2]]}]~Join~{(0*\[Epsilon]list[[1]])};
Return[{p,dp}]]

qdq[\[Epsilon]list_,\[Beta]_,dv_]:=pdp[\[Epsilon]list,\[Beta],dv]


(* ::Subsection::Closed:: *)
(*Deriving inverse equations for highest order terms db, ds, dqi, dpi,..*)


(*Deriving inverse equations for highest order terms db, ds, dqi, dpi,..*)
(* dipq2[...]=inverseEquations[...]*)

inverseEquations::usage="propagation eqns for dipq2 (db,ds,dq1 or dp1,...)";
inverseEquations[surfacecurvatures_,Lambda_,eps_,lambda1_,lambda2_]:=Module[{i12,eqns,M2K,M2Kb},
i12=i12assign[lambda1,lambda2,eps,Lambda];
M2K=(Inverse[M2[eps,i12,Lambda,lambda1,lambda2,Length[Lambda]]]//.tobspq[eps]//.removeUV[eps]//Simplify).surfacecurvatures;
M2Kb=(Inverse[M2[eps,i12,Lambda,lambda1,lambda2,Length[Lambda]]].Kb[eps,i12,Lambda,lambda1,lambda2,Length[Lambda]]//.tobspq[eps]//.removeUV[eps]//Simplify);
(*eqns=(Simplify[Inverse[M2[eps,i12,Lambda,lambda1,lambda2,Length[Lambda]]]].(surfacecurvatures-Kb[eps,i12,Lambda,lambda1,lambda2,Length[Lambda]])//.tobspq[eps]//.removeUV[eps]);*)
eqns=M2K-M2Kb;
Return[eqns];]

dipq2::usage="(db,ds,dq)..., A modified version of di in the paper, for ease of integration. Such that M2.dipq2= M.(\!\(\*FractionBox[\(db\), \(\[Beta]\)]\),\!\(\*FractionBox[\(ds\), \(\[Alpha]\)]\),\!\(\*FractionBox[\(dq\), \(\[Beta]\)]\)) ..";
dipq2[lambda1_,lambda2_,ep_,Lambda_]:=Module[{i12,dp,dq,d},
i12=i12assign[lambda1,lambda2,ep,Lambda];
dp=stringlist["dp",Length[ep]];
dq=stringlist["dq",Length[ep]];
d=Join[{db,ds},Table[If[i12[[i]]==1,dq[[i]],dp[[i]]],{i,Length[ep]}]];
Return[d]];

duvFordipq[lambda1_,lambda2_,ep_,Lambda_]:=Module[{i12,d},
i12=i12assign[lambda1,lambda2,ep,Lambda];
dp=stringlist["dp",Length[ep]];
dq=stringlist["dq",Length[ep]];
d=Join[{dv,du},Table[If[i12[[i]]==1,dv ,du],{i,Length[ep]}]];
Return[d]];


(* ::Subsection::Closed:: *)
(*Internal functions for deriving inverse equations*)


(*Internal functions for deriving inverse equations *)
M[ep_,i12_,Lambda_,lambda1_,lambda2_,target_]:=
Module[{eps,m},
eps=toUV[ep];
m=Table[
Join[{1/lambda2[eps,Lambda[[surface]]]^2,-1/lambda1[eps,Lambda[[surface]]]^2},
Table[If[i12[[i]]==1,-D[Log[lambda1[eps,Lambda[[surface]]]],eps[[i]]]/lambda2[eps,Lambda[[surface]]]^2,-D[Log[lambda2[eps,Lambda[[surface]]]],eps[[i]]]/lambda1[eps,Lambda[[surface]]]^2],{i,Length[eps]}]],{surface,target}]//Simplify;
Return[m]];

M2::usage="A modified version of M in the paper, such that M2.(db,ds,dq)=M.(\!\(\*FractionBox[\(db\), \(\[Beta]\)]\),\!\(\*FractionBox[\(ds\), \(\[Alpha]\)]\),\!\(\*FractionBox[\(dq\), \(\[Beta]\)]\)) ..";
M2[ep_,i12_,Lambda_,lambda1_,lambda2_,target_]:=
Module[{eps,m},
eps=toUV[ep];
m=Table[Join[{1/(\[Beta][u,v]*lambda2[eps,Lambda[[surface]]]^2),-1/(\[Alpha][u,v]*lambda1[eps,Lambda[[surface]]]^2)},
Table[If[i12[[i]]==1,-D[Log[lambda1[eps,Lambda[[surface]]]],eps[[i]]]/(\[Beta][u,v]*lambda2[eps,Lambda[[surface]]]^2),-D[Log[lambda2[eps,Lambda[[surface]]]],eps[[i]]]/(\[Alpha][u,v]*lambda1[eps,Lambda[[surface]]]^2)],{i,Length[eps]}]],{surface,target}]//Simplify;
Return[m]];


di\[Epsilon][ep_,i12_]:=Module[{eps,d},eps=toUV[ep];
d=Join[{1/\[Beta][u,v] D[-D[\[Alpha][u,v],v]/(\[Alpha][u,v]\[Beta][u,v]),v],
1/\[Alpha][u,v] D[D[\[Beta][u,v],u]/(\[Alpha][u,v]\[Beta][u,v]),u]},
Table[
If[i12[[i]]==1,
	D[D[eps[[i]],v]/\[Beta][u,v],v]/\[Beta][u,v],
	D[D[eps[[i]],u]/\[Alpha][u,v],u]/\[Alpha][u,v]],
{i,Length[eps]}]];
Return[d]];

dipq[ep_,i12_]:=Module[{eps,p,q,pp,qq,d},
pp=toUV[stringlist["p",Length[ep]]];
qq=toUV[stringlist["q",Length[ep]]];
eps=toUV[ep];pp=toUV[p];
d=Join[{1/\[Beta][u,v] D[-D[\[Alpha][u,v],v]/(\[Alpha][u,v]\[Beta][u,v]),v],1/\[Alpha][u,v] D[D[\[Beta][u,v],u]/(\[Alpha][u,v]\[Beta][u,v]),u]},Table[If[i12[[i]]==1,D[qq[[i]],v]/\[Beta][u,v],D[pp[[i]],u]/\[Alpha][u,v]],{i,Length[eps]}]];
Return[d]];

gaussiancurvature[ep_,Lambda_,lambda1_,lambda2_,target_]:=
Module[{eps,curvature},
eps=toUV[ep];
curvature=Table[1/2 R[DiagonalMatrix[{\[Alpha][u,v]^2lambda1[eps,Lambda[[surface]]]^2,\[Beta][u,v]^2lambda2[eps,Lambda[[surface]]]^2}],{u,v}],{surface,target}];
Return[curvature]]

Kb[eps_,i12_,Lambda_,lambda1_,lambda2_,target_]:=Simplify[gaussiancurvature[eps,Lambda,lambda1,lambda2,target]-M[eps,i12,Lambda,lambda1,lambda2,target].di\[Epsilon][eps,i12]]



(* ::Section:: *)
(*Common functions*)


(* ::Subsection:: *)
(*Geometrical functions*)


(*Geometrical functions *)
\[CapitalGamma][g_,x_,i_,j_,k_](*\[CapitalGamma]^i_jk*):=
FullSimplify[(1/2)Sum[Inverse[g][[i,n]] (D[g[[n,j]],x[[k]]]+D[g[[n,k]],x[[j]]]-D[g[[j,k]],x[[n]]]),{n,Length[x]}]];
Riemann[g_,x_,i_,j_,k_,l_](*R^i_jkl*):=
Simplify[D[\[CapitalGamma][g,x,i,j,l],x[[k]]]-D[\[CapitalGamma][g,x,i,j,k],x[[l]]]+Sum[\[CapitalGamma][g,x,i,m,k]\[CapitalGamma][g,x,m,j,l]-\[CapitalGamma][g,x,i,m,l]\[CapitalGamma][g,x,m,j,k],{m,Length[x]}]];

Ricci[g_,x_,i_,j_](*R_ij*):=Simplify[Sum[Riemann[g,x,k,i,k,j],{k,Length[x]}]];

R[g_,x_](*R*):=Sum[Inverse[g][[i,j]]Ricci[g,x,i,j],{i,Length[x]},{j,Length[x]}];

\[CapitalGamma][g_,x_,i_,j_,k_,x0_](*\[CapitalGamma]^i_jk(x_0)*):=
ReplaceAll[\[CapitalGamma][g,x,i,j,k],MapThread[#1->#2&,{x,x0}]];
(*geodesic curvature \[Kappa]g*)
\[Kappa]g[\[Gamma]_,t0_,g_,x_](*\[Kappa]_g(\!\(\(\[Gamma]\((t)\)\)
\*SubscriptBox[\(|\), \(t0\)]\);g,x)*):=FullSimplify[ReplaceAll[Sqrt[Det[g[\[Gamma]]]]/(D[\[Gamma],t].g[\[Gamma]].D[\[Gamma],t])^(3/2)
( \[CapitalGamma][g[x],x,2,1,1,\[Gamma]]D[\[Gamma],t][[1]]^3
+(2\[CapitalGamma][g[x],x,2,1,2,\[Gamma]]-\[CapitalGamma][g[x],x,1,1,1,\[Gamma]])D[\[Gamma],t][[1]]^2 D[\[Gamma],t][[2]]
+(\[CapitalGamma][g[x],x,2,2,2,\[Gamma]]-2\[CapitalGamma][g[x],x,1,1,2,\[Gamma]])D[\[Gamma],t][[2]]^2 D[\[Gamma],t][[1]]
-\[CapitalGamma][g[x],x,1,2,2,\[Gamma]]D[\[Gamma],t][[2]]^3
+D[\[Gamma],{t,2}][[2]] D[\[Gamma],t][[1]]-D[\[Gamma],{t,2}][[1]] D[\[Gamma],t][[2]]),t->t0]];


(*fpn= find Subscript[n, \[Perpendicular]] given n *)
fnp[n_,g_]:=1/Sqrt[Det[g]] {-(g.n)[[2]],(g.n)[[1]]};


(* ::Subsection::Closed:: *)
(*Replacement functions*)


(* Replacement functions *)
toUV[x_]:=ToExpression[ToString[#]<>"[u,v]"]&/@x;

tobspq[ep_]:=Module[{eps,pp,qq,result},
eps=toUV[ep];pp=toUV[stringlist["p",Length[ep]]];qq=toUV[stringlist["q",Length[ep]]];result=(D[{MapThread[#1->#2&,{D[eps,u],\[Alpha][u,v]pp}],D[\[Beta][u,v],u]->\[Alpha][u,v]\[Beta][u,v] s[u,v]},{u,#}]&/@{0,1})~Join~(D[{MapThread[#1->#2&,{D[eps,v],\[Beta][u,v]qq}],D[\[Alpha][u,v],v]->-\[Alpha][u,v]\[Beta][u,v] b[u,v]},{v,#}]&/@{0,1})//Flatten;Return[result]];

removeUV[eps_]:=Module[{epsuv,ppuv,qquv,pp,qq,dp,dq,result},
pp=stringlist["p",Length[eps]];
qq=stringlist["q",Length[eps]];
dp=stringlist["dp",Length[eps]];
dq=stringlist["dq",Length[eps]];
epsuv=toUV[eps];
ppuv=toUV[pp];
qquv=toUV[qq];
result=
Flatten[{{\[Alpha][u,v]->\[Alpha],\[Beta][u,v]->\[Beta],b[u,v]->b,s[u,v]->s,D[b[u,v],v]->db,D[s[u,v],u]->ds},MapThread[#1->#2&,{epsuv,eps}],MapThread[#1->#2&,{ppuv,pp}],MapThread[#1->#2&,{qquv,qq}],MapThread[#1->#2&,{D[qquv,v],dq}],MapThread[#1->#2&,{D[ppuv,u],dp}]}];
Return[result];];


(* ::Subsection::Closed:: *)
(*List making and naming functions*)


(*List making and naming functions*)
(* used to create lists such as {{x1,y1},{x2,y2},...}, data*)
datastruct[eps_]:=0{{u,v}, xy[Length[eps]+2],nxy[Length[eps]+2],\[Alpha],\[Beta],eps,b,s,eps,eps}
datatable[eps_,maxu_,maxv_,du_,dv_]:=Table[datastruct[eps],{i,Ceiling[maxu/du]},{j,Ceiling[maxv/dv]}]

xy[number_]:=ToExpression["{x"<>ToString[#]<>",y"<>ToString[#]<>"}"]&/@Range[number]
nxy[number_]:=ToExpression["{nx"<>ToString[#]<>",ny"<>ToString[#]<>"}"]&/@Range[number]
stringlist[string_,length_]:=ToExpression[string<>ToString[#]]&/@Range[length]


(* ::Subsection::Closed:: *)
(*Organizing functions*)


(* Organizing functions *)
i12assign[lambda1_,lambda2_,eps_,Lambda_]:=Table[If[#,1,2]&[ Nand@@Table[0===FullSimplify[D[lambda1[eps,Lambda[[i]]],eps[[j]]]],{i,Length[Lambda]}]],{j,Length[eps]}];

q\[Lambda][lambda1_,eps_,Lambda_]:=Table[Nand@@Table[0===FullSimplify[D[lambda1[eps,Lambda[[i]]],eps[[j]]]],{i,Length[Lambda]}],{j,Length[Lambda]-2}];

CheckSolverValidity::usage="Is the uniaxial sheet's inverse design described by initial value or boundary problems";
CheckSolverValidity[lambda1_,lambda2_,eps_,Lambda_]:=If[AnyTrue[MapThread[#1&&#2&,{q\[Lambda][lambda1,eps,Lambda],q\[Lambda][lambda2,eps,Lambda]}],#&],"Elliptic equations - search for a different solver",
Print["Go ahead"];
Table[If[q\[Lambda][lambda1,eps,Lambda][[i]],Print[ToString[eps[[i]]]<>" and q"<>ToString[i]<>" initial conditions along u"],Print[ToString[eps[[i]]]<>" and p"<>ToString[i]<>"initial conditions along v"]],{i,Length[eps]}];];


(* ::Subsection:: *)
(*General functions*)


(* General functions *)
\[Lambda]times[f_,lambdai_,eps_,Lambda_]:=If[Length[f[[1]]]>1,MapThread[#1.#2&,{IdentityMatrix[2]lambdai[eps,#]&/@Lambda,f}],MapThread[#1*#2&,{lambdai[eps,#]&/@Lambda,f}]];
normalizevec[v_,g_]:=v/Sqrt[v.g.v]

(* embedding to metric: (x,y,z)[xi,eta] \[Rule] g[xi,eta] *)
exprToFunction[expr_,vars_]:=Function[Evaluate[expr/.Thread[vars->Array[Slot,Length[vars]]]]]
embeddingToMetric[embedding_,coordinates_]:=exprToFunction[Table[D[embedding,i].D[embedding,j],{i,coordinates},{j,coordinates}],coordinates];


(* ::Section:: *)
(*Plotting*)


tubeEmbed::usage="Plot u-curves as tubes from data for a given embedding {x,y,z}[\[Theta],\[Phi]]";
tubeEmbed[data_,umaxtable_,vmax_,du_,dv_,vjump_,surfacenum_,embedding_,variable_]:=Table[
Graphics3D[{CapForm["Butt"],GrayLevel[.8],Specularity[White,1000],
Tube[exprToFunction[embedding,variable]@@@data[[;;Floor[umaxtable[[vv]]/du],
vv,2,surfacenum]]]},Boxed->False,Lighting->"Neutral"],{vv,1,vmax/dv,vjump}];

tubeEmbedColor[data_,umaxtable_,vmax_,du_,dv_,vjump_,surfacenum_,embedding_,variable_]:=Table[
Graphics3D[{CapForm["Butt"],
Tube[exprToFunction[embedding,variable]@@@data[[;;Floor[umaxtable[[vv]]/du],
vv,2,surfacenum]]]},Boxed->False],{vv,1,vmax/dv,vjump}];

