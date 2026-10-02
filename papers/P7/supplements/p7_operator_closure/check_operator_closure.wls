(* Exact closeout controls, not an all-rank proof by finite extrapolation. *)
ClearAll["Global`*"];
checks = {};
need[name_String, value_] := If[TrueQ[value], AppendTo[checks, name],
  Print["P7_OPERATOR_CHECK_FAILED:" <> name]; Exit[1]];
eq[a_, b_] := TrueQ[FullSimplify[a == b]];
zeroMatrix[a_] := And @@ (TrueQ[Simplify[#] == 0] & /@ Flatten[a]);

rot = {{Sqrt[t], -Sqrt[1-t]}, {Sqrt[1-t], Sqrt[t]}};
need["full_normal_symbol_orthogonal", FullSimplify[
  Transpose[rot].rot == IdentityMatrix[2], 0 < t < 1]];
hom = {{Sqrt[t], -s Sqrt[1-t]}, {Sqrt[1-t], Sqrt[t]}};
need["homotopy_determinant", eq[Det[hom], t+s (1-t)]];
need["homotopy_positive", FullSimplify[t+s (1-t)>0, 0<t<1 && 0<=s<=1]];
tjr = r (r+1)/(j (j+1));
ce2 = 2 (j+1)/j; cc2 = 2 (j-1)/j;
need["sparse_a_column", eq[r (r+1)/(2 (j+1)^2) ce2, tjr]];
need["sparse_b_column", eq[(j-r)(j+r+1)/(2 (j+1)^2) ce2, 1-tjr]];
need["sparse_c_column", eq[r (r+1)/(2 (j^2-1)) cc2, tjr]];
need["full_raw_block_det_square", eq[1/(ce2 cc2), j^2/(4 (j^2-1))]];
need["source_rank", eq[(j+1)^2+j^2-1, 2j(j+1)]];
need["target_rank", eq[j(j+2)+j^2, 2j(j+1)]];
need["endpoint_multiplicity", eq[1+2j+1, 2j+2]];
need["intermediate_multiplicity", eq[Sum[2r+1,{r,1,j-1}],j^2-1]];
sparse2=(2/(2(2+1)))^(2*2+2) (1^2*2^2/(4*3^2*(2^2-1)))^3;
full2=(2/(2*3))^9 (2/(2*1))^3;
need["sparse_j2", sparse2==1/14348907];
need["full_j2", full2==1/19683];
need["conflation_detector", sparse2 != full2];
need["j1_raw_vs_clifford", (1/4)^4==1/256 && 4^4 (1/4)^4==1];

weightSum[n_] := n (n+1)(n+2)/3;
need["general_weight_sum", eq[Sum[(n-2k)^2,{k,0,n}],weightSum[n]]];
density=Expand[(j weightSum[j] (p1-2 e)/4+
 (j+1) weightSum[j-1] (p1+2 e)/4)/2-j(j+1)p1/24];
need["general_density", eq[density,j(j+1)((j^2+j-1)p1-(2j+1)e)/12]];
need["j1_density", eq[density/.{j->1,p1->3 tau,e->chi},(tau-chi)/2]];
need["chirality_detector", !eq[-density/.j->1,(p1-3e)/6]];
Do[need["central_descent_"<>ToString[n],
  And@@(EvenQ[Total[#]]& /@ {{n,n},{n,n-2},{n+1,n-1},{n-1,n-1}}) &&
  OddQ[n+n-1]], {n,1,8}];

(* Direct Hermitian spin matrices: independent of the polynomial kernel route. *)
spin[n_Integer] := Module[{a,z},
 a=Table[If[i==k-1,Sqrt[(k-1)(n-k+2)],0],{i,n+1},{k,n+1}];
 z=DiagonalMatrix[Table[n/2-k,{k,0,n}]];
 {(a+Transpose[a])/2,(a-Transpose[a])/(2 I),z}];
Do[Module[{ss,ll,rr,dd,tt,dim,rhs},
 ss=KroneckerProduct[#,IdentityMatrix[p+1],IdentityMatrix[q+1]]& /@ spin[1];
 ll=KroneckerProduct[IdentityMatrix[2],#,IdentityMatrix[q+1]]& /@ spin[p];
 rr=KroneckerProduct[IdentityMatrix[2],IdentityMatrix[p+1],#]& /@ spin[q];
 dd=2 Sum[ss[[k]].(rr[[k]]-ll[[k]]),{k,3}];
 tt=ss+ll+rr; dim=2(p+1)(q+1);
 rhs=(p(p+2)/2+q(q+2)/2+3/4) IdentityMatrix[dim]-Sum[tt[[k]].tt[[k]],{k,3}];
 need["matrix_square_"<>ToString[p]<>"_"<>ToString[q],zeroMatrix[dd.dd-rhs]];
 need["matrix_selfadjoint_"<>ToString[p]<>"_"<>ToString[q],zeroMatrix[dd-ConjugateTranspose[dd]]];
],{p,0,4},{q,0,4}];
cas[x_]:=x(x+2)/4;
need["positive_extremal_branch",eq[(cas[n+2r]-cas[n])/r,n+r+1]];
need["finite_signed_branch",eq[(cas[2r-2-m]-cas[m])/r,r-1-m]];
need["zero_branch",eq[2(cas[r-1]+cas[r-1])+1-r^2,0]];

Print[ExportString[<|"schema"->"paper7.operator_closure.wolfram.v1",
 "result"->"P7_OPERATOR_CLOSURE_EXACT_CONTROLS_PASS",
 "wolfram_version"->$Version,"checks"->checks,"count"->Length[checks],
 "sparse_j2"->ToString[sparse2,InputForm],"full_j2"->ToString[full2,InputForm],
 "scope"->"symbol, weights and finite angular controls; analytic proof in manuscript"|>,"RawJSON"]];
Exit[0];
