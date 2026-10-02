#!/usr/bin/env wolframscript
(* Exact four-sample arrow and fixed-arrow checks E1-E7; see ../README.md. *)
u[th_, ph_] := 9/20 + Cos[th]/4 + Cos[2 th + ph]/5;
bern[q_, 1] := q;
bern[q_, 0] := 1 - q;
qsLaw[ph_, taus_List] := Association @ Table[
    x -> Integrate[TrigReduce[Expand[Product[bern[u[th + taus[[k]], ph], x[[k]]], {k, Length[taus]}]]],
        {th, 0, 2 Pi}]/(2 Pi),
    {x, Tuples[{0, 1}, Length[taus]]}];
grid = {0, Pi/3, 2 Pi/3, Pi};
P4 = qsLaw[Pi/3, grid];
R4 = Association @ Table[x -> P4[Reverse[x]], {x, Tuples[{0, 1}, 4]}];
P3 = qsLaw[Pi/3, {0, Pi/3, Pi}];
Q3 = qsLaw[Pi/3, {0, 2 Pi/3, Pi}];
R3 = Association @ Table[x -> Q3[Reverse[x]], {x, Tuples[{0, 1}, 3]}];
marg[L_, x_] := Sum[L[{x[[1]], x[[2]], y, x[[3]]}], {y, 0, 1}];

e1a = AllTrue[Values[P4], MatchQ[#, _Rational | _Integer] && # > 0 &] && Total[Values[P4]] == 1;
e1b = AllTrue[Tuples[{0, 1}, 3], marg[P4, #] == P3[#] &];
e1c = AllTrue[Tuples[{0, 1}, 3], marg[R4, #] == R3[#] &];
P40 = qsLaw[0, grid];
e1d = AllTrue[Tuples[{0, 1}, 4], P40[#] == P40[Reverse[#]] &];
e2 = {P4[{0, 0, 0, 0}], R4[{0, 0, 0, 0}], P4[{0, 0, 1, 0}], R4[{0, 0, 1, 0}]} == {50681, 50681, 42799, 45799}/640000;
fibres = GatherBy[Tuples[{0, 1}, 4], #[[{1, 2, 4}]] &];
ratioSets = Map[Union[Map[P4[#]/R4[#] &, #]] &, fibres];
e3 = AnyTrue[ratioSets, Length[#] > 1 &];
c3 = Sum[P3[x] Log[P3[x]/R3[x]], {x, Tuples[{0, 1}, 3]}];
c4 = Sum[P4[x] Log[P4[x]/R4[x]], {x, Tuples[{0, 1}, 4]}];
sealed = 0.000738887716627336479919968717778`33;
c3n = N[c3, 50];
c4n = N[c4, 50];
e4rel = Abs[c3n - sealed]/sealed;
e4 = e4rel <= 10^-20;
e5pos = c4n - c3n > 0;

d = (2 - Sqrt[3])/8;
gridMin = Min[Table[N[u[th, ph]], {th, 0., 2 Pi, 2 Pi/2000}, {ph, 0., Pi/3, Pi/120}]];
gridMax = Max[Table[N[u[th, ph]], {th, 0., 2 Pi, 2 Pi/2000}, {ph, 0., Pi/3, Pi/120}]];
nmin = NMinimize[{u[th, ph], 0 <= ph <= Pi/3, 0 <= th <= 2 Pi}, {th, ph}, WorkingPrecision -> 30];
nmax = NMaximize[{u[th, ph], 0 <= ph <= Pi/3, 0 <= th <= 2 Pi}, {th, ph}, WorkingPrecision -> 30];
e6 = gridMin >= N[d] && gridMax <= 9/10 && First[nmin] >= N[d] && First[nmax] <= 9/10;
cad = TimeConstrained[
    FindInstance[c^2 + s^2 == 1 && cc^2 + ss^2 == 1 && cc >= 1/2 && ss >= 0 &&
        9/20 + c/4 + ((2 c^2 - 1) cc - 2 c s ss)/5 < d, {c, s, cc, ss}, Reals],
    60, $TimedOut];

e7sym = Simplify[(10/dd) (13/20)^2 (2 Pi/T)^2 (T/2) - 169 Pi^2/(20 dd T)] === 0;
run2 = Import[FileNameJoin[{DirectoryName[$InputFileName], "..", "results", "EP1B_RUN2.json"}], "RawJSON"];
toNum[v_] := If[StringQ[v], ToExpression[v], v];
e7rows = Map[{#["T"], N[169 Pi^2/(20 d #["T"]), 20], toNum[#["L_sigma"]]} &, run2["BA"]];
e7 = e7sym && AllTrue[e7rows, #[[2]] >= #[[3]] &];

Print["EP1ADD1 version ", $VersionNumber];
Print["E1 ", {e1a, e1b, e1c, e1d}];
Print["E2 ", e2];
Print["E3 ", e3, " nonconstant fibres ", Select[Transpose[{fibres[[All, 1, {1, 2, 4}]], ratioSets}], Length[#[[2]]] > 1 &]];
Print["P4 ", Normal[P4]];
Print["P3 ", Normal[P3]];
Print["R3 ", Normal[R3]];
Print["c3 ", c3n];
Print["c4 ", c4n];
Print["c4-c3 ", N[c4 - c3, 30], " positive ", e5pos];
Print["E4 relerr ", N[e4rel, 5], " pass ", e4];
Print["E6 gridMin ", gridMin, " gridMax ", gridMax, " delta ", N[d, 20], " nmin ", nmin, " nmax ", nmax, " pass ", e6];
Print["E6 CAD violation search ", cad];
Print["E7 ", e7sym, " rows ", e7rows, " pass ", e7];
Print["EP1ADD1 ALL ", And[e1a, e1b, e1c, e1d, e2, e3, e4, e5pos, e6, e7]];
