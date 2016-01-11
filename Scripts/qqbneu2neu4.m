(* ::Package:: *)

(*  	qqbneu1neu1.m
		generates the Fortran code for
		q q-bar -> chi chi in the MSSM
		last modified Oct 2014

Note: the QED contributions are not taken into account. To plug
the QED part back in, comment out the parts in DiagramSelect that
eliminate a V[1|3]...
*)


Clear["Global`*"]
SetDirectory[NotebookDirectory[]];
<< FeynArts` \.1d
<< FormCalc`
ClearProcess[]

time1 = SessionTime[]


CalcProcess = "uubar";

(*Process definitions*)
If[CalcProcess === "uubar",
	process = {F[3, {1}], -F[3, {1}]} -> {F[11,{2}], F[11,{4}]};
	name = CalcProcess;
	gen1 = 1;
	,
	(*else*)
	If[CalcProcess === "ddbar",
		process = {F[4, {1}], -F[4, {1}]} -> {F[11,{2}], F[11,{4}]};
		name = CalcProcess;
		gen1 = 1;
		,
		(*else*)
		Print["No valid input."]
	]
]


(*Neglect Masses (URL)*)
Neglect[ME] = Neglect[ME2] = 0;
(*Neglect[MQU] = Neglect[MQD] = 0;*)
Neglect[MU] = Neglect[MU2] = 0;
Neglect[MC] = Neglect[MC2] = 0;
(*Neglect[MT] = Neglect[MT2] = 0;*)
Neglect[MD] = Neglect[MD2] = 0;
Neglect[MS] = Neglect[MS2] = 0;
(*Neglect[MB] = Neglect[MB2] = 0;*)

(*Diagonale CKM Matrix*)
CKM = IndexDelta;
CKMC = IndexDelta;


(*Options*)
SetOptions[InsertFields, Model -> "MSSMCT",
           Restrictions -> {NoLightFHCoupling}(*No Fermion-Higgs coupling*),
		   ExcludeParticles -> {S[1|2|3|4|5|6|11|12], F[1|2], V[1|3]} 
			(*Exclude Top, Higgs, Neutrinos, massive Leptons, Sneutrinos, Sleptons, Photon, W-Boson*)];

SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}, AutoEdit -> False];

(*Reduce tensor to scalar integrals and choose regularisazation scheme*)
(*D = dimensional regularization (default),*)
(*4 = constrained differential renormalization,*)
(*0 = keeps the whole amplitude D-dimensional*)
SetOptions[CalcFeynAmp,Dimension->D];
(*Note: There is currently a bug in FormCalc which does not allow to compile*)
(*the generated code with PaVeReduce\[Rule]True set.*)
(*One has to replace "Derivative(1)(IGram)(MS2)" with "(-1/(MS2**2))"*)

(*Save the Diagrams*)
$PaintSE = MkDir["Diagrams"];
DoPaint[diags_, type_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, name <> "." <> type <> ".eps"], #]&)]

(*SUSY Counterterm*)
GSY = GS;
ELY = EL(*+dZe1y*);


Print["Born"]

tops = CreateTopologies[0, 2 -> 2];
ins = InsertFields[tops, process];
DoPaint[ins, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]];
born = born//.{Alfa2->0}


Print["Counter Terms"]

top = CreateCTTopologies[1, 2 -> 2, ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs}]; (*Exclude Tadpole- and Self-Energy CT on external legs*)
ins = InsertFields[top, process];

(*All Counter Terms*)
DoPaint[ins, "counterAll"];
counterAll = CreateFeynAmp[ins];

(*Self Energy Counter Terms*)
insSelf = DiagramSelect[ins,(FieldMemberQ[FieldPoints[##],FieldPoint[_][_, _]]) &];
insSelf = DiagramDelete[insSelf, 1];(*Delete CT for EW Corrections*)
DoPaint[insSelf, "counterSelf"];
counterSelf = CreateFeynAmp[insSelf];

(*Vertex Correction Counter Terms*)
insVert = DiagramSelect[ins,FieldMemberQ[FieldPoints[##],FieldPoint[1][_, _, _]] &];
insVert = DiagramDelete[insVert, 1];
DoPaint[insVert, "counterVert"];
counterVert = CreateFeynAmp[insVert];


Print["Self Energies"]

top = CreateTopologies[1, 2 -> 2, SelfEnergiesOnly];
ins = InsertFields[top, process];
ins = DiagramSelect[ins,(FieldMemberQ[FieldPoints[##],FieldPoint[_][_,_,V[5,{_}]]] || FieldMemberQ[LoopFields[##],F[15,{_}]]
						|| FieldMemberQ[LoopFields[##], S[13,{_,_,_}]] || FieldMemberQ[LoopFields[##], S[14,{_,_,_}]] || FieldMemberQ[LoopFields[##], V[5,{_}]]) &]; (*Diagram contains g, sg*)


ins = DiagramSelect[ins,(FreeQ[FieldPoints[##],FieldPoint[_][_,_,V[2]]]) &]; (*Diagram does not contain Z*)

ins = DiagramSelect[ins,FreeQ[LoopFields[##], V[2]] &]; (*Loop does not contain Z*)
ins = DiagramSelect[ins,FreeQ[FieldPoints[##],FieldPoint[_][_,_,V[2],V[2]]] &]; (*Diagram does not contain 2Z*)
ins = DiagramSelect[ins,FreeQ[FieldPoints[##],FieldPoint[_][F[gen1,{_,_}],-F[gen1,{_,_}],V[2]]] &]; (*Diagram does not contain 2Z*)

DoPaint[ins, "self"];
self = CalcFeynAmp[CreateFeynAmp[ins], counterSelf];
self = self//.{Alfa2->0}


Print["Vertices"]

top = CreateTopologies[1, 2 -> 2, TrianglesOnly];
ins = InsertFields[top, process];

(*All Vertices*)
(*DoPaint[ins, "vertAll"];
vertAll = CreateFeynAmp[ins];*)

ins = DiagramSelect[ins,(MemberQ[LoopFields[##], V[5,{_}]] || MemberQ[LoopFields[##], F[15,{_}]]) &];(*Loop contains a g or sg*)

DoPaint[ins, "vert"];
vert = CalcFeynAmp[CreateFeynAmp[ins], counterVert];
vert = vert//.{Alfa2->0}


Print["Boxes"]

top = CreateTopologies[1, 2 -> 2, BoxesOnly];
ins = InsertFields[top, process];

ins = DiagramSelect[ins,(MemberQ[LoopFields[##], V[5,{_}]] || MemberQ[LoopFields[##], F[15,{_}]] || 
				MemberQ[FieldPoints[##],FieldPoint[_][_,_,V[5,{_}]]] ) &];(*Loop contains a g or sg or Loop couples via gluon*)

DoPaint[ins, "box"];
box = CalcFeynAmp[CreateFeynAmp[ins]];
box = box//.{Alfa2->0}


amps = {born, self, vert, box};
{born, self, vert, box} = Abbreviate[amps, 6, Preprocess -> OnSize[100, Simplify, 500, DenCollect]];

col = ColourME[All, born];

abbr = OptimizeAbbr[Abbr[]]
subexpr = OptimizeAbbr[Subexpr[]]//.{Alfa2->0}


(*Calculate the Renormalization constants*)
CalcRenConstRestrict = 
 ExcludeFieldPoints -> {
	FieldPoint[_][S[1|2|3|4|5|6], _, _],(*Exclude Higgs*)
	FieldPoint[_][S[1|2|3|4|5|6], _, _, _],
	FieldPoint[_][V[1|2|3], _, _],(*Exclude Photon, Z, W*)
	FieldPoint[_][V[1|2|3], _, _, _],
	FieldPoint[_][F[11], _, _],(*Exclude Neutralino*)
	FieldPoint[_][F[11], _, _, _],
	FieldPoint[_][F[12], _, _],(*Exclude Chargino*)
	FieldPoint[_][F[12], _, _, _]
};

(*Optionen von InsertFields \[UDoubleDot]berschreiben*)
SetOptions[InsertFields, Model -> "MSSMCT",
		   Restrictions -> {NoLightFHCoupling, CalcRenConstRestrict},
		   ExcludeParticles -> {S[1|2|3|4|5|6|11|12], F[1|2], V[1|2|3], F[11|12], U[1|2|3|4|5]} 
];

dir = SetupCodeDir[name <> ".fortran", Drivers -> name <> ".drivers"];
WriteSquaredME[born, {self, vert, box}, col, abbr, subexpr, dir];

(*neue Syntax durch ge\[ADoubleDot]nderte Datei "FormCalc84.m" und "FormCalc.m" im FormCalc Installationsordner.*)
(*Nun ist es m\[ODoubleDot]glich die RenConst. zu berechnen ohne, dass sofort der Fortran Code erzeugt wird.*)
(*Auf die berechneten RenConst. k\[ODoubleDot]nnen weitere Operationen angewendet werden. Danach werden*)
(*die RenConst. mit der Funktion "WriteRenConst[]" in Fortran Code umgewandelt.*)
(*
renConst = CalcRenConst[amps];(*ge\[ADoubleDot]nderte Funktion, berechnet die RenConst, wandelt sie aber nicht sofort in Fortran Code um.*)
renConst = renConst//.{Alfa -> 0}//FullSimplify (*consider just strong corrections*)
Export["renConst.wdx",renConst,"WDX"];
WriteRenConst[renConst,dir];(*Neue Funktion! Hinzugef\[UDoubleDot]gt in "FormCalc84.m" und "FormCalc.m"*)
*)


Print["time used: ", SessionTime[] - time1]
Exit[];
