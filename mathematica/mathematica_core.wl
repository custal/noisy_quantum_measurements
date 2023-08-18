(* ::Package:: *)

(*Trine*)

phi1Trine={{1},{0}};
phi2Trine={{1/2},{Sqrt[3]/2}};
phi3Trine={{1/2},{-Sqrt[3]/2}};

(*Defining the matrices M0,M1,and M2*)
povmsTrine={2/3*Matrixify[phi1Trine],
2/3*Matrixify[phi2Trine],
2/3*Matrixify[phi3Trine]};

outcomesTrine = PovmPermutations[povmsTrine];

(*Tetra*)

phi1Tetra={1,0};
phi2Tetra={Rationalize[1/Sqrt[3]],Rationalize[Sqrt[2/3]]};
phi3Tetra={Rationalize[1/Sqrt[3]],E^(I*2*Pi/3)*Rationalize[Sqrt[2/3]]};
phi4Tetra={Rationalize[1/Sqrt[3]],E^(-I*2*Pi/3)*Rationalize[Sqrt[2/3]]};

povmsTetra = { 1/2*Matrixify[phi1Tetra],
		1/2*Matrixify[phi2Tetra],
		1/2*Matrixify[phi3Tetra],
		1/2*Matrixify[phi4Tetra]};

outcomesTetra = PovmPermutations[povmsTetra];


MeasureAbsolute[x_, y_] := ((1 - x)^2 + y^2)/(1 - x + y)

OptimalBasisValues[eigVals_] := 
 Module[{x, 
	   y, \[Lambda]min = Min[eigVals], \[Lambda]max = Max[eigVals]},
	  {x, y} = Piecewise[{
	     {{\[Lambda]max, (Sqrt[2] - 1)*(1 - \[Lambda]max)}, \[Lambda]max <=
	        1 - (Sqrt[2] + 1)*\[Lambda]min}, 
	     {{\[Lambda]max, \[Lambda]min}, 
	      1 - (Sqrt[2] + 1)*\[Lambda]min < \[Lambda]max <= 
	       1 - (Sqrt[2] - 1)*\[Lambda]min},
	     {{1 - (Sqrt[2] - 1)*\[Lambda]min, \[Lambda]min}, 
	      1 - (Sqrt[2] - 1)*\[Lambda]min < \[Lambda]max}
	     }];
	  {x, y}
	  ]
  
PovmPermutations[povms_] := 
 Module[{outcomes = <||>}, 
	  Do[outcomes[{i, j}] = KroneckerProduct[povms[[i]], povms[[j]]];
	   , {i, Length[povms]}, {j, Length[povms]}]; 
	  outcomes]

Matrixify[v_] := Outer[Times, v, ConjugateTranspose[v]]

GenerateProtocols[outcomes_] := 
 Module[{permutations, protocols},
  
	  permutations = Keys[outcomes];
	  protocols = <|{} -> Null,
	    permutations -> Null|>;
	  
	  Do[
	   protocols[selection] = Null;
	   protocols[Complement[permutations, selection]] = Null;
	   , {selection, Subsets[permutations, Floor[Length[outcomes]/2]]}];
	  protocols]
  
SumPovms[protocol_, povmOutcomes_] := Module[{sum},
(*Use the first povm outcome assuming all povms have the same dimension*)
	  sum = ConstantArray[0, Dimensions[Values[povmOutcomes][[1]]]];
	  
	  Do[sum += povmOutcomes[outcome]
	   , {outcome, protocol}];
	  sum]

CalculateMeasures[protocol_, povmOutcomes_] := 
 Module[{Q, eigenValues, eigenVectors, x, y, measureValue, result},
	  
	  Q = SumPovms[protocol, povmOutcomes];
	  {eigenValues, eigenVectors} = Chop[N[Eigensystem[Q]]];
	  {x, y} = OptimalBasisValues[eigenValues];
	  measureValue = MeasureAbsolute[x, y];
	 
	  result = <|"protocol" -> protocol, "eigenValues" -> eigenValues, 
	    "eigenVectors" -> eigenVectors, "measureValue" -> measureValue, 
	    "measureFunction" -> "MeasureAbsolute", "xStarYStar" -> {x, y}|>;
	  result
	  ]
  
CalculateAllProtocols[outcomes_] := 
 Module[{protocols, protocol, result, count,protocolNum},
	  protocols = GenerateProtocols[outcomes];
	  
	  count=0;
	  protocolNum=Length[protocols];
	  Monitor[Do[
	    result = CalculateMeasures[protocol, outcomes];
	    protocols[protocol] = result;
	    count+=1;
	    , {protocol, Keys[protocols]}], 
	    Row[{ProgressIndicator[count,{1,protocolNum}],StringForm["``\[LessEqual]``",count,protocolNum]}," "]];
	  Print[StringForm["Calculated `` results!",protocolNum]];
	  protocols]
  
FindLowestMeasure[protocols_]:=Module[{lowestMeasure=Infinity, lowestMeasureProtocols},
	
	Do[
		If[protocol["measureValue"] <= lowestMeasure, 
		If[protocol["measureValue"] < lowestMeasure, Module[{}, lowestMeasureProtocols={protocol};
													 lowestMeasure=protocol["measureValue"];], 
			lowestMeasureProtocols=Append[lowestMeasureProtocols, protocol]]]
		, {protocol, protocols}];
	lowestMeasureProtocols]

BraKetExpectation[braket_, operator_]:=ConjugateTranspose[braket] . operator . braket

PartialTrace[\[Rho]_, system_, dimensions_] := Module[{operatorList, basis, operatorBasis, operator, systemDims, sum, traceOperator},
	operatorList = IdentityMatrix /@ dimensions;
	systemDims = dimensions[[system]];
	
	sum = ConstantArray[0, Dimensions[\[Rho]]/systemDims];
	Do[operatorBasis = operatorList;
		operatorBasis[[system]] = ConjugateTranspose[{UnitVector[systemDims, basis]}] (*This is given as a row vector so need to conjugate transpose*);
		traceOperator = KroneckerProduct @@ operatorBasis;
		sum += BraKetExpectation[traceOperator, \[Rho]],
	{basis, Range[dimensions[[system]]]}];
	sum]

EntropyOfEntanglement[state_, trace_:{1}, dimensions_:{2,2}] := Module[{\[Rho], system, \[Lambda]},
	\[Rho] = Chop[Matrixify[Normalize[state]]];
	Do[\[Rho]=Chop[PartialTrace[\[Rho], system, dimensions]], {system, trace}];
	\[Lambda] = DeleteCases[Chop[Eigenvalues[\[Rho]]], 0]; (*Only sum non-zero eigen values*)
	-Chop[Total[\[Lambda]*Log[\[Lambda]]]]]

PlotProtocol[protocol_, dimensions_] := Module[{grid, outcome}, 
	grid = ConstantArray[0, dimensions];
	
	Do[
		grid[[outcome[[1]],outcome[[2]]]]=1,
		{outcome, protocol}];
		MatrixPlot[grid, ColorFunction->"Monochrome"]]
	
ProbabilityGrid[state_, outcomes_] := Module[{grid, protocol},

	grid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];

	Do[
		grid[[protocol[[1]],protocol[[2]]]] = Chop[BraKetExpectation[Normalize[state], SumPovms[{protocol}, outcomes]]], {protocol, Keys[outcomes]}];
	grid
	]
		
ProbabilityRatioGrid[state1_, state2_, outcomes_] := Module[{grid, protocol, measureOperator, probability1, probability2},
	grid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];

	Do[
		measureOperator = SumPovms[{protocol}, outcomes];
		probability1 = Chop[BraKetExpectation[Normalize[state1], measureOperator]];
		probability2 = Chop[BraKetExpectation[Normalize[state2], measureOperator]];
		grid[[protocol[[1]],protocol[[2]]]] = probability1/probability2, {protocol, Keys[outcomes]}];

	grid
	]
		
PlotProtocolVerbose[protocol_, outcomes_, ratio_:False] := Module[{mat, protocolGrid, probabilityGrid},
(*Probability grid is calcualted assuming we are in the boring case where the |Subscript[\[Psi], 0]> = |Subscript[\[Lambda], max]> and |Subscript[\[Psi], 1]> = |Subscript[\[Lambda], min]>*)

	protocolGrid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];
	If[
		ratio, 
		probabilityGrid = ProbabilityRatioGrid[protocol["eigenVectors"][[1]], protocol["eigenVectors"][[4]], outcomes],
		Module[{},probabilityGrid = ProbabilityGrid[protocol["eigenVectors"][[1]], outcomes];
		probabilityGrid = probabilityGrid/Max[probabilityGrid]]
	];
	
	Do[
		protocolGrid[[measure[[1]], measure[[2]]]]=1,
		{measure, protocol["protocol"]}];

	MatrixPlot[probabilityGrid, Mesh -> All,ColorFunction->"LightTemperatureMap",ColorFunctionScaling->False,
	      Epilog -> MapIndexed[Text[Style[Round[#, .01], Large], {#2[[2]], Length[protocolGrid]+1-#2[[1]]} - .5] &, protocolGrid, {2}]]]
	      
KeyToList[protocols_, keyFunc_]:= Module[{list},
	list = {};
	Do[list=Append[list,keyFunc[out]],{out,protocols}];
	list]
	
SplitAndSplice[state_]:=KroneckerProduct[PartialTrace[state, 2, {2,2}], PartialTrace[state, 1, {2,2}]]/Tr[state]

NearTo[A_, B_, round_:-10]:=(Round[#1, 10^(round)] & /@ A)== (Round[#1, 10^(round)] & /@ B)

CheckSeparability[protocol_] := Module[{Q, Qsplice, vect, productStateTruths, \[Rho], \[Rho]splice},
	Q= Chop[N[SumPovms[protocol["protocol"], outcomesTetra]]];
	Qsplice = SplitAndSplice[Q];
	productStateTruths = {NearTo[Q, Qsplice]};

	Do[
	\[Rho] = Chop[N[Matrixify[Normalize[vect]]]];
	\[Rho]splice = Chop[SplitAndSplice[\[Rho]]];
	AppendTo[productStateTruths, NearTo[\[Rho], \[Rho]splice, -9]], {vect, protocol["eigenVectors"]}
];

	productStateTruths
]
